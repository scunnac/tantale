

.tale_parts_from_file <- function(fasta) {
  if (grepl("TALE_Protein_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readAAStringSet(fasta)
  if (grepl("TALE_DNA_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readDNAStringSet(fasta)
  if (length(taleStrings) == 0L) {
    logger::log_warn("No part sequence found in: {fasta}. Returning an empty tibble.")
    tbl <- tibble::tibble(arrayIDs = character(),
                   domainType = character(),
                   positionInArray = character(),
                   positionInCrd = character(),
                   string = character(),
                   sourceDirectory = character())
    return(tbl)
  } 
  tbl <- tibble::tibble(arrayID = gsub("(.*): .*", "\\1", names(taleStrings)),
                        domainType = gsub(".*: (.*?)[ ]?[0-9]{0,}$", "\\1", names(taleStrings)),
                        positionInCrd = gsub(".*: repeat[ ]([0-9]{0,})$", "\\1", names(taleStrings)) %>%
                          as.integer() %>%
                          suppressWarnings(),
                        string = as.character(taleStrings) %>% as.vector(),
                        sourceDirectory = dirname(fasta)
  )
  missingTerm <- setdiff(c("N-terminus", "C-terminus"), unique(tbl$domainType))
  if (length(missingTerm) != 0L) {
    logger::log_warn("Array {unique(tbl$arrayID)} is missing a {missingTerm} domain in {fasta}")
    warning()
    missingTerm <- tibble::tibble(arrayID = unique(tbl$arrayID),
                                  domainType = missingTerm,
                                  positionInCrd = NA,
                                  string = NA,
                                  sourceDirectory = dirname(fasta))
    tbl <- dplyr::bind_rows(tbl, missingTerm)
    logger::skip_formatter(as.character(knitr::kable(missingTerm))) %>%
      logger::log_debug()
  }
  tbl %<>% 
    dplyr::rowwise() %>%
    dplyr::mutate(positionInArray = switch(domainType,
                                                 `N-terminus` = 1,
                                                 `repeat` = positionInCrd + 1,
                                                 `C-terminus` = nrow(tbl)
                                                 ))
  return(tbl)
}




.rvds_from_annotale_file <- function(fasta) {
  if (!grepl("TALE_RVDs.fasta", basename(fasta))) {
    logger::log_error("The provided file does not seem to be an AnnoTALE RVDs file: {fasta}")
    stop()
  } else {
    rvdTble <- .split_list(fasta) %>%
      lapply(function(x) tibble::tibble(string = x,
                                        positionInCrd = 1:length(x))
             ) %>%
      dplyr::bind_rows(.id = "arrayID")
    rvdTble <- rvdTble %>% dplyr::mutate(sourceDirectory = dirname(fasta),
                                         domainType = "repeat")
  }
  return(rvdTble)
}



#' Read TALE parts from a tell_tales output directory
#'
#' Implementation behind \code{\link{tales_from_telltale}} and the deprecated
#' \code{\link{tales_from_telltale}}. Internal so that package code can call it without
#' tripping the deprecation warning.
#' @noRd
.tale_parts <- function(telltale_dir) {
  # Get info from telltale output dir
  # !!!! arrayID are assumed to be unique !!!!
  protPartsFiles <- list.files(telltale_dir, "TALE_Protein_parts.fasta", recursive = T, full.names = T)
  dnaPartsFiles <- list.files(telltale_dir, "TALE_DNA_parts.fasta", recursive = T, full.names = T)
  if (telltale_dir %>% dirname() %>% unique() %>% length() != 1L) {
    log_error("The provided path most likely does not correspond to a SINGLE tell_tales output directory.")
  }
  # Fetch info from annotale/telltale files with .tale_parts_from_file
  taleProtString <- lapply(protPartsFiles, .tale_parts_from_file) %>% dplyr::bind_rows()
  taleDnaString <- lapply(dnaPartsFiles, .tale_parts_from_file) %>% dplyr::bind_rows()
  #stopifnot(nrow(taleProtString) == nrow(taleDnaString))
  # Join info in a table with one domain per row
  tale_parts <- dplyr::full_join(taleDnaString %>% dplyr::rename(dnaSeq = string),
                                taleProtString %>% dplyr::rename(aaSeq = string),
                                by = c("arrayID", "domainType", "positionInArray", "positionInCrd", "sourceDirectory"),
                                relationship = "one-to-one") %>%
    dplyr::mutate(aaSeq = gsub("[*]", "", aaSeq))

  # Get RVDs
  # NOTE: could be easier to get the RVDs directly from AnnoTALE output with
  # .rvds_from_annotale_file() but I currently feel that it is good to
  # be aware of disagreements between AnnoTALE diagnostic on terminal domains presence in AA seqs
  # and nhmmer diagnostic on terminal domains CDS presence on DNA.
  rvds <- .split_list(list.files(path = telltale_dir,
                                        pattern = "rvdSequences.fas",
                                        recursive = F,
                                        full.names = T)
                             ) %>%
    lapply(function(x) tibble::tibble(rvd = x, positionInArray = 1:length(x) )) %>%
    dplyr::bind_rows(.id = "arrayID")  
  anchorCodes <- c("NTERM", "CTERM", "XXXXX")
  
  
  # Some checks on the consistency between parts and rvd sequences
  # if nhmmer did not report on a C-Term CDS, the corresponding domain
  # "CTERM" tag will not be written in the rvd slot of the table.
  arraysConsistency <- dplyr::full_join(tale_parts %>% dplyr::count(arrayID, name = "AnnoTALELength"),
                                        rvds %>% dplyr::count(arrayID, name = "rvdFileLength"),
                                        by = dplyr::join_by(arrayID)) %>%
    dplyr::mutate(sameLength = AnnoTALELength == rvdFileLength)
  
  if (any(is.na(arraysConsistency$sameLength))) {
    logger::log_warn("There are mismatches in array IDs between rvd seq file and AnnoTALE parts files:")
    logger::skip_formatter(as.character(knitr::kable(arraysConsistency %>% dplyr::filter(is.na(sameLength))))) %>%
      logger::log_warn()
    warning()
  } else if (!all(arraysConsistency$sameLength, na.rm = TRUE)) {
    logger::log_error("Array lengths are inconsistent between rvd ",
                      "seq file and AnnoTALE parts files:")
    logger::skip_formatter(as.character(knitr::kable(arraysConsistency %>% dplyr::filter(!sameLength)))) %>%
      logger::log_error()
    stop()
  }
  
  # Include RVDs in the talParts tibble
  tale_parts <- dplyr::left_join(tale_parts, 
                                rvds,
                                by = c("arrayID", "positionInArray"),
                                unmatched = "drop", relationship = "one-to-one")
  # Include seqnames in the talParts tibble
  tale_parts %<>% dplyr::left_join(
    readr::read_tsv(list.files(telltale_dir, "hitsReport.tsv", recursive = T, full.names = T),
                    show_col_types = FALSE) %>%
      dplyr::select(arrayID, seqnames) %>%
      dplyr::distinct(),
    by = "arrayID", relationship = "many-to-one"
  )
  # Check talparts
  partsWithMissingAaSeq <- tale_parts %>% dplyr::filter(is.na(aaSeq)) %>% dplyr::pull(arrayID) %>% unique()
  partsWithMissingDnaSeq <- tale_parts %>% dplyr::filter(is.na(dnaSeq)) %>% dplyr::pull(arrayID) %>% unique()
  partsWithMissingRvdSeq <- tale_parts %>% dplyr::filter(is.na(rvd)) %>% dplyr::pull(arrayID) %>% unique()
  if (any(sapply(list(partsWithMissingAaSeq, partsWithMissingDnaSeq, partsWithMissingRvdSeq), length) != 0L)) {
    logger::log_warn("Be aware that the output tale_parts tibble has records with missing sequences:")
    logger::log_warn("{unique(c(partsWithMissingAaSeq, partsWithMissingDnaSeq, partsWithMissingRvdSeq))}")
    tale_parts %>%
      #dplyr::select(arrayID, domainType, positionInArray, sourceDirectory) %>%
      dplyr::filter(arrayID %in% c(partsWithMissingAaSeq, partsWithMissingDnaSeq, partsWithMissingRvdSeq)) %>%
      knitr::kable() %>% as.character() %>%
      logger::skip_formatter() %>% logger::log_debug()
    warning()
  }
  return(tale_parts)
}

identSubMat <- matrix(data = rep(0, times = length(Biostrings::AA_PROTEINOGENIC) ^ 2),
                      nrow = length(Biostrings::AA_PROTEINOGENIC),
                      ncol = length(Biostrings::AA_PROTEINOGENIC),
                      dimnames = list(Biostrings::AA_PROTEINOGENIC, Biostrings::AA_PROTEINOGENIC))
diag(identSubMat) <- 1


.pairwise_align_biostrings <- function(part_aa_set, ncores = 1) {
  bpparam <- BiocParallel::MulticoreParam(ncores, progressbar = TRUE)
  
  if (anyDuplicated(names(part_aa_set))) {
    stop("Parts in the provided input have duplicated names. Cannot proceeed...")
  }
  
  pair_align_scores <- BiocParallel::bplapply(
    X = 1:length(part_aa_set),
    FUN = function(i) {
      singleSubAln <- pwalign::pairwiseAlignment(pattern = part_aa_set,
                                                    subject = part_aa_set[i],
                                                    substitutionMatrix = identSubMat, #"BLOSUM62",
                                                    gapOpening = 1, gapExtension = 0.5,
                                                    type = "global", scoreOnly = FALSE)
      tibble::tibble(subj = names(part_aa_set[i]),
                     pattern = names(pwalign::alignedPattern(singleSubAln)),
                     score = BiocGenerics::score(singleSubAln),
                     nedit = Biostrings::nmismatch(singleSubAln))
    },
    BPPARAM = bpparam) %>%
    dplyr::bind_rows()
  
  pair_align_scores %<>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      maxLength = max(nchar(part_aa_set[subj]), nchar(part_aa_set[pattern])),
      # This is an approximate equivalent of how Alvaro computed dissimilarity in distal
      Dissim = 100 - 100 * (maxLength - score) / maxLength,
      Dissim = ifelse(Dissim < 0, 100, 100 - Dissim),
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-maxLength)
  
  # Check the pair_align_scores tibble
  .check_pair_align_tbl(pair_align_scores = pair_align_scores, part_aa_set = part_aa_set)
  
  return(pair_align_scores)
}


.pairwise_align_mmseq2 <- function(part_aa_set, ncores = 1, conda_bin = "auto") {
  outdir <- tempfile(pattern = "distalPairwiseAlign2")
  dir.exists(outdir) || dir.create(outdir, recursive = TRUE)
  partAaStringSetFile <- file.path(outdir, "taleAsParts.fsa")
  mmseq2DbPath <- file.path(outdir, 'mmseq2DB')
  prefDbPath <- file.path(outdir, 'resultDB_pref')
  alnDbPath <- file.path(outdir, 'resultDB_aln')
  alnTabFile <- file.path(outdir, 'alnRes.tab')
  
  df <- expand.grid(names(part_aa_set),
                    names(part_aa_set),
                    stringsAsFactors = FALSE
                    ) %>%
    tibble::as_tibble()
  colnames(df) <- c("query", "target")
  if(anyDuplicated(df) != 0) stop("The provided sequences must have unique names.")
  
  Biostrings::writeXStringSet(part_aa_set, filepath = partAaStringSetFile)
  
  mmseq2createdb <- glue::glue("mmseqs createdb {partAaStringSetFile} {mmseq2DbPath}")

  mmseq2prefilter <- glue::glue("mmseqs prefilter {mmseq2DbPath} {mmseq2DbPath} {prefDbPath}",
                               "-v 3 --threads {max(floor(ncores/2), 1)} --max-seqs 1000 -s 7.5 --add-self-matches 1",
                               "--cov-mode 0", .sep = " ")
  
  mmseq2align <- glue::glue("mmseqs align {mmseq2DbPath} {mmseq2DbPath} {prefDbPath} {alnDbPath}",
                               "-v 3 --threads {ncores} --add-self-matches 1 --min-seq-id 0",
                               "--cov-mode 0 --gap-open aa:11,nucl:5 --gap-extend aa:1,nucl:2",
                               "-a 1 --alignment-mode 3 --alignment-output-mode 0 --seq-id-mode 1",
                               .sep = " ")
  
  mmseq2convertalis <- glue::glue("mmseqs convertalis {mmseq2DbPath} {mmseq2DbPath} {alnDbPath} {alnTabFile}",
                               "--format-mode 4 -v 3",
                               "--format-output query,target,evalue,raw,pident,nident,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qcov,tcov",
                               .sep = " ")
  
  if (!as.logical(.create_tantale_env(conda_bin = conda_bin))) {
    logger::log_debug("Invoking mmseq2 using the following command:\n {stringr::str_wrap(mmseq2Cmd, 80)}")
    res <- .run_in_conda(env_name = "tantale",
                            conda_bin = conda_bin,
                            command = mmseq2createdb)
    res <- .run_in_conda(env_name = "tantale",
                            conda_bin = conda_bin,
                            command = mmseq2prefilter)
    res <- .run_in_conda(env_name = "tantale",
                            conda_bin = conda_bin,
                            command = mmseq2align)
    res <- .run_in_conda(env_name = "tantale",
                            conda_bin = conda_bin,
                            command = mmseq2convertalis)
  } else {
    stop("Could not create the tantale conda environment on your machine to run mmseq2...")
  }
  

  pair_align_scores <- readr::read_tsv(alnTabFile, show_col_types = FALSE) %>%
    dplyr::mutate(target = as.character(target), query = as.character(query)) %>%
    dplyr::group_by(target, query) %>%
    dplyr::slice_max(raw, n = 1, with_ties = FALSE)
  pair_align_scores <- dplyr::left_join(df, pair_align_scores) %>%
    dplyr::rename(subj = target, pattern = query) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Dissim = ifelse(is.na(pident), 100, 100 - pident*min(qcov,tcov))) %>%
    dplyr::ungroup()

  # Check the pair_align_scores tibble
  .check_pair_align_tbl(pair_align_scores = pair_align_scores, part_aa_set = part_aa_set)
  unlink(outdir, recursive = TRUE)
  return(pair_align_scores)
}




.pairwise_align_decipher <- function(part_aa_set, ncores = 1) {
  if (anyDuplicated(names(part_aa_set))) {
    stop("Parts in the provided input have duplicated names. Cannot proceeed...")
  }
  msa <- DECIPHER::AlignSeqs(myXStringSet = part_aa_set, normPower = 0,
                             processors = ncores, verbose = FALSE)
  # When no guide tree is supplied, DECIPHER::StaggerAlignment() builds one
  # internally using DistanceMatrix(..., correction = "TN93+F"), a nucleotide
  # substitution model that errors out on an AAStringSet. Build the guide
  # tree ourselves with a correction-free (protein-safe) distance matrix.
  staggerTree <- if (length(msa) >= 3) {
    distForTree <- DECIPHER::DistanceMatrix(msa, processors = ncores, verbose = FALSE)
    suppressWarnings(DECIPHER::Treeline(myDistMatrix = distForTree, method = "NJ",
                                        processors = ncores, verbose = FALSE))
  } else {
    NULL
  }
  msa <- DECIPHER::StaggerAlignment(msa, tree = staggerTree, fullLength = TRUE,
                                    processors = ncores, verbose = FALSE)
  distMat <- DECIPHER::DistanceMatrix(msa, method = "longest",
                                      includeTerminalGaps = TRUE,
                                      processors = ncores, verbose = FALSE)
  pair_align_scores <- reshape2::melt(as.matrix(distMat)) %>% tibble::as_tibble()
  colnames(pair_align_scores) <- c("pattern", "subj", "Dissim")
  pair_align_scores %<>% dplyr::mutate(pattern = as.character(pattern), subj = as.character(subj))
  pair_align_scores %<>% dplyr::mutate(Dissim = Dissim*100) %>%
    dplyr::ungroup()
  
  # Check the pair_align_scores tibble
  .check_pair_align_tbl(pair_align_scores = pair_align_scores, part_aa_set = part_aa_set)
  
  return(pair_align_scores)
}





.check_pair_align_tbl <- function(pair_align_scores, part_aa_set) {
  partCombinCounts <- pair_align_scores %>% dplyr::select(pattern, subj) %>%
    dplyr::count(pattern, subj) %>%
    dplyr::pull(n)
  if (!all(partCombinCounts == 1L)) {
    stop("Some alignment pairs have more than one record...",)
  }
  if (length(names(part_aa_set))^2 != nrow(pair_align_scores)) {
    stop("Some parts pairs are absent from the pairwise parts distance table")
  }
}

#' Report on potential 'pseudo TALEs' in a tale_parts object
#' @description
#' NOT TESTED!!!!
#' This displays a compact but information rich view of the TALEs stored in a
#' tale_parts object.
#' 
#' @param tale_parts a table of TALE parts as returned by the
#' \code{\link{tales_from_telltale}} function or
#' \code{\link{tales_compare}}
#' @param sanitize If \code{FALSE}, will return all the arrays with at least one 
#' part with a missing sequence. If \code{TRUE}, will return all the arrays that have
#' no part with a missing sequence.
#' 
#'
#' @return a tale_parts object
#' @export
diagnose_tale_parts <- function(tale_parts, sanitize = FALSE) {
  # Check talparts
  partsWithMissingAaSeq <- tale_parts %>% dplyr::filter(is.na(aaSeq)) %>%
    dplyr::select(arrayID, sourceDirectory) %>%
    dplyr::distinct()
  partsWithMissingDnaSeq <- tale_parts %>% dplyr::filter(is.na(dnaSeq)) %>%
    dplyr::select(arrayID, sourceDirectory) %>%
    dplyr::distinct()
  partsWithMissingRvdSeq <- tale_parts %>% dplyr::filter(is.na(rvd)) %>%
    dplyr::select(arrayID, sourceDirectory) %>%
    dplyr::distinct()
  problems <- dplyr::bind_rows(partsWithMissingRvdSeq,
                               partsWithMissingDnaSeq,
                               partsWithMissingAaSeq
                               ) %>%
    dplyr::distinct()
  pseudoTales <- dplyr::left_join(problems, tale_parts,
                                  relationship = "one-to-many",
                                  by = dplyr::join_by(arrayID, sourceDirectory)
                                  ) %>%
    dplyr::arrange(sourceDirectory, arrayID, positionInArray)
  if (nrow(problems) != 0L) {
    logger::log_warn("Be aware that the output tale_parts tibble has records with missing sequences")
    warning()
  }
  if (!sanitize) {
    pseudoTales %>% return()
  } else {
    logger::log_info("Returning TALE arrays with no empty sequence parts")
    dplyr::setdiff(tale_parts, pseudoTales) %>% return()
  }
}




#' Visualize TALE content in a tale_parts object
#' @description
#' NOT TESTED!!!!
#' This displays a compact but information rich view of the TALEs stored in a
#' tale_parts object.
#' 
#' @param tale_parts a table of TALE parts as returned by the
#' \code{\link{tales_from_telltale}} function or
#' \code{\link{tales_compare}}
#'
#' @return The ggplot object
#' @export
plot_tale_composition <- function(tale_parts) {
  partsForPlots <- tale_parts %>%
    mutate(label = if_else(domainType == "repeat", rvd, ""),
           aaSeqLength = factor(nchar(aaSeq))
    )
  p  <- partsForPlots %>% ggplot(mapping = aes(fill = aaSeqLength,
                                               color = domainType,
                                               label = label,
                                               y = arrayID,
                                               x = positionInArray),
                                 color = isNaAaSeq) +
    scale_color_viridis_d(option = "rocket") +
    scale_fill_discrete() +
    scale_x_continuous(breaks = 1:50, minor_breaks = NULL) +
    geom_point(shape = 21, size = 5, stroke = 0.9) +
    ggnewscale::new_scale_color() +
    ggnewscale::new_scale_fill() +
    geom_text(size = 2.1, color = "white") + 
    facet_grid(seqnames~ ., scales = "free_y", space = "free") +
    labs(title = "Overview of TALE composition by genome") +
    theme_light()
  print(p)
  return(p)
}


# tale_parts <- readRDS("/home/cunnac/TEMP/talePartsForDistalr.rds")
# h_cut = 10
# ncores = 1
# aln_method = "DECIPHER"
# conda_bin = "/home/cunnac/bin/miniconda3/condabin/conda"



#' Compute TALE and repeat relatedness
#'
#' Quantifies how TALE arrays, and the individual repeat units they are built
#' from, relate to one another. An R re-implementation of the original DisTAL
#' Perl program: it still uses the ARLEM binary for the repeat-array alignment
#' step, but performs the rest with R support and parallelization, which makes
#' it much faster (the exact speedup depends on \code{aln_method}).
#'
#' Two products are irreducible and expensive — the pairwise protein alignment
#' between repeat units, and ARLEM on the coded arrays. Everything else the
#' former \code{tales_compare()} returned was a projection of its inputs, so this
#' function returns only what cannot be recomputed cheaply.
#'
#' This is where \code{dom_code} is minted, over the whole set of parts
#' supplied, and where the resulting objects are stamped with a namespace
#' identifying that set — see \code{\link{tales_namespace}}. Passing a subset
#' later is safe; re-running on a different part set mints different codes,
#' and the differing namespace is what stops the two being joined by mistake.
#'
#' @param x A \code{\link{tales}} object whose parts carry amino acid
#'   sequences.
#' @param ncores Number of cores for the pairwise alignment step.
#' @param aln_method Approach for pairwise similarities between part amino acid
#'   sequences: \code{"DECIPHER"} (default), \code{"Biostrings"} or
#'   \code{"mmseq2"}.
#' @param conda_bin Path to a Conda binary, if \code{reticulate} cannot find it.
#'
#' @return A list of three objects, all describing the same run:
#' \itemize{
#'   \item \code{tales}: the input, with a \code{dom_code} column added.
#'   \item \code{domain_distances}: a \code{\link{domain_distances}} between repeat units,
#'     keyed by \code{dom_code}.
#'   \item \code{tale_distances}: a \code{\link{tale_distances}} between whole arrays,
#'     keyed by \code{array_id}.
#' }
#'
#' @references
#' Pérez-Quintero A.L. et al. (2015). QueTAL: a suite of tools to classify and
#' compare TAL effectors functionally and phylogenetically.
#' \emph{Frontiers in Plant Science} \strong{6}, 545.
#' \doi{10.3389/fpls.2015.00545}
#'
#' Abouelhoda M.I., Giegerich R., Behzadi B., Steyaert J.-M. (2009). Alignment
#' of minisatellite maps based on run-length encoding scheme.
#' \emph{Journal of Bioinformatics and Computational Biology} \strong{7}(2),
#' 287--308. \doi{10.1142/S0219720009004060}
#'
#' @seealso \code{\link{tales_group}} to cluster arrays from the returned
#'   \code{tale_distances}.
#' @export
tales_compare <- function(x, ncores = 1, aln_method = "DECIPHER",
                              conda_bin = "auto") {
  if (!is_tales(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  if (!"aa_seq" %in% names(x)) {
    cli::cli_abort(
      c("{.arg x} must carry an {.field aa_seq} column.",
        "i" = "Repeat similarity is computed from part amino acid sequences."),
      class = c("tantale_error_compare_no_aa", "tantale_error")
    )
  }
  if ("dom_code" %in% names(x)) {
    cli::cli_warn(
      c("{.arg x} already carries a {.field dom_code} column; it will be re-minted.",
        "i" = "Similarity tables from the earlier run are keyed by the old codes and must not be reused with this result.",
        "i" = "Their differing {.fn tales_namespace} is what will catch such a mix-up."),
      class = c("tantale_warning_relatedness_remint", "tantale_warning")
    )
    x <- x[setdiff(names(x), "dom_code")]
  }

  namespace <- .tales_dom_code_namespace(x$aa_seq)

  # The core still speaks the legacy column vocabulary; translate either side.
  legacy <- .tales_to_legacy(x)
  core <- .tales_compare_core(tale_parts = legacy, ncores = ncores,
                                  aln_method = aln_method, conda_bin = conda_bin)

  list(
    tales = tales(core$tale_parts, dom_code_namespace = namespace),
    domain_distances = domain_distances(
      core$dissim_long %>% dplyr::rename(RepU1 = subj, RepU2 = pattern),
      dom_code_namespace = namespace
    ),
    # Keyed by array_id, not dom_code, so deliberately unstamped: array ids are
    # meaningful names that do not silently collide across runs the way
    # cur_group_id() codes do (class-design.md §3.5).
    tale_distances = tale_distances(core$tal_sim)
  )
}

#' Rename a tales object's columns back to the legacy vocabulary
#'
#' Temporary bridge, the mirror of \code{.tales_rename_legacy()}. Removable once
#' the column sweep of restructuring-notes.md §9.2 lands.
#' @noRd
.tales_to_legacy <- function(x) {
  x <- tibble::as_tibble(x)
  hit <- intersect(names(x), unname(TALES_LEGACY_NAMES))
  if (length(hit) == 0L) return(x)
  back <- stats::setNames(names(TALES_LEGACY_NAMES), unname(TALES_LEGACY_NAMES))
  names(x)[match(hit, names(x))] <- unname(back[hit])
  x
}


#' The expensive part of the relatedness computation
#'
#' Shared by \code{\link{tales_compare}} and the deprecated
#' \code{\link{tales_compare}}. Speaks the legacy column vocabulary and returns raw
#' pieces; classing, stamping and assembly happen in the callers. Deliberately
#' does no clustering: that was a stored field with no consumers, recomputed by
#' its only would-be user at a different cut height (restructuring-notes.md §1).
#' @noRd
.tales_compare_core <- function(tale_parts, ncores = 1,
                                    aln_method = "DECIPHER", conda_bin = "auto") {
  
  #### Reality checks ####
  
  ## Make sure we are dealing only with parts that have defined protein sequences.
  if (any(is.na(tale_parts$aaSeq) | tale_parts$aaSeq == "")) {
    logger::log_error("It seems that some of the provided TALE parts miss an amino acid sequence. Cannot proceed!")
    tale_parts %>% dplyr::filter(is.na(aaSeq)) %>% 
      knitr::kable() %>% as.character() %>%
      logger::skip_formatter() %>% logger::log_error()
    stop()
  }
  if (any(is.na(tale_parts$dnaSeq) | tale_parts$dnaSeq == "")) {
    logger::log_warn("It seems that some of the provided TALE parts miss the DNA sequence!")
  } 
  
  ## Make sure that arrayID - position combinations are unique
  # in case someone would not have made arrayIDs unique before
  # mixing tale predictions from several genomes...
  arayPosCombinCounts <- tale_parts %>%
    dplyr::group_by(arrayID, positionInArray) %>%
    dplyr::count() %>%
    dplyr::pull(n)
  if (!all(arayPosCombinCounts == 1L)) {
    logger::log_error("Your tale arrays identifers are probably not unique.",
         "\n",
         "Make sure that there is only one part per position per arrayID.")
    stop()
  }
  
  # Assign domain codes
  tale_parts %<>% dplyr::group_by(aaSeq) %>%
    dplyr::mutate(domCode = dplyr::cur_group_id() %>% unlist() %>% as.character()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(domCode = dplyr::if_else(is.na(aaSeq), as.character(NA), domCode))
  
  
  #### Assemble repeat code strings and write in a file for arlem ####
  logger::log_info("Assemble repeat code TALE strings and write in a file for ARLEM")
  
  codesSeqsfile <- tempfile(fileext = ".fasta")
  repeatStrings <- tale_parts %>%
    dplyr::group_by(arrayID) %>%
    dplyr::arrange(positionInArray) %>%
    dplyr::summarise(repeatString = paste(domCode, collapse = " "),
                     posString = paste(positionInArray, collapse = " "))
  codesSeqSet <- Biostrings::BStringSet(repeatStrings$repeatString)
  names(codesSeqSet) <- repeatStrings$arrayID
  
  # # Must use seqinr because Biostrings wraps sequences in fasta file which messes up Arlem...
  # codesSeqLst <- as.list(repeatStrings$repeatString)
  # names(codesSeqLst) <- repeatStrings$arrayID
  # seqinr::write.fasta(codesSeqLst, names = names(codesSeqLst),
  #                     file.out = codesSeqsfile, as.string = TRUE, nbchar = 10000)
  
  Biostrings::writeXStringSet(x = codesSeqSet,
                              filepath = codesSeqsfile,
                              format = "fasta", width = 20000L)
  
  
  
  #### Compute systematic pairwise dissimilarities (distances) between 'repeat' units. ####
  # Get unique domains sequences
  taleAaParts <- Biostrings::AAStringSet(tale_parts$aaSeq)
  names(taleAaParts) <- tale_parts$domCode
  uniqueTaleAaParts <- unique(taleAaParts)
  stopifnot(!anyDuplicated(names(uniqueTaleAaParts)))
  stopifnot(!anyDuplicated(names(unique(taleAaParts))))
  
  # Get pairwise repeat aa sequence dissimilarity scores in a long tibble
  logger::log_info("Computing a distance matrix between TALE parts amino acid sequences ",
                   "using: {aln_method}")
  if (aln_method == "mmseq2") {
    dissimLong <- .pairwise_align_mmseq2(part_aa_set = uniqueTaleAaParts, ncores = ncores,
                                        conda_bin = conda_bin)
    #saveRDS(dissimLong, file = "/home/cunnac/TEMP/dissimLong")
  } else if (aln_method == "Biostrings") {
    dissimLong <- .pairwise_align_biostrings(part_aa_set = uniqueTaleAaParts, ncores = ncores)
  } else if (aln_method == "DECIPHER") {
    dissimLong <- .pairwise_align_decipher(part_aa_set = uniqueTaleAaParts, ncores = ncores)
  } else {
    logger::log_errors() && stop("'aln_method' parameter must be either 'Biostrings', 'mmseq2' or 'DECIPHER'")
  }
  dissimLong %<>% dplyr::mutate(Sim = 100 - Dissim)
  # Convert Distance (dissimilarity) measures to Similarity with a four-parameter logistic function
  # pair_align_scores %<>% dplyr::mutate(Sim = 100/(1+exp(-1*-0.9*(Dissim-3))))
  
  
  # Convert to square matrix
  dissimMat <- reshape2::acast(dissimLong, formula = subj ~ pattern, value.var = "Dissim")
  stopifnot(nrow(dissimMat) == ncol(dissimMat))
  # reorder row and colnames because I suspect arlem expect them in increasing order
  dissimMat <- dissimMat[rownames(dissimMat) %>% as.numeric() %>% order(),
                         colnames(dissimMat) %>% as.numeric() %>% order()]
  
  #### Generate an ARLEM cost matrix ####
  if (TRUE) {
    method <- "minkowski"
    logger::log_info("Generate an ARLEM cost matrix which meets triangle inequality criteria by computing ",
                     "the {method} distance between pairwise distance vectors.")
    dissimMat <- as.matrix(stats::dist(dissimMat, method = method, p = 3.5, diag = TRUE, upper = TRUE))
    dissimMat <- dissimMat/max(dissimMat) * 100
  }
  # if (!fossil::tri.ineq(dissimMat)) {
  #   logger::log_error("TALE domains dissimilarity (distance) matrix does not respect the triangle inequality",
  #                     "Arlem will fail. Aborting...")
  #   stop()
  # }
  
  #### Prepare Arlem cfile with systematic pairwise distances between 'repeat' units. ####
  # Get parameters for arlem
  TypeNo <- glue::glue("# Type no. ", nrow(dissimMat))
  Types <- glue::glue("# Types ", paste(1:nrow(dissimMat), collapse = " "))
  # Convert mat to character, beware of the ceiling in conversion...
  dissimMat <- matrix(ceiling(dissimMat) %>% format(),
                      ncol = ncol(dissimMat),
                      dimnames = list(rownames(dissimMat), colnames(dissimMat))
                      )
  # 'erase' lower triangle and diag
  dissimMat[lower.tri(dissimMat)] <- ""
  diag(dissimMat) <- ""
  
  # Convert mat rows to strings of space separated values
  dissimMatLines <- apply(dissimMat, 1, function(row) {
    string <- paste(row, collapse = " ")
    gsub("^[ ]+", "", string)
  })
  # remove empty last line
  dissimMatLines <- dissimMatLines[1:(length(dissimMatLines) - 1)]
  # Add the arlem stuff
  dissimMatLines <- c(TypeNo, Types,
                      "# Indel align 10", "# Indel hist 10", "# Dup 10",
                      "# matrix",
                      dissimMatLines)
  
  # write lines in a temp cfile
  cfile <- tempfile()
  writeLines(cfile, text = dissimMatLines)
  
  #### run arlem with a system call and parse std output ####
  arlemPath <- system.file("tools", "arlem", "arlem", package = "tantale", mustWork = T)
  arlemCmd <- glue::glue("{arlemPath} -f {codesSeqsfile} -cfile {cfile} -align -insert -showalign")
  logger::log_info("Running ARLEM version 1.0 : ")
  logger::log_info("Copyright by Mohamed I. Abouelhoda")
  logger::log_info("Plz. cite Abouelhoda, Giegerich, Behzadi, and Steyaert")
  arlemRawRes <- system(arlemCmd, intern = TRUE)
  logger::log_debug(logger::skip_formatter(arlemRawRes))
  arlemSelfRes <- grep("Processed Seq[.]:", arlemRawRes, value = TRUE) 
  arlemSelfScores <- gsub("Processed Seq[.]: ([0-9]{1,}) Score: ([0-9]{1,}),.*", "\\1|\\2",
                          substring(arlemSelfRes, 1, 35)) %>%
    strsplit(split = "\\|") %>%
    lapply(function(s) t(as.matrix(as.numeric(s)))) %>%
    do.call(rbind, .) %>% tibble::as_tibble(.name_repair = "minimal")
  arlemRes <- grep("Score of aligning Seq:",
                   arlemRawRes, value = TRUE)
  arlemScores <- gsub("Score of aligning Seq:([0-9]+), Seq:([0-9]+) =([0-9]+)", "\\1|\\2|\\3",
                      arlemRes)
  arlemScores <- strsplit(arlemScores, split = "\\|")
  arlemScores <- lapply(arlemScores, function(s) {t(as.matrix(as.numeric(s)))}) %>%
    do.call(rbind, .) %>%
    tibble::as_tibble(.name_repair = "minimal")
  colnames(arlemScores) <- c("TAL1", "TAL2", "arlemScore")
  # Shaping into matrix to have scores in both directions (fill diag and triangle)
  arlemScoresMat <- reshape2::acast(arlemScores, formula = TAL1 ~ TAL2, value.var = "arlemScore")
  arlemScoresMat <- cbind("0" = NA, arlemScoresMat)
  arlemScoresMat <- rbind(arlemScoresMat, NA)
  rownames(arlemScoresMat)[length(codesSeqSet)] <- length(codesSeqSet) - 1
  arlemScores <- stats::as.dist(t(arlemScoresMat), diag = TRUE, upper = TRUE) %>% as.matrix() %>%
    reshape2::melt(value.name = "arlemScore") %>%
    tibble::as_tibble()
  colnames(arlemScores) <- c("TAL1", "TAL2", "arlemScore")
  
  #### Compute normalized arlem scores and include arrayIDs rather than arlem index ####
  arrayLengths <- tale_parts %>% dplyr::group_by(arrayID) %>% dplyr::count()
  
  normArlemScoresTble <- arlemScores %>%
    dplyr::mutate(
      TAL1 = names(codesSeqSet)[TAL1 + 1],
      TAL2 = names(codesSeqSet)[TAL2 + 1]
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      maxLength = max(arrayLengths$n[arrayLengths$arrayID == TAL1],
                      arrayLengths$n[arrayLengths$arrayID == TAL2]),
      normArlemScore = arlemScore/maxLength,
      Sim = 100 - normArlemScore
    ) %>%
    dplyr::ungroup()
  
  #### Check features of the Arlem results table
  arraysCount <- codesSeqSet %>% length()
  if (nrow(normArlemScoresTble) != arraysCount^2) {
    allCombs <- expand.grid(names(codesSeqSet), names(codesSeqSet), stringsAsFactors = FALSE) %>% tibble::as_tibble()
    colnames(allCombs) <- c("TAL1", "TAL2")
    absentCombs <- dplyr::left_join(allCombs, normArlemScoresTble) %>%
      dplyr::filter(is.na(Sim))
    cat(knitr::kable(absentCombs), sep = "\n")
    logger::log_error("The TALE similarity table does not have the expected number",
    "of comparisons (number of missing comps: {nrow(absentCombs)})...")
    stop()
  }
  
  
  #### return the raw pieces; callers class and assemble them ####
  logger::log_info("Finished computing TALE and repeat relatedness.")
  list(
    tale_parts = tale_parts,
    dissim_long = dissimLong,
    tal_sim = normArlemScoresTble,
    coded_seq_set = codesSeqSet
  )
}




