

.tale_parts_from_file <- function(fasta) {
  if (grepl("TALE_Protein_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readAAStringSet(fasta)
  if (grepl("TALE_DNA_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readDNAStringSet(fasta)
  if (length(taleStrings) == 0L) {
    cli::cli_warn("No part sequence found in: {fasta}. Returning an empty tibble.")
    tbl <- tibble::tibble(array_id = character(),
                   domain_type = character(),
                   position_in_array = character(),
                   position_in_crd = character(),
                   string = character(),
                   source_directory = character())
    return(tbl)
  } 
  tbl <- tibble::tibble(array_id = gsub("(.*): .*", "\\1", names(taleStrings)),
                        domain_type = gsub(".*: (.*?)[ ]?[0-9]{0,}$", "\\1", names(taleStrings)),
                        position_in_crd = gsub(".*: repeat[ ]([0-9]{0,})$", "\\1", names(taleStrings)) %>%
                          as.integer() %>%
                          suppressWarnings(),
                        string = as.character(taleStrings) %>% as.vector(),
                        source_directory = dirname(fasta)
  )
  missingTerm <- setdiff(c("N-terminus", "C-terminus"), unique(tbl$domain_type))
  if (length(missingTerm) != 0L) {
    cli::cli_warn("Array {unique(tbl$array_id)} is missing a {missingTerm} domain in {fasta}")
    missingTerm <- tibble::tibble(array_id = unique(tbl$array_id),
                                  domain_type = missingTerm,
                                  position_in_crd = NA,
                                  string = NA,
                                  source_directory = dirname(fasta))
    tbl <- dplyr::bind_rows(tbl, missingTerm)
  }
  tbl %<>% 
    dplyr::rowwise() %>%
    dplyr::mutate(position_in_array = switch(domain_type,
                                                 `N-terminus` = 1,
                                                 `repeat` = position_in_crd + 1,
                                                 `C-terminus` = nrow(tbl)
                                                 ))
  return(tbl)
}




.rvds_from_annotale_file <- function(fasta) {
  if (!grepl("TALE_RVDs.fasta", basename(fasta))) {
    cli::cli_abort("The provided file does not seem to be an AnnoTALE RVDs file: {fasta}", class = c("tantale_error"))
  } else {
    rvdTble <- .split_list(fasta) %>%
      lapply(function(x) tibble::tibble(string = x,
                                        position_in_crd = 1:length(x))
             ) %>%
      dplyr::bind_rows(.id = "array_id")
    rvdTble <- rvdTble %>% dplyr::mutate(source_directory = dirname(fasta),
                                         domain_type = "repeat")
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
  # !!!! array_id are assumed to be unique !!!!
  protPartsFiles <- list.files(telltale_dir, "TALE_Protein_parts.fasta", recursive = T, full.names = T)
  dnaPartsFiles <- list.files(telltale_dir, "TALE_DNA_parts.fasta", recursive = T, full.names = T)
  if (telltale_dir %>% dirname() %>% unique() %>% length() != 1L) {
    cli::cli_warn("The provided path most likely does not correspond to a SINGLE tell_tales output directory.")
  }
  # Fetch info from annotale/telltale files with .tale_parts_from_file
  taleProtString <- lapply(protPartsFiles, .tale_parts_from_file) %>% dplyr::bind_rows()
  taleDnaString <- lapply(dnaPartsFiles, .tale_parts_from_file) %>% dplyr::bind_rows()
  #stopifnot(nrow(taleProtString) == nrow(taleDnaString))
  # Join info in a table with one domain per row
  tale_parts <- dplyr::full_join(taleDnaString %>% dplyr::rename(dna_seq = string),
                                taleProtString %>% dplyr::rename(aa_seq = string),
                                by = c("array_id", "domain_type", "position_in_array", "position_in_crd", "source_directory"),
                                relationship = "one-to-one") %>%
    dplyr::mutate(aa_seq = gsub("[*]", "", aa_seq))

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
    lapply(function(x) tibble::tibble(rvd = x, position_in_array = 1:length(x) )) %>%
    dplyr::bind_rows(.id = "array_id")  
  anchorCodes <- tales_anchor_codes()
  
  
  # Some checks on the consistency between parts and rvd sequences
  # if nhmmer did not report on a C-Term CDS, the corresponding domain
  # "CTERM" tag will not be written in the rvd slot of the table.
  arraysConsistency <- dplyr::full_join(tale_parts %>% dplyr::count(array_id, name = "AnnoTALELength"),
                                        rvds %>% dplyr::count(array_id, name = "rvdFileLength"),
                                        by = dplyr::join_by(array_id)) %>%
    dplyr::mutate(sameLength = AnnoTALELength == rvdFileLength)
  
  if (any(is.na(arraysConsistency$sameLength))) {
    cli::cli_warn("There are mismatches in array IDs between rvd seq file and AnnoTALE parts files:")
    cli::cli_warn("There are mismatches in array IDs between rvd seq file and AnnoTALE parts files.")
  } else if (!all(arraysConsistency$sameLength, na.rm = TRUE)) {
    cli::cli_abort(
      c("Array lengths are inconsistent between the rvd seq file and the AnnoTALE parts files.",
        "i" = "Affected arrays: {.val {arraysConsistency$array_id[!arraysConsistency$sameLength]}}"),
      class = c("tantale_error_parts_inconsistent", "tantale_error"))
  }
  
  # Include RVDs in the talParts tibble
  tale_parts <- dplyr::left_join(tale_parts, 
                                rvds,
                                by = c("array_id", "position_in_array"),
                                unmatched = "drop", relationship = "one-to-one")
  # Include seqnames in the talParts tibble
  tale_parts %<>% dplyr::left_join(
    readr::read_tsv(list.files(telltale_dir, "hitsReport.tsv", recursive = T, full.names = T),
                    show_col_types = FALSE) %>%
      dplyr::select(array_id, seqnames) %>%
      dplyr::distinct(),
    by = "array_id", relationship = "many-to-one"
  )
  # Check talparts
  partsWithMissingAaSeq <- tale_parts %>% dplyr::filter(is.na(aa_seq)) %>% dplyr::pull(array_id) %>% unique()
  partsWithMissingDnaSeq <- tale_parts %>% dplyr::filter(is.na(dna_seq)) %>% dplyr::pull(array_id) %>% unique()
  partsWithMissingRvdSeq <- tale_parts %>% dplyr::filter(is.na(rvd)) %>% dplyr::pull(array_id) %>% unique()
  if (any(sapply(list(partsWithMissingAaSeq, partsWithMissingDnaSeq, partsWithMissingRvdSeq), length) != 0L)) {
    cli::cli_warn(c("The returned tale_parts has records with missing sequences.",
                    "i" = "Affected array{?s}: {.val {unique(c(partsWithMissingAaSeq, partsWithMissingDnaSeq, partsWithMissingRvdSeq))}}"))
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





#' Visualise the domain composition of a set of TALE arrays
#'
#' @description
#' A compact, information-rich view of the arrays in a \code{tales} object: one
#' point per part, positioned by its place in the array, coloured by domain type
#' and filled by amino-acid length, with the RVD printed on each repeat.
#'
#' @param position Which coordinate to lay the parts out on. \code{"array"}
#'   (default) uses \code{position_in_array}, so each array starts at 1 and runs
#'   contiguously. \code{"alignment"} uses \code{alignment_position}, which
#'   requires an aligned object (or one demoted from a \code{\link{tales_msa}},
#'   which keeps the column): gaps then appear as empty columns and shared
#'   features line up. Aberrant repeats, for instance, are visible as a column
#'   in the aligned layout and scattered in the unaligned one.
#' @param x A \code{\link{tales}} object, as returned by
#'   \code{\link{tales_from_telltale}} or in the \code{tales} element of
#'   \code{\link{tales_compare}}'s output. A legacy \code{tale_parts} data
#'   frame is accepted and converted.
#' @return The ggplot object, invisibly printed as a side effect.
#' @export
#' @family TALE plots
plot_tales_composition <- function(x, position = c("array", "alignment")) {
  position <- match.arg(position)
  if (!is_tales(x)) x <- tales(x)
  .tales_require(x, "plot_tales_composition")
  if (identical(position, "alignment") && !"alignment_position" %in% names(x)) {
    cli::cli_abort(
      c("{.code position = \"alignment\"} needs the {.field alignment_position} column.",
        "i" = "Align first with {.fn tales_align}, or use {.code position = \"array\"}."),
      class = c("tantale_error_projection_column", "tantale_error")
    )
  }
  partsForPlots <- x %>%
    dplyr::mutate(label = dplyr::if_else(domain_type == "repeat", rvd, ""),
                  aa_length = factor(nchar(aa_seq)),
                  .x = if (identical(position, "alignment")) .data$alignment_position
                       else .data$position_in_array)

  p <- partsForPlots %>%
    ggplot2::ggplot(mapping = ggplot2::aes(fill = aa_length,
                                           color = domain_type,
                                           label = label,
                                           y = array_id,
                                           x = .x)) +
    ggplot2::scale_color_viridis_d(option = "rocket") +
    ggplot2::scale_fill_discrete() +
    ggplot2::scale_x_continuous(
      name = if (identical(position, "alignment")) "Position in alignment" else "Position in array",
      breaks = 1:100, minor_breaks = NULL) +
    ggplot2::geom_point(shape = 21, size = 5, stroke = 0.9) +
    ggnewscale::new_scale_color() +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_text(size = 2.1, color = "white") +
    ggplot2::labs(title = "Overview of TALE composition by genome") +
    ggplot2::theme_light()

  # seqnames groups arrays by source contig. It is optional in a tales, so the
  # facet is added only when it is there -- a fasta-derived object has none.
  if ("seqnames" %in% names(x)) {
    p <- p + ggplot2::facet_grid(seqnames ~ ., scales = "free_y", space = "free")
  }

  print(p)
  invisible(p)
}


#' Plot the domain composition of a tales object
#'
#' @description
#' The \code{\link[=plot]{plot}} method for \code{\link{tales}}, delegating to
#' \code{\link{plot_tales_composition}}. A \code{\link{tales_msa}} dispatches
#' to \code{\link{plot.tales_msa}} instead, being the more specific class.
#'
#' @param x A \code{\link{tales}} object.
#' @param ... Passed to \code{\link{plot_tales_composition}}.
#' @return The ggplot object, invisibly.
#' @method plot tales
#' @export
#' @family TALE plots
plot.tales <- function(x, ...) {
  plot_tales_composition(x, ...)
}


# tale_parts <- readRDS("/home/cunnac/TEMP/talePartsForDistalr.rds")
# h_cut = 10
# ncores = 1
# aln_method = "DECIPHER"
# conda_bin = "/home/cunnac/bin/miniconda3/condabin/conda"



#' Derive aa_seq from dna_seq by translation
#'
#' \code{tales_compare()} needs protein sequences, but an object may carry only
#' the DNA. TALE part coding sequences are in frame, so translating them
#' recovers \code{aa_seq} exactly.
#'
#' Two details are easy to get silently wrong. \code{no.init.codon = TRUE} is
#' required: TALE repeats begin on \code{CTG}/\code{TTG}, which are alternative
#' start codons, so the default forces the first residue to \code{M} -- that
#' alone accounted for 865 of 955 mismatches on the reference fixture. And the
#' C-terminal parts carry a trailing stop codon, which is stripped.
#'
#' With both handled, translation reproduces the stored \code{aa_seq} for all
#' 955 parts of the reference fixture.
#'
#' @param dna A character vector of in-frame coding sequences.
#' @return A character vector of amino-acid sequences.
#' @keywords internal
.translate_parts <- function(dna) {
  bad <- is.na(dna) | !nzchar(dna)
  out <- rep(NA_character_, length(dna))
  if (all(bad)) return(out)
  ok <- !bad
  if (any(nchar(dna[ok]) %% 3 != 0)) {
    cli::cli_abort(
      c("Cannot translate {.field dna_seq}: some sequences are not a whole number of codons.",
        "i" = "TALE part coding sequences are expected to be in frame."),
      class = c("tantale_error_translate_frame", "tantale_error")
    )
  }
  aa <- suppressWarnings(as.character(Biostrings::translate(
    Biostrings::DNAStringSet(dna[ok]),
    no.init.codon = TRUE, if.fuzzy.codon = "solve"
  )))
  out[ok] <- sub("[*]$", "", aa)
  out
}


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
#' @family pairwise distances
tales_compare <- function(x, ncores = 1, aln_method = "DECIPHER",
                              conda_bin = "auto") {
  if (!is_tales(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  # aa_seq is what the comparison needs, but dna_seq satisfies it by
  # translation -- the same "one of these two" shape the column contract
  # already uses for rvd / dom_code.
  if (!"aa_seq" %in% names(x) && "dna_seq" %in% names(x)) {
    x$aa_seq <- .translate_parts(x$dna_seq)
    cli::cli_warn(
      c("Derived {.field aa_seq} by translating {.field dna_seq}.",
        "i" = "The translated column is kept in the returned {.cls tales}."),
      class = c("tantale_warning_translated_aa", "tantale_warning")
    )
  }
  if (!"aa_seq" %in% names(x)) {
    cli::cli_abort(
      c("{.arg x} must carry an {.field aa_seq} or {.field dna_seq} column.",
        "i" = "Repeat similarity is computed from part amino acid sequences,",
        "i" = "which can be translated from {.field dna_seq} if needed."),
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

  core <- .tales_compare_core(tale_parts = tibble::as_tibble(x), ncores = ncores,
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


#' The expensive part of the relatedness computation
#'
#' Called by \code{\link{tales_compare}}. Takes a plain tibble in the canonical
#' column vocabulary and returns raw pieces; classing, stamping and assembly
#' happen in the caller. Deliberately
#' does no clustering: that was a stored field with no consumers, recomputed by
#' its only would-be user at a different cut height (restructuring-notes.md §1).
#' @noRd
.tales_compare_core <- function(tale_parts, ncores = 1,
                                aln_method = "DECIPHER", conda_bin = "auto") {
  
  #### Reality checks ####
  
  ## Make sure we are dealing only with parts that have defined protein sequences.
  if (any(is.na(tale_parts$aa_seq) | tale_parts$aa_seq == "")) {
    # must match the guard above, or an empty-string part lists nothing
    badArrays <- unique(tale_parts$array_id[is.na(tale_parts$aa_seq) | tale_parts$aa_seq == ""])
    cli::cli_abort(
      c("Some of the provided TALE parts have no amino acid sequence.",
        "i" = "Affected array{?s}: {.val {badArrays}}"),
      class = c("tantale_error_parts_no_aa", "tantale_error"))
  }
  if (any(is.na(tale_parts$dna_seq) | tale_parts$dna_seq == "")) {
    cli::cli_warn("It seems that some of the provided TALE parts miss the DNA sequence!")
  } 
  
  ## Make sure that array_id - position combinations are unique
  # in case someone would not have made array_ids unique before
  # mixing tale predictions from several genomes...
  arayPosCombinCounts <- tale_parts %>%
    dplyr::group_by(array_id, position_in_array) %>%
    dplyr::count() %>%
    dplyr::pull(n)
  if (!all(arayPosCombinCounts == 1L)) {
    cli::cli_abort(paste0("Your tale arrays identifers are probably not unique.",
         "\n",
         "Make sure that there is only one part per position per array_id."), class = c("tantale_error"))
  }
  
  # Assign domain codes
  tale_parts %<>% dplyr::group_by(aa_seq) %>%
    dplyr::mutate(dom_code = dplyr::cur_group_id() %>% unlist() %>% as.character()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dom_code = dplyr::if_else(is.na(aa_seq), as.character(NA), dom_code))
  
  
  #### Assemble repeat code strings and write in a file for arlem ####
  cli::cli_inform("Assemble repeat code TALE strings and write in a file for ARLEM")
  
  codesSeqsfile <- tempfile(fileext = ".fasta")
  repeatStrings <- tale_parts %>%
    dplyr::group_by(array_id) %>%
    dplyr::arrange(position_in_array) %>%
    dplyr::summarise(repeatString = paste(dom_code, collapse = " "),
                     posString = paste(position_in_array, collapse = " "))
  codesSeqSet <- Biostrings::BStringSet(repeatStrings$repeatString)
  names(codesSeqSet) <- repeatStrings$array_id
  
  # # Must use seqinr because Biostrings wraps sequences in fasta file which messes up Arlem...
  # codesSeqLst <- as.list(repeatStrings$repeatString)
  # names(codesSeqLst) <- repeatStrings$array_id
  # seqinr::write.fasta(codesSeqLst, names = names(codesSeqLst),
  #                     file.out = codesSeqsfile, as.string = TRUE, nbchar = 10000)
  
  Biostrings::writeXStringSet(x = codesSeqSet,
                              filepath = codesSeqsfile,
                              format = "fasta", width = 20000L)
  
  
  
  #### Compute systematic pairwise dissimilarities (distances) between 'repeat' units. ####
  # Get unique domains sequences
  taleAaParts <- Biostrings::AAStringSet(tale_parts$aa_seq)
  names(taleAaParts) <- tale_parts$dom_code
  uniqueTaleAaParts <- unique(taleAaParts)
  stopifnot(!anyDuplicated(names(uniqueTaleAaParts)))
  stopifnot(!anyDuplicated(names(unique(taleAaParts))))
  
  # Get pairwise repeat aa sequence dissimilarity scores in a long tibble
  cli::cli_inform(paste0("Computing a distance matrix between TALE parts amino acid sequences ",
                   "using: {aln_method}"))
  if (aln_method == "mmseq2") {
    dissimLong <- .pairwise_align_mmseq2(part_aa_set = uniqueTaleAaParts, ncores = ncores,
                                        conda_bin = conda_bin)
    #saveRDS(dissimLong, file = "/home/cunnac/TEMP/dissimLong")
  } else if (aln_method == "Biostrings") {
    dissimLong <- .pairwise_align_biostrings(part_aa_set = uniqueTaleAaParts, ncores = ncores)
  } else if (aln_method == "DECIPHER") {
    dissimLong <- .pairwise_align_decipher(part_aa_set = uniqueTaleAaParts, ncores = ncores)
  } else {
    cli::cli_abort("{.arg aln_method} must be one of {.val Biostrings}, {.val mmseq2} or {.val DECIPHER}, not {.val {aln_method}}.",
                   class = c("tantale_error_aln_method", "tantale_error"))
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
    cli::cli_inform(paste0("Generate an ARLEM cost matrix which meets triangle inequality criteria by computing ",
                     "the {method} distance between pairwise distance vectors."))
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
  cli::cli_inform("Running ARLEM version 1.0 : ")
  cli::cli_inform("Copyright by Mohamed I. Abouelhoda")
  cli::cli_inform(paste0("Plz. cite Abouelhoda, Giegerich, Behzadi, and Steyaert"))
  arlemRawRes <- system(arlemCmd, intern = TRUE)
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
  
  #### Compute normalized arlem scores and include array_ids rather than arlem index ####
  arrayLengths <- tale_parts %>% dplyr::group_by(array_id) %>% dplyr::count()
  
  normArlemScoresTble <- arlemScores %>%
    dplyr::mutate(
      TAL1 = names(codesSeqSet)[TAL1 + 1],
      TAL2 = names(codesSeqSet)[TAL2 + 1]
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      maxLength = max(arrayLengths$n[arrayLengths$array_id == TAL1],
                      arrayLengths$n[arrayLengths$array_id == TAL2]),
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
      cli::cli_abort(
        c("The TALE similarity table does not have the expected number of comparisons.",
          "x" = "Expected {arraysCount^2}, got {nrow(normArlemScoresTble)}; {nrow(absentCombs)} missing.",
          "i" = "First missing pair{?s}: {.val {paste(utils::head(absentCombs$TAL1, 3), utils::head(absentCombs$TAL2, 3), sep = \"/\")}}"),
        class = c("tantale_error_arlem_incomplete", "tantale_error"))
  }
  
  
  #### return the raw pieces; callers class and assemble them ####
  cli::cli_inform("Finished computing TALE and repeat relatedness.")
  list(
    tale_parts = tale_parts,
    dissim_long = dissimLong,
    tal_sim = normArlemScoresTble,
    coded_seq_set = codesSeqSet
  )
}




