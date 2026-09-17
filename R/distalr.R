

#### The three steps of tales_compare() ####
#
# tales_compare() used to be one call around a ~200 line internal. The three
# things it does are separable and each independently useful, so they are
# exported separately and the wrapper composes them (restructuring-notes.md
# 8.5b).
#
# The steps are a chain, not a menu:
#
#   x  <- tales_assign_domain_codes(x)            # 1
#   dd <- tales_domain_distances(x)               # 2
#   td <- tales_tale_distances(x, dd)             # 3, consumes dd
#
# Step 3 consuming step 2 is the part worth exposing, because it states
# something about the biology that the monolith concealed: two TALEs are
# compared by aligning their repeat arrays, where the cost of substituting
# one repeat for another is how different those repeats are as proteins. The
# repeat-level comparison is not a by-product of the TALE-level one, it is
# its input.


#' Assign a domain code to every distinct part sequence
#'
#' @description
#' Step 1 of [tales_compare()]. Gives each distinct `aa_seq` in `x` an
#' integer code, recorded in a `dom_code` column, so that the parts can be
#' compared once each rather than once per occurrence.
#'
#' A code names a distinct **domain** sequence, not a repeat: the N- and
#' C-termini are parts like the repeats are, and they get codes too.
#'
#' @details
#' TALEs reuse domains heavily, within an array and between arrays, so the
#' number of distinct sequences is far smaller than the number of parts. That
#' ratio is what makes [tales_domain_distances()] affordable -- it compares
#' distinct domains, not parts.
#'
#' @section The codes are only meaningful within one call:
#' Codes are assigned with [dplyr::cur_group_id()] over the distinct `aa_seq`
#' values **present in `x`**. Add an array, remove one, or reorder the
#' sequences, and the same protein can get a different number. They are
#' positions in this table's own vocabulary, not identifiers of anything.
#'
#' So **comparing codes between two calls is an error**, and a similarity
#' table keyed by one call's codes must never be used with another call's.
#' This is not a caution to remember: it is enforced. Every object minted
#' here is stamped with a `dom_code_namespace` derived from the sequences
#' that produced it, and the classes refuse to join across namespaces. Read
#' the stamp with [tales_namespace()].
#'
#' @param x A [tales] object carrying an `aa_seq` column.
#' @return `x` with a `dom_code` column and a `dom_code_namespace` stamp.
#' @seealso [tales_compare()], which runs all three steps;
#'   [tales_domain_codes()], which reads the code-to-sequence table back out
#'   of an object that already has codes.
#' @export
#' @family pairwise distances
tales_assign_domain_codes <- function(x) {
  if (!is_tales(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  if (!"aa_seq" %in% names(x)) {
    cli::cli_abort(
      c("{.arg x} must carry an {.field aa_seq} column.",
        "i" = "Domain codes name distinct part amino acid sequences."),
      class = c("tantale_error_compare_no_aa", "tantale_error"))
  }
  .assert_parts_have_aa(x)
  
  namespace <- .tales_dom_code_namespace(x$aa_seq)
  tales(.assign_dom_codes(tibble::as_tibble(x)),
        dom_code_namespace = namespace)
}


#' Pairwise distances between distinct TALE domains
#'
#' @description
#' Step 2 of [tales_compare()]. Aligns every distinct domain sequence against
#' every other and returns their pairwise dissimilarity.
#'
#' @details
#' This is the expensive step, and the one worth having on its own: the
#' repeat-level distances answer questions about repeat diversity that need
#' no TALE-level alignment at all, and computing them does not require
#' running ARLEM.
#'
#' Distances are between **distinct domains**, keyed by `dom_code`, so the
#' cost goes with the number of distinct sequences rather than the number of
#' parts. See [summary()][summary.tales] for that ratio on a given object.
#'
#' @param x A [tales] object carrying `dom_code`, as returned by
#'   [tales_assign_domain_codes()].
#' @param aln_method One of `"DECIPHER"` (the default), `"Biostrings"` or
#'   `"mmseq2"`. `"mmseq2"` runs in the conda environment.
#' @param ncores Number of cores for the pairwise alignment.
#' @param conda_bin Passed to `reticulate`, for `aln_method = "mmseq2"`.
#' @return A [domain_distances] object, stamped with `x`'s namespace.
#' @seealso [tales_tale_distances()], which consumes this.
#' @export
#' @family pairwise distances
tales_domain_distances <- function(x, aln_method = "DECIPHER", ncores = 1,
                                   conda_bin = "auto") {
  .assert_coded_tales(x, "tales_domain_distances")
  parts <- tibble::as_tibble(x)
  
  aa <- Biostrings::AAStringSet(parts$aa_seq)
  names(aa) <- parts$dom_code
  unique_aa <- unique(aa)
  stopifnot(!anyDuplicated(names(unique_aa)))
  
  cli::cli_inform("Computing a distance matrix between TALE parts amino acid sequences using: {aln_method}")
  dissim <- switch(
    aln_method,
    mmseq2 = .pairwise_align_mmseq2(part_aa_set = unique_aa, ncores = ncores,
                                    conda_bin = conda_bin),
    Biostrings = .pairwise_align_biostrings(part_aa_set = unique_aa, ncores = ncores),
    DECIPHER = .pairwise_align_decipher(part_aa_set = unique_aa, ncores = ncores),
    cli::cli_abort(
      "{.arg aln_method} must be one of {.val Biostrings}, {.val mmseq2} or {.val DECIPHER}, not {.val {aln_method}}.",
      class = c("tantale_error_aln_method", "tantale_error"))
  )
  dissim <- dplyr::mutate(dissim, sim = 100 - dissim)
  
  domain_distances(dissim, dom_code_namespace = tales_namespace(x))
}


#' Pairwise distances between TALE arrays
#'
#' @description
#' Step 3 of [tales_compare()]. Aligns the arrays against each other as
#' strings of domain codes, using the domain distances from step 2 as the
#' cost of substituting one domain for another.
#'
#' @details
#' The two arguments are not independent, and that is the point. ARLEM aligns
#' each array's sequence of `dom_code`s; what it costs to align one domain
#' against a different one is taken from `domain_distances`, so the TALE-level
#' comparison is built on the domain-level one rather than computed beside it.
#'
#' The domain distances are first passed through a Minkowski distance
#' (`p = 3.5`) between their rows and rescaled to 0-100. That step is not
#' cosmetic: ARLEM needs a cost matrix satisfying the triangle inequality,
#' and raw pairwise alignment dissimilarities do not.
#'
#' @section Both arguments must come from the same call:
#' `domain_distances` is keyed by `dom_code`, and those codes mean what they
#' mean only within the call that minted them
#' ([tales_assign_domain_codes()]). Passing distances from one run with codes
#' from another silently compares the wrong domains, which is why both
#' objects carry a namespace stamp and this function refuses when they
#' disagree.
#'
#' @param x A [tales] object carrying `dom_code`.
#' @param domain_distances A [domain_distances] object for the same `x`, as
#'   returned by [tales_domain_distances()].
#' @return A [tale_distances] object, keyed by `array_id`.
#' @seealso [tales_compare()], which runs all three steps.
#' @export
#' @family pairwise distances
tales_tale_distances <- function(x, domain_distances) {
  .assert_coded_tales(x, "tales_tale_distances")
  .assert_same_namespace(x, domain_distances)
  
  parts <- tibble::as_tibble(x)
  coded <- .coded_seq_set(parts)
  cost <- .arlem_cost_file(domain_distances)
  scores <- .run_arlem(coded, cost)
  tale_distances(.normalise_arlem_scores(scores, parts, coded))
}


#### shared internals ####

#' Every part must have an amino acid sequence to be compared
#' @noRd
.assert_parts_have_aa <- function(x) {
  missing <- is.na(x$aa_seq) | x$aa_seq == ""
  if (any(missing)) {
    cli::cli_abort(
      c("Some of the provided TALE parts have no amino acid sequence.",
        "i" = "Affected array{?s}: {.val {unique(x$array_id[missing])}}"),
      class = c("tantale_error_parts_no_aa", "tantale_error"))
  }
  if ("dna_seq" %in% names(x) && any(is.na(x$dna_seq) | x$dna_seq == "")) {
    cli::cli_warn("It seems that some of the provided TALE parts miss the DNA sequence!")
  }
  invisible(NULL)
}

#' The dom_code assignment itself, on a plain tibble
#' @noRd
.assign_dom_codes <- function(parts) {
  parts %>%
    dplyr::group_by(aa_seq) %>%
    dplyr::mutate(dom_code = as.character(unlist(dplyr::cur_group_id()))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dom_code = dplyr::if_else(is.na(aa_seq), NA_character_, dom_code))
}

#' Steps 2 and 3 both require codes to already exist
#' @noRd
.assert_coded_tales <- function(x, fn) {
  if (!is_tales(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  if (!"dom_code" %in% names(x)) {
    cli::cli_abort(
      c("{.fn {fn}} needs a {.field dom_code} column.",
        "i" = "Assign codes first with {.fn tales_assign_domain_codes}."),
      class = c("tantale_error_projection_column", "tantale_error"))
  }
  invisible(NULL)
}

#' Codes from one run must not be used with distances from another
#' @noRd
.assert_same_namespace <- function(x, dd) {
  a <- tales_namespace(x)
  b <- tales_namespace(dd)
  if (!is.null(a) && !is.null(b) && nzchar(a) && nzchar(b) && !identical(a, b)) {
    cli::cli_abort(
      c("{.arg x} and {.arg domain_distances} come from different runs.",
        "x" = "Namespaces {.val {a}} and {.val {b}}.",
        "i" = "Domain codes are only meaningful within the call that minted them.",
        "i" = "Recompute the distances from this {.arg x} with {.fn tales_domain_distances}."),
      class = c("tantale_error_namespace_mismatch", "tantale_error"))
  }
  invisible(NULL)
}

#' Each array as a space-separated string of its domain codes
#'
#' Not tales_coded_strings(): ARLEM is given these through a file whose
#' record order defines the integer index it reports results by, so the
#' construction and that ordering have to stay together.
#' @noRd
.coded_seq_set <- function(parts) {
  strings <- parts %>%
    dplyr::group_by(array_id) %>%
    dplyr::arrange(position_in_array) %>%
    dplyr::summarise(repeatString = paste(dom_code, collapse = " "),
                     posString = paste(position_in_array, collapse = " "))
  out <- Biostrings::BStringSet(strings$repeatString)
  names(out) <- strings$array_id
  out
}

#' Write ARLEM's substitution cost matrix
#'
#' The Minkowski pass is what makes the matrix usable: ARLEM requires the
#' triangle inequality and raw pairwise dissimilarities do not satisfy it.
#' @return The path of the cfile.
#' @noRd
.arlem_cost_file <- function(dd) {
  cli::cli_inform("Generate an ARLEM cost matrix which meets triangle inequality criteria by computing the minkowski distance between pairwise distance vectors.")
  mat <- reshape2::acast(tibble::as_tibble(dd), formula = id1 ~ id2,
                         value.var = "dissim")
  stopifnot(nrow(mat) == ncol(mat))
  # ARLEM's types are 1..n in order, so the rows must be in numeric code
  # order for type i to mean domain code i.
  mat <- mat[order(as.numeric(rownames(mat))), order(as.numeric(colnames(mat)))]
  
  mat <- as.matrix(stats::dist(mat, method = "minkowski", p = 3.5,
                               diag = TRUE, upper = TRUE))
  mat <- mat / max(mat) * 100
  
  header <- c(glue::glue("# Type no. ", nrow(mat)),
              glue::glue("# Types ", paste(seq_len(nrow(mat)), collapse = " ")),
              "# Indel align 10", "# Indel hist 10", "# Dup 10", "# matrix")
  
  mat <- matrix(format(ceiling(mat)), ncol = ncol(mat),
                dimnames = list(rownames(mat), colnames(mat)))
  mat[lower.tri(mat)] <- ""
  diag(mat) <- ""
  lines <- apply(mat, 1, function(row) gsub("^[ ]+", "", paste(row, collapse = " ")))
  lines <- lines[seq_len(length(lines) - 1L)]   # drop the empty last row
  
  cfile <- tempfile()
  writeLines(cfile, text = c(header, lines))
  cfile
}

#' Run ARLEM over the coded strings and parse its stdout
#' @noRd
.run_arlem <- function(coded, cfile) {
  seqfile <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(coded, filepath = seqfile, format = "fasta",
                              width = 20000L)
  arlem <- system.file("tools", "arlem", "arlem", package = "tantale",
                       mustWork = TRUE)
  cmd <- glue::glue("{shQuote(arlem)} -f {shQuote(seqfile)} -cfile {shQuote(cfile)} -align -insert -showalign")
  cli::cli_inform("Running ARLEM version 1.0 : ")
  cli::cli_inform("Copyright by Mohamed I. Abouelhoda")
  cli::cli_inform("Plz. cite Abouelhoda, Giegerich, Behzadi, and Steyaert")
  raw <- system(cmd, intern = TRUE)
  
  scores <- grep("Score of aligning Seq:", raw, value = TRUE)
  scores <- gsub("Score of aligning Seq:([0-9]+), Seq:([0-9]+) =([0-9]+)",
                 "\\1|\\2|\\3", scores)
  scores <- strsplit(scores, split = "\\|")
  scores <- tibble::as_tibble(
    do.call(rbind, lapply(scores, function(s) t(as.matrix(as.numeric(s))))),
    .name_repair = "minimal")
  colnames(scores) <- c("id1", "id2", "arlem_score")
  
  # ARLEM reports one direction only; mirror it so every ordered pair is
  # present, which is what the completeness check downstream expects.
  mat <- reshape2::acast(scores, formula = id1 ~ id2, value.var = "arlem_score")
  mat <- cbind("0" = NA, mat)
  mat <- rbind(mat, NA)
  rownames(mat)[length(coded)] <- length(coded) - 1
  out <- tibble::as_tibble(reshape2::melt(
    as.matrix(stats::as.dist(t(mat), diag = TRUE, upper = TRUE)),
    value.name = "arlem_score"))
  colnames(out) <- c("id1", "id2", "arlem_score")
  out
}

#' Turn ARLEM's indices into array ids and normalise by array length
#' @noRd
.normalise_arlem_scores <- function(scores, parts, coded) {
  lengths <- dplyr::count(dplyr::group_by(parts, array_id))
  
  out <- scores %>%
    dplyr::mutate(id1 = names(coded)[id1 + 1],
                  id2 = names(coded)[id2 + 1]) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      max_length = max(lengths$n[lengths$array_id == id1],
                       lengths$n[lengths$array_id == id2]),
      norm_arlem_score = arlem_score / max_length) %>%
    dplyr::ungroup()
  
  n <- length(coded)
  if (nrow(out) != n^2) {
    all_combs <- tibble::as_tibble(expand.grid(names(coded), names(coded),
                                               stringsAsFactors = FALSE))
    colnames(all_combs) <- c("id1", "id2")
    absent <- dplyr::filter(dplyr::left_join(all_combs, out,
                                             by = c("id1", "id2")),
                            is.na(norm_arlem_score))
    cli::cli_abort(
      c("The TALE similarity table does not have the expected number of comparisons.",
        "x" = "Expected {n^2}, got {nrow(out)}; {nrow(absent)} missing.",
        "i" = "First missing pair{?s}: {.val {paste(utils::head(absent$id1, 3), utils::head(absent$id2, 3), sep = '/')}}"),
      class = c("tantale_error_arlem_incomplete", "tantale_error"))
  }
  out
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
      tibble::tibble(id1 = names(part_aa_set[i]),
                     id2 = names(pwalign::alignedPattern(singleSubAln)),
                     score = BiocGenerics::score(singleSubAln),
                     nedit = Biostrings::nmismatch(singleSubAln))
    },
    BPPARAM = bpparam) %>%
    dplyr::bind_rows()
  
  pair_align_scores %<>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      max_length = max(Biostrings::nchar(part_aa_set[id1]), Biostrings::nchar(part_aa_set[id2])),
      # This is an approximate equivalent of how Alvaro computed dissimilarity in distal
      dissim = 100 - 100 * (max_length - score) / max_length,
      dissim = ifelse(dissim < 0, 100, 100 - dissim),
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-max_length)
  
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
  
  mmseqs <- shQuote(.tantale_bin("mmseqs", conda_bin = conda_bin))

  mmseq2createdb <- glue::glue("{mmseqs} createdb {shQuote(partAaStringSetFile)} {shQuote(mmseq2DbPath)}")

  mmseq2prefilter <- glue::glue("{mmseqs} prefilter {shQuote(mmseq2DbPath)} {shQuote(mmseq2DbPath)} {shQuote(prefDbPath)}",
                               "-v 3 --threads {max(floor(ncores/2), 1)} --max-seqs 1000 -s 7.5 --add-self-matches 1",
                               "--cov-mode 0", .sep = " ")
  
  mmseq2align <- glue::glue("{mmseqs} align {shQuote(mmseq2DbPath)} {shQuote(mmseq2DbPath)} {shQuote(prefDbPath)} {shQuote(alnDbPath)}",
                               "-v 3 --threads {ncores} --add-self-matches 1 --min-seq-id 0",
                               "--cov-mode 0 --gap-open aa:11,nucl:5 --gap-extend aa:1,nucl:2",
                               "-a 1 --alignment-mode 3 --alignment-output-mode 0 --seq-id-mode 1",
                               .sep = " ")
  
  mmseq2convertalis <- glue::glue("{mmseqs} convertalis {shQuote(mmseq2DbPath)} {shQuote(mmseq2DbPath)} {shQuote(alnDbPath)} {shQuote(alnTabFile)}",
                               "--format-mode 4 -v 3",
                               "--format-output query,target,evalue,raw,pident,nident,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qcov,tcov",
                               .sep = " ")
  
  if (!as.logical(.create_tantale_env(conda_bin = conda_bin))) {
    # Each stage is checked. Previously all four statuses were assigned to
    # `res` and none was tested, so a failed prefilter surfaced only as a
    # confusing error from convertalis -- or not at all.
    .tantale_exec(mmseq2createdb,    what = "mmseqs createdb")
    .tantale_exec(mmseq2prefilter,   what = "mmseqs prefilter")
    .tantale_exec(mmseq2align,       what = "mmseqs align")
    .tantale_exec(mmseq2convertalis, what = "mmseqs convertalis")
  } else {
    stop("Could not create the tantale conda environment on your machine to run mmseq2...")
  }
  

  pair_align_scores <- readr::read_tsv(alnTabFile, show_col_types = FALSE) %>%
    dplyr::mutate(target = as.character(target), query = as.character(query)) %>%
    dplyr::group_by(target, query) %>%
    dplyr::slice_max(raw, n = 1, with_ties = FALSE)
  pair_align_scores <- dplyr::left_join(df, pair_align_scores) %>%
    dplyr::rename(id1 = target, id2 = query) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(dissim = ifelse(is.na(pident), 100, 100 - pident*min(qcov,tcov))) %>%
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
  colnames(pair_align_scores) <- c("id2", "id1", "dissim")
  pair_align_scores %<>% dplyr::mutate(id2 = as.character(id2), id1 = as.character(id1))
  pair_align_scores %<>% dplyr::mutate(dissim = dissim*100) %>%
    dplyr::ungroup()
  
  # Check the pair_align_scores tibble
  .check_pair_align_tbl(pair_align_scores = pair_align_scores, part_aa_set = part_aa_set)
  
  return(pair_align_scores)
}





.check_pair_align_tbl <- function(pair_align_scores, part_aa_set) {
  partCombinCounts <- pair_align_scores %>% dplyr::select(id1, id2) %>%
    dplyr::count(id1, id2) %>%
    dplyr::pull(n)
  if (!all(partCombinCounts == 1L)) {
    stop("Some alignment pairs have more than one record...",)
  }
  if (length(names(part_aa_set))^2 != nrow(pair_align_scores)) {
    stop("Some parts pairs are absent from the pairwise parts distance table")
  }
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

  # The three steps, composed. Each is exported and separately useful; see
  # R/tales_compare_steps.R and restructuring-notes.md 8.5b. Step 3 takes
  # step 2's output because the TALE alignment is scored *using* the domain
  # distances, which is a fact about the method worth making visible.
  coded <- tales_assign_domain_codes(x)
  dd <- tales_domain_distances(coded, aln_method = aln_method,
                               ncores = ncores, conda_bin = conda_bin)
  td <- tales_tale_distances(coded, dd)

  cli::cli_inform("Finished computing TALE and repeat relatedness.")
  list(
    tales = coded,
    domain_distances = dd,
    # Keyed by array_id, not dom_code, so deliberately unstamped: array ids are
    # meaningful names that do not silently collide across runs the way
    # cur_group_id() codes do (class-design.md §3.5).
    tale_distances = td
  )
}






