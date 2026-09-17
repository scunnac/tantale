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
