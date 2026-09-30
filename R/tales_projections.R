##### Projections of a 'tales' object #####
# Views that were formerly stored as fields of distalr()'s returned list
# ('coded.repeats.str', 'repeats.code'). Both are pure functions of the parts
# table, so storing them invited the copy and the source drifting apart; see
# restructuring-notes.md §1.
#
# Both need `dom_code`, which is minted by tales_compare_distal().


#' Domain-coded strings, one per TALE array
#'
#' Renders each array as a separated string of its `dom_code`s in part order:
#' the encoding the array alignments consume (MAFFT's text mode, and
#' [tales_tale_distances()]), where each distinct
#' domain sequence is one "residue".
#'
#' @details
#' This is the sibling of [tales_rvd_strings()] and takes the same two
#' arguments, but **both defaults differ**, because the two projections feed
#' different consumers:
#'
#' \tabular{lll}{
#'   \tab [tales_rvd_strings()] \tab [tales_coded_strings()] \cr
#'   `sep` \tab `"-"`, the AnnoTALE convention \tab `" "`, what MAFFT
#'     `--text` splits on \cr
#'   filter \tab `rvd_only = TRUE` \tab `repeats_only = FALSE` \cr
#' }
#'
#' The separator is free to choose here in a way it is not for RVDs: a
#' `dom_code` is a bare integer rendered as text, so `"1 2 3"` and `"1-2-3"`
#' are equally unambiguous. It defaults to a space because that is what the
#' documented consumers of this encoding expect.
#'
#' The termini are kept by default, where [tales_rvd_strings()] drops them.
#' Target prediction concerns the repeat domain only, so dropping them there
#' is right; alignment is the consumer here, and the two termini are the most
#' reliable anchors an alignment of TALE arrays has. Set `repeats_only = TRUE`
#' to compare bare repeat arrays.
#'
#' @param x A [tales] object carrying a `dom_code` column.
#' @param sep Separator between codes. Defaults to `" "`; see Details.
#' @param repeats_only Drop the terminus parts, keeping only repeats.
#'   `FALSE` by default; see Details. Needs a `domain_type` column.
#' @return A [Biostrings::BStringSet], named by `array_id`.
#' @seealso [tales_domain_codes()] for the code-to-sequence lookup,
#'   [tales_rvd_strings()] for the sibling projection.
#' @export
#' @family tales projections
#' @examples
#' x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
#'                                      package = "tantale"))
#' xa <- tales_assign_domain_codes(x)
#' tales_coded_strings(xa)[1]
#' tales_coded_strings(xa, sep = "-", repeats_only = TRUE)[1]
tales_coded_strings <- function(x, sep = " ", repeats_only = FALSE) {
  .tales_assert_dom_code(x, "tales_coded_strings")
  if (isTRUE(repeats_only)) {
    if (!"domain_type" %in% names(x)) {
      cli::cli_abort(
        c("{.code repeats_only = TRUE} needs a {.field domain_type} column.",
          "i" = "Use {.code repeats_only = FALSE} to render every part."),
        class = c("tantale_error_projection_column", "tantale_error"))
    }
    x <- x[x$domain_type == "repeat", ]
    if (nrow(x) == 0L) {
      cli::cli_abort("No parts left to render.",
                     class = c("tantale_error_projection_empty", "tantale_error"))
    }
  }
  ord <- x[order(x$array_id, x$position_in_array), ]
  strings <- vapply(split(ord$dom_code, ord$array_id),
                    paste, character(1), collapse = sep)
  out <- Biostrings::BStringSet(unname(strings))
  names(out) <- names(strings)
  out
}

#' The domain code lookup table
#'
#' One row per distinct domain sequence: its \code{dom_code}, the amino acid
#' sequence it stands for, and the RVD carried by that domain.
#'
#' @param x A \code{\link{tales}} object carrying a \code{dom_code} column.
#' @return A \code{\link[tibble]{tibble}} with \code{dom_code} and
#'   \code{aa_seq} columns, plus \code{rvd} when \code{x} carries one, one
#'   row per code.
#' @seealso \code{\link{tales_coded_strings}}
#' @export
#' @family tales projections
#' @examples
#' x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
#'                                      package = "tantale"))
#' xa <- tales_assign_domain_codes(x)
#' head(tales_domain_codes(xa))
tales_domain_codes <- function(x) {
  .tales_assert_dom_code(x, "tales_domain_codes")
  .tales_require(x, "tales_domain_codes")
  # rvd is included when present but not required: the dom_code <-> aa_seq
  # correspondence is the substance here, and it is a hard invariant of the
  # class. An object carrying dom_code and aa_seq but no rvd is perfectly
  # valid, and blocking it bought nothing.
  cols <- intersect(c("dom_code", "aa_seq", "rvd"), names(x))
  out <- unique(tibble::as_tibble(x)[cols])
  out[order(out$dom_code), ]
}

#' RVD strings, one per TALE array
#'
#' Renders each array as a separated string of its RVDs in part order — the
#' form target-prediction tools consume.
#'
#' @param x A \code{\link{tales}} object carrying an \code{rvd} column.
#' @param sep Separator between RVDs. Defaults to \code{"-"}, the convention
#'   used by AnnoTALE and by this package's own sample files.
#' @param rvd_only Drop the terminus parts, whose \code{rvd} holds an anchor
#'   code rather than a real RVD (see \code{\link{tales_anchor_codes}}).
#'   \code{TRUE} by default, since target prediction concerns the central
#'   repeat domain only.
#' @return A \code{\link[Biostrings]{BStringSet}}, named by \code{array_id}.
#' @seealso \code{\link{tales_coded_strings}}
#' @export
#' @family tales projections
#' @examples
#' x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
#'                                      package = "tantale"))
#' tales_rvd_strings(x)[1]
#' tales_rvd_strings(x, rvd_only = FALSE)[1] # keeps NTERM/CTERM markers
tales_rvd_strings <- function(x, sep = "-", rvd_only = TRUE) {
  if (!is_tales(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  if (!"rvd" %in% names(x)) {
    cli::cli_abort("{.fn tales_rvd_strings} needs an {.field rvd} column.",
                   class = c("tantale_error_projection_column", "tantale_error"))
  }
  if (isTRUE(rvd_only)) {
    x <- x[!x$rvd %in% tales_anchor_codes(), ]
  }
  if (nrow(x) == 0L) {
    cli::cli_abort("No parts left to render.",
                   class = c("tantale_error_projection_empty", "tantale_error"))
  }
  ord <- x[order(x$array_id, x$position_in_array), ]
  strings <- vapply(split(ord$rvd, ord$array_id),
                    paste, character(1), collapse = sep)
  out <- Biostrings::BStringSet(unname(strings))
  names(out) <- names(strings)
  out
}

#' Guard the precondition both projections share
#'
#' `dom_code` must be present and complete. The completeness half is not
#' pedantry: `paste(dom_code, collapse = " ")` would silently emit the literal
#' string "NA" for a missing code, producing a corrupt alignment residue rather
#' than an error (restructuring-notes.md §1).
#' @noRd
.tales_assert_dom_code <- function(x, fn) {
  if (!is_tales(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  if (!"dom_code" %in% names(x)) {
    cli::cli_abort(
      c("{.fn {fn}} needs a {.field dom_code} column.",
        "i" = "{.fn tales_compare_distal} mints one."),
      class = c("tantale_error_projection_column", "tantale_error")
    )
  }
  if (anyNA(x$dom_code)) {
    cli::cli_abort(
      c("{.field dom_code} must not contain {.val NA}.",
        "i" = "A missing code would be pasted into a coded string as the literal text {.val NA}."),
      class = c("tantale_error_projection_na", "tantale_error")
    )
  }
  invisible(NULL)
}


#' Whole-array protein sequence, one per TALE array
#'
#' Reassembles each array's parts (N-terminus, repeats and C-terminus, in
#' order) into one full-length amino acid sequence, by pasting `aa_seq`
#' together with no separator: the result is a real protein sequence, where
#' [tales_rvd_strings()] joins RVD tokens with hyphens.
#'
#' @details
#' Parts are ordered by `position_in_array`. An alignment never reorders an
#' array's parts, only inserts gaps between them, and a `tales_msa`'s gaps are
#' never rows: a gap is a column with no row for that array. So `position_in_array`
#' gives the same part order `alignment_position` would, without needing to
#' skip anything, and works unchanged whether `x` is a bare `tales` or a
#' `tales_msa`.
#'
#' @param x A [tales] object carrying an `aa_seq` column.
#' @return An [Biostrings::AAStringSet], named by `array_id`.
#' @seealso [tales_get_dna_seq()] for the DNA sequence sibling.
#' @export
#' @family tales projections
#' @examples
#' x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
#'                                      package = "tantale"))
#' tales_get_protein_seq(x)[1]
tales_get_protein_seq <- function(x) {
  seqs <- .tales_assemble_seq(x, "aa_seq", "tales_get_protein_seq")
  out <- Biostrings::AAStringSet(unname(seqs))
  names(out) <- names(seqs)
  out
}

#' Whole-array DNA sequence, one per TALE array
#'
#' The DNA sibling of [tales_get_protein_seq()]: reassembles each array's
#' parts, in order, into one full-length nucleotide sequence, by pasting
#' `dna_seq` together directly. See [tales_get_protein_seq()]'s Details for
#' why `position_in_array` is the ordering key even for a `tales_msa`.
#'
#' @param x A [tales] object carrying a `dna_seq` column.
#' @return A [Biostrings::DNAStringSet], named by `array_id`.
#' @seealso [tales_get_protein_seq()] for the protein sequence sibling.
#' @export
#' @family tales projections
#' @examples
#' x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
#'                                      package = "tantale"))
#' tales_get_dna_seq(x)[1]
tales_get_dna_seq <- function(x) {
  seqs <- .tales_assemble_seq(x, "dna_seq", "tales_get_dna_seq")
  out <- Biostrings::DNAStringSet(unname(seqs))
  names(out) <- names(seqs)
  out
}

#' Shared assembly logic for the whole-array sequence converters
#'
#' `col` is deliberately not concatenated with a separator (unlike
#' [tales_rvd_strings()]/[tales_coded_strings()]): those render token
#' strings for tools that split on a separator, this renders a real
#' biological sequence.
#' @return A named character vector, one string per `array_id`.
#' @noRd
.tales_assemble_seq <- function(x, col, fn) {
  if (!is_tales(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  if (!col %in% names(x)) {
    cli::cli_abort("{.fn {fn}} needs a {.field {col}} column.",
                   class = c("tantale_error_projection_column", "tantale_error"))
  }
  if (anyNA(x[[col]])) {
    cli::cli_abort(
      c("{.field {col}} must not contain {.val NA}.",
        "i" = "A missing sequence would be pasted into the assembled sequence as the literal text {.val NA}.",
        "i" = "{.fn tales_anomalies} lists the arrays concerned; {.code tales(x, sanitize = TRUE)} removes them."),
      class = c("tantale_error_projection_na", "tantale_error")
    )
  }
  ord <- x[order(x$array_id, x$position_in_array), ]
  vapply(split(ord[[col]], ord$array_id), paste, character(1), collapse = "")
}
