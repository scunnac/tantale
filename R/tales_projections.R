##### Projections of a 'tales' object #####
# Views that were formerly stored as fields of distalr()'s returned list
# ('coded.repeats.str', 'repeats.code'). Both are pure functions of the parts
# table, so storing them invited the copy and the source drifting apart; see
# restructuring-notes.md §1.
#
# Both need `dom_code`, which is minted by tales_relatedness().


#' Domain-coded strings, one per TALE array
#'
#' Renders each array as a space-separated string of its \code{dom_code}s in
#' part order — the encoding ARLEM and MAFFT's text mode consume, where each
#' distinct domain sequence is one "residue".
#'
#' @param x A \code{\link{tales}} object carrying a \code{dom_code} column.
#' @return A \code{\link[Biostrings]{BStringSet}}, named by \code{array_id}.
#' @seealso \code{\link{tales_domain_codes}} for the code-to-sequence lookup.
#' @export
tales_coded_strings <- function(x) {
  .tales_assert_dom_code(x, "tales_coded_strings")
  ord <- x[order(x$array_id, x$position_in_array), ]
  strings <- vapply(split(ord$dom_code, ord$array_id),
                    paste, character(1), collapse = " ")
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
#' @return A \code{\link[tibble]{tibble}} with \code{dom_code}, \code{aa_seq}
#'   and \code{rvd} columns, one row per code.
#' @seealso \code{\link{tales_coded_strings}}
#' @export
tales_domain_codes <- function(x) {
  .tales_assert_dom_code(x, "tales_domain_codes")
  needed <- c("aa_seq", "rvd")
  missing <- setdiff(needed, names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      "{.fn tales_domain_codes} needs the column{?s} {.field {missing}}.",
      class = c("tantale_error_projection_column", "tantale_error")
    )
  }
  out <- unique(tibble::as_tibble(x)[c("dom_code", "aa_seq", "rvd")])
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
        "i" = "{.fn tales_relatedness} mints one."),
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
