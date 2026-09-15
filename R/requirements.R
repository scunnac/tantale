#### What each consumer needs of a tales object ####
#
# The column contract in validate_tales() answers "is this a valid tales?".
# This table answers a different question: "what does *this operation* need of
# one?". The two are not the same, and the gap between them is real -- aa_seq is
# optional in the contract, because a tales built from RVD strings has none,
# yet tales_compare() cannot work without it (or dna_seq to translate).
#
# Before this existed, every function invented its own guard. Two of them
# checked columns they never read, one reported an "output tibble" it was about
# to not produce, and plot_tales_composition() checked nothing at all and simply
# broke. Stating the requirements in one place is what makes them uniform.
#
# all_of  : every column must be present
# any_of  : at least one must be present (the "or" the contract already uses)
# optional: used when present, never required -- documented so it is a promise

TALES_REQUIREMENTS <- list(
  tales_compare = list(
    any_of = c("aa_seq", "dna_seq"),
    note   = "dna_seq is translated to aa_seq when aa_seq is absent"
  ),
  tales_rvd_strings      = list(all_of = "rvd"),
  tales_coded_strings    = list(all_of = "dom_code"),
  tales_domain_codes     = list(all_of = c("dom_code", "aa_seq"),
                                optional = "rvd"),
  tales_align            = list(any_of = c("rvd", "dom_code"),
                                note = "whichever residue_col names"),
  plot.tales             = list(all_of = c("rvd", "aa_seq", "domain_type"),
                                optional = c("seqnames", "alignment_position"),
                                note = "seqnames adds a facet; alignment_position enables position = \"alignment\""),
  tale_parts_to_rvd      = list(all_of = "rvd"),
  repeat_to_rvd_map_distalr = list(all_of = c("dom_code", "rvd"))
)


#' Check a tales against the documented requirements of a consumer
#'
#' Reads \code{TALES_REQUIREMENTS} so every function reports a missing column
#' the same way, rather than each inventing its own guard.
#'
#' @param x A \code{tales} object.
#' @param fn Name of the consumer, as a key of \code{TALES_REQUIREMENTS}.
#' @return \code{x}, invisibly. Aborts if a requirement is unmet.
#' @keywords internal
.tales_require <- function(x, fn) {
  req <- TALES_REQUIREMENTS[[fn]]
  if (is.null(req)) {
    cli::cli_abort("No requirements recorded for {.fn {fn}}.",
                   class = c("tantale_error_requirements", "tantale_error"))
  }
  missing <- setdiff(req$all_of, names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      "{.fn {fn}} needs the column{?s} {.field {missing}}.",
      class = c("tantale_error_projection_column", "tantale_error")
    )
  }
  if (!is.null(req$any_of) && !any(req$any_of %in% names(x))) {
    cli::cli_abort(
      c("{.fn {fn}} needs {.field {req$any_of}}.",
        "i" = "At least one of those columns must be present.",
        if (!is.null(req$note)) c("i" = req$note)),
      class = c("tantale_error_projection_column", "tantale_error")
    )
  }
  invisible(x)
}


#' Column requirements of the tales consumers
#'
#' @description
#' What each function needs of a \code{\link{tales}} object. This is distinct
#' from the class's own column contract, which says what makes an object
#' *valid*: a valid \code{tales} may still lack what a given operation needs.
#'
#' @details
#' \describe{
#'   \item{\code{tales_compare()}}{\code{aa_seq} or \code{dna_seq} (translated)}
#'   \item{\code{tales_align()}}{whichever of \code{rvd} / \code{dom_code} \code{residue_col} names}
#'   \item{\code{tales_rvd_strings()}}{\code{rvd}}
#'   \item{\code{tales_coded_strings()}}{\code{dom_code}}
#'   \item{\code{tales_domain_codes()}}{\code{dom_code} and \code{aa_seq}; \code{rvd} included when present}
#'   \item{\code{plot()} on a tales}{\code{rvd}, \code{aa_seq}, \code{domain_type}; \code{seqnames} adds a facet, \code{alignment_position} enables the aligned layout}
#'   \item{\code{tale_parts_to_rvd()}}{\code{rvd}}
#'   \item{\code{repeat_to_rvd_map_distalr()}}{\code{dom_code}, \code{rvd}}
#' }
#'
#' @return A tibble of function, requirement kind, and columns.
#' @export
#' @family tales objects
tales_requirements <- function() {
  rows <- lapply(names(TALES_REQUIREMENTS), function(fn) {
    req <- TALES_REQUIREMENTS[[fn]]
    kinds <- c("all_of", "any_of", "optional")
    present <- kinds[vapply(kinds, function(k) !is.null(req[[k]]), logical(1))]
    tibble::tibble(
      fn = fn,
      requirement = present,
      columns = vapply(present, function(k) paste(req[[k]], collapse = ", "), character(1))
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$fn, out$requirement), ]
}
