##### The 'tales_msa' class #####
# A 'tales' that additionally carries an alignment coordinate. Gaps are
# *implicit*: a gap at (array, alignment_position) is the absence of a row, so
# every row remains a real part and every 'tales' invariant holds unchanged.
# See dev/class-design.md §4.


#### Constructors ####

#' Low-level constructor for a tales_msa object
#'
#' Attaches the class without validating. Use \code{\link{tales_msa}} unless the
#' invariants are already established.
#'
#' @param x A \code{tales} object carrying an \code{alignment_position} column.
#' @param alignment_width Integer width of the alignment. Stored as an
#'   attribute rather than derived, because subsetting arrays can empty the last
#'   column and would silently shrink a derived value.
#' @return A \code{tales_msa} object.
#' @keywords internal
new_tales_msa <- function(x, alignment_width = NULL) {
  stopifnot(is.data.frame(x))
  if (!is_tales(x)) x <- new_tales(x)
  if (!is.null(alignment_width)) {
    attr(x, "alignment_width") <- as.integer(alignment_width)
  }
  class(x) <- unique(c("tales_msa", class(x)))
  x
}

#' Is this a tales_msa object?
#' @param x An object.
#' @return A logical scalar.
#' @export
is_tales_msa <- function(x) inherits(x, "tales_msa")

#' Width of a TALE alignment
#'
#' The number of columns in the alignment, including those that are all gaps in
#' the object at hand. Stored rather than derived: subsetting arrays can empty
#' the last column, which would silently shrink \code{max(alignment_position)}.
#'
#' @param x A \code{tales_msa} object.
#' @return An integer scalar, or \code{NULL} if unset.
#' @export
tales_width <- function(x) {
  attr(x, "alignment_width", exact = TRUE)
}

#' Create a tales_msa object
#'
#' A \code{\link{tales}} object plus an \code{alignment_position} column. Gaps
#' are implicit: a gap is simply the absence of a row at that
#' (\code{array_id}, \code{alignment_position}).
#'
#' Normally produced by \code{\link{tales_align}} rather than called directly.
#'
#' @param x A data frame with the \code{\link{tales}} columns plus
#'   \code{alignment_position}.
#' @param alignment_width Integer alignment width; defaults to
#'   \code{max(alignment_position)}.
#' @param dom_code_namespace Optional scalar string, see
#'   \code{\link{tales_namespace}}.
#' @return A validated \code{tales_msa} object.
#' @export
tales_msa <- function(x, alignment_width = NULL, dom_code_namespace = NULL) {
  x <- tales(x, dom_code_namespace = dom_code_namespace)
  if ("alignment_position" %in% names(x) && is.numeric(x$alignment_position)) {
    x$alignment_position <- as.integer(x$alignment_position)
  }
  if (is.null(alignment_width) && "alignment_position" %in% names(x) &&
      nrow(x) > 0L) {
    alignment_width <- max(x$alignment_position, na.rm = TRUE)
  }
  validate_tales_msa(new_tales_msa(x, alignment_width = alignment_width))
}


#### Validator ####

#' Validate a tales_msa object
#'
#' Checks every \code{\link{validate_tales}} invariant, then those specific to
#' an alignment. As for \code{tales}, only properties closed under row
#' subsetting are checked here; grid completeness is a precondition of the
#' functions that need it.
#'
#' @param x A \code{tales_msa} object.
#' @return \code{x}, invisibly, if valid; otherwise an error.
#' @export
validate_tales_msa <- function(x) {
  validate_tales(x)

  if (!"alignment_position" %in% names(x)) {
    cli::cli_abort(
      "A {.cls tales_msa} object requires an {.field alignment_position} column.",
      class = c("tantale_error_tales_missing_column", "tantale_error")
    )
  }
  .tales_check_type(x, "alignment_position", is.integer, "an integer vector")

  if (nrow(x) == 0L) return(invisible(x))

  ## 11. positive, never NA
  if (anyNA(x$alignment_position) || any(x$alignment_position < 1L)) {
    cli::cli_abort(
      "{.field alignment_position} must be a positive integer without {.val NA}.",
      class = c("tantale_error_msa_position", "tantale_error")
    )
  }

  ## 12. unique within an array
  if (any(duplicated(data.frame(a = x$array_id, p = x$alignment_position)))) {
    cli::cli_abort(
      "{.field alignment_position} must be unique within an array.",
      class = c("tantale_error_msa_duplicate", "tantale_error")
    )
  }

  ## 13. an alignment may insert gaps but never reorder parts
  .tales_msa_check_order(x)

  ## 14. within the declared width
  w <- tales_width(x)
  if (!is.null(w) && any(x$alignment_position > w)) {
    cli::cli_abort(
      c("{.field alignment_position} must not exceed the alignment width ({w}).",
        "x" = "Found up to {max(x$alignment_position)}."),
      class = c("tantale_error_msa_width", "tantale_error")
    )
  }

  invisible(x)
}

#' Check that the alignment preserves part order
#'
#' Within an array, ranking by \code{alignment_position} must agree with ranking
#' by \code{position_in_array}: an alignment inserts gaps, it never reorders.
#' @noRd
.tales_msa_check_order <- function(x) {
  o <- x[order(x$array_id, x$position_in_array), c("array_id", "alignment_position")]
  bad <- unique(o$array_id[
    stats::ave(o$alignment_position, o$array_id, FUN = function(z) c(0L, diff(z))) < 0L
  ])
  if (length(bad) > 0L) {
    cli::cli_abort(
      c("An alignment may insert gaps but must not reorder parts.",
        "x" = "{.field alignment_position} disagrees with {.field position_in_array} in {.val {utils::head(bad, 5)}}."),
      class = c("tantale_error_msa_order", "tantale_error")
    )
  }
  invisible(NULL)
}

#' @noRd
.tales_msa_contract_holds <- function(x) {
  "alignment_position" %in% names(x) && is.integer(x$alignment_position)
}


#### Views ####

#' Render a TALE alignment as a matrix
#'
#' Materialises the rectangular form: one row per array, one column per
#' alignment position, cells holding the requested layer. This is the only place
#' the gapped matrix is built — the long object is the canonical storage.
#'
#' @param x A \code{tales_msa} object.
#' @param value Name of the column to fill cells with. Defaults to the first
#'   available of \code{rvd}, \code{dom_code}.
#' @param gap Value to use for gaps. Defaults to \code{NA}.
#' @param ... Ignored.
#' @return A character matrix with arrays as rows and alignment positions as
#'   columns.
#' @method as.matrix tales_msa
#' @export
as.matrix.tales_msa <- function(x, value = NULL, gap = NA, ...) {
  if (is.null(value)) {
    value <- intersect(TALES_RESIDUE_COLS, names(x))[1]
  }
  if (!value %in% names(x)) {
    cli::cli_abort(
      c("This alignment has no {.field {value}} layer.",
        "i" = "Available: {.field {setdiff(names(x), c('array_id', 'alignment_position'))}}"),
      class = c("tantale_error_msa_layer", "tantale_error")
    )
  }
  arrays <- unique(x$array_id)
  width <- tales_width(x) %||% max(x$alignment_position, na.rm = TRUE)

  m <- matrix(gap, nrow = length(arrays), ncol = width,
              dimnames = list(arrays, as.character(seq_len(width))))
  idx <- cbind(match(x$array_id, arrays), x$alignment_position)
  m[idx] <- as.character(x[[value]])
  m
}

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#### Alignment ####

#' Align the repeat arrays of a tales object
#'
#' Aligns TALE arrays on one of their residue layers with MAFFT's text mode,
#' returning the alignment as a \code{\link{tales_msa}} — the input object with
#' an \code{alignment_position} column added, so every other layer
#' (\code{rvd}, \code{dom_code}, \code{aa_seq}, ...) remains available.
#'
#' The mapping back from MAFFT's output is **positional**: the k-th non-gap cell
#' of an aligned row is the k-th part fed in. That is well defined only because
#' this function builds MAFFT's input from \code{x} itself, which is why
#' \code{\link{tales_assert_complete}} is enforced first.
#'
#' @param x A \code{\link{tales}} object holding complete arrays.
#' @param residue_col Which layer to align on: \code{"rvd"} (default) or
#'   \code{"dom_code"}. Given explicitly rather than guessed from the values.
#' @param repeat_sims Optional similarity table passed to MAFFT as a scoring
#'   matrix, as accepted by \code{\link{build_repeat_msa}}.
#' @param ... Passed to \code{\link{build_repeat_msa}} (e.g. \code{mafft_opts}).
#' @return A \code{tales_msa} object.
#' @export
tales_align <- function(x, residue_col = c("rvd", "dom_code"),
                        repeat_sims = NULL, ...) {
  residue_col <- match.arg(residue_col)
  if (!is_tales(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  if (!residue_col %in% names(x)) {
    cli::cli_abort(
      c("{.arg x} has no {.field {residue_col}} column to align on.",
        "i" = "Available residue column{?s}: {.field {intersect(TALES_RESIDUE_COLS, names(x))}}"),
      class = c("tantale_error_msa_layer", "tantale_error")
    )
  }
  if (nrow(x) == 0L) {
    cli::cli_abort("{.arg x} has no parts to align.",
                   class = c("tantale_error_msa_empty", "tantale_error"))
  }
  tales_assert_complete(x, arg = "x")

  ## One space-separated string per array, in position_in_array order.
  ## A space is safe for both layers: neither RVDs nor repeat codes contain one.
  ord <- x[order(x$array_id, x$position_in_array), ]
  seqs <- split(as.character(ord[[residue_col]]), ord$array_id)
  seqs <- lapply(seqs, paste, collapse = " ")

  m <- build_repeat_msa(input_seqs = seqs, sep = " ", repeat_sims = repeat_sims,
                        gap_symbol = NA, ...)

  if (!setequal(rownames(m), unique(x$array_id))) {
    cli::cli_abort(
      "MAFFT returned a different set of arrays than was submitted.",
      class = c("tantale_error_msa_backmap", "tantale_error")
    )
  }

  ## Positional back-mapping. MAFFT is run with --reorder by default, so rows
  ## must be matched by name, never by position.
  mapping <- lapply(rownames(m), function(a) {
    nz <- which(!is.na(m[a, ]))
    parts <- sort(x$position_in_array[x$array_id == a])
    if (length(nz) != length(parts)) {
      cli::cli_abort(
        c("Alignment of {.val {a}} returned {length(nz)} residue{?s} for {length(parts)} part{?s}.",
          "i" = "The back-mapping is positional and requires them to agree."),
        class = c("tantale_error_msa_backmap", "tantale_error")
      )
    }
    tibble::tibble(array_id = a, position_in_array = parts,
                   alignment_position = as.integer(nz))
  })
  mapping <- dplyr::bind_rows(mapping)

  out <- dplyr::left_join(x, mapping, by = c("array_id", "position_in_array"))
  tales_msa(out, alignment_width = ncol(m),
            dom_code_namespace = tales_namespace(x))
}
