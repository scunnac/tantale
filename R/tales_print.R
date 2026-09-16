#### Rendering tales and tales_msa objects ####
#
# Both classes are tibbles underneath, so without these they render as
# "# A tibble: 955 x 10" -- which says nothing about what the object is. The
# things a reader actually needs are the ones the class knows and a tibble
# does not: how many arrays, as opposed to how many parts; which residue
# layers are present, since that decides what tales_align() and plot() can do;
# and the dom_code namespace stamp, whose whole purpose is to catch tables
# from different runs being mixed and which is otherwise invisible until
# something fails.
#
# The sequence preview follows Biostrings: the first few and the last few,
# with the middle elided. An alignment of TALEs is a thing you recognise by
# looking at it, and a count of rows is not.
#
# format() builds the lines and print() emits them, rather than print() doing
# both at once. Both classes inherit a format() method from tibble, so without
# this the two would disagree about what the object looks like: print() would
# show the view below and format() the tibble one.


#' Which layer to preview
#' @noRd
.tales_preview_layer <- function(x) {
  available <- intersect(TALES_RESIDUE_COLS, names(x))
  if (!length(available)) return(NULL)
  # dom_code first: where both exist the repeat codes are what distinguishes
  # arrays that happen to share an RVD sequence
  if ("dom_code" %in% available) "dom_code" else available[1]
}


#' Truncate a string to the console width, marking that it was cut
#' @noRd
.tales_fit <- function(s, room) {
  room <- max(room, 12L)
  ifelse(nchar(s) > room, paste0(substr(s, 1, room - 4L), " ..."), s)
}


#' The header lines common to both classes
#' @noRd
.tales_header <- function(x, summary_line) {
  layers <- intersect(TALES_RESIDUE_COLS, names(x))
  bits <- character()
  if (length(layers)) bits <- c(bits, paste0("layers: ", paste(layers, collapse = ", ")))
  ns <- tales_namespace(x)
  if (!is.null(ns) && nzchar(ns)) bits <- c(bits, paste0("namespace: ", ns))
  other <- setdiff(names(x), c("array_id", "position_in_array", "alignment_position", layers))
  if (length(other)) bits <- c(bits, paste0(length(other), " other column",
                                            if (length(other) > 1) "s" else ""))
  c(summary_line,
    if (length(bits)) paste0("  ", paste(bits, collapse = "   |   ")))
}


#' The preview table: the first n arrays, an ellipsis, the last n
#'
#' @param ids Array names, in the order to show them.
#' @param values One rendered sequence per array, already padded where the
#'   caller wants columns to line up.
#' @noRd
.tales_preview <- function(ids, values, layer, n) {
  if (!length(ids)) return(character())
  idWidth <- min(max(nchar(ids)), 28L)
  room <- max(getOption("width", 80L) - idWidth - 6L, 20L)

  block <- function(i) sprintf("  %-*s  %s", idWidth, ids[i], .tales_fit(values[i], room))
  elide <- length(ids) > 2L * n
  head <- if (elide) seq_len(n) else seq_along(ids)
  tail <- if (elide) seq.int(length(ids) - n + 1L, length(ids)) else integer()

  c(paste0("  ", strrep(" ", idWidth), "  ", layer),
    block(head),
    if (length(tail)) c(sprintf("  %-*s  ...", idWidth, "..."), block(tail)))
}


#' Render a tales object as lines of text
#'
#' @description
#' The representation \code{\link{print.tales}} emits, returned as a character
#' vector rather than written out -- so an object's description can go
#' somewhere other than the console: a log, an error message, a report.
#'
#' It shows what the class knows and a tibble would not: the number of arrays
#' as distinct from the number of parts, which residue layers the object
#' carries, the \code{dom_code} namespace when it is stamped, and a preview of
#' the first and last arrays as sequences.
#'
#' The preview uses \code{dom_code} when present and \code{rvd} otherwise.
#' Repeat codes are the more discriminating of the two: two arrays can share
#' an RVD sequence while being built from different repeats.
#'
#' This replaces the \code{format()} a \code{tales} would otherwise inherit
#' from tibble. Call \code{format(tibble::as_tibble(x))} for that one.
#'
#' @param x A \code{\link{tales}} object.
#' @param n Number of arrays to show at each end. Everything is shown when the
#'   object holds \code{2 * n} arrays or fewer.
#' @param ... Unused, present for compatibility with the \code{format} generic.
#' @return A character vector, one element per line.
#' @method format tales
#' @export
#' @family tales objects
format.tales <- function(x, n = 2L, ...) {
  arrays <- unique(x$array_id)
  out <- .tales_header(x, cli::format_inline(
    "{.cls {class(x)[1]}} {length(arrays)} array{?s}, {nrow(x)} part{?s}"))

  layer <- .tales_preview_layer(x)
  if (is.null(layer) || !length(arrays)) return(out)

  ord <- x[order(match(x$array_id, arrays), x$position_in_array), ]
  seqs <- vapply(split(as.character(ord[[layer]]), factor(ord$array_id, levels = arrays)),
                 paste, character(1), collapse = " ")
  c(out, .tales_preview(arrays, seqs, layer, n))
}


#' Render a tales_msa object as lines of text
#'
#' @description
#' As \code{\link{format.tales}}, plus the alignment width, and a preview that
#' shows aligned rows rather than bare sequences: gaps are drawn and every
#' cell padded to a common width, so the columns line up down the page.
#'
#' @param x A \code{\link{tales_msa}} object.
#' @param n Number of arrays to show at each end.
#' @param gap What to draw where an array has no residue at a position.
#' @param ... Unused, present for compatibility with the \code{format} generic.
#' @return A character vector, one element per line.
#' @method format tales_msa
#' @export
#' @family tales objects
format.tales_msa <- function(x, n = 2L, gap = "-", ...) {
  arrays <- unique(x$array_id)
  width <- tales_width(x)
  positions <- width %||% 0L
  out <- .tales_header(x, cli::format_inline(
    "{.cls {class(x)[1]}} {length(arrays)} array{?s}, {positions} alignment position{?s}"))

  layer <- .tales_preview_layer(x)
  if (is.null(layer) || !length(arrays) || is.null(width)) return(out)

  m <- as.matrix(x, value = layer)
  # pad every cell to the same width so the columns align down the page,
  # which is the point of looking at an alignment at all
  cell <- max(nchar(as.character(m)), nchar(gap), na.rm = TRUE)
  drawn <- apply(m, 1, function(r) {
    paste(formatC(ifelse(is.na(r), gap, r), width = cell), collapse = " ")
  })
  c(out, .tales_preview(arrays, drawn[match(arrays, names(drawn))], layer, n))
}


#' Print a tales object
#'
#' @description
#' Emits what \code{\link{format.tales}} builds; see there for what is shown
#' and why.
#'
#' @inheritParams format.tales
#' @return \code{x}, invisibly.
#' @method print tales
#' @export
#' @family tales objects
print.tales <- function(x, n = 2L, ...) {
  cat(format(x, n = n, ...), sep = "\n")
  invisible(x)
}


#' Print a tales_msa object
#'
#' @description
#' Emits what \code{\link{format.tales_msa}} builds; see there for what is
#' shown and why.
#'
#' @inheritParams format.tales_msa
#' @return \code{x}, invisibly.
#' @method print tales_msa
#' @export
#' @family tales objects
print.tales_msa <- function(x, n = 2L, gap = "-", ...) {
  cat(format(x, n = n, gap = gap, ...), sep = "\n")
  invisible(x)
}
