#### Printing tales and tales_msa objects ####
#
# Both classes are tibbles underneath, so without these they print as
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


#' Which layer to preview, and a label for it
#' @noRd
.tales_preview_layer <- function(x) {
  available <- intersect(TALES_RESIDUE_COLS, names(x))
  if (!length(available)) return(NULL)
  # dom_code first: it is the more compact of the two, and where both exist
  # the repeat codes are what distinguishes arrays that share an RVD sequence
  if ("dom_code" %in% available) "dom_code" else available[1]
}


#' Truncate a string to the console width, marking that it was cut
#' @noRd
.tales_fit <- function(s, room) {
  room <- max(room, 12L)
  ifelse(nchar(s) > room, paste0(substr(s, 1, room - 4L), " ..."), s)
}


#' The rows to show: the first n and the last n, or everything if it fits
#' @noRd
.tales_preview_rows <- function(ids, n) {
  if (length(ids) <= 2L * n) return(list(head = seq_along(ids), tail = integer()))
  list(head = seq_len(n), tail = seq.int(length(ids) - n + 1L, length(ids)))
}


#' Print one block of the preview table
#' @noRd
.tales_print_block <- function(ids, values, idWidth, room) {
  for (i in seq_along(ids)) {
    cat(sprintf("  %-*s  %s\n", idWidth, ids[i], .tales_fit(values[i], room)))
  }
}


#' Print a tales object
#'
#' @description
#' Shows what the class knows and a tibble would not: the number of arrays as
#' distinct from the number of parts, which residue layers the object carries,
#' the \code{dom_code} namespace when it is stamped, and a preview of the
#' first and last arrays as sequences.
#'
#' The preview uses \code{dom_code} when present and \code{rvd} otherwise.
#' Repeat codes are the more discriminating of the two: two arrays can share
#' an RVD sequence while being built from different repeats.
#'
#' Use \code{tibble::as_tibble(x)} to see the underlying table instead.
#'
#' @param x A \code{\link{tales}} object.
#' @param n Number of arrays to show at each end. Everything is shown when the
#'   object holds \code{2 * n} arrays or fewer.
#' @param ... Unused, present for compatibility with the \code{print} generic.
#' @return \code{x}, invisibly.
#' @method print tales
#' @export
#' @family tales objects
print.tales <- function(x, n = 2L, ...) {
  arrays <- unique(x$array_id)
  # format_inline(), not cli_text(): a print method must write to stdout, and
  # cli_text() emits on the message connection
  cat(cli::format_inline(
    "{.cls {class(x)[1]}} {length(arrays)} array{?s}, {nrow(x)} part{?s}"), "\n", sep = "")

  layers <- intersect(TALES_RESIDUE_COLS, names(x))
  other <- setdiff(names(x), c("array_id", "position_in_array", layers))
  bits <- character()
  if (length(layers)) bits <- c(bits, paste0("layers: ", paste(layers, collapse = ", ")))
  ns <- tales_namespace(x)
  if (!is.null(ns) && nzchar(ns)) bits <- c(bits, paste0("namespace: ", ns))
  if (length(other)) bits <- c(bits, paste0(length(other), " other column",
                                            if (length(other) > 1) "s" else ""))
  if (length(bits)) cat("  ", paste(bits, collapse = "   |   "), "\n", sep = "")

  layer <- .tales_preview_layer(x)
  if (is.null(layer) || !length(arrays)) return(invisible(x))

  ord <- x[order(match(x$array_id, arrays), x$position_in_array), ]
  seqs <- vapply(split(as.character(ord[[layer]]), factor(ord$array_id, levels = arrays)),
                 paste, character(1), collapse = " ")

  idWidth <- min(max(nchar(arrays)), 28L)
  room <- max(getOption("width", 80L) - idWidth - 6L, 20L)
  sel <- .tales_preview_rows(arrays, n)

  cat("  ", strrep(" ", idWidth), "  ", layer, "\n", sep = "")
  .tales_print_block(arrays[sel$head], seqs[sel$head], idWidth, room)
  if (length(sel$tail)) {
    cat(sprintf("  %-*s  ...\n", idWidth, "..."))
    .tales_print_block(arrays[sel$tail], seqs[sel$tail], idWidth, room)
  }
  invisible(x)
}


#' Print a tales_msa object
#'
#' @description
#' As \code{\link{print.tales}}, plus the alignment width, and a preview that
#' shows aligned rows rather than bare sequences: gaps are drawn, so the
#' columns line up and the shape of the alignment is visible.
#'
#' @param x A \code{\link{tales_msa}} object.
#' @param n Number of arrays to show at each end.
#' @param gap What to draw where an array has no residue at a position.
#' @param ... Unused, present for compatibility with the \code{print} generic.
#' @return \code{x}, invisibly.
#' @method print tales_msa
#' @export
#' @family tales objects
print.tales_msa <- function(x, n = 2L, gap = "-", ...) {
  arrays <- unique(x$array_id)
  width <- tales_width(x)
  shown <- width %||% 0L
  cat(cli::format_inline(
    "{.cls {class(x)[1]}} {length(arrays)} array{?s}, {shown} alignment position{?s}"),
    "\n", sep = "")

  layers <- intersect(TALES_RESIDUE_COLS, names(x))
  bits <- character()
  if (length(layers)) bits <- c(bits, paste0("layers: ", paste(layers, collapse = ", ")))
  ns <- tales_namespace(x)
  if (!is.null(ns) && nzchar(ns)) bits <- c(bits, paste0("namespace: ", ns))
  if (length(bits)) cat("  ", paste(bits, collapse = "   |   "), "\n", sep = "")

  layer <- .tales_preview_layer(x)
  if (is.null(layer) || !length(arrays) || is.null(width)) return(invisible(x))

  m <- as.matrix(x, value = layer)
  # pad every cell to the same width so the columns align down the page,
  # which is the point of looking at an alignment at all
  cell <- max(nchar(as.character(m)), nchar(gap), na.rm = TRUE)
  drawn <- apply(m, 1, function(r) {
    paste(formatC(ifelse(is.na(r), gap, r), width = cell), collapse = " ")
  })
  drawn <- drawn[match(arrays, names(drawn))]

  idWidth <- min(max(nchar(arrays)), 28L)
  room <- max(getOption("width", 80L) - idWidth - 6L, 20L)
  sel <- .tales_preview_rows(arrays, n)

  cat("  ", strrep(" ", idWidth), "  ", layer, "\n", sep = "")
  .tales_print_block(arrays[sel$head], drawn[sel$head], idWidth, room)
  if (length(sel$tail)) {
    cat(sprintf("  %-*s  ...\n", idWidth, "..."))
    .tales_print_block(arrays[sel$tail], drawn[sel$tail], idWidth, room)
  }
  invisible(x)
}
