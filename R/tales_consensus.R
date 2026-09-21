#### Consensus of a TALE alignment ####
#
# What was left of msa.R once the plotting moved to tales_plot.R. Both are
# exported and both take a plain matrix rather than a tales_msa, so they are
# usable on any alignment; summary.tales_msa() and plot.tales_msa() are the
# in-package callers.


#' Compute a consensus from a TALE msa
#' @description Pick the most frequent element in each column of the alignment
#'   matrix.
#'
#' @details A column has a consensus only when one element is strictly more
#'   common than every other. Where two or more are tied for most frequent --
#'   as happens whenever each array carries a different repeat at that
#'   position -- the result is \code{NA}, because there is no agreement to
#'   report. \code{NA} is likewise returned when the most common thing at a
#'   position is a gap.
#'
#' @param align A multiple Tal sequences alignment in the form of a
#'   matrix.
#' @return A vector of consensus elements in each column of \code{align}.
#' 
#' @export
#' @family TALE alignment
#' @examples
#' # column 1: HD is a clear majority. column 2: a three-way tie, no consensus.
#' aln <- matrix(c("HD", "HD", "NI",
#'                "NG", "NI", "HD"),
#'              nrow = 3, dimnames = list(c("A1", "A2", "A3"), NULL))
#' aln
#' tales_consensus(aln)
tales_consensus <- function(align) {
  sapply(1:ncol(align), function(x) {
  allElements <- align[,x]
  candidates <- sort(unique(allElements), na.last = TRUE)
  freq <- sapply(candidates, function(p) S4Vectors::countMatches(p, allElements))
  # No consensus unless one element is strictly more common than every other.
  # A column in which each array carries a different repeat has no majority,
  # and reporting one of the tied values would invent agreement that is not
  # there.
  if (sum(freq == max(freq)) > 1L) return(NA_character_)
  candidates[which.max(freq)]
})
}

#' Do elements in a TALE msa match the consensus?
#' @description Compute a logical matrix corresponding to the input \code{align}
#' input with \code{TRUE} where an element matches the consensus at that
#' position and \code{FALSE} where it does not. Columns with no consensus --
#' see \code{\link{tales_consensus}} -- are \code{NA} throughout, since
#' there is nothing there to match.
#'
#' @param align A multiple Tal sequences alignment in the form of a
#'   matrix.
#' @param long Set to \code{TRUE} (default) to return a long tibble, or
#'   \code{FALSE} to return a logical matrix with the same shape as
#'   \code{align}.
#' @return A multiple Tal sequences alignment in the form of a
#'   matrix filled with logical values if \code{long} is \code{FALSE} and
#'   a long tibble representing the original alignment otherwise (default).
#' 
#' @export
#' @family TALE alignment
#' @examples
#' aln <- matrix(c("HD", "HD", "NI",
#'                "NG", "NI", "HD"),
#'              nrow = 3, dimnames = list(c("A1", "A2", "A3"), NULL))
#' tales_consensus_match(aln, long = FALSE)
#' tales_consensus_match(aln)
tales_consensus_match <- function(align, long = TRUE) {
  consensus <- tales_consensus(align)
  # A logical matrix of its own, rather than overwriting the character one:
  # assigning TRUE into a character matrix stores "TRUE", so the result used
  # to be strings, and sum()/which()/! on it did the wrong thing.
  out <- matrix(FALSE, nrow = nrow(align), ncol = ncol(align),
                dimnames = dimnames(align))
  for (k in seq_len(ncol(align))) {
    rept <- consensus[k]
    if (is.na(rept)) {
      # The column has no consensus, so "does this match it?" has no answer.
      out[, k] <- NA
    } else {
      # A gap is not a match.
      out[, k] <- !is.na(align[, k]) & toupper(align[, k]) == toupper(rept)
    }
  }
  align <- out
  if (!long) return(align)
  matchConsensusLong <- .matrix_to_long(align) %>%
    dplyr::as_tibble()
  colnames(matchConsensusLong) <- c("array_id", "alignment_position", "tales_consensus_match")
  return(matchConsensusLong)
}


#### tales_msa-native counterparts ####
#
# tales_consensus()/tales_consensus_match() take a matrix, the natural shape
# when this file was written. A tales_msa is a long tibble now, and
# plot.tales_msa() already holds one -- routing it through as.matrix() and
# back just to get a consensus was pure overhead (ledger section on
# plot.tales_msa()'s refactor). These two remove that overhead: same
# semantics, a tales_msa in, a long tibble out. Not exported yet, but written
# to the same standard as their matrix-based counterparts above, since they
# are candidates to eventually replace them as the public API -- see the
# ledger for the open question of whether tales_consensus() itself should
# move to this shape.
#
# Deliberately NOT a generic long-table function: it takes x, a real
# tales_msa, not a caller-assembled long table. A tales_msa's own gap
# semantics are an *absent* row, never an explicit NA one (unlike a matrix,
# which always has every cell) -- a first version of this took a generic
# table instead and pushed the job of building a "complete" grid onto the
# caller, which is exactly backwards: get it wrong once here, and every
# caller is safe, rather than every caller needing to get it right.
#
# Both check x's type themselves, unlike as.matrix.tales_msa() (which
# leans on S3 dispatch for that): these are plain functions, not methods,
# so nothing stops a bad x from reaching them if they are ever exported.

#' Guard the precondition both tales_msa-native consensus functions share
#' @noRd
.assert_tales_msa_layer <- function(x, value_col) {
  if (!is_tales_msa(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales_msa} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  if (!value_col %in% names(x)) {
    cli::cli_abort(
      c("This alignment has no {.field {value_col}} layer.",
        "i" = "Available: {.field {setdiff(names(x), c('array_id', 'alignment_position'))}}"),
      class = c("tantale_error_msa_layer", "tantale_error")
    )
  }
  invisible(NULL)
}

#' Compute a consensus from a tales_msa
#'
#' The \code{\link{tales_consensus}}-equivalent for a \code{\link{tales_msa}}
#' object instead of a matrix. Written so a caller that already holds the
#' alignment as a \code{tales_msa} -- \code{\link{plot.tales_msa}} chief
#' among them -- never needs to build a matrix just to ask "what does most
#' arrays carry at this position".
#'
#' @details
#' Same rule as \code{\link{tales_consensus}}, applied per alignment position
#' instead of per matrix column: a position has a consensus only when one
#' value is strictly more common than every other. A tie -- including a tie
#' with the gap itself -- reports \code{NA}, and so does a position where a
#' gap is outright the most common thing.
#'
#' A \code{tales_msa} has no row at all for a gap (unlike a matrix, which
#' always has a cell, \code{NA} or not), so a gap's count at a position is
#' recovered as \code{(number of arrays) - (number of real rows at that
#' position)}, not read off any column -- getting this arithmetic right here
#' is the whole reason this takes \code{x} directly rather than a
#' caller-assembled long table (see the file header).
#'
#' @param x A \code{\link{tales_msa}} object.
#' @param value_col Name of the residue column to take the consensus of
#'   (e.g. \code{"dom_code"} or \code{"rvd"}).
#' @return A tibble with one row per alignment position, columns
#'   \code{alignment_position} and \code{value_col}, the latter holding the
#'   consensus (\code{NA} where there is none).
#' @seealso \code{\link{tales_consensus}}, the matrix-based original this
#'   mirrors exactly; \code{\link{.tales_consensus_match_msa}}, its
#'   match-flag counterpart.
#' @noRd
.tales_consensus_long <- function(x, value_col) {
  .assert_tales_msa_layer(x, value_col)
  arrays <- unique(x$array_id)
  n_arrays <- length(arrays)
  n_positions <- tales_width(x) %||% max(x$alignment_position, na.rm = TRUE)
  values_by_position <- split(x[[value_col]], x$alignment_position)

  consensus <- vapply(as.character(seq_len(n_positions)), function(p) {
    present <- values_by_position[[p]]
    if (is.null(present)) present <- character(0)
    n_gap <- n_arrays - length(present)

    candidates <- sort(unique(present))
    freq <- vapply(candidates, function(v) sum(present == v), integer(1))
    if (n_gap > 0L) {
      candidates <- c(candidates, NA_character_)
      freq <- c(freq, n_gap)
    }
    if (!length(freq) || sum(freq == max(freq)) > 1L) return(NA_character_)
    candidates[which.max(freq)]
  }, character(1), USE.NAMES = FALSE)

  tibble::tibble(alignment_position = seq_len(n_positions), !!value_col := consensus)
}


#' Do elements of a tales_msa match its consensus?
#'
#' The \code{\link{tales_consensus_match}}-equivalent for a
#' \code{\link{tales_msa}} object, using \code{\link{.tales_consensus_long}}
#' rather than \code{\link{tales_consensus}} to get the consensus it compares
#' against. Same semantics: a gap never matches, and a position with no
#' consensus reports \code{NA} for every array at that position, not
#' \code{FALSE} -- which means every array needs a row at every position in
#' the result, including the ones \code{x} itself has no row for (a gap);
#' this materialises that full grid, x does not need to already have it.
#'
#' @param x A \code{\link{tales_msa}} object.
#' @param value_col As \code{\link{.tales_consensus_long}}.
#' @param consensus A precomputed result of \code{\link{.tales_consensus_long}},
#'   to avoid recomputing it when the caller already has one. Computed from
#'   \code{x} when not supplied.
#' @return A tibble with \code{array_id}, \code{alignment_position} and a
#'   logical \code{match} column, one row per array per position.
#' @seealso \code{\link{.tales_consensus_long}}
#' @noRd
.tales_consensus_match_long <- function(x, value_col, consensus = NULL) {
  .assert_tales_msa_layer(x, value_col)
  arrays <- unique(x$array_id)
  n_positions <- tales_width(x) %||% max(x$alignment_position, na.rm = TRUE)
  if (is.null(consensus)) {
    consensus <- .tales_consensus_long(x, value_col)
  }

  grid <- tidyr::expand_grid(array_id = arrays,
                             alignment_position = seq_len(n_positions))
  grid <- dplyr::left_join(grid, x[c("array_id", "alignment_position", value_col)],
                           by = c("array_id", "alignment_position"))

  refByPosition <- stats::setNames(consensus[[value_col]],
                                   as.character(consensus$alignment_position))
  rept <- unname(refByPosition[as.character(grid$alignment_position)])
  grid$match <- ifelse(
    is.na(rept), NA,
    !is.na(grid[[value_col]]) & toupper(grid[[value_col]]) == toupper(rept)
  )
  grid[c("array_id", "alignment_position", "match")]
}
