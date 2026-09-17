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
  matchConsensusLong <- align %>% reshape2::melt() %>%
    dplyr::as_tibble()
  colnames(matchConsensusLong) <- c("array_id", "position_in_array", "tales_consensus_match")
  return(matchConsensusLong)
}
