



#' Split TALE sequence strings into vectors
#'
#' Implementation behind \code{\link{as_tales}} and the deprecated
#' \code{\link{as_tales}}. Internal so that package code can call it without
#' tripping the deprecation warning.
#' @noRd
.split_list <- function(strings, sep = "-") {
  if (is.list(strings) &&
      any(sapply(strings, length) > 1)) {
      cli::cli_abort("The value provided for strings seems already to be splitted.", class = c("tantale_error"))
  } else if (length(strings) == 1 && is.character(strings)) {
    stopifnot(fs::file_exists(strings))
    seqs <- as.character(Biostrings::readBStringSet(strings), use.names = TRUE)
  } else if (class(strings) %in% c("AAStringSet", "BStringSet")) {
    seqs <- as.character(strings, use.names = TRUE)
  } else if (length(strings) >= 1 && is.list(strings)) {
    seqs <- strings
  } else {
    cli::cli_abort("Something is wrong with the value provided for strings.", class = c("tantale_error"))
  }
  
  seqsAsVectors <- stringr::str_split(seqs, pattern = glue::glue("[{sep}]"))
  seqsAsVectors <- lapply(seqsAsVectors, function(x) { # Remove last residue if it is empty string
    if ( x[length(x)] == "") {
      cli::cli_warn("Dropped an empty last element from a vectorized sequence.",
                    class = "tantale_warning_empty_element")
      x[-length(x)]
    } else x
  }
  )
  names(seqsAsVectors) <- names(seqs)
  return(seqsAsVectors)
}



#### Repeat-code / RVD alignment matrix conversion ####
# Not called by production code (a tales_msa carries both layers at once, so
# plot() names them directly -- fill = "dom_code", label = "rvd"). Kept
# because it builds fixture data in test_plot_tales_msa.R. The inverse
# direction, .rvd_to_repeat_align(), is in inst/legacy/conversion_retired.R
# (§37).

#' Substitute Distal repeat IDs for RVDs in a TALE alignment matrix
#'
#' @param repeat_align A multiple TALE repeat sequences alignment in the form
#'   of a matrix, as returned by \code{\link{tales_align}}.
#' @param rvd_map A data frame mapping each repeat code (\code{repeatID})
#'   to its RVD (\code{RVD}).
#' @return A TALE alignment matrix made up of RVD sequences.
#' @noRd
.repeat_to_rvd_align <- function(repeat_align, rvd_map) {
  states <- unique(as.vector(repeat_align))
  ##### TODO: check that all values in states are present in the rvd_map df ####
  # If not, error
  rvd_align <- t(
    apply(repeat_align, 1,
          function(repeatSeq){
            rvdSeq <- rvd_map$RVD[match(repeatSeq, rvd_map$repeatID)]
          }
    )
  )
  rvd_align <- matrix(rvd_align, nrow = nrow(repeat_align)) # in case of 1-row matrix
  rownames(rvd_align) <- rownames(repeat_align)
  colnames(rvd_align) <- colnames(repeat_align)
  return(rvd_align)
}
