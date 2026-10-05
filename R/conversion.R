



#' Split TALE sequence strings into vectors
#'
#' Implementation behind \code{\link{as_tales}} and the deprecated
#' \code{\link{as_tales}}. Internal so that package code can call it without
#' tripping the deprecation warning.
#' @noRd
.split_list <- function(strings, sep = "-") {
  if (is.list(strings) &&
      any(lengths(strings) > 1)) {
      cli::cli_abort("{.arg strings} already holds split sequences.",
                     class = c("tantale_error_bad_argument", "tantale_error"))
  } else if (length(strings) == 1 && is.character(strings)) {
    stopifnot(fs::file_exists(strings))
    seqs <- as.character(Biostrings::readBStringSet(strings), use.names = TRUE)
  } else if (inherits(strings, c("AAStringSet", "BStringSet"))) {
    seqs <- as.character(strings, use.names = TRUE)
  } else if (length(strings) >= 1 && is.list(strings)) {
    seqs <- strings
  } else {
    cli::cli_abort("{.arg strings} must be a file path, an {.cls AAStringSet} or {.cls BStringSet}, or a list of strings.",
                   class = c("tantale_error_bad_argument", "tantale_error"))
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
