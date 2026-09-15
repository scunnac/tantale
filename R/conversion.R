



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
      warning("Last element in 'vectorized' sequence is empty. It was removed from output.")
      x[-length(x)]
    } else x
  }
  )
  names(seqsAsVectors) <- names(seqs)
  return(seqsAsVectors)
}




.format_repeat_dist_mat <- function(dist_mat_file) {
  # Distal-1.2 repeat distance matrix is 'almost' symetrical but does not contains the diagonal
  # Top triangle of a symetrical :   Symetrical :
  # 1234                              1234
  #  234                              2234
  #   34                              3334
  #    4                              4444
  # It is like this:
  # 234
  # 34
  # 4
  # So, we need to do a few transformations:
  # Loading file content as a matrix
  distalRepeatDist <- as.matrix(
    read.table(dist_mat_file,
               sep = " ",
               fill = TRUE,
               blank.lines.skip = TRUE,
               header = FALSE,
               check.names = FALSE,
               comment.char = "#"
    )
  )
  # Getting ride of the last columns that appears because the lines in the file have a final space
  distalRepeatDist <- distalRepeatDist[, -ncol(distalRepeatDist)]
  # Adding a first column
  distalRepeatDist <- cbind(V0 = NA, distalRepeatDist)
  # Adding a last line
  distalRepeatDist <- rbind(distalRepeatDist, NA)
  
  # Shifting values to the right in rows with NAs
  distalRepeatDist <- t(apply(distalRepeatDist, 1, function(x) {
    row <- x
    c(row[is.na(row)], row[!is.na(row)])
  }))
  # Filling diagonal with 0 values
  diag(distalRepeatDist) <- 0
  # Filling NA values with diagonal symetric values
  newmat <- distalRepeatDist
  for (i in 1:nrow(distalRepeatDist)) {
    for (j in 1:ncol(distalRepeatDist)) {
      if (is.na(distalRepeatDist[i, j])) {
        distalRepeatDist[i, j] <- distalRepeatDist[j, i]
      } else {
        next()
      }
    }
  }
  # Names
  stopifnot(exprs= all.equal(nrow(distalRepeatDist), ncol(distalRepeatDist)))
  rownames(distalRepeatDist) <- 0:(nrow(distalRepeatDist) - 1)
  colnames(distalRepeatDist) <- 0:(ncol(distalRepeatDist) - 1)
  # Similarity rather than dissimilarity (that is what Alvaro does in his scripts)
  distalRepeatSim <- 100 - distalRepeatDist
  # str(distalRepeatSim)
  # image(t(apply(distalRepeatSim, 2, rev)))
  
  # Output in 'long format'
  distalRepeatSimTable <- reshape2::melt(distalRepeatSim,
                                         as.is = TRUE,
                                         varnames = c("id1", "id2"),
                                         value.name = "sim")
  # str(distalRepeatSimTable)
  return(distalRepeatSimTable)
}





#' Generate a mapping between Distal repeat IDs and their cognate RVD.
#'
#' Uses Distal repeat sequences and RVD sequences from a set of TALEs to return
#' the association between repeat ID and RVD.
#'
#' Care must be taken that TALEs in the two sets of sequences have the same name.
#' In addition, the function tries hard to make sure that the two sets of sequences are identical in every ways but the individual 'values' they contain.
#' It is therefore notably important to make sure that the sequences are consistent in whether they include N-term and C-term domains IDs/Tags or not.
#'
#' @param repeat_vecs Expects a list of Distal repeat IDs character
#'   vectors. Each \strong{named} element corresponding to a TALE.
#' @param repeat_vecs Expects a list of Distal RVDs character
#'   vectors. Each \strong{named} element corresponding to a TALE.
#' @param rvd_vecs A named list of RVD vectors, one per array, parallel to
#'   \code{repeat_vecs}.
#' @return A two columns repeatID - RVD data frame.
#' @export
#' @family tales projections
repeat_to_rvd_map <- function(repeat_vecs, rvd_vecs) {
  # Making sure, these objects are indentical in every ways but the actual values of the vectors
  stopifnot(setequal(names(repeat_vecs), names(rvd_vecs)))
  lrep <- sapply(repeat_vecs, length)
  lrep <- lrep[order(names(lrep))]
  lrvd <- sapply(rvd_vecs, length)
  lrvd <- lrvd[order(names(lrvd))]
  stopifnot(names(lrvd) == names(lrep))
  stopifnot(apply(cbind(lrep, lrvd), 1, function(x) x[1] == x[2]))
  
  # Function to melt the lists
  .l2df <- function(l) {
    dplyr::bind_rows(lapply(l, function(v) data.frame(idx = 1:length(v),
                                                      repeats = v,
                                                      stringsAsFactors = FALSE)
    ),
    .id = "Name")
  }
  talesRepeatDf <- .l2df(repeat_vecs)
  talesRvdDf <- .l2df(rvd_vecs)
  stopifnot(nrow(talesRepeatDf) == nrow(talesRvdDf))
  # Merging to have the repeat ID vs RVDs
  repeat2rvd <- dplyr::full_join(talesRepeatDf, talesRvdDf, by = c("idx" = "idx", "Name" = "Name"))
  colnames(repeat2rvd) <- c("names", "idx", "repeatID", "RVD")
  # Check that for each repeat there is only one corresponding RVD (the converse is NOT true)
  repeat2rvd <- repeat2rvd %>% dplyr::group_by(repeatID, RVD) %>%
    dplyr::summarise(count = dplyr::n())
  check <- repeat2rvd %>%
    dplyr::arrange(repeatID) %>%
    dplyr::group_by(repeatID) %>%
    dplyr::summarise(count = dplyr::n())
  stopifnot(sum(check$count) == length(unique(repeat2rvd$repeatID)))
  # Simplifiy the df
  repeat2rvd <- repeat2rvd %>% dplyr::select(repeatID, RVD) %>% dplyr::ungroup()
  
  return(repeat2rvd)
}


#' Generate a mapping between Distal repeat IDs and their cognate RVD.
#'
#' Uses Distal repeat sequences and RVD sequences from a set of TALEs 
#' analyzed with the \code{\link{tales_compare}} function to return
#' the association between repeat ID and RVD.
#'
#' @param tale_parts The tale_parts object in a \code{\link{tales_compare}} output.
#' @return A two columns repeatID - RVD data frame.
#' @export
#' @family tales projections
repeat_to_rvd_map_distalr <- function(tale_parts) {
  .tales_require(tale_parts, "repeat_to_rvd_map_distalr")
  tale_parts %>% 
    dplyr::select(dom_code, rvd) %>%
    dplyr::distinct() %>%
    dplyr::rename(repeatID = dom_code,  RVD = rvd) %>%
    dplyr::arrange(repeatID)
}
























#' Generates a RVD sequences set from a tale_parts object
#'
#' Uses a tale_parts object in a \code{\link{tales_compare}} output
#' to return a \code{\link[Biostrings]{BStringSet}} of RVD sequences.
#' RVDs are separated by the character specified in the \code{sep} parameter.
#' 
#'
#'
#'
#' Uses Distal repeat sequences and RVD sequences from a set of TALEs 
#' analyzed with the \code{\link{tales_compare}} function to return
#' the association between repeat ID and RVD.
#'
#' @param tale_parts The tale_parts object in a \code{\link{tales_compare}} output.
#' @param sep Used as a RVD separatator
#' @param rvd_only Retrun only RVDs and ommit N- and C- terminal domains 
#' @return A two columns repeatID - RVD data frame.
#' @export
#' @family tales projections
tale_parts_to_rvd <- function(tale_parts, sep = "-", rvd_only = FALSE) {
  
  if(rvd_only) {
    # tales_anchor_codes() rather than a retyped list: the set has THREE
    # members, and the hardcoded pair here silently kept "XXXXX" -- a
    # terminus detected but not identified -- in a repeats-only string.
    tale_parts %<>% dplyr::filter(!rvd %in% tales_anchor_codes())
  }
  
  rvdStrings <- tale_parts %>%
    dplyr::group_by(array_id) %>%
    dplyr::arrange(position_in_array) %>%
    dplyr::summarise(
      rvdString = paste(rvd, collapse = "-"),
      posString = paste(position_in_array, collapse = sep)
    )
  rvdStringsSet <- Biostrings::BStringSet(rvdStrings$rvdString)
  names(rvdStringsSet) <- rvdStrings$array_id
  return(rvdStringsSet)
}




