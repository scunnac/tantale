



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




.format_domain_distances_mat <- function(dist_mat_file) {
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
  distalRepeatSimTable <- as.data.frame(as.table(distalRepeatSim), stringsAsFactors = FALSE)
  colnames(distalRepeatSimTable) <- c("id1", "id2", "sim")
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
#' @param rvd_vecs A named list of RVD vectors, one per array, parallel to
#'   \code{repeat_vecs}.
#' @return A two columns repeatID - RVD data frame.
#' @export
#' @family tales projections
#' @examples
#' repeat_vecs <- list(A1 = c("12", "45", "12"), A2 = c("45", "78"))
#' rvd_vecs <- list(A1 = c("HD", "NI", "HD"), A2 = c("NI", "NG"))
#' repeat_to_rvd_map(repeat_vecs, rvd_vecs)
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
#' analyzed with the \code{\link{tales_compare_distal}} function to return
#' the association between repeat ID and RVD.
#'
#' @param tale_parts The tale_parts object in a \code{\link{tales_compare_distal}} output.
#' @return A two columns repeatID - RVD data frame.
#' @export
#' @family tales projections
#' @examples
#' parts <- data.frame(dom_code = c(1, 2, 1, 3), rvd = c("HD", "NI", "HD", "NG"))
#' repeat_to_rvd_map_distalr(parts)
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
#' Uses a tale_parts object in a \code{\link{tales_compare_distal}} output
#' to return a \code{\link[Biostrings]{BStringSet}} of RVD sequences, one per
#' array, ordered by \code{position_in_array} and joined with \code{sep}.
#'
#' @param tale_parts The tale_parts object in a \code{\link{tales_compare_distal}} output.
#' @param sep Separator joining consecutive RVDs within an array's string.
#' @param rvd_only Return only RVDs and omit N- and C-terminal domains.
#' @return A \code{\link[Biostrings]{BStringSet}}, one element per array,
#'   named by \code{array_id}.
#' @export
#' @family tales projections
#' @examples
#' parts <- data.frame(
#'   array_id = c("A1", "A1", "A1", "A2", "A2"),
#'   position_in_array = c(1, 2, 3, 1, 2),
#'   rvd = c("NTERM", "HD", "CTERM", "NTERM", "NI")
#' )
#' tale_parts_to_rvd(parts)
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
    dplyr::summarise(rvdString = paste(rvd, collapse = sep))
  rvdStringsSet <- Biostrings::BStringSet(rvdStrings$rvdString)
  names(rvdStringsSet) <- rvdStrings$array_id
  return(rvdStringsSet)
}


#### Repeat-code / RVD alignment matrix conversion ####
# Neither direction is called by current production code any more (a
# tales_msa carries both layers at once, so plot() names them directly --
# fill = "dom_code", label = "rvd" -- rather than converting between two
# separate matrices). Kept here, not in inst/legacy/, because both
# directions are still real test dependencies: repeat_to_rvd_align() builds
# fixture data in test_plot_tales_msa.R, and .rvd_to_repeat_align() is
# tested directly, for its own error conditions, in test_error_conditions.R.

#' Substitute Distal repeat IDs for RVDs in a TALE alignment matrix
#'
#' @param repeat_align A multiple TALE repeat sequences alignment in the form
#'   of a matrix, as returned by \code{\link{tales_align}}.
#' @param rvd_map The return value of \code{\link{repeat_to_rvd_map}} or
#'   \code{\link{repeat_to_rvd_map_distalr}} (the latter if the alignment
#'   came from \code{\link{tales_compare_distal}}).
#' @return A TALE alignment matrix made up of RVD sequences.
#' @noRd
repeat_to_rvd_align <- function(repeat_align, rvd_map) {
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

#' View an RVD alignment in terms of repeat codes
#'
#' The inverse direction of \code{repeat_to_rvd_align()}: given an alignment
#' computed on RVDs, substitute each non-gap cell with the corresponding repeat
#' code. The two alphabets differ greatly -- a few dozen RVDs against hundreds
#' of mostly-singleton repeat codes -- so aligning on one and viewing as the
#' other is a genuinely different result from aligning on the other directly.
#'
#' The mapping is **positional**: the k-th non-gap cell of a row is taken to be
#' the k-th element of that row's repeat vector. That is only well defined when
#' the two agree in length, which is now checked.
#'
#' @param rvd_msa_by_group A character matrix of aligned RVDs, rows named by array.
#' @param repeat_vecs A named list of repeat-code vectors, one per row of
#'   \code{rvd_msa_by_group}.
#' @return A character matrix with the dimensions and dimnames of \code{rvd_msa_by_group}.
#' @keywords internal
.rvd_to_repeat_align <- function(rvd_msa_by_group, repeat_vecs) {
  missingRows <- setdiff(rownames(rvd_msa_by_group), names(repeat_vecs))
  if (length(missingRows) > 0L) {
    cli::cli_abort(
      c("Every aligned row needs a matching entry in {.arg repeat_vecs}.",
        "x" = "Missing: {.val {missingRows}}"),
      class = c("tantale_error_rvd_repeat_missing", "tantale_error"))
  }
  nonGap <- rowSums(!is.na(rvd_msa_by_group))
  lens <- lengths(repeat_vecs[rownames(rvd_msa_by_group)])
  bad <- which(nonGap != lens)
  if (length(bad) > 0L) {
    cli::cli_abort(
      c("The back-mapping is positional, so each row must have as many non-gap \\
         cells as it has repeat codes.",
        "x" = "Mismatched row{?s}: {.val {rownames(rvd_msa_by_group)[bad]}}",
        "i" = "non-gap cells {nonGap[bad]} vs {lens[bad]} repeat codes"),
      class = c("tantale_error_rvd_repeat_length", "tantale_error"))
  }
  repSeqs <- lapply(rownames(rvd_msa_by_group), function(r) {
    rvdSeq <- rvd_msa_by_group[r,]
    repSeq <- repeat_vecs[[r]]

    n = 1
    for (i in 1:length(rvdSeq)) {
      if (is.na(rvdSeq[i])) {
        next()
      } else {
        rvdSeq[i] <- repSeq[n]
        n <- n + 1
      }
    }
    rvdSeq <- matrix(rvdSeq, nrow = 1)
    return(rvdSeq)
  })
  repeatMsaByGroup <- do.call(rbind, repSeqs)
  repeatMsaByGroup <- matrix(repeatMsaByGroup, nrow = nrow(rvd_msa_by_group))
  rownames(repeatMsaByGroup) <- rownames(rvd_msa_by_group)
  colnames(repeatMsaByGroup) <- colnames(rvd_msa_by_group)
  return(repeatMsaByGroup)
}


