



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
                                         varnames = c("RepU1", "RepU2"),
                                         value.name = "Sim")
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
  if (!any("domCode" %in% colnames(tale_parts))) {
    cli::cli_abort("The provided object does not contain a 'domCode' column. Are you using a tale_parts object from distalr()", class = c("tantale_error"))
  }
  if (nrow(diagnose_tale_parts(tale_parts)) != 0L) {
    cli::cli_abort("The provided object does not seem to be sanitized. Have you used a tale_parts object with no empty sequences?", class = c("tantale_error"))
  }
  tale_parts %>% 
    dplyr::select(domCode, rvd) %>%
    dplyr::distinct() %>%
    dplyr::rename(repeatID = domCode,  RVD = rvd) %>%
    dplyr::arrange(repeatID)
}






#' Substitute Distal repeat IDs for RVDs in a TALE alignment matrix.
#'
#'
#' @param repeat_align A multiple TALE repeat sequences alignment in the form of
#'   a matrix as returned by
#'   \code{\link{tales_align}}.
#' @param rvd_map The return value of the
#'   \code{\link[tantale:repeat_to_rvd_map]{repeat_to_rvd_map}} function or the 
#'   \code{\link[tantale:repeat_to_rvd_map_distalr]{repeat_to_rvd_map_distalr}} function
#'   if you used the \code{\link{tales_compare}} function.
#'   
#'
#' @return A TALE alignment matrix made up of RVD sequences.
#' @export
#' @family tales projections
repeat_to_rvd_align <- function(repeat_align , rvd_map) {
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



.repeat_to_sim_align <- function(repeat_align, repeat_sim, ref_tag = NULL) {
  # A function that substitute the repeatIDs with the aa similarity relative to a
  # reference repeat for each column. The ref repeat is the one from a TALE that
  # is defined as a reference in the alignment. This function takes as input, the
  # repeat alignment and the df output by `.format_repeat_dist_mat()` This
  # function outputs the modified alignment matrix
  
  refRowIdx <- match(.pick_ref_name(repeat_align, ref_tag = ref_tag), rownames(repeat_align))
  simAlign <- apply(repeat_align, 2,
                    function(column) {
                      refState <- column[refRowIdx]
                      relevantSims <- subset(repeat_sim, subset = RepU1 == refState)
                      sim <- relevantSims$Sim[match(column, relevantSims$RepU2, nomatch = NA)]
                      if (is.na(refState)) sim[!is.na(column)] <- 0 # if reference repeat is NA, set the aligned repeat sim = 0
                      return(sim)
                    }
  )
  simAlign <- matrix(simAlign, nrow = nrow(repeat_align)) # in case of 1-row matrix
  rownames(simAlign) <- rownames(repeat_align)
  colnames(simAlign) <- colnames(repeat_align)
  return(simAlign)
}



#' Convert repeat alignment to clusterID alignment
#'
#' @param repeat_sim A long, three columns data frame with pairwise similarity scores between repeats as available in the \code{domain_distances} element of the object returned by the \code{\link{tales_compare}} function.
#' @param repeat_align a multiple Tal repeat sequences alignment in the form of a matrix as returned by \code{\link{tales_align}}.
#' @param h_cut a numeric value indicating the height at which to cut the hclust tree of repeats. Interpreted on a distance scale (0 = identical).
#' @return a matrix with exactly the same dimension as the input \code{repeat_sim} but containing clusterID instead of
#' repeatID.
#' @noRd
.repeat_to_cluster_align <- function(repeat_sim, repeat_align, h_cut = 10) {
  # as.dist() expects a DISTANCE. Feeding it the similarity built an inverted
  # dendrogram, so the clusters were wrong (ledger 6). Invert first.
  repeat_dissim <- 100 - as.matrix(reshape2::acast(repeat_sim, RepU1 ~ RepU2, value.var = "Sim"))
  dist_clust <- hclust(as.dist(repeat_dissim))
  dist_cut <- as.data.frame(cbind(RepID = dist_clust$labels, Rep_clust = cutree(dist_clust, h = h_cut)))
  clustIDAlign <- apply(repeat_align, 2,
                        function(column){
                          as.numeric(dist_cut$Rep_clust[match(column, dist_cut$RepID)])
                        })
  clustIDAlign <- matrix(clustIDAlign, nrow = nrow(repeat_align)) # in case of 1-row matrix
  rownames(clustIDAlign) <- rownames(repeat_align)
  colnames(clustIDAlign) <- colnames(repeat_align)
  return(clustIDAlign)
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



#' Recode an RVD alignment as similarity to a reference row
#'
#' Substitutes each RVD with a score expressing how similar its DNA-binding
#' preference is to the RVD of a reference TALE, column by column. This is the
#' RVD-level counterpart of \code{.repeat_to_sim_align()}, which works on
#' protein sequence similarity instead: the two come apart, since repeats can be
#' sequence-divergent yet share an RVD, or near-identical yet differ at
#' positions 12-13.
#'
#' Currently unwired: no \code{fill_type} in either plotting function requests
#' an RVD-level layer. It is the only consumer of the internal
#' \code{rvdSimDf} dataset.
#'
#' @param rvd_align A character matrix of aligned RVDs.
#' @param rvd_sims A data frame of pairwise RVD similarity with columns
#'   \code{rvd1}, \code{rvd2} and \code{Cor}. Defaults to the package's
#'   internal \code{rvdSimDf}.
#' @param ref_tag Pattern selecting the reference row; see
#'   \code{.pick_ref_name()}.
#' @return A numeric matrix with the dimensions and dimnames of
#'   \code{rvd_align}.
#' @keywords internal
.rvd_to_match_align <- function(rvd_align, rvd_sims = rvdSimDf, ref_tag = NULL) {
  refRowIdx <- match(.pick_ref_name(rvd_align, ref_tag = ref_tag), rownames(rvd_align))
  simAlign <- apply(rvd_align, 2,
                    function(column) {
                      refState <- column[refRowIdx]
                      relevantSims <- subset(rvd_sims, subset = rvd1 == refState)
                      relevantSims$Cor[match(column, relevantSims$rvd2, nomatch = NA)]
                    }
  )
  simAlign <- matrix(simAlign, nrow = nrow(rvd_align)) # in case of 1-row matrix
  rownames(simAlign) <- rownames(rvd_align)
  colnames(simAlign) <- colnames(rvd_align)
  return(simAlign)
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
  if (nrow(diagnose_tale_parts(tale_parts)) != 0L) {
    cli::cli_abort("The provided object does not seem to be sanitized. Have you used a tale_parts object with no empty sequences?", class = c("tantale_error"))
  }
  
  if(rvd_only) {
    # tales_anchor_codes() rather than a retyped list: the set has THREE
    # members, and the hardcoded pair here silently kept "XXXXX" -- a
    # terminus detected but not identified -- in a repeats-only string.
    tale_parts %<>% dplyr::filter(!rvd %in% tales_anchor_codes())
  }
  
  rvdStrings <- tale_parts %>%
    dplyr::group_by(arrayID) %>%
    dplyr::arrange(positionInArray) %>%
    dplyr::summarise(
      rvdString = paste(rvd, collapse = "-"),
      posString = paste(positionInArray, collapse = sep)
    )
  rvdStringsSet <- Biostrings::BStringSet(rvdStrings$rvdString)
  names(rvdStringsSet) <- rvdStrings$arrayID
  return(rvdStringsSet)
}




