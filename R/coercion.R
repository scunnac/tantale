


#'
#' Split strings of TALE sequences (`sep`-separated rvd or distal repeat IDs)
#'
#' @description Load the content of fasta file containing TALE sequences (either RVD or Distal repeat code)
#' and return a list of vectors each one composed of the individual elements of the sequence.
#'
#' @param atomicStrings Either, the path to a fasta file, an AAStringSet or "BStringSet"
#' or a list. In all cases, each element of these objects is a string of a
#' tale sequence (`sep`-separated rvd or distal repeat IDs)
#' @param sep Separator of the elements of the sequence
#'
#' @return A list of named vectors representing the 'splited' sequence.
#'
#' @export
toListOfSplitedStr <- function(atomicStrings, sep = "-") {
  if (is.list(atomicStrings) &&
      any(sapply(atomicStrings, length) > 1)) {
      logger::log_error("The value provided for atomicStrings seems already to be splitted.")
      stop()
  } else if (length(atomicStrings) == 1 && is.character(atomicStrings)) {
    stopifnot(fs::file_exists(atomicStrings))
    seqs <- as.character(Biostrings::readBStringSet(atomicStrings), use.names = TRUE)
  } else if (class(atomicStrings) %in% c("AAStringSet", "BStringSet")) {
    seqs <- as.character(atomicStrings, use.names = TRUE)
  } else if (length(atomicStrings) >= 1 && is.list(atomicStrings)) {
    seqs <- atomicStrings
  } else {
    logger::log_error("Something is wrong with the value provided for atomicStrings.")
    stop()
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




formatDistalRepeatDistMat <- function(distalRepeatDistMatFile) {
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
    read.table(distalRepeatDistMatFile,
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
#' @param talesRepeatVectors Expects a list of Distal repeat IDs character
#'   vectors. Each \strong{named} element corresponding to a TALE.
#' @param talesRepeatVectors Expects a list of Distal RVDs character
#'   vectors. Each \strong{named} element corresponding to a TALE.
#' @return A two columns repeatID - RVD data frame.
#' @export
getRepeat2RvdMapping <- function(talesRepeatVectors, talesRvdVectors) {
  # Making sure, these objects are indentical in every ways but the actual values of the vectors
  stopifnot(setequal(names(talesRepeatVectors), names(talesRvdVectors)))
  lrep <- sapply(talesRepeatVectors, length)
  lrep <- lrep[order(names(lrep))]
  lrvd <- sapply(talesRvdVectors, length)
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
  talesRepeatDf <- .l2df(talesRepeatVectors)
  talesRvdDf <- .l2df(talesRvdVectors)
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





convertRvd2RepeatAlign <- function(rvdMsaByGroup, repeatVectors) {
  repSeqs <- lapply(rownames(rvdMsaByGroup), function(r) {
    rvdSeq <- rvdMsaByGroup[r,]
    repSeq <- repeatVectors[[r]]
    
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
  repeatMsaByGroup <- matrix(repeatMsaByGroup, nrow = nrow(rvdMsaByGroup))
  rownames(repeatMsaByGroup) <- rownames(rvdMsaByGroup)
  colnames(repeatMsaByGroup) <- colnames(rvdMsaByGroup)
  return(repeatMsaByGroup)
}


#' Substitute Distal repeat IDs for RVDs in a TALE alignment matrix.
#'
#'
#' @param repeatAlign A multiple TALE repeat sequences alignment in the form of
#'   a matrix as returned by
#'   \code{\link[tantale:buildRepeatMsa]{buildRepeatMsa}} or as one of the
#'   elements of the \code{SeqOfRepsAlignments} or the \code{SeqOfRvdAlignments}
#'   slot in the return object of the \code{\link{buildDisTalGroups}} function.
#' @param repeat2RvdMapping The return value of the
#'   \code{\link[tantale:getRepeat2RvdMapping]{getRepeat2RvdMapping}} function.
#'
#' @return A TALE alignment matrix made up of RVD sequences.
#' @export
convertRepeat2RvdAlign <-  function(repeatAlign , repeat2RvdMapping) {
  states <- unique(as.vector(repeatAlign))
  # TODO: check that all values in states are present in the repeat2RvdMapping df.
  # If not, error
  rvdAlign <- t(
    apply(repeatAlign, 1,
          function(repeatSeq){
            rvdSeq <- repeat2RvdMapping$RVD[match(repeatSeq, repeat2RvdMapping$repeatID)]
          }
    )
  )
  rvdAlign <- matrix(rvdAlign, nrow = nrow(repeatAlign)) # in case of 1-row matrix
  rownames(rvdAlign) <- rownames(repeatAlign)
  colnames(rvdAlign) <- colnames(repeatAlign)
  return(rvdAlign)
}



convertRepeat2SimAlign <-  function(repeatAlign, repeatSim, refTag = NULL) {
  # A function that substitute the repeatIDs with the aa similarity relative to a
  # reference repeat for each column. The ref repeat is the one from a TALE that
  # is defined as a reference in the alignment. This function takes as input, the
  # repeat alignment and the df output by `formatDistalRepeatDistMat()` This
  # function outputs the modified alignment matrix
  
  refRowIdx <- match(pickRefName(repeatAlign, refTag = refTag), rownames(repeatAlign))
  simAlign <- apply(repeatAlign, 2,
                    function(column) {
                      refState <- column[refRowIdx]
                      relevantSims <- subset(repeatSim, subset = RepU1 == refState)
                      sim <- relevantSims$Sim[match(column, relevantSims$RepU2, nomatch = NA)]
                      if (is.na(refState)) sim[!is.na(column)] <- 0 # if reference repeat is NA, set the aligned repeat sim = 0
                      return(sim)
                    }
  )
  simAlign <- matrix(simAlign, nrow = nrow(repeatAlign)) # in case of 1-row matrix
  rownames(simAlign) <- rownames(repeatAlign)
  colnames(simAlign) <- colnames(repeatAlign)
  return(simAlign)
}




convertRvd2MatchAlign <-  function(rvdAlign, rvdSims = tantale::rvdSimDf, refTag = NULL) {
  # A function that substitute the RVDs with a 'RVD match score' relative to a
  # reference rvd for each column. The ref repeat is the one from a TALE that is
  # defined as a reference in the alignment. This function take as input, the
  # repeat alignment and the tantale::rvdSimDf This function output the modified
  # alignment matrix
  refRowIdx <- match(pickRefName(rvdAlign, refTag = refTag), rownames(rvdAlign))
  simAlign <- apply(rvdAlign, 2,
                    function(column) {
                      refState <- column[refRowIdx]
                      relevantSims <- subset(rvdSims, subset = rvd1 == refState)
                      relevantSims$Cor[match(column, relevantSims$rvd2, nomatch = NA)]
                    }
  )
  simAlign <- matrix(simAlign, nrow = nrow(rvdAlign)) # in case of 1-row matrix
  rownames(simAlign) <- rownames(rvdAlign)
  colnames(simAlign) <- colnames(rvdAlign)
  return(simAlign)
}

#' Convert repeat alignment to clusterID alignment
#'
#' @param repeatSim A long, three columns data frame with pairwise similarity scores between repeats as available in the \code{repeat.similarity slot} of the object returned by the \code{\link[tantale:runDistal]{runDistal}} function. \strong{(CORRECT???!!!)}
#' @param repeatAlign a multiple Tal repeat sequences alignment in the form of a matrix as returned by \code{\link[tantale:buildRepeatMsa]{buildRepeatMsa}} or as one of the elements of the \code{SeqOfRepsAlignments} slot in the return object of the \code{\link{buildDisTalGroups}} function.
#' @param h.cut a numeric value indicating the position where to cut the hclust tree of repeats.
#'
#' @return a matrix with exactly the same dimension as the input \code{repeatSim} but containing clusterID instead of repeatID.
convertRepeat2ClusterIDAlign <- function(repeatSim, repeatAlign, h.cut = 10) {
  repeatSim <-  as.matrix(reshape2::acast(repeatSim, RepU1 ~ RepU2, value.var="Sim"))
  dist_clust <- hclust(as.dist(repeatSim))
  dist_cut <- as.data.frame(cbind(RepID = dist_clust$labels, Rep_clust = cutree(dist_clust, h = h.cut)))
  clustIDAlign <- apply(repeatAlign, 2,
                        function(column){
                          as.numeric(dist_cut$Rep_clust[match(column, dist_cut$RepID)])
                        })
  clustIDAlign <- matrix(clustIDAlign, nrow = nrow(repeatAlign)) # in case of 1-row matrix
  rownames(clustIDAlign) <- rownames(repeatAlign)
  colnames(clustIDAlign) <- colnames(repeatAlign)
  return(clustIDAlign)
}

#' Generate a mapping between Distal repeat IDs and their cognate RVD.
#'
#' Uses Distal repeat sequences and RVD sequences from a set of TALEs 
#' analyzed with the \code{\link[tantale:distalr]{distalr}} function to return
#' the association between repeat ID and RVD.
#'
#' @param distalrTaleParts The taleParts object in a \code{\link[tantale:distalr]{distalr}} output.
#' @return A two columns repeatID - RVD data frame.
#' @export
getRepeat2RvdMappingFromDistalr <- function(distalrTaleParts) {
  if (!any("domCode" %in% colnames(distalrTaleParts))) {
    logger::log_error("The provided object does not contain a 'domCode' column. Are you using a taleParts object from distalr()")
    stop()
  }
  if (nrow(diagnoseTaleParts(distalrTaleParts)) != 0L) {
    logger::log_error("The provided object does not seem to be sanitized. Have you used a taleParts object with no empty sequences?")
    stop()
  }
  distalrTaleParts %>% 
    dplyr::select(domCode, rvd) %>%
    dplyr::distinct() %>%
    dplyr::rename(repeatID = domCode,  RVD = rvd) %>%
    dplyr::arrange(repeatID)
}

#' Generates a RVD sequences set from a taleParts object
#'
#' Uses a taleParts object in a \code{\link[tantale:distalr]{distalr}} output
#' to return a \code{\link[Biostrings::BStringSet]{BStringSet}} of RVD sequences.
#' RVDs are separated by the character specified in the \code{sep} parameter.
#' 
#'
#'
#'
#' Uses Distal repeat sequences and RVD sequences from a set of TALEs 
#' analyzed with the \code{\link[tantale:distalr]{distalr}} function to return
#' the association between repeat ID and RVD.
#'
#' @param distalrTaleParts The taleParts object in a \code{\link[tantale:distalr]{distalr}} output.
#' @param sep Used as a RVD separatator
#' @return A two columns repeatID - RVD data frame.
#' @export
taleParts2RvdStringSet <- function(taleParts, sep = "-") {
  if (nrow(diagnoseTaleParts(taleParts)) != 0L) {
    logger::log_error("The provided object does not seem to be sanitized. Have you used a taleParts object with no empty sequences?")
    stop()
  }
  rvdStrings <- taleParts %>%
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




