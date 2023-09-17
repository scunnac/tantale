

##### General utility functions ####

pickRefName <- function(align, refTag = NULL) {
  # How do we select the reference TALE in an alignement?
  #   - the reference could be defined by name or by a string match in the name (eg a strain ID)
  #   - the reference could by default be defined as the longest tal and picked by
  #     ordering their names in case of ties...

  # Find refTag in seq names if provided and output the corresponding unique match
  if (!is.null(refTag)) {
    match <- grepl(refTag, rownames(align))
    if (sum(match) != 1) {
      warning("Cannot identify a single unambiguous sequence to define as a reference using the string in refTag.\n",
              "Using the default method for reference selection.")
      refTag <- NULL
    } else {
      refName <- rownames(align)[match]
    }
  }
  # If no refTag is provided, pick the longest seq(s) and if there are ties, pick the first one alphabetically
  if (is.null(refTag)) {
    strippedAlignLengths <- apply(align, 1, function(seq) length(seq[!is.na(seq)]))
    longest <- rownames(align)[strippedAlignLengths == max(strippedAlignLengths)]
    ifelse(length(longest) == 1, refName <- longest, refName <- sort(longest)[1])
  }
  return(refName)
}

#' Compute a consensus from a TALE msa
#' @description Pick the most frequent element in each column of the alignment matrix.
#'
#' @param align A multiple Tal sequences alignment in the form of a
#'   matrix.
#' @return A vector of consensus elements in each column of \code{align}.
#' 
#' @export
taleAlignConsensus <- function(align) {
  sapply(1:ncol(align), function(x) {
  allElements <- align[,x]
  freq <- sapply(unique(allElements), function(p) S4Vectors::countMatches(p, allElements))
  unique(allElements)[which.max(freq)]
})
}

#' Do elements in a TALE msa match the consensus?
#' @description Compute a logical matrix corresponding to the input \code{align}
#' input with \code{TRUE} if an element match the consensus element at that position
#' or \code{FALSE} otherwise.
#'
#' @param align A multiple Tal sequences alignment in the form of a
#'   matrix.
#' @return A multiple Tal sequences alignment in the form of a
#'   matrix filled with logical values if \code{} is \code{FALSE} and
#'   a long tibble representing the original alignment otherwise (default).
#' 
#' @export
matchConsensus <- function(align, returnLong = TRUE) {
  consensus <- taleAlignConsensus(align)
  align <- align
  for (k in 1:ncol(align)){
    rept <- consensus[k]
    if (is.na(rept)) {
      align[,k] <- FALSE
    } else {
      align[,k] <- ifelse(toupper(align[,k]) == toupper(rept), TRUE, FALSE)
    }
  }
  if (!returnLong) return(align)
  matchConsensusLong <- align %>% reshape2::melt() %>%
    dplyr::as_tibble()
  colnames(matchConsensusLong) <- c("arrayID", "positionInArray", "matchConsensus")
  return(matchConsensusLong)
}


##### Tale domains sequences multiple alignment ####


#' Perform multiple alignment of TALE repeat or RVD sequences
#'
#' @description Perform multiple alignment of TALE repeat or RVD sequences using the \code{--text} mode of the \href{https://mafft.cbrc.jp/alignment/software/}{MAFFT} Multiple alignment program.
#'
#' By default, uses the simple scoring matrix defined in \href{https://mafft.cbrc.jp/alignment/software/textcomparison.html}{the text mode of MAFFT}. Users can optionally provide a custom scoring matrix.
#'
#' @param inputSeqs Any object accepted as input by the
#'  \code{\link[tantale:toListOfSplitedStr]{toListOfSplitedStr}} function, such as the path to a fasta file containing the TALE sequences to be aligned or the \code{coded.repeats.str} slot of the object returned by the \code{\link[tantale:runDistal]{runDistal}} function. Can also be the return value of the \code{\link[tantale:taleParts2RvdStringSet]{taleParts2RvdStringSet}} function if one wants to align RVD sequences.
#'
#' @param sep Passed to \code{toListOfSplitedStr()} to split the TALEs strings in input.
#' @param distalRepeatSims A long, three columns data frame with pairwise similarity scores between repeats as available in the \code{repeat.similarity slot} of the object returned by the \code{\link[tantale:runDistal]{runDistal}} function.
#' @param mafftOpts A character string containing additional options for the MAFFT command. This is notably useful to tweak the Gap opening and gap extension penalties.
#' @param mafftPath Path to a MAFFT installation directory. By default uses the MAFFT version included in tantale.
#' @param gapSymbol Specify a alternative symbol for gaps in the alignments.
#'
#' @return A character matrix representing the multiple alignment.
#' @export
buildRepeatMsa <- function(inputSeqs, sep = " ", distalRepeatSims = NULL,
                           mafftOpts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
                           mafftPath = system.file("tools", "mafft-linux64",package = "tantale", mustWork = TRUE),
                           gapSymbol = NA) {
  # A bunch of tempfiles
  simMatHexFile <- tempfile(pattern = "simMatHexFile")
  simMatAsciiFile <- tempfile(pattern = "simMatAsciiFile")
  hexFile <- tempfile(pattern = "hexFile")
  asciFile <- tempfile(pattern = "asciFile")
  mafftAsciiOutFile <- tempfile(pattern = "mafftAsciiOutFile")
  mafftHexOutFile <- tempfile(pattern = "mafftHexOutFile")

  # This is the data frame that will enable conversion of RVDs to Hexadecimal codes
  asciitable = data.frame(hex = as.raw(1:255),
                          printable =rawToChar(as.raw(1:255),multiple=TRUE),
                          stringsAsFactors = FALSE)
  mafftExcludedHex <- as.raw(c(0x0, 0x3E, 0x3D, 0x3C, 0x2D, 0x20, 0x0d, 0x0a))
  asciitableForMafft <- asciitable[! asciitable$hex %in% mafftExcludedHex, ]
  #Encoding(asciitableForMafft$printable) <- "bytes"
  
  # Load repeat/RVD sequences
  seqsAsVectors <- suppressWarnings(toListOfSplitedStr(inputSeqs, sep = sep))
  residues <- unique(unlist(seqsAsVectors))
  
  # Deals with cases where the nomber of sequences is < 2
  if (length(seqsAsVectors) == 0L) {
    logger::log_warn("The provided object in inputSeqs is empty. Returning an empty matrix")
    warning()
    return(matrix())
  }
  if (length(seqsAsVectors) == 1L) {
    logger::log_info("The provided object in inputSeqs has only one sequence. Returning it as a matrix.")
    msaOfResiduesAsMatrix <- as.matrix(as.data.frame(seqsAsVectors))
    msaOfResiduesAsMatrix <- matrix(msaOfResiduesAsMatrix, nrow = 1)
    rownames(msaOfResiduesAsMatrix) <- colnames(as.data.frame(seqsAsVectors))
    colnames(msaOfResiduesAsMatrix) <- 1:length(msaOfResiduesAsMatrix)
    return(msaOfResiduesAsMatrix)
  }
  
  # Determine the type of 'elements' (rvd or repeat) contained in the sequences
  frequentRvds <- c("NN", "NG", "HD", "NI", "N*", "NS")
  if(! any(residues %in% frequentRvds)) {
    logger::log_info("Will be assuming sequences contain repeat unit codes because ",
                     "none of the RVDs obtained from input sequences matches ",
                     "a list of 'frequent RVDs': {paste(frequentRvds, collapse = ' ')}")
    repeatType <- "repeatUnit"
  } else {
    logger::log_info("Input sequences are detected as RVD sequences.")
    repeatType <- "rvds"
  }
  if( length(residues) > nrow(asciitableForMafft) ) {
    logger::log_error("Number of unique resisues (RVDs or repeat units) must be =< 248.")
    logger::log_error("Currently, your set of sequences contains {length(residues)} unique residues...")
    stop()
  }
  

  
  # Coding residues in hexadecimal representations and concatenating them for mafft --text
  seqsOfHex <- sapply(seqsAsVectors, function(x) {
    idxs <- match(x, residues)
    paste(asciitableForMafft$hex[idxs], collapse = " ")
    }
  )
  seqsOfHex <- Biostrings::BStringSet(seqsOfHex)
  names(seqsOfHex) <- names(seqsAsVectors)
  # Write to a temp file
  Biostrings::writeXStringSet(seqsOfHex, filepath = hexFile)


  # If provided recode also the distance matrix
  if(is.null(distalRepeatSims) || repeatType == "rvds") {
    maffMatOpt <- ""
  } else if (!is.null(distalRepeatSims)) {
    logger::log_info("The provided similarity matrix file will be used to compute msa.")
    if (length(distalRepeatSims) > 1 &&
        (is.data.frame(distalRepeatSims) | tibble::is_tibble(distalRepeatSims))
    ) {
      repeatSims <- distalRepeatSims[,c("RepU1", "RepU2", "Sim")]
    } else if (length(distalRepeatSims) == 1 && is.character(distalRepeatSims)) {
      repeatSims <- formatDistalRepeatDistMat(distalRepeatSims)
    } else {
      logger::log_error("Somthing is wrong with the value provided for distalRepeatSims. It must be either")
      logger::log_error("the path to a '*_Repeatmatrix.mat' file produced by Distal or")
      logger::log_error("table like object with three columns, usually produced by the")
      logger::log_error("formatDistalRepeatDistMat() function")
      stop()
    }

    stopifnot(all(residues %in% unique(repeatSims$RepU1)))
    repeatSims <- subset(repeatSims, RepU1 %in% residues & RepU2 %in% residues)
    stopifnot(all.equal(nrow(repeatSims), length(residues)^2))
    repeatSims$RepU1 <- asciitableForMafft$hex[match(repeatSims$RepU1, residues)]
    repeatSims$RepU2 <- asciitableForMafft$hex[match(repeatSims$RepU2, residues)]
    colnames(repeatSims) <-  NULL
    write.table(repeatSims, file = simMatHexFile, row.names = FALSE, fileEncoding = "ASCII")
    maffMatOpt <- glue::glue("--textmatrix {simMatHexFile}")
  }

  # Running mafft msa
  logger::log_info("Now running MAFFT (Copyright 2002-2007 Kazutaka Katoh) on TALE array sequences.")
  asciiConverstionCmd <- glue::glue("{mafftPath}/mafftdir/libexec/hex2maffttext {hexFile} > {asciFile}")
  mafftCmd <-  glue::glue("{mafftPath}/mafft.bat {maffMatOpt} --text {mafftOpts} {asciFile} > {mafftAsciiOutFile}")
  MsaConversionToHexCmd <- glue::glue("{mafftPath}/mafftdir/libexec/maffttext2hex {mafftAsciiOutFile} > {mafftHexOutFile}")
  res <- system(command = paste(asciiConverstionCmd, mafftCmd, MsaConversionToHexCmd, sep = "; "),
         ignore.stdout = FALSE, ignore.stderr = FALSE, intern = FALSE)

  # Getting msa output and converting back to alignment of residues
  msaOfHex <- Biostrings::readBStringSet(mafftHexOutFile)
  if (length(msaOfHex) == 0L) {
    logger::log_error("MAFFT failled to complete sucessfully...")
    res %>% logger::skip_formatter() %>% logger::log_error()
    stop()
  } else {
    logger::log_debug("MAFFT completed sucessfully!! Yeah!")
  }
  #cat(as.character(msaOfHex), sep = "\n")
  msaAsHexVectors <- stringr::str_split(as.character(msaOfHex), pattern = " ")
  msaAsHexVectors <- lapply(msaAsHexVectors, function(x) x[-length(x)]) # Remove last "" element
  #cat(knitr::kable(t(matrix(msaAsHexVectors))), sep = " ")
  msaOfResiduesAsMatrix <- t(
    sapply(msaAsHexVectors, function(x) {
      idxs <- match(x, asciitableForMafft$hex, nomatch = NA)
      #cat("idx in asciiTable: ", idxs, "\n")
      seqOfResidues <- residues[idxs]
      seqOfResidues[is.na(seqOfResidues)] <- gapSymbol
      #cat("Seq of residues: ", seqOfResidues, "\n")
      seqOfResidues
    }
    )
  )
  rownames(msaOfResiduesAsMatrix) <- names(msaOfHex)
  colnames(msaOfResiduesAsMatrix) <- 1:ncol(msaOfResiduesAsMatrix)
  return(msaOfResiduesAsMatrix)
}





##### Tale domains msa plotting ####
#' Plotting a multiple alignment of TALE sequences
#' @description Plot in frame of \code{\link[gplots:heatmap.2]{heatmap.2}} for Tals alignment.
#' 
#' @details
#'  "repeat.similarity" plot shows RVD alignment of Tals, a hierarchical dendrogram
#'  reflecting overall similarities between TALEs, similarity between repeats alignment by the color
#'    of cells, and (if rvdSim is provided) similarity between RVDs alignment in the 
#'    color of the RVD labels.
#'    
#'  "repeat.clusters" plot shows repeat alignment of Tals with cells filled with colors representing
#'  the repeat clustering group and a hierarchical dendrogram reflecting overall similarities
#'  between TALEs.
#'   
#'  "repeat.clusters.with.rvd" plots repeat alignment of Tals with rvd labeled.
#'  If plotting from the outputs of \code{buildDistalGroups},
#'  you supply repeatClustID/Similarity alignment to \code{forMatrix}, with
#'  \code{talsim}, \code{forCellNote} - repeat/rvd alignment, and refgrep optionally.
#'  
#'  But if you don't have these alignments, you can provide \strong{repeat alignment}
#'  to \code{forMatrix} with \code{repeatSim}, the repeat similarity data frame,
#'  the function will convert it into repeatClustID/Similarity alignment depending
#'  on the plot type. In case of \emph{repeat.clusters}, you may want to adjust
#'  the param \code{repeat.clust.h.cut} to decrease/increase the number of repeat clusters.
#' 
#' 
#' 
#' 
#' @param talsim a \emph{three columns Tals similarity table} as obtained
#'  with \code{\link[tantale:runDistal]{runDistal}} in the 'tal.similarity' slot of the returned object.
#' @param repeatAlign a multiple Tal repeat sequences alignment in the
#'  form of a matrix as returned by \code{\link[tantale:buildRepeatMsa]{buildRepeatMsa}}
#'  or as one of the elements of the \code{SeqOfRepsAlignments} slot in the return object
#'  of the \code{\link{buildDisTalGroups}} function.
#' @param repeatSim A long, three columns data frame with pairwise similarity 
#' scores between repeats as available in the \code{repeat.similarity slot}
#' of the object returned by the \code{\link[tantale:runDistal]{runDistal}} function.
#' \strong{(CORRECT???!!!)}
#' @param plot.type Either \code{"repeat.similarity"}, \code{"repeat.clusters"} ,
#'  \code{"repeat.clusters.with.rvd"}. Defines the type of plot that will be produced
#'   by the function. See below for details.
#' @param repeat.clust.h.cut height for tree cutting when plot type
#'  in "repeat.clusters".
#' @param rvdAlign (optional) when the rvds need to be labeled in the
#'  plot (plot.type = "repeat.similarity" or "repeat.clusters.with.rvd",
#'  a multiple Tal repeat sequences alignment in the form of a matrix as
#'  returned by \code{\link[tantale:buildRepeatMsa]{buildRepeatMsa}} or as one
#'  of the elements of the \code{SeqOfRepsAlignments} slot in the return object
#'  of the \code{\link{buildDisTalGroups}} function. 
#' @param refgrep regular expression pattern that will be used to search Tal names
#' to select the reference in the alignment.
#' @param consensusSeq (logical) whether to display the consensus sequence when 
#' the plot type is "repeat.clusters.with.rvd".
#' @param noteColSet In case rvdSim = NULL, vector of 2 colors for rvd alignment,
#' the first color is for matched rvds, and the second color is for mismatched ones.
#' In the other case, more colors should be supplied.
#' @param save.path file path to save the plot. If save.path is NULL, the heatmap 
#' will be printed. If save.path is specified, the image file will be created with
#' the format based on file extension.
#' @param ... any other arguments of \code{\link[gplots:heatmap.2]{heatmap.2}}
#' 
#' @return the return value of \code{\link[gplots:heatmap.2]{heatmap.2}}
#' 
#' @export
heatmap_msa <- function(talsim, repeatAlign, rvdAlign = NULL, repeatSim, repeat.clust.h.cut = 90, refgrep = NULL, consensusSeq = FALSE, noteColSet = NULL, plot.type, save.path, ...) {
  
  
  if (startsWith(plot.type, "repeat.clusters")) {
    forMatrix <- convertRepeat2ClusterIDAlign(repeatAlign = repeatAlign, repeatSim = repeatSim, h.cut = repeat.clust.h.cut)
  } else if (plot.type == "repeat.similarity") {
    forMatrix <- convertRepeat2SimAlign(repeatAlign = repeatAlign, repeatSim = repeatSim, refTag = refgrep)
  } else if (!hasArg(plot.type) || is.null(plot.type) || is.na(plot.type)) {
    stop("Missing plot.type")
  } else{
    stop(glue::glue("\"{plot.type}\" plot is not available."))
  }
  
  
  if (plot.type == "repeat.clusters") {
    forCellNote <- repeatAlign
  } else {
    if (is.null(rvdAlign)) stop("\"{plot.type}\" plot requires rvdAlign!")
    forCellNote <- rvdAlign
  }
  
  forCellNote <- forCellNote[rownames(forMatrix),]
  
  
  
  # for 'Rowv' = dend
  
  ###!!! THIS FAILS IF repeatAlign has a single sequence
  
  talsim <- talsim[talsim$TAL1 %in% rownames(forCellNote), ]
  talsim <- talsim[talsim$TAL2 %in% rownames(forCellNote), ]
  talsim <- as.matrix(reshape2::acast(talsim, TAL1 ~ TAL2, value.var = "Sim")) # melt then unmelt ...
  taldist <- 100 -talsim
  taldist <- taldist[rownames(forCellNote), ]
  taldist <- taldist[, rownames(forCellNote)]
  
  clust <- hclust(as.dist(taldist))
  forRowv <- as.dendrogram(clust)
  
  
  
  if (startsWith(plot.type, "repeat.clusters")) { # in case of repcode plotting
    # define 'col'
    # forCol <- function(x) scales::hue_pal(l = 55)(n=100)[0:x]
    forCol <- function(x) viridis::inferno(n=100, end = .9)[0:x]
    forBreaks <- 0:max(forMatrix, na.rm = T)
    
    # define 'key'
    forKeyxlab <- NA
    forKey <- FALSE
    forNoteCex <- 1
    
    # define 'notecol'
    if (is.null(noteColSet)) {
      noteColSet <- list(matched = "white", mismatched = "#01FFFF")
    } else {
      noteColSet <- list(matched = noteColSet[1], mismatched = noteColSet[2])
    }
    if (endsWith(plot.type, "with.rvd")) {
      # consensus rvds
      # rvdsAlignedStrings <- apply(forCellNote, 1, function(x) paste(x, collapse = "-")) %>% BStringSet()
      # consensusRVD <- Biostrings::consensusString(rvdsAlignedStrings, ambiguityMap = "+")
      
      consensusRVD <- sapply(1:ncol(forCellNote), function(x) {
        allRVDs <- forCellNote[,x]
        freq <- sapply(unique(allRVDs), function(p) countMatches(p, allRVDs))
        unique(allRVDs)[which.max(freq)]
      })
      
      # notecol
      rvdcol <- forCellNote
      # if (is.null(refgrep)) {
      #   reftalID <- which.max(apply(forCellNote, 1, function(x) length(x[!is.na(x)])))
      # } else {
      #   reftalID <- which(grepl(refgrep, rownames(forCellNote)))
      # }
      # rownames(forMatrix)[reftalID] <- paste0(rownames(forMatrix)[reftalID], "_")
      # for (k in 1:ncol(rvdcol)){
      #   rvd <- as.character(forCellNote[reftalID, k])
      #   if (is.na(rvd)) {
      #     rvdcol[,k] <- noteColSet$mismatched
      #   } else {
      #     rvdcol[,k] <- ifelse(toupper(rvdcol[,k]) == toupper(rvd), noteColSet$matched, noteColSet$mismatched)
      #   }
      # }
      for (k in 1:ncol(rvdcol)){
        rvd <- consensusRVD[k]
        if (is.na(rvd)) {
          rvdcol[,k] <- noteColSet$mismatched
        } else {
          rvdcol[,k] <- ifelse(toupper(rvdcol[,k]) == toupper(rvd), noteColSet$matched, noteColSet$mismatched)
        }
      }
      rvdcol <- rvdcol[order.dendrogram(forRowv), ]
      forNoteCol <- t(as.matrix(rvdcol))
      
      
      
    } else if (endsWith(plot.type, "repeat.clusters")) {
      forNoteCol <- "white"
    } else {
      stop(glue::glue("\"{plot.type}\" plot is not available."))
    }
  }
  else if (plot.type == "repeat.similarity") { # in case of rvd plotting
    # define 'col'
    # forCol <- colorRampPalette(c("dodgerblue4", "dodgerblue3", "dodgerblue", "deepskyblue", "white"))
    forCol <- colorRampPalette(c("#421727", "#6e2742", "#9a365c", "#b03e69", "#ffffff"))
    forBreaks <- 0:100
    
    rvdSim <- NULL
    # define 'notecol'
    if (!is.null(rvdSim) && is.matrix(rvdSim) && is.numeric(rvdSim)) { # in case rvd similarity matrix is provided
      if (is.null(noteColSet)) {
        noteColSet <- viridis::viridis(n=200, direction = -1)
      }
      col_range <- function(x) noteColSet[as.integer(x)]
      rvdSim <- rvdSim[rownames(forCellNote), ]
      rvdcol <- rvdSim
      for (i in 1:nrow(rvdSim)) {
        for (j in 1:ncol(rvdSim)) {
          if (is.na(rvdcol[i,j])) {
            rvdcol[i,j] <- "grey"
          } else {
            rvdcol[i,j] <- col_range(rvdSim[i,j] * 100 + 100)
          }
        }
      }
      
    } else { # default
      if (is.null(noteColSet)) {
        noteColSet <- list(matched = "black", mismatched = "red")
      } else {
        noteColSet <- list(matched = noteColSet[1], mismatched = noteColSet[2])
      }
      
      if (is.null(refgrep)) {
        sim_len <- apply(forMatrix, 1, function(x) {length(grep("100", x))})
        reftal <- names(which.max(sim_len))
        reftalID <- which(rownames(forCellNote) == reftal)
      } else {
        reftalID <- which(grepl(refgrep, rownames(forCellNote)))
      }
      rownames(forMatrix)[reftalID] <- paste0(rownames(forMatrix)[reftalID], "_#")
      rvdcol <- forCellNote
      for (k in 1:ncol(rvdcol)){
        rvd <- as.character(forCellNote[reftalID, k])
        if (is.na(rvd)) {
          rvdcol[,k] <- noteColSet$mismatched
        } else {
          rvdcol[,k] <- ifelse(toupper(rvdcol[,k]) == toupper(rvd), noteColSet$matched, noteColSet$mismatched)
        }
      }
    }
    rvdcol <- rvdcol[order.dendrogram(forRowv), ]
    forNoteCol <- t(as.matrix(rvdcol))
    
    # define 'key'
    forKeyxlab <- "AA similarity"
    forKey <-  TRUE
    forNoteCex <- 1.2
  }
  
  
  # adjust size, position, ... of plot's elements
  wid_left <- 1
  wid_right <- 0.125 * ncol(forMatrix)
  if (max(nchar(forCellNote), na.rm = T) >= 5) wid_right <- wid_right * 2
  hei_top <- ifelse(isFALSE(forKey), .5, .75)
  hei_bottom <- 0.125 * (nrow(forMatrix) + 1.5)
  
  # add "notecol" legend
  extra.key <- function(x = NULL, check.plot.type = plot.type, hei = hei_top * 2.54, wid = wid_left * 2.54, colSet = noteColSet) {
    if (check.plot.type == "repeat.similarity" || endsWith(check.plot.type, "with.rvd")) {
      if (!is.null(x) && is.matrix(x) && is.numeric(x)) {
        par(mai = c(hei*.4, 0, hei*.2, wid*.1), mgp = c(2, 1, 0))
        image(z = matrix(seq(-1, 1, by = .01), ncol = 1), col = colSet, yaxt = "n", xaxt = "n", xlab = "RVD similarity")
        axis(1, at = seq(0, 1, by = .25), labels = c("-1", NA, "0", NA, "1"))
      } else {
        par(mai = c(hei*.2, wid*.25, hei*.2,  wid*.25), mgp = c(2, 1, 0))
        image(z = matrix(c(0, 1), ncol = 2), col = "grey50", yaxt = "n", xaxt = "n")
        abline(h = 0.5, col = "grey", lwd = 1.5)
        text(0, 1, labels = ifelse(check.plot.type == "repeat.similarity", "reference", "consensus"), col = colSet$matched, font = 2)
        text(0, 0, labels = "other", col = colSet$mismatched, font = 2)
        mtext(side = 1, at = 0, text = "RVD alignment", cex = .75, col = "black", padj = 0.5)
        if (check.plot.type == "repeat.clusters.with.rvd" && consensusSeq) {
          par(mar = c(0,0,0,0))
          image(z = matrix(1:length(consensusRVD), ncol = 1), col = "grey50", bg = "grey", yaxt = "n", xaxt = "n")
          for (i in 1:length(consensusRVD)) {
            abline(v = (i-.5)/(length(consensusRVD)-1), col = "grey")
            text((i-1)/(length(consensusRVD)-1), 0, labels = consensusRVD[i], font = 2, col = colSet$matched, cex = 1.2)
          }
          par(xpd = NA)
          lab <- axis(side = 4, at = 0, labels = "", las = 2, line = -.5, tick = 0)
          text(x = par("usr")[2] + 1.5 * strwidth("M"), adj = c(0,NA),
               y = lab, labels = "#Consensus", col = "grey50", cex = 1.2, font = 2)
        }
        
      }
    } else {
      return(NULL)
    }
    
  }
  
  # save plot to "save.path" if provided, the format of image depends on extension of "save.path"
  if (hasArg(save.path)) {
    img_format <- gsub(".*\\.", "", basename(save.path))
    img_size <-  list(save.path, width = (wid_left + wid_right + wid_left) * 2, height = (hei_top + hei_bottom) * 2)
    if (img_format %in% c("bmp", "jpeg", "png", "tiff")) {
      img_size <- c(img_size, units = "in", res = 1440)
    }
    do.call(img_format, img_size)
  }
  
  default_arg_list <- list(x = forMatrix, #
                           
                           # dendrogram control
                           Rowv = forRowv, ##
                           Colv = FALSE,
                           dendrogram = "row",
                           
                           # colors
                           col = forCol, #
                           breaks = forBreaks, ##
                           
                           # cell labeling
                           cellnote  = forCellNote, #
                           notecex = forNoteCex,
                           notecol = forNoteCol, ##
                           na.color = "grey",
                           extrafun = extra.key, ## key for notecol
                           
                           # Row/Column Labeling
                           margins = c(3, 0),
                           cexRow = 1,
                           cexCol = 1,
                           labCol = c(1:ncol(forMatrix)), ##
                           adjCol = c(NA, .5),
                           offsetCol = 0,
                           offsetRow = 0,
                           srtCol = 0,
                           
                           # block sepration
                           colsep = 0:ncol(forMatrix), ##
                           rowsep = 0:nrow(forMatrix), ##
                           sepcolor = "#bbbbbb",
                           sepwidth = c(0.005,0.005),
                           
                           # color key + density info
                           density.info = "none",
                           key = forKey, ##
                           key.title = NA,
                           key.par = list(mai = c(hei_top*.6*2.54, wid_left*.1*2.54, hei_top*.1*2.54, 0),
                                          mgp = c(2, 1, 0)
                           ),
                           key.xlab = forKeyxlab, ##
                           
                           # plot labels
                           xlab = "Domain",
                           # ylab = "TAL ID",
                           
                           # plot layout
                           lmat = rbind(c(4, 3, 5), c(4,6,7), c(2, 1, 0)),
                           lhei = c(hei_top,
                                    ifelse(consensusSeq, hei_bottom/nrow(forMatrix), 0), hei_bottom), ##
                           lwid = c(wid_left, wid_right, wid_left), ##
                           
                           # trace
                           trace = "none")
  
  custom_arg_list <- as.list(substitute(list(...)))[-1L]
  
  default_arg_list <- default_arg_list[!names(default_arg_list) %in% names(custom_arg_list)]
  
  heatmap_plot <- do.call(gplots::heatmap.2, c(default_arg_list, custom_arg_list))
  
  if (hasArg(save.path)) {
    dev.off()
  }
  return(invisible(heatmap_plot))
}



#' 'Nice' plotting a multiple alignment of TALE sequences
#' @description Plot TALEs msa in the ggplot2 framework.
#'
#' @details This function as a similar purpose as
#' \code{\link[tantale:heatmap_msa]{heatmap_msa}} but has been implemented with
#' \code{\link[ggplot2:ggplot]{ggplot}}. It is more versatile (takes single row
#' matrices of alignment) but a bit slower.
#'
#' The type of plot that you will get will depend on the provided information in
#' the form of parameter values See the tantale website for detailed usage
#' cases.
#'
#' The only mandatory argument is either \code{repeatAlign} \strong{or}
#' \code{rvdAlign}.
#'
#' The plot is printed and returned for further modifications is necessary.
#'
#'
#' @param talsim a \emph{three columns Tals similarity table} as obtained with
#'   \code{\link[tantale:runDistal]{runDistal}} in the 'tal.similarity' slot of
#'   the returned object.
#' @param repeatAlign A multiple Tal repeat sequences alignment in the form of a
#'   matrix as returned by \code{\link[tantale:buildRepeatMsa]{buildRepeatMsa}}
#'   or as one of the elements of the \code{SeqOfRepsAlignments} slot in the
#'   return object of the \code{\link{buildDisTalGroups}} function.
#' @param repeatSim A long, three columns data frame with pairwise similarity
#'   scores between repeats as available in the \code{repeat.similarity slot} of
#'   the object returned by the \code{\link[tantale:runDistal]{runDistal}} or the 
#'   the \code{\link[tantale:distalr]{distalr}} functions.
#' @param repeat.clust.h.cut height for tree cutting when defining domain/repeat
#'   clusters.
#' @param rvdAlign A multiple Tal RVD sequences alignment in the form of a
#'   matrix as returned by \code{\link[tantale:convertRepeat2RvdAlign]{convertRepeat2RvdAlign}}, 
#'    \code{\link[tantale:buildRepeatMsa]{buildRepeatMsa}} 
#'   or as one of the elements of the \code{SeqOfRepsAlignments} slot in the
#'   return object of the \code{\link{buildDisTalGroups}} function.
#' @param refgrep Regular expression pattern that will be used to search TALE
#'   names to select the reference in the alignment.
#' @param consensusSeq (logical) Whether to display the consensus sequence
#'  **NOT IMPLEMENTED YET**
#' @param fillType Either "repeatClust" or "repeatSim". If both options are
#'   possible because the necessary information is there (at least a
#'   \code{repeatSim} value), this argument will decide what type of 'box color
#'   filling' is employed and it is either based on the cluster where the repeat
#'   falls after clustering all the repeat in the alignment or it is based on
#'   the amino acid similarity between a repeat at a position and the repeat of
#'   the 'reference' TALE at this position.
#' @return An \code{\link[aplot:insert_left]{aplot}} object.
#' 
#' @export
ggplotTalesMsa <- function(repeatAlign,
                           talsim = NULL,
                           rvdAlign = NULL,
                           repeatSim = NULL,
                           repeat.clust.h.cut = 90,
                           refgrep = NULL,
                           consensusSeq = FALSE,
                           fillType = "repeatClust" #"repeatSim"
) {
  
  # Arguments checking
  if (is.null(rvdAlign) & is.null(repeatAlign)) {
    logger::log_error("You must provide at least either a value for `repeatAlign` or for `rvdAlign`")
    stop()
  }
  if (!is.null(rvdAlign)) {
    countOfTales <- nrow(rvdAlign)
    arrayNames <- rownames(rvdAlign)
  }
  if (!is.null(repeatAlign)) {
    countOfTales <-  nrow(repeatAlign)
    arrayNames <- rownames(repeatAlign)
  }
  if (!is.null(repeatAlign) & is.null(nrow(repeatAlign))) {
    logger::log_error("Check the provided input repeatAlign matrix.",
                      "It may conain a single sequence that was coerced to vector rather than remaining a matrix...",
                      .sep = " ")
    stop()
  }
  if (!is.null(rvdAlign) & is.null(nrow(rvdAlign))) {
    logger::log_error("Check the provided input rvdAlign matrix.",
                      "It may conain a single sequence that was coerced to vector rather than remaining a matrix...",
                      .sep = " ")
    stop()
  }
  if (countOfTales < 1) {
    logger::log_error("The provided input repeatAlign matrix has less than one sequence. Cannot proceed...")
    stop()
  }
  
  
  # Getting repeat align
  if (!is.null(repeatAlign)) {
    repeatAlignLong <- repeatAlign %>% reshape2::melt() %>%
    dplyr::as_tibble()
  colnames(repeatAlignLong) <- c("arrayID", "positionInArray", "domCode")
  repeatAlignLong %<>% dplyr::mutate(arrayID = as.character(arrayID))
  #### TODO: pad domCode for it to be 4 characters long ####
  repeatMatchConsensusLong <- matchConsensus(repeatAlign)
  colnames(repeatMatchConsensusLong) <- c("arrayID", "positionInArray", "matchConsensusRepeat")
  repeatAlignLong %<>% dplyr::left_join(repeatMatchConsensusLong,
                                        by = dplyr::join_by(arrayID, positionInArray))
  }
  
  
  # Getting rvd align if available
  if (!is.null(rvdAlign)) {
    rvdAlignLong <- rvdAlign %>% reshape2::melt() %>%
      dplyr::as_tibble()
    colnames(rvdAlignLong) <- c("arrayID", "positionInArray", "rvd")
    rvdAlignLong %<>% dplyr::mutate(rvd = gsub("NTERM", "N-", rvd),
                                       rvd = gsub("CTERM", "-C", rvd)
    )
    
    # Tale rvd text color if possible
    # consensus RVD sequence
    # Coloring of RVDs in alignment depending on whether they match the consensus at
    # the position
    consensusRVD <- taleAlignConsensus(rvdAlign)
    rvdConsensusSeqLong <- tibble::tibble(arrayID = "Consensus",
                                          positionInArray = seq_along(consensusRVD),
                                          rvd = consensusRVD,
                                          matchConsensusRvd = TRUE,
                                          domCode = NA,
                                          repeatClusterId = NA,
                                          repeatSimVsRef = NA
    )
    rvdMatchConsensusLong <- matchConsensus(rvdAlign)
    colnames(rvdMatchConsensusLong) <- c("arrayID", "positionInArray", "matchConsensusRvd")
    # Join with rvd tible
    rvdAlignLong %<>% dplyr::left_join(rvdMatchConsensusLong,
                                          by = dplyr::join_by(arrayID, positionInArray))
  }
  
  # Assign main alignment object in long format
  if (!is.null(repeatAlign) & !is.null(rvdAlign)) {
    repeatAlignLong %<>% dplyr::inner_join(rvdAlignLong,
                                          by = dplyr::join_by(arrayID, positionInArray),
                                          unmatched = "error",
                                          relationship = "one-to-one")
  } else if (!is.null(repeatAlign) & is.null(rvdAlign)) {
    repeatAlignLong <- repeatAlignLong
  } else if (is.null(repeatAlign)) {
    repeatAlignLong <- rvdAlignLong
  } else {
    stop("something wrong with parameters values")
  }

  # joining repeat cluster if possible
  # joining repeat similarity relative to ref
  if (!is.null(repeatSim) & !is.null(repeatAlign)) {
    repeatClusterAlignLong <- convertRepeat2ClusterIDAlign(repeatAlign = repeatAlign,
                                                           repeatSim = repeatSim,
                                                           h.cut = repeat.clust.h.cut) %>%
      reshape2::melt() %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(value = as.character(value))
    colnames(repeatClusterAlignLong) <- c("arrayID", "positionInArray", "repeatClusterId")
    
    refTaleId <- pickRefName(align = repeatAlign, refTag = refgrep)
    repeatSimAlignLong <- convertRepeat2SimAlign(repeatAlign = repeatAlign,
                                                 repeatSim = repeatSim,
                                                 refTag = refgrep) %>%
      reshape2::melt() %>%
      dplyr::as_tibble()
    colnames(repeatSimAlignLong) <- c("arrayID", "positionInArray", "repeatSimVsRef")
    # Join with main tible
    repeatAlignLong %<>%
      dplyr::left_join(repeatClusterAlignLong,
                       by = dplyr::join_by(arrayID, positionInArray)) %>%
      dplyr::left_join(repeatSimAlignLong,
                       by = dplyr::join_by(arrayID, positionInArray))
  }
  
  
  # Building TALE tree if possible
  if (!is.null(talsim) & countOfTales > 1) {
    talsimForDendo <- talsim[talsim$TAL1 %in% arrayNames, ]
    talsimForDendo <- talsimForDendo[talsimForDendo$TAL2 %in% arrayNames, ]
    talsimForDendo <- as.matrix(reshape2::acast(talsimForDendo, TAL1 ~ TAL2, value.var = "Sim"))
    taldist <- 100 - talsimForDendo
    taldist <- taldist[arrayNames, ]
    taldist <- taldist[, arrayNames]
    taleshclust <- stats::hclust(as.dist(taldist))
  }
  
  
  # Add a symbol to designate the reference if necessary
  if (exists("refTaleId")) { # in the tibble
    repeatAlignLong$arrayID[repeatAlignLong$arrayID == refTaleId] <-  paste0(
      repeatAlignLong$arrayID[repeatAlignLong$arrayID == refTaleId],
      "_#"
    )
  }
  if (exists("refTaleId") & exists("taleshclust")) { # in the tree
    taleshclust$labels[taleshclust$labels == refTaleId] <- paste0(
      taleshclust$labels[taleshclust$labels == refTaleId],
      "_#"
    )
  }
  
  #### TODO: Bind a 'consensus' tibbe or a consensus plot if requested ####
  
  
  # Create base plot
  bp <- repeatAlignLong %>% ggplot2::ggplot(mapping = ggplot2::aes(
    x = positionInArray, y = arrayID)
  ) +
    ggplot2::scale_x_discrete(
      name = "Position in array",
      limits = factor(1:max(repeatAlignLong$positionInArray))
    ) +
    ggplot2::scale_y_discrete(name = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
  
  # COLORS in plots
  repeatClusterFillPaletteFunct <- colorRampPalette(c("#421727", "#6e2742", "#9a365c", "#b03e69", "azure2"))
  repeatClusterFillScale <- ggplot2::scale_fill_manual(name = "Repeats cluster",
                                                       drop = TRUE,
                                                       na.translate = FALSE,
                                                       palette = repeatClusterFillPaletteFunct,
                                                       guide = NULL)
  # repeatSimFillScale <- ggplot2::scale_fill_gradient(name = "Similarity relative to reference",
  #                                                    limits = c(70, 100),
  #                                                    low = "red", high = "lightgrey")
  repeatSimFillScale <- ggplot2::scale_fill_distiller(name = "Similarity relative to reference",
                                                      direction = -1)
  # labelConsensusColorScale <- ggplot2::scale_color_manual(name = "Match consensus?",
  #                                                         values = c(`TRUE` = "black",
  #                                                                    `FALSE` = "red")
  # )
  labelConsensusColorScale <- ggplot2::scale_color_manual(name = "Match consensus?",
                                                          values = c(`FALSE` = "deeppink2",
                                                                     `TRUE` = "cyan3")
  )
  # Add aesthetics as requested AND possible
  
  #### TODO: add the consensus in the plot ####
  
  if (!is.null(repeatSim) & !is.null(rvdAlign)) {
    if (fillType == "repeatSim") {
      p <- bp +
        repeatSimFillScale +
        labelConsensusColorScale +
        ggplot2::geom_label(mapping = ggplot2::aes(fill = repeatSimVsRef,
                                                   label = rvd,
                                                   color = matchConsensusRvd),
                            label.size = NA,
                            family = "mono",
                            size = 3, fontface = "bold",
                            na.rm = TRUE
        )
    } else if (fillType == "repeatClust") {
      p <- bp +
        repeatClusterFillScale +
        labelConsensusColorScale +
        ggplot2::geom_label(mapping = ggplot2::aes(fill = repeatClusterId,
                                                   label = rvd,
                                                   color = matchConsensusRvd),
                            label.size = NA,
                            family = "mono",
                            size = 3, fontface = "bold",
                            na.rm = TRUE
        )
    } else {
      logger::log_error("the fillType value must be either 'repeatSim' or 'repeatClust'")
      stop()
    }
  } else if (!is.null(repeatSim) & is.null(rvdAlign)) {
    if (fillType == "repeatSim") {
      p <- bp +
        repeatSimFillScale +
        labelConsensusColorScale +
        ggplot2::geom_label(mapping = ggplot2::aes(fill = repeatSimVsRef,
                                                   label = domCode,
                                                   color = matchConsensusRepeat),
                            label.size = NA,
                            family = "mono",
                            size = 3, fontface = "bold",
                            na.rm = TRUE
        )
    } else if (fillType == "repeatClust") {
      p <- bp +
        repeatClusterFillScale +
        labelConsensusColorScale +
        ggplot2::geom_label(mapping = ggplot2::aes(fill = repeatClusterId,
                                                   label = domCode,
                                                   color = matchConsensusRepeat),
                            label.size = NA,
                            family = "mono",
                            size = 3, fontface = "bold",
                            na.rm = TRUE
        )
    } else {
      logger::log_error("the fillType value must be either 'repeatSim' or 'repeatClust'")
      stop()
    }
  } else if (is.null(repeatSim) & !is.null(rvdAlign)) {
    p <- bp + 
      labelConsensusColorScale +
      ggplot2::geom_label(mapping = ggplot2::aes(label = rvd,
                                                 color = matchConsensusRvd),
                          fill = "grey80",
                          label.size = NA,
                          family = "mono",
                          size = 3, fontface = "bold",
                          na.rm = TRUE
      )
  } else if (is.null(repeatSim) & is.null(rvdAlign)) {
    p <- bp +
      labelConsensusColorScale +
      ggplot2::geom_label(mapping = ggplot2::aes(label = domCode,
                                                 color = matchConsensusRepeat),
                          fill = "grey80",
                          label.size = NA,
                          family = "mono",
                          size = 3, fontface = "bold",
                          na.rm = TRUE
      )
  } else {
    logger::log_error("Cannot ouput a plot based on the suppplied combination of parameter values...")
    stop()
  }
  # Merge tree and align
  if (exists("taleshclust")) {
    t <- ggtree::ggtree(ape::as.phylo(taleshclust))
    finalPlot <- p %>% aplot::insert_left(t, width = .08)
  } else {
    finalPlot <- p
  }
  print(finalPlot)
  return(finalPlot)
}


