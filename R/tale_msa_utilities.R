

##### General utility functions ####


pickRefName <- function(align, refTag = NULL) {
  # How do we select the reference TALE in an alignement?
  #   - the reference could be defined by name or by a string match in the name (eg a strain ID)
  #   - the reference could by default be defined as the longest tal and picked by
  #     ordering their names in case of ties...

  # Find refTag in seq names if provided and output the corresponding unique match
  if (! is.null(refTag)) {
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
    strippedAlignLengths <- apply(align, 1, function(seq) length(seq[! is.na(seq)]))
    longest <- rownames(align)[strippedAlignLengths == max(strippedAlignLengths)]
    ifelse(length(longest) == 1, refName <- longest, refName <- sort(longest)[1])
  }
  return(refName)
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


##### Tale classification ####

#' Grouping TALEs
#'
#' @description
#'
#'   Classifying Tal groups by hierchical clustering or by k-medoids clustering based on their similarity.
#'
#' @param method one of two methods: "hclust" (see \code{\link[stats:cutree]{cutree}}) and "k-medoids" (see \code{\link[cluster:pam]{pam}}). 
#' @param taleSim a \emph{three columns Tals similarity table} as obtained with \code{\link[tantale:runDistal]{runDistal}} in the 'tal.similarity' slot of the returned object.
#' @param plotTree logical indicating whether to plot hclust tree or not. If the method is "k-medoids", no tree will be plotted (but instead, a plot of silhoutte value).
#' @param k_test integer vector of 2 indicating the range of k to test, only available when method = "k-medoids". Note that the minimum value for k is 2.
#' @param k integer indicating number of groups you want Tals to be classified. Or only in case that method is "k-medoids", k = "auto" to automatically pick the optimum k or k = NULL to interactively pick it. Do not always trust the automatic picking, it is better to choose k interactively or test with different values.
#' @return a data frame containing name of tals from taleSim and their classified groups.
#'
#' @export
groupTales <- function(taleSim, plotTree = FALSE, k = NULL, k_test = NULL, method = "k-medoids") {
  distMat <- 100 - reshape2::acast(taleSim, formula = TAL1 ~ TAL2, value.var = "Sim")

  if (method == "k-medoids") {
    if (is.null(k_test) || !is.numeric(k_test)) stop("invalid k values!")
    if (plotTree) {
      plotTree <- FALSE
      message("tale tree will not be plotted!")
    }
    allPam <- lapply(k_test, function(kpam) {
      set.seed(7)
      kmeanClust <- cluster::pam(as.dist(distMat), kpam)
      return(as.list(kmeanClust))
      })

    silhVals <- sapply(allPam, function(a) a$silinfo$avg.width)

    if (is.null(k)) {
      plot(k_test, pch = 19, col = "cornflowerblue", silhVals, xlab = "number of groups", ylab = "average silhouette values")
      cat("Choose a number of groups:\t")
      numGroups <- as.numeric(readLines(con = stdin(), 1))
    } else if (k == "auto") {
      # numGroups <- k_test[which.max(silhVals)]
      ## I tried applying the Kneedle algorithm to find the elbow point of the curve
      ## the algorithm is in this paper: https://raghavan.usc.edu//papers/kneedle-simplex11.pdf
      ## but my knowledge in linear algebra is all gone, so I have just applied the first step 
      ## however, the result is already quite good so far for the silhoutte curve as I tested.
      ## I will go back to this if interested...
      find_elbow <- function(v, k) {
        stopifnot(length(v) == length(k))
        n <- length(v)
        a <- (v[n] - v[1])/(k[n] - k[1])
        b <- -1
        c <- (v[1]*k[n] - v[n]*k[1])/(k[n]-k[1])
        d <- sapply(1:n, function(i) abs(a*i + b*v[i] + c)/sqrt(a*a + b*b))
        return(k[which.max(d)])
      }
      numGroups <- find_elbow(silhVals, k_test)
      point_col <- sapply(k_test,
                          function(v) ifelse(v != numGroups, "cornflowerblue", "red"),
                          simplify = T)
      plot(k_test, pch = 19, col = point_col, silhVals, xlab = "number of groups", ylab = "average silhouette values")
      message(paste("The number of groups is automatically decided based on the silhoutte value:", numGroups))
    } else if (length(k) == 1 && is.numeric(k)) {
      numGroups <- k
      message(paste("Number of groups is decided based on the provided value of k:", numGroups))
      point_col <- sapply(k_test,
                          function(v) ifelse(v != numGroups, "cornflowerblue", "red"),
                          simplify = T)
      plot(k_test, pch = 19, col = point_col, silhVals, xlab = "number of groups", ylab = "average silhouette values")
    } else {
      stop("invalid k value!")
    }
    group <- allPam[[which(k_test == numGroups)]]$clustering
    taleGroups <- data.frame(name = names(group), group = group, row.names = NULL)
  } else if (method == "hclust") {
    numGroups <- k
    if (!is.null(k)) message("'k' will be ignored")
    if (!is.numeric(numGroups) || length(numGroups) != 1) stop("'k' must be specified as a number!")
    if (method == "euclidean") {
      taleTree <- hclust(
        d = dist(
          distMat,
          method = "euclidean"
          ),
        method = "ward.D"
        )
    } else if (method == "distal") {
      taleTree <- hclust(as.dist(distMat), method = "ward.D")
    }
    hi <- max(taleTree$height)
    lo <- min(taleTree$height)
    repeat {
      if (lo >= hi) stop(glue::glue("Cannot determine {numGroups} groups!"))
      cutOff <- mean(c(lo, hi))
      treeCuts <- cutree(taleTree, h = cutOff)
      if (max(treeCuts) < numGroups) {
        hi <- cutOff
      } else if (max(treeCuts) > numGroups) {
        lo <- cutOff
      } else break()
    }
    # treeCuts <- cutree(taleTree, h = cutOff)
    taleGroups <- data.frame(name = names(treeCuts), group = treeCuts, row.names = NULL)

    # plot(taleTree, xlab = "TALEs", main = )
    # abline(h = cutOff, lty = 2)
    g <- split(names(treeCuts), treeCuts)
    p <- taleTree  %>% ggtree::ggtree()
    clades <- sapply(g, function(n) tidytree::MRCA(p, n))
    p <- tidytree::groupClade(p, clades, group_name='subtree') + ggtree::aes(color=subtree)
    p <- p + ggtree::layout_dendrogram() +
      #ggtree::geom_tippoint(size=5, shape=21) +
      ggtree::geom_tiplab(ggtree::aes(label=label),
                          angle= 90, hjust=1,
                          offset = -0.4,
                          align = TRUE, color='black') +
      # ggplot2::scale_color_brewer("Groups", palette="BrBG") + # allowed maximum for palette BrBG is 11
      viridis::scale_color_viridis(discrete = T, option = "C", breaks = 1:numGroups) +
      ggplot2::geom_vline(xintercept = -(cutOff/2), linetype = 2) +
      ggtree::geom_text(x = (cutOff/2 - max(taleTree$height)/50),
                        y = 4, label = paste("half of 'cutOff' value: ", sprintf("%.2f", cutOff/2)),
                        color = "darkgrey", fontface = "plain") +
      ggplot2::xlab("Height/2") +
      ggtree::theme_dendrogram(plot.margin = ggplot2::margin(6,6,220,6))
    }





  if(plotTree == TRUE) {print(p)}
  return(taleGroups)


}



##### Tale domains sequences multiple alignment ####


#' Perform multiple alignment of TALE repeat or RVD sequences
#'
#' @description Perform multiple alignment of TALE repeat or RVD sequences using the \code{--text} mode of the \href{https://mafft.cbrc.jp/alignment/software/}{MAFFT} Multiple alignment program.
#'
#' By default, uses the simple scoring matrix defined in \href{https://mafft.cbrc.jp/alignment/software/textcomparison.html}{the text mode of MAFFT}. Users can optionally provide a custom scoring matrix.
#'
#' @param inputSeqs Either the path to a fasta file containing the TALE sequences to be aligned or a named list of vectors, each made up of the individual residues that will be aligned, as returned by the \code{\link[tantale:fa2liststr]{fa2liststr}} function or available in the \code{coded.repeats.str} slot of the object returned by the \code{\link[tantale:runDistal]{runDistal}} function.  Can also be the return value of the \code{\link[tantale:convertRepeat2RvdAlign]{convertRepeat2RvdAlign}} function if one wants to align RVD sequences.
#'
#' @param distalRepeatSims A long, three columns data frame with pairwise similarity scores between repeats as available in the \code{repeat.similarity slot} of the object returned by the \code{\link[tantale:runDistal]{runDistal}} function.
#' @param mafftOpts A character string containing additional options for the MAFFT command. This is notably useful to tweak the Gap opening and gap extension penalties.
#' @param mafftPath Path to a MAFFT installation directory. By default uses the MAFFT version included in tantale.
#' @param gapSymbol Specify a alternative symbol for gaps in the alignments.
#'
#' @return A character matrix representing the multiple alignment.
#' @export
buildRepeatMsa <- function(inputSeqs, distalRepeatSims = NULL,
                           mafftOpts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5", # "--globalpair --weighti 1 --maxiterate 1000 --reorder --op 0 --ep 5 --thread 8"
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
  if (length(inputSeqs) > 1 && is.list(inputSeqs)) {
    seqsAsVectors <- inputSeqs
    # WOULD BE BETTER to accept a Biostrings XStringSet object as imput rather than a list...
    # Just need to write a toListOfStr function based on fa2liststr and that accept both fasta file or XStringSet
  } else if (length(inputSeqs) == 1 && is.character(inputSeqs)) {
    seqsAsVectors <- fa2liststr(inputSeqs)
  } else {
    stop("##  Something is wrong with the value provided for inputSeqs. It must be either\n",
         "##  the path to a fasta file containing strings of tale sequences (space/dash-separated\n",
         "##  rvd or distal repeat IDs) or a list of vectors of individual rvd or distal repeat IDs.")
  }
  residues <- unique(unlist(seqsAsVectors))
  # Determine the type of 'elements' (rvd or repeat) contained in the sequences
  frequentRvds <- c("NN", "NG", "HD", "NI", "N*", "NS")
  if(! any(residues %in% frequentRvds)) {
    cat(glue::glue("## [{date()}] Will be assuming sequences contain repeat unit codes because
                    ##    none of the RVDs obtained from input sequences matches
                    ##    a list of 'frequent RVDs': {paste(frequentRvds, collapse = ' ')}"),
        "\n"
    )
    repeatType <- "repeatUnit"
  } else {
    cat(glue::glue("[{date()}] Input sequences are detected as RVD sequences."),
        "\n"
    )
    repeatType <- "rvds"
  }
  if( length(residues) > nrow(asciitableForMafft) ) {
    stop(glue::glue("## [{date()}] Number of unique resisues (RVDs or repeat units) must be =< 248.
                    ## Currently, your set of sequences contains {length(residues)} unique residues...")
    )
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
    cat(glue::glue("##  [{date()}] The provided similarity matrix file will be used to compute msa."), "\n")
    if (length(distalRepeatSims) > 1 &&
        (is.data.frame(distalRepeatSims) | tibble::is_tibble(distalRepeatSims))
    ) {
      repeatSims <- distalRepeatSims
    } else if (length(distalRepeatSims) == 1 && is.character(distalRepeatSims)) {
      repeatSims <- formatDistalRepeatDistMat(distalRepeatSims)
    } else {
      stop("##  Somthing is wrong with the value provided for distalRepeatSims. It must be either\n",
           "##  the path to a '*_Repeatmatrix.mat' file produced by Distal or\n",
           "##  a table like object with three columns, usually produced by the\n",
           "##  formatDistalRepeatDistMat() function.")
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
  asciiConverstionCmd <- glue::glue("{mafftPath}/mafftdir/libexec/hex2maffttext {hexFile} > {asciFile}")
  mafftCmd <-  glue::glue("{mafftPath}/mafft.bat {maffMatOpt} --text {mafftOpts} {asciFile} > {mafftAsciiOutFile}")
  MsaConversionToHexCmd <- glue::glue("{mafftPath}/mafftdir/libexec/maffttext2hex {mafftAsciiOutFile} > {mafftHexOutFile}")
  system(command = paste(asciiConverstionCmd, mafftCmd, MsaConversionToHexCmd, sep = "; "),
         ignore.stdout = FALSE, ignore.stderr = FALSE, intern = FALSE)

  # Getting msa output and converting back to alignment of residues
  msaOfHex <- Biostrings::readBStringSet(mafftHexOutFile)
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
#' Plotting of multiple alinged Tal sequences
#' @description Plot in frame of \code{\link[gplots:heatmap.2]{heatmap.2}} for Tals alignment.
#' @param talsim a \emph{three columns Tals similarity table} as obtained with \code{\link[tantale:runDistal]{runDistal}} in the 'tal.similarity' slot of the returned object.
#' @param repeatAlign a multiple Tal repeat sequences alignment in the form of a matrix as returned by \code{\link[tantale:buildRepeatMsa]{buildRepeatMsa}} or as one of the elements of the \code{SeqOfRepsAlignments} slot in the return object of the \code{\link{buildDisTalGroups}} function.
#' @param repeatSim A long, three columns data frame with pairwise similarity scores between repeats as available in the \code{repeat.similarity slot} of the object returned by the \code{\link[tantale:runDistal]{runDistal}} function. \strong{(CORRECT???!!!)}
#' @param plot.type Either \code{"repeat.similarity"}, \code{"repeat.clusters"} , \code{"repeat.clusters.with.rvd"}. Defines the type of plot that will be produced by the function. See below for details.
#' @param repeat.clust.h.cut height for tree cutting when plot type in "repeat.clusters".
#' @param rvdAlign (optional) when the rvds need to be labeled in the plot (plot.type = "repeat.similarity" or "repeat.clusters.with.rvd", a multiple Tal repeat sequences alignment in the form of a matrix as returned by \code{\link[tantale:buildRepeatMsa]{buildRepeatMsa}} or as one of the elements of the \code{SeqOfRepsAlignments} slot in the return object of the \code{\link{buildDisTalGroups}} function. 
#' @param refgrep regular expression pattern that will be used to search Tal names to select the reference in the alignment.
#' @param consensusSeq (logical) whether to display the consensus sequence when the plot type is "repeat.clusters.with.rvd".
#' @param noteColSet In case rvdSim = NULL, vector of 2 colors for rvd alignment, the first color is for matched rvds, and the second color is for mismatched ones. In the other case, more colors should be supplied.
#' @param save.path file path to save the plot. If save.path is NULL, the heatmap will be printed. If save.path is specified, the image file will be created with the format based on file extension.
#' @param ... any other arguments of \code{\link[gplots:heatmap.2]{heatmap.2}}
#' @details
#'  "repeat.similarity" plot shows RVD alignment of Tals, hierarchical relationship between them in the dendrogram, similarity between repeats alignment by the color of cells, and (if rvdSim is provided) similarity between RVDs alignment in the color of cellnotes.
#'
#'  "repeat.clusters" plot shows repeat alignment of Tals, hierarchical relationship between them in the dendrogram, and clustering groups of repeats by the color of cells.
#'  "repeat.clusters.with.rvd" plots repeat alignment of Tals with rvd labeled.
#'  If plotting from the outputs of \code{buildDistalGroups}, you supply repeatClustID/Similarity alignment to \code{forMatrix}, with \code{talsim}, \code{forCellNote} - repeat/rvd alignment, and refgrep optionally.
#'  But if you don't have these alignments, you can provide \strong{repeat alignment} to \code{forMatrix} with \code{repeatSim}, the repeat similarity data frame, the function will convert it into repeatClustID/Similarity alignment depending on the plot type. In case of \emph{repeat.clusters}, you may want to adjust the param \code{repeat.clust.h.cut} to decrease/increase the number of repeat clusters.
#'  
#' @return the return value of \code{\link[gplots:heatmap.2]{heatmap.2}}
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
  talsim <- talsim[talsim$TAL1 %in% rownames(forCellNote), ]
  talsim <- talsim[talsim$TAL2 %in% rownames(forCellNote), ]
  talsim <- as.matrix(reshape2::acast(talsim, TAL1 ~ TAL2, value.var="Sim")) # melt then unmelt ...
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
    forKeyxlab <- "AA identity"
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


