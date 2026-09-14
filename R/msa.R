

##### General utility functions ####

.pick_ref_name <- function(align, ref_tag = NULL) {
  # How do we select the reference TALE in an alignement?
  #   - the reference could be defined by name or by a string match in the name (eg a strain ID)
  #   - the reference could by default be defined as the longest tal and picked by
  #     ordering their names in case of ties...

  # Find ref_tag in seq names if provided and output the corresponding unique match
  if (!is.null(ref_tag)) {
    match <- grepl(ref_tag, rownames(align))
    if (sum(match) != 1) {
      warning("Cannot identify a single unambiguous sequence to define as a reference using the string in ref_tag.\n",
              "Using the default method for reference selection.")
      ref_tag <- NULL
    } else {
      refName <- rownames(align)[match]
    }
  }
  # If no ref_tag is provided, pick the longest seq(s) and if there are ties, pick the first one alphabetically
  if (is.null(ref_tag)) {
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
#' @family TALE alignment
tales_consensus <- function(align) {
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
  align <- align
  for (k in 1:ncol(align)){
    rept <- consensus[k]
    if (is.na(rept)) {
      align[,k] <- FALSE
    } else {
      align[,k] <- ifelse(toupper(align[,k]) == toupper(rept), TRUE, FALSE)
    }
  }
  if (!long) return(align)
  matchConsensusLong <- align %>% reshape2::melt() %>%
    dplyr::as_tibble()
  colnames(matchConsensusLong) <- c("arrayID", "positionInArray", "tales_consensus_match")
  return(matchConsensusLong)
}


##### Tale domains sequences multiple alignment ####



#' Align TALE sequences with MAFFT text mode
#'
#' Implementation behind \code{\link{tales_align}} and the deprecated
#' \code{\link{tales_align}}. Internal so that package code can call it
#' without tripping the deprecation warning.
#' @noRd
.build_repeat_msa <- function(input_seqs, sep = " ", repeat_sims = NULL,
                           mafft_opts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
                           mafft_path = system.file("tools", "mafft-linux64",package = "tantale", mustWork = TRUE),
                           gap_symbol = NA) {
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
  seqsAsVectors <- suppressWarnings(.split_list(input_seqs, sep = sep))
  residues <- unique(unlist(seqsAsVectors))
  
  # Deals with cases where the nomber of sequences is < 2
  if (length(seqsAsVectors) == 0L) {
    cli::cli_warn("The provided object in input_seqs is empty. Returning an empty matrix")
    return(matrix())
  }
  if (length(seqsAsVectors) == 1L) {
    cli::cli_inform("The provided object in input_seqs has only one sequence. Returning it as a matrix.")
    msaOfResiduesAsMatrix <- as.matrix(as.data.frame(seqsAsVectors))
    msaOfResiduesAsMatrix <- matrix(msaOfResiduesAsMatrix, nrow = 1)
    rownames(msaOfResiduesAsMatrix) <- colnames(as.data.frame(seqsAsVectors))
    colnames(msaOfResiduesAsMatrix) <- 1:length(msaOfResiduesAsMatrix)
    return(msaOfResiduesAsMatrix)
  }
  
  # Determine the type of 'elements' (rvd or repeat) contained in the sequences
  frequentRvds <- c("NN", "NG", "HD", "NI", "N*", "NS")
  if(! any(residues %in% frequentRvds)) {
    cli::cli_inform(paste0("Will be assuming sequences contain repeat unit codes because ",
                     "none of the RVDs obtained from input sequences matches ",
                     "a list of 'frequent RVDs': {paste(frequentRvds, collapse = ' ')}"))
    repeatType <- "repeatUnit"
  } else {
    cli::cli_inform("Input sequences are detected as RVD sequences.")
    repeatType <- "rvds"
  }
  if( length(residues) > nrow(asciitableForMafft) ) {
    cli::cli_warn("Number of unique resisues (RVDs or repeat units) must be =< 248.")
    cli::cli_abort(paste0("Currently, your set of sequences contains {length(residues)} unique residues..."), class = c("tantale_error"))
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
  if(is.null(repeat_sims) || repeatType == "rvds") {
    maffMatOpt <- ""
  } else if (!is.null(repeat_sims)) {
    cli::cli_inform("The provided similarity matrix file will be used to compute msa.")
    if (length(repeat_sims) > 1 &&
        (is.data.frame(repeat_sims) | tibble::is_tibble(repeat_sims))
    ) {
      repeatSims <- repeat_sims[,c("RepU1", "RepU2", "Sim")]
    } else if (length(repeat_sims) == 1 && is.character(repeat_sims)) {
      repeatSims <- .format_repeat_dist_mat(repeat_sims)
    } else {
      cli::cli_warn("Somthing is wrong with the value provided for repeat_sims. It must be either")
      cli::cli_warn("the path to a '*_Repeatmatrix.mat' file produced by Distal or")
      cli::cli_warn(paste0("table like object with three columns, usually produced by the"))
      cli::cli_abort(".format_repeat_dist_mat() function", class = c("tantale_error"))
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
  cli::cli_inform("Now running MAFFT (Copyright 2002-2007 Kazutaka Katoh) on TALE array sequences.")
  asciiConverstionCmd <- glue::glue("{mafft_path}/mafftdir/libexec/hex2maffttext {hexFile} > {asciFile}")
  mafftCmd <-  glue::glue("{mafft_path}/mafft.bat {maffMatOpt} --text {mafft_opts} {asciFile} > {mafftAsciiOutFile}")
  MsaConversionToHexCmd <- glue::glue("{mafft_path}/mafftdir/libexec/maffttext2hex {mafftAsciiOutFile} > {mafftHexOutFile}")
  res <- system(command = paste(asciiConverstionCmd, mafftCmd, MsaConversionToHexCmd, sep = "; "),
         ignore.stdout = FALSE, ignore.stderr = FALSE, intern = FALSE)

  # Getting msa output and converting back to alignment of residues
  msaOfHex <- Biostrings::readBStringSet(mafftHexOutFile)
  if (length(msaOfHex) == 0L) {
    cli::cli_warn("MAFFT failled to complete sucessfully...")
    cli::cli_abort("MAFFT exit status: {res}", class = c("tantale_error"))
  } else {
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
      seqOfResidues[is.na(seqOfResidues)] <- gap_symbol
      #cat("Seq of residues: ", seqOfResidues, "\n")
      seqOfResidues
    }
    )
  )
  rownames(msaOfResiduesAsMatrix) <- names(msaOfHex)
  colnames(msaOfResiduesAsMatrix) <- 1:ncol(msaOfResiduesAsMatrix)
  return(msaOfResiduesAsMatrix)
}








#' 'Nice' plotting a multiple alignment of TALE sequences
#' @description Plot TALEs msa in the ggplot2 framework.
#'
#' @details This function as a similar purpose as
#' the retired \code{msa_heatmap()} (see \code{inst/legacy/}) but is implemented with
#' \code{\link[ggplot2:ggplot]{ggplot}}. It is more versatile (takes single row
#' matrices of alignment) but a bit slower.
#'
#' The type of plot that you will get will depend on the provided information in
#' the form of parameter values See the tantale website for detailed usage
#' cases.
#'
#' The only mandatory argument is either \code{repeat_align} \strong{or}
#' \code{rvd_align}.
#'
#' The plot is printed and returned for further modifications is necessary.
#'
#'
#' @param tal_sim a \emph{three columns Tals similarity table} as obtained with
#'   \code{\link{tales_compare}} in the \code{tale_distances} element of
#'   the returned object.
#' @param repeat_align A multiple Tal repeat sequences alignment in the form of a
#'   matrix as returned by \code{\link{tales_align}}.
#' @param repeat_sim A long, three columns data frame with pairwise similarity
#'   scores between repeats as available in the \code{domain_distances} element of
#'   the object returned by the \code{\link{tales_compare}} function.
#' @param h_cut height for tree cutting when defining domain/repeat
#'   clusters.
#' @param rvd_align A multiple Tal RVD sequences alignment in the form of a
#'   matrix as returned by \code{\link[tantale:repeat_to_rvd_align]{repeat_to_rvd_align}}
#'   or \code{\link{tales_align}}.
#' @param ref_pattern Regular expression pattern that will be used to search TALE
#'   names to select the reference in the alignment.
#' @param consensus (logical) Whether to add a consensus row above the
#'   alignment. The consensus is the most frequent element in each column, taken
#'   from \code{rvd_align} when supplied and from \code{repeat_align}
#'   otherwise, so it always matches whatever the cells are labelled with. It is
#'   drawn as its own panel rather than an extra alignment row, because
#'   \code{aplot} reorders the alignment's y axis onto the tree's leaves and
#'   would drop a row the tree has no leaf for.
#' @param fill_type Either "repeat_clust" or "repeat_sim". If both options are
#'   possible because the necessary information is there (at least a
#'   \code{repeat_sim} value), this argument will decide what type of 'box color
#'   filling' is employed and it is either based on the cluster where the repeat
#'   falls after clustering all the repeat in the alignment or it is based on
#'   the amino acid similarity between a repeat at a position and the repeat of
#'   the 'reference' TALE at this position.
#' @return An \code{\link[aplot:insert_left]{aplot}} object.
#' 
#' @export
#' @family TALE plots
#' Build the one-row consensus panel used by plot_tales_msa()
#'
#' Returns a standalone ggplot holding a single "Consensus" row, styled to match
#' the main alignment so the two read as one figure when composed with aplot.
#'
#' This has to be a separate panel rather than an extra row of the alignment:
#' \code{aplot::insert_left()} reorders the main plot's y axis onto the tree's
#' leaves, and a y level with no matching leaf is silently dropped -- the
#' consensus row simply disappears.
#'
#' @param align The alignment matrix to take the consensus of.
#' @param n_positions Width of the alignment, so the x scale matches the main plot.
#' @param pad Whether to pad labels to three characters, as the main plot does
#'   for \code{domCode}.
#' @return A ggplot.
#' @noRd
.consensus_panel <- function(align, n_positions, pad = FALSE) {
  cons <- tales_consensus(align)
  cons <- gsub("NTERM", "N-", cons)
  cons <- gsub("CTERM", "-C", cons)
  if (isTRUE(pad)) cons <- stringr::str_pad(cons, 3, "left")
  df <- tibble::tibble(positionInArray = seq_along(cons),
                       arrayID = "Consensus",
                       label = cons)
  ggplot2::ggplot(df, mapping = ggplot2::aes(x = positionInArray, y = arrayID)) +
    ggplot2::geom_label(mapping = ggplot2::aes(label = label),
                        fill = "grey92", color = "grey15",
                        label.size = NA, family = "mono",
                        size = 3, fontface = "bold", na.rm = TRUE) +
    ggplot2::scale_x_discrete(name = NULL, limits = factor(1:n_positions)) +
    ggplot2::scale_y_discrete(name = NULL,
                              expand = ggplot2::expansion(mult = c(0.15, 0.15))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   # keep the vertical rules so columns stay traceable between
                   # this panel and the alignment below it
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
}


plot_tales_msa <- function(repeat_align,
                           tal_sim = NULL,
                           rvd_align = NULL,
                           repeat_sim = NULL,
                           h_cut = 10,
                           ref_pattern = NULL,
                           consensus = FALSE,
                           fill_type = "repeat_clust" #"repeat_sim"
) {
  
  # Arguments checking
  if (is.null(rvd_align) & is.null(repeat_align)) {
    cli::cli_abort("You must provide at least either a value for `repeat_align` or for `rvd_align`", class = c("tantale_error"))
  }
  if (!is.null(rvd_align)) {
    countOfTales <- nrow(rvd_align)
    arrayNames <- rownames(rvd_align)
  }
  if (!is.null(repeat_align)) {
    countOfTales <-  nrow(repeat_align)
    arrayNames <- rownames(repeat_align)
  }
  if (!is.null(repeat_align) & is.null(nrow(repeat_align))) {
    cli::cli_abort(paste0("Check the provided input repeat_align matrix.",
                      "It may conain a single sequence that was coerced to vector rather than remaining a matrix...",
                      .sep = " "), class = c("tantale_error"))
  }
  if (!is.null(rvd_align) & is.null(nrow(rvd_align))) {
    cli::cli_abort(paste0("Check the provided input rvd_align matrix.",
                      "It may conain a single sequence that was coerced to vector rather than remaining a matrix...",
                      .sep = " "), class = c("tantale_error"))
  }
  if (countOfTales < 1) {
    cli::cli_abort("The provided input repeat_align matrix has less than one sequence. Cannot proceed...", class = c("tantale_error"))
  }
  
  
  # Getting repeat align
  if (!is.null(repeat_align)) {
    repeatAlignLong <- repeat_align %>% reshape2::melt() %>%
    dplyr::as_tibble()
  colnames(repeatAlignLong) <- c("arrayID", "positionInArray", "domCode")
  repeatAlignLong %<>% dplyr::mutate(arrayID = as.character(arrayID),
                                     domCode = stringr::str_pad(domCode, 3, "left"))
  repeatMatchConsensusLong <- tales_consensus_match(repeat_align)
  colnames(repeatMatchConsensusLong) <- c("arrayID", "positionInArray", "matchConsensusRepeat")
  repeatAlignLong %<>% dplyr::left_join(repeatMatchConsensusLong,
                                        by = dplyr::join_by(arrayID, positionInArray))
  }
  
  
  # Getting rvd align if available
  if (!is.null(rvd_align)) {
    rvdAlignLong <- rvd_align %>% reshape2::melt() %>%
      dplyr::as_tibble()
    colnames(rvdAlignLong) <- c("arrayID", "positionInArray", "rvd")
    rvdAlignLong %<>% dplyr::mutate(rvd = gsub("NTERM", "N-", rvd),
                                       rvd = gsub("CTERM", "-C", rvd)
    )
    
    # Tale rvd text color if possible
    # consensus RVD sequence
    # Coloring of RVDs in alignment depending on whether they match the consensus at
    # the position
    consensusRVD <- tales_consensus(rvd_align)
    rvdConsensusSeqLong <- tibble::tibble(arrayID = "Consensus",
                                          positionInArray = seq_along(consensusRVD),
                                          rvd = consensusRVD,
                                          matchConsensusRvd = TRUE,
                                          domCode = NA,
                                          repeatClusterId = NA,
                                          repeatSimVsRef = NA
    )
    rvdMatchConsensusLong <- tales_consensus_match(rvd_align)
    colnames(rvdMatchConsensusLong) <- c("arrayID", "positionInArray", "matchConsensusRvd")
    # Join with rvd tible
    rvdAlignLong %<>% dplyr::left_join(rvdMatchConsensusLong,
                                          by = dplyr::join_by(arrayID, positionInArray))
  }
  
  # Assign main alignment object in long format
  if (!is.null(repeat_align) & !is.null(rvd_align)) {
    repeatAlignLong %<>% dplyr::inner_join(rvdAlignLong,
                                          by = dplyr::join_by(arrayID, positionInArray),
                                          unmatched = "error",
                                          relationship = "one-to-one")
  } else if (!is.null(repeat_align) & is.null(rvd_align)) {
    repeatAlignLong <- repeatAlignLong
  } else if (is.null(repeat_align)) {
    repeatAlignLong <- rvdAlignLong
  } else {
    stop("something wrong with parameters values")
  }

  # joining repeat cluster if possible
  # joining repeat similarity relative to ref
  if (!is.null(repeat_sim) & !is.null(repeat_align)) {
    repeatClusterAlignLong <- .repeat_to_cluster_align(repeat_align = repeat_align,
                                                           repeat_sim = repeat_sim,
                                                           h_cut = h_cut) %>%
      reshape2::melt() %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(value = as.character(value))
    colnames(repeatClusterAlignLong) <- c("arrayID", "positionInArray", "repeatClusterId")
    
    refTaleId <- .pick_ref_name(align = repeat_align, ref_tag = ref_pattern)
    repeatSimAlignLong <- .repeat_to_sim_align(repeat_align = repeat_align,
                                                 repeat_sim = repeat_sim,
                                                 ref_tag = ref_pattern) %>%
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
  if (!is.null(tal_sim) & countOfTales > 1) {
    talsimForDendo <- tal_sim[tal_sim$TAL1 %in% arrayNames, ]
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
  # scale_fill_manual() takes `values`, not `palette`: the name collided with
  # the `palette` discrete_scale() supplies internally, so every call to this
  # function aborted with "formal argument 'palette' matched by multiple actual
  # arguments" regardless of fill_type. discrete_scale() is the scale that
  # actually accepts a palette *function*, which is what is wanted here.
  repeatClusterFillScale <- ggplot2::discrete_scale(aesthetics = "fill",
                                                    name = "Repeats cluster",
                                                    palette = repeatClusterFillPaletteFunct,
                                                    drop = TRUE,
                                                    na.translate = FALSE,
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
  
  if (!is.null(repeat_sim) & !is.null(rvd_align)) {
    if (fill_type == "repeat_sim") {
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
    } else if (fill_type == "repeat_clust") {
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
      cli::cli_abort("the fill_type value must be either 'repeat_sim' or 'repeatClust'", class = c("tantale_error"))
    }
  } else if (!is.null(repeat_sim) & is.null(rvd_align)) {
    if (fill_type == "repeat_sim") {
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
    } else if (fill_type == "repeat_clust") {
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
      cli::cli_abort("the fill_type value must be either 'repeat_sim' or 'repeatClust'", class = c("tantale_error"))
    }
  } else if (is.null(repeat_sim) & !is.null(rvd_align)) {
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
  } else if (is.null(repeat_sim) & is.null(rvd_align)) {
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
    cli::cli_abort("Cannot ouput a plot based on the suppplied combination of parameter values...", class = c("tantale_error"))
  }
  # Merge tree and align
  if (exists("taleshclust")) {
    t <- ggtree::ggtree(ape::as.phylo(taleshclust))
    finalPlot <- p %>% aplot::insert_left(t, width = .08)
  } else {
    finalPlot <- p
  }

  # Consensus as its own panel on top. See .consensus_panel() for why it cannot
  # simply be another row of the alignment.
  if (isTRUE(consensus)) {
    consensusAlign <- if (!is.null(rvd_align)) rvd_align else repeat_align
    # aplot's height is a *ratio* of the main plot, so a fixed value would grow
    # with the number of arrays -- several rows tall for a large group. Scale it
    # so the consensus stays about one alignment row high whatever the count.
    consensusHeight <- max(0.08, min(0.45, 1.0 / countOfTales))
    finalPlot <- aplot::insert_top(
      finalPlot,
      .consensus_panel(consensusAlign,
                       n_positions = max(repeatAlignLong$positionInArray),
                       pad = is.null(rvd_align)),
      height = consensusHeight
    )
  }

  print(finalPlot)
  return(finalPlot)
}


