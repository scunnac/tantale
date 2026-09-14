

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
  colnames(matchConsensusLong) <- c("array_id", "position_in_array", "tales_consensus_match")
  return(matchConsensusLong)
}


##### Tale domains sequences multiple alignment ####



#' Align TALE sequences with MAFFT text mode
#'
#' Implementation behind \code{\link{tales_align}} and the deprecated
#' \code{\link{tales_align}}. Internal so that package code can call it
#' without tripping the deprecation warning.
#' Normalise a pairwise table to what MAFFT's --textmatrix wants
#'
#' MAFFT scores matches, so it wants a *similarity*: higher means more alike.
#' Accepts either the canonical \code{\link{pairwise_distances}} vocabulary
#' (\code{id1}, \code{id2}, \code{dissim}) or the legacy one
#' (\code{RepU1}, \code{RepU2}, \code{Sim}), inverting the distance on the
#' way in.
#'
#' @param x A data frame of pairwise scores.
#' @return A three-column data frame named \code{id1}, \code{id2}, \code{sim}.
#' @keywords internal
.as_mafft_score_table <- function(x) {
  nms <- names(x)
  if (all(c("id1", "id2") %in% nms)) {
    if ("dissim" %in% nms) {
      out <- data.frame(id1 = x$id1, id2 = x$id2, sim = 100 - x$dissim,
                        stringsAsFactors = FALSE)
    } else if ("sim" %in% nms) {
      out <- data.frame(id1 = x$id1, id2 = x$id2, sim = x$sim,
                        stringsAsFactors = FALSE)
    } else {
      cli::cli_abort(
        c("A pairwise table must carry {.field dissim} or {.field sim}.",
          "i" = "Got: {.field {nms}}"),
        class = c("tantale_error_msa_sim_table", "tantale_error")
      )
    }
    return(out)
  }
  if (all(c("RepU1", "RepU2", "Sim") %in% nms)) {
    out <- as.data.frame(x[, c("RepU1", "RepU2", "Sim")])
    names(out) <- c("id1", "id2", "sim")
    return(out)
  }
  cli::cli_abort(
    c("Cannot read {.arg repeat_sims}.",
      "i" = "Expected {.field id1}/{.field id2}/{.field dissim}, or the legacy {.field RepU1}/{.field RepU2}/{.field Sim}.",
      "x" = "Got: {.field {nms}}"),
    class = c("tantale_error_msa_sim_table", "tantale_error")
  )
}


#' RVD similarity matrix for aligning RVD sequences
#'
#' @description
#' The repeat-level similarity matrix is keyed by \code{dom_code}, which is
#' meaningless for an RVD sequence, so RVD alignments have had no scoring matrix
#' at all -- MAFFT treated \code{NI} and \code{NN} as no more alike than
#' \code{NI} and \code{HD}. This supplies the missing one, from the Spearman
#' correlation of each RVD's A/C/G/T preference profile
#' (\code{rvdSimDf}, derived from TALVEZ's \code{mat1}).
#'
#' @details
#' \code{XX} -- a terminus detected but not identified -- has a uniform base
#' profile, so its correlation with everything is undefined. Those cells are
#' filled as **neutral** (0): we know nothing about it, so it should neither
#' attract nor repel.
#'
#' The one exception is \code{XX} against itself, which is set to the maximum.
#' That is forced, not chosen: MAFFT produces unusable output when the diagonal
#' is not high, and a low diagonal would anyway assert that an \code{XX} must
#' *not* align with an \code{XX}, which is a stronger claim than ignorance.
#'
#' Scale is irrelevant -- MAFFT normalises the matrix, so a linear rescale or
#' offset leaves the alignment unchanged. Only relative structure
#' matters, so the correlations are used as they are.
#'
#' @param residues Character vector of the RVDs present in the alignment.
#' @return A data frame of \code{id1}, \code{id2}, \code{sim} covering
#'   every ordered pair of \code{residues}.
#' @keywords internal
.rvd_score_table <- function(residues) {
  residues <- unique(as.character(residues))
  grid <- expand.grid(id1 = residues, id2 = residues,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  key <- paste(rvdSimDf$rvd1, rvdSimDf$rvd2)
  grid$sim <- rvdSimDf$Cor[match(paste(grid$id1, grid$id2), key)]
  # unknown pairings, and anything involving XX, are neutral
  grid$sim[is.na(grid$sim)] <- 0
  # ... except a symbol against itself, which must stay high for MAFFT
  grid$sim[grid$id1 == grid$id2] <- 1
  grid
}


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
  # "rvd" opts in to the built-in RVD matrix; NULL and FALSE mean no matrix.
  # Opt-in rather than on-by-default: the matrix demonstrably makes alignments
  # more compact (fewer gaps, narrower) but there is no evidence it makes them
  # biologically better, and defaulting it on would silently change every
  # existing RVD alignment.
  if (is.null(repeat_sims) || isFALSE(repeat_sims)) {
    maffMatOpt <- ""
  } else if (identical(repeat_sims, "rvd")) {
    if (repeatType != "rvds") {
      cli::cli_abort(
        c('{.code repeat_sims = "rvd"} only applies to an RVD alignment.',
          "i" = "These sequences look like repeat unit codes."),
        class = c("tantale_error_msa_sim_table", "tantale_error")
      )
    }
    cli::cli_inform("Scoring the RVD alignment with the built-in RVD similarity matrix.")
    repeatSims <- .rvd_score_table(residues)
    repeatSims$id1 <- asciitableForMafft$hex[match(repeatSims$id1, residues)]
    repeatSims$id2 <- asciitableForMafft$hex[match(repeatSims$id2, residues)]
    colnames(repeatSims) <- NULL
    write.table(repeatSims, file = simMatHexFile, row.names = FALSE, fileEncoding = "ASCII")
    maffMatOpt <- glue::glue("--textmatrix {simMatHexFile}")
  } else if (repeatType == "rvds") {
    cli::cli_abort(
      c("A repeat similarity table cannot score an RVD alignment.",
        "i" = "It is keyed by {.field dom_code}, which has no meaning for RVDs.",
        "i" = 'Use {.code repeat_sims = "rvd"} for the built-in RVD matrix, or leave it {.code NULL} for none.'),
      class = c("tantale_error_msa_sim_table", "tantale_error")
    )
  } else {
    cli::cli_inform("The provided similarity matrix file will be used to compute msa.")
    if (length(repeat_sims) > 1 &&
        (is.data.frame(repeat_sims) | tibble::is_tibble(repeat_sims))
    ) {
      repeatSims <- .as_mafft_score_table(repeat_sims)
    } else if (length(repeat_sims) == 1 && is.character(repeat_sims)) {
      repeatSims <- .format_repeat_dist_mat(repeat_sims)
    } else {
      cli::cli_warn("Somthing is wrong with the value provided for repeat_sims. It must be either")
      cli::cli_warn("the path to a '*_Repeatmatrix.mat' file produced by Distal or")
      cli::cli_warn(paste0("table like object with three columns, usually produced by the"))
      cli::cli_abort(".format_repeat_dist_mat() function", class = c("tantale_error"))
    }

    stopifnot(all(residues %in% unique(repeatSims$id1)))
    repeatSims <- subset(repeatSims, id1 %in% residues & id2 %in% residues)
    stopifnot(all.equal(nrow(repeatSims), length(residues)^2))
    repeatSims$id1 <- asciitableForMafft$hex[match(repeatSims$id1, residues)]
    repeatSims$id2 <- asciitableForMafft$hex[match(repeatSims$id2, residues)]
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
#'   for \code{dom_code}.
#' @return A ggplot.
#' @noRd
.consensus_panel <- function(align, n_positions, pad = FALSE) {
  cons <- tales_consensus(align)
  cons <- gsub("NTERM", "N-", cons)
  cons <- gsub("CTERM", "-C", cons)
  if (isTRUE(pad)) cons <- stringr::str_pad(cons, 3, "left")
  df <- tibble::tibble(position_in_array = seq_along(cons),
                       array_id = "Consensus",
                       label = cons)
  ggplot2::ggplot(df, mapping = ggplot2::aes(x = position_in_array, y = array_id)) +
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


#' 'Nice' plotting a multiple alignment of TALE sequences
#' @description Plot TALEs msa in the ggplot2 framework.
#'
#' @details
#' Three things are decided independently, and it helps to read the figure that
#' way: what each cell *says*, what colour that text is, and what colour the
#' block behind it is.
#'
#' \strong{Cell text} is the RVD when \code{rvd_align} is supplied, and the
#' repeat code otherwise. Termini are relabelled \code{N-} and \code{-C}; an
#' unidentified terminus keeps its \code{XXXXX} code, which is deliberately not
#' mistakable for an RVD. Repeat codes are padded to three characters so columns
#' line up.
#'
#' \strong{Text colour} always answers one question: does this element match
#' the consensus of its column? Cyan for yes, pink for no. The consensus is the
#' most frequent element in the column (\code{\link{tales_consensus}}), taken
#' over RVDs when \code{rvd_align} is supplied and over repeat codes otherwise
#' -- so the text colour and the text itself always describe the same layer.
#'
#' \strong{Block fill} is what \code{fill_type} selects, and it is the only
#' part that can be unavailable:
#'
#' \tabular{lll}{
#'   \strong{fill_type} \tab \strong{shows} \tab \strong{needs} \cr
#'   \code{"repeat_clust"} \tab which cluster the repeat falls in, cut at \code{h_cut} \tab \code{repeat_sim} \cr
#'   \code{"repeat_sim"} \tab protein-sequence similarity to the reference, 0-100 \tab \code{repeat_sim} \cr
#'   \code{"rvd_sim"} \tab how alike the RVD's DNA-binding preference is to the reference's, -1 to 1 \tab \code{rvd_align} \cr
#' }
#'
#' With no \code{repeat_sim} and no \code{rvd_sim}, every block is flat grey:
#' the text still carries the consensus comparison, but there is nothing to
#' colour blocks by.
#'
#' A cell with no value for the chosen layer keeps its text and loses its
#' colour. In \code{"rvd_sim"} that is the termini, which have no DNA-binding
#' preference and so no position on a specificity scale.
#'
#' \strong{The reference} matters for both similarity fills. \code{ref_pattern}
#' is matched against the array names and must identify exactly one, otherwise
#' the default is used with a warning; by default it is the array with the most
#' non-gap parts, ties broken alphabetically. The reference row is marked with a
#' trailing \code{_#}.
#'
#' \strong{Two panels may be attached.} Supplying \code{tal_sim} with more
#' than one array adds a dendrogram panel on the left; \code{consensus = TRUE}
#' adds a consensus panel on top. When either is present the return value is an
#' \code{aplot} composition rather than a single ggplot, so modify the
#' alignment through its \code{plotlist} element rather than adding layers to
#' the result directly.
#'
#' The only mandatory argument is either \code{repeat_align} \strong{or}
#' \code{rvd_align}; what the figure shows depends on which of the optional
#' inputs you supply, as described below.
#'
#' The plot is printed and returned for further modifications is necessary.
#'
#'
#' @param tal_sim Pairwise distances between whole TALEs, as the
#'   \code{tale_distances} element of a \code{\link{tales_compare}} result.
#'   Used to build the tree panel that orders the alignment rows.
#' @param repeat_align A multiple Tal repeat sequences alignment in the form of a
#'   matrix as returned by \code{\link{tales_align}}.
#' @param repeat_sim Pairwise distances between repeat units, as the
#'   \code{domain_distances} element of a \code{\link{tales_compare}} result.
#'   Used to group repeats into clusters, and to score each repeat against the
#'   reference TALE\'s repeat at the same alignment column.
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
#'   otherwise, so it always matches whatever the cells are labelled with. It
#'   is drawn as its own panel above the alignment.
#' @param fill_type One of \code{"repeat_clust"}, \code{"repeat_sim"} or
#'   \code{"rvd_sim"}. The first two colour cells by repeat cluster or by
#'   protein-sequence similarity to the reference. \code{"rvd_sim"} colours
#'   them instead by how alike each RVD's *DNA-binding preference* is to the
#'   reference TALE's RVD at that position, on a diverging scale over
#'   \code{[-1, 1]}; it needs \code{rvd_align} but not \code{repeat_sim}.
#'   The repeat- and RVD-level views genuinely differ: \code{HD} and \code{ND}
#'   are distinct repeats with identical specificity, while repeats differing
#'   only at positions 12-13 are near-identical proteins targeting different
#'   bases.
#'   Legacy note: "repeat_clust" or "repeat_sim". If both options are
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

  # Both tables are addressed as id1/id2/dissim below. pairwise_distances()
  # also accepts the older TAL1/RepU1/Sim spellings, so either is allowed in.
  if (!is.null(tal_sim))    tal_sim    <- tibble::as_tibble(pairwise_distances(tal_sim))
  if (!is.null(repeat_sim)) repeat_sim <- tibble::as_tibble(pairwise_distances(repeat_sim))
  
  
  # Getting repeat align
  if (!is.null(repeat_align)) {
    repeatAlignLong <- repeat_align %>% reshape2::melt() %>%
    dplyr::as_tibble()
  colnames(repeatAlignLong) <- c("array_id", "position_in_array", "dom_code")
  repeatAlignLong %<>% dplyr::mutate(array_id = as.character(array_id),
                                     dom_code = stringr::str_pad(dom_code, 3, "left"))
  repeatMatchConsensusLong <- tales_consensus_match(repeat_align)
  colnames(repeatMatchConsensusLong) <- c("array_id", "position_in_array", "matchConsensusRepeat")
  repeatAlignLong %<>% dplyr::left_join(repeatMatchConsensusLong,
                                        by = dplyr::join_by(array_id, position_in_array))
  }
  
  
  # Getting rvd align if available
  if (!is.null(rvd_align)) {
    rvdAlignLong <- rvd_align %>% reshape2::melt() %>%
      dplyr::as_tibble()
    colnames(rvdAlignLong) <- c("array_id", "position_in_array", "rvd")
    rvdAlignLong %<>% dplyr::mutate(rvd = gsub("NTERM", "N-", rvd),
                                       rvd = gsub("CTERM", "-C", rvd)
    )
    
    # Tale rvd text color if possible
    # consensus RVD sequence
    # Coloring of RVDs in alignment depending on whether they match the consensus at
    # the position
    consensusRVD <- tales_consensus(rvd_align)
    rvdConsensusSeqLong <- tibble::tibble(array_id = "Consensus",
                                          position_in_array = seq_along(consensusRVD),
                                          rvd = consensusRVD,
                                          matchConsensusRvd = TRUE,
                                          dom_code = NA,
                                          repeatClusterId = NA,
                                          repeatSimVsRef = NA
    )
    rvdMatchConsensusLong <- tales_consensus_match(rvd_align)
    colnames(rvdMatchConsensusLong) <- c("array_id", "position_in_array", "matchConsensusRvd")
    # Join with rvd tible
    rvdAlignLong %<>% dplyr::left_join(rvdMatchConsensusLong,
                                          by = dplyr::join_by(array_id, position_in_array))
  }
  
  # Assign main alignment object in long format
  if (!is.null(repeat_align) & !is.null(rvd_align)) {
    repeatAlignLong %<>% dplyr::inner_join(rvdAlignLong,
                                          by = dplyr::join_by(array_id, position_in_array),
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
    colnames(repeatClusterAlignLong) <- c("array_id", "position_in_array", "repeatClusterId")
    
    refTaleId <- .pick_ref_name(align = repeat_align, ref_tag = ref_pattern)
    repeatSimAlignLong <- .repeat_to_sim_align(repeat_align = repeat_align,
                                                 repeat_sim = repeat_sim,
                                                 ref_tag = ref_pattern) %>%
      reshape2::melt() %>%
      dplyr::as_tibble()
    colnames(repeatSimAlignLong) <- c("array_id", "position_in_array", "repeatSimVsRef")
    # Join with main tible
    repeatAlignLong %<>%
      dplyr::left_join(repeatClusterAlignLong,
                       by = dplyr::join_by(array_id, position_in_array)) %>%
      dplyr::left_join(repeatSimAlignLong,
                       by = dplyr::join_by(array_id, position_in_array))
  }
  

  # joining RVD similarity relative to the reference
  # This is the RVD-level counterpart of repeatSimVsRef: that one scores protein
  # sequence similarity, this one scores how alike two RVDs' DNA-binding
  # preferences are. They come apart -- HD and ND are different repeats with
  # identical specificity, while repeats differing only at positions 12-13 are
  # nearly identical proteins targeting different bases.
  if (!is.null(rvd_align)) {
    if (!exists("refTaleId")) {
      refTaleId <- .pick_ref_name(align = rvd_align, ref_tag = ref_pattern)
    }
    rvdSimAlignLong <- .rvd_to_match_align(rvd_align = rvd_align,
                                           ref_tag = ref_pattern) %>%
      reshape2::melt() %>%
      dplyr::as_tibble()
    colnames(rvdSimAlignLong) <- c("array_id", "position_in_array", "rvdSimVsRef")
    rvdSimAlignLong %<>% dplyr::mutate(array_id = as.character(array_id))
    repeatAlignLong %<>% dplyr::left_join(rvdSimAlignLong,
                                          by = dplyr::join_by(array_id, position_in_array))
  }
  
  # Building TALE tree if possible
  if (!is.null(tal_sim) & countOfTales > 1) {
    talsimForDendo <- tal_sim[tal_sim$id1 %in% arrayNames, ]
    talsimForDendo <- talsimForDendo[talsimForDendo$id2 %in% arrayNames, ]
    taldist <- as.matrix(reshape2::acast(talsimForDendo, id1 ~ id2, value.var = "dissim"))
    taldist <- taldist[arrayNames, ]
    taldist <- taldist[, arrayNames]
    taleshclust <- stats::hclust(as.dist(taldist))
  }
  
  
  # Add a symbol to designate the reference if necessary
  if (exists("refTaleId")) { # in the tibble
    repeatAlignLong$array_id[repeatAlignLong$array_id == refTaleId] <-  paste0(
      repeatAlignLong$array_id[repeatAlignLong$array_id == refTaleId],
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
    x = position_in_array, y = array_id)
  ) +
    ggplot2::scale_x_discrete(
      name = "Position in array",
      limits = factor(1:max(repeatAlignLong$position_in_array))
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
  # The RVD score is a correlation on [-1, 1], so it wants a diverging scale
  # centred on zero rather than the sequential one used for repeat similarity.
  rvdSimFillScale <- ggplot2::scale_fill_gradient2(
    name = "RVD specificity vs reference",
    limits = c(-1, 1), midpoint = 0,
    low = "#B2182B", mid = "grey92", high = "#2166AC", na.value = "grey80")
  # labelConsensusColorScale <- ggplot2::scale_color_manual(name = "Match consensus?",
  #                                                         values = c(`TRUE` = "black",
  #                                                                    `FALSE` = "red")
  # )
  labelConsensusColorScale <- ggplot2::scale_color_manual(name = "Match consensus?",
                                                          values = c(`FALSE` = "deeppink2",
                                                                     `TRUE` = "cyan3")
  )
  # Add aesthetics as requested AND possible
  
  if (identical(fill_type, "rvd_sim")) {
    if (is.null(rvd_align)) {
      cli::cli_abort(
        c('{.code fill_type = "rvd_sim"} needs an {.arg rvd_align}.',
          "i" = "It scores each RVD against the reference TALE's RVD at that position."),
        class = c("tantale_error_msa_layer", "tantale_error")
      )
    }
    p <- bp +
      rvdSimFillScale +
      labelConsensusColorScale +
      ggplot2::geom_label(mapping = ggplot2::aes(fill = rvdSimVsRef,
                                                 label = rvd,
                                                 color = matchConsensusRvd),
                          label.size = NA,
                          family = "mono",
                          size = 3, fontface = "bold",
                          na.rm = TRUE
      )
  } else if (!is.null(repeat_sim) & !is.null(rvd_align)) {
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
      cli::cli_abort("{.arg fill_type} must be {.val repeat_clust}, {.val repeat_sim} or {.val rvd_sim}.", class = c("tantale_error_msa_layer", "tantale_error"))
    }
  } else if (!is.null(repeat_sim) & is.null(rvd_align)) {
    if (fill_type == "repeat_sim") {
      p <- bp +
        repeatSimFillScale +
        labelConsensusColorScale +
        ggplot2::geom_label(mapping = ggplot2::aes(fill = repeatSimVsRef,
                                                   label = dom_code,
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
                                                   label = dom_code,
                                                   color = matchConsensusRepeat),
                            label.size = NA,
                            family = "mono",
                            size = 3, fontface = "bold",
                            na.rm = TRUE
        )
    } else {
      cli::cli_abort("{.arg fill_type} must be {.val repeat_clust}, {.val repeat_sim} or {.val rvd_sim}.", class = c("tantale_error_msa_layer", "tantale_error"))
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
      ggplot2::geom_label(mapping = ggplot2::aes(label = dom_code,
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
                       n_positions = max(repeatAlignLong$position_in_array),
                       pad = is.null(rvd_align)),
      height = consensusHeight
    )
  }

  print(finalPlot)
  return(finalPlot)
}


