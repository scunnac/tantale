
# load(file.path("/home/cunnac/TEMP/test_tantale/mining.RData"))
# saveRDS(distalr_deci_output, file = testthat::test_path("data_for_tests", "sampleDistalrOutput.rds"))
# saveRDS(grp, file = testthat::test_path("data_for_tests", "sampleDistalrGroups.rds"))
# saveRDS(repeatMsaByGroup_withSim, file = testthat::test_path("data_for_tests", "repeatMsaByGroup.rds"))

load_all()
distalrOut <- readRDS(file = testthat::test_path("data_for_tests", "sampleDistalrOutput.rds"))
#taleGroups <- readRDS(file = testthat::test_path("data_for_tests", "sampleDistalrGroups.rds"))
repeatMsaByGroup <- readRDS(file = testthat::test_path("data_for_tests", "repeatMsaByGroup.rds"))


# all(taleGroups$name %in% distalrOut$tal.similarity$TAL1) &
#   all(taleGroups$name %in% distalrOut$tal.similarity$TAL2) 

#'  "repeat.similarity" plot shows RVD alignment of Tals, hierarchical relationship
#'   between them in the dendrogram, similarity between repeats alignment by the color
#'    of cells, and (if rvdSim is provided) similarity between RVDs alignment in the 
#'    color of cellnotes.#'
#'    
#'  "repeat.clusters" plot shows repeat alignment of Tals, hierarchical relationship
#'   between them in the dendrogram, and clustering groups of repeats by the color of cells.
#'   
#'  "repeat.clusters.with.rvd" plots repeat alignment of Tals with rvd labeled.
#'  If plotting from the outputs of \code{buildDistalGroups},
#'  you supply repeatClustID/Similarity alignment to \code{forMatrix}, with
#'  \code{talsim}, \code{forCellNote} - repeat/rvd alignment, and refgrep optionally.
ggplotTalesMsa <- function(repeatAlign,
                           talsim = NULL,
                           rvdAlign = NULL,
                           repeatSim = NULL,
                           repeat.clust.h.cut = 90,
                           refgrep = NULL,
                           consensusSeq = FALSE,
                           noteColSet = NULL,
                           fillType = c("repeat clusters", "repeat similarities")) {
  
}

fillType = "repeatSim" # "repeatClust"
talsim = distalrOut$tal.similarity
repeatAlign = repeatMsaByGroup[[6]]
repeatSim = distalrOut$repeat.similarity
rvdAlign = convertRepeat2RvdAlign(repeatAlign = repeatAlign,
                                  repeat2RvdMapping = getRepeat2RvdMappingFromDistalr(distalrOut$taleParts))

repeat.clust.h.cut = 90
refgrep = NULL
consensusSeq = FALSE
noteColSet = NULL




# Getting repeat align
countOfTales <- nrow(repeatAlign)
if (countOfTales < 1) {
  logger::log_error("The provided input repeatAlign matrix has less than one sequence. Cannot proceed...")
  stop()
}
repeatAlignLong <- repeatAlign %>% reshape2::melt() %>%
  dplyr::as_tibble()
colnames(repeatAlignLong) <- c("arrayID", "positionInArray", "domCode")
repeatAlignLong %<>% dplyr::mutate(arrayID = as.character(arrayID))

# Getting rvd align if available and joining
if(!is.null(rvdAlign)) {
  rvdAlignLong <- rvdAlign %>% reshape2::melt() %>%
    dplyr::as_tibble()
  colnames(rvdAlignLong) <- c("arrayID", "positionInArray", "rvd")
  # Join with main tible
  repeatAlignLong %<>% dplyr::left_join(rvdAlignLong,
                                 by = dplyr::join_by(arrayID, positionInArray))
  repeatAlignLong %<>% dplyr::mutate(rvd = gsub("NTERM", "N-", rvd),
                                     rvd = gsub("CTERM", "-C", rvd)
                                     )
}

# Tale rvd text color if possible and requested
# consensus RVD sequence
# Coloring of RVDs in alignment depending on whether they match the consensus at
# the position
if (!is.null(rvdAlign)) {
  consensusRVD <- sapply(1:ncol(rvdAlign), function(x) {
    allRVDs <- rvdAlign[,x]
    freq <- sapply(unique(allRVDs), function(p) S4Vectors::countMatches(p, allRVDs))
    unique(allRVDs)[which.max(freq)]
  })
  rvdConsensusSeqLong <- tibble::tibble(arrayID = "Consensus",
                                positionInArray = seq_along(consensusRVD),
                                rvd = consensusRVD)
  rvdcol <- rvdAlign
  for (k in 1:ncol(rvdcol)){
    rvd <- consensusRVD[k]
    if (is.na(rvd)) {
      rvdcol[,k] <- FALSE
    } else {
      rvdcol[,k] <- ifelse(toupper(rvdcol[,k]) == toupper(rvd), TRUE, FALSE)
    }
  }
  rvdMatchConsensusLong <- rvdcol %>% reshape2::melt() %>%
    dplyr::as_tibble()
  colnames(rvdMatchConsensusLong) <- c("arrayID", "positionInArray", "matchConsensusRvd")
  # Join with main tible
  repeatAlignLong %<>% dplyr::left_join(rvdMatchConsensusLong,
                                 by = dplyr::join_by(arrayID, positionInArray))
}

# joining repeat cluster if possible and requested
# joining repeat similarity relative to ref if possible and requested
if (!is.null(repeatSim)) {
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
# Building TALE tree
if (!is.null(talsim) & countOfTales > 1) {
  talsimForDendo <- talsim[talsim$TAL1 %in% rownames(repeatAlign), ]
  talsimForDendo <- talsimForDendo[talsimForDendo$TAL2 %in% rownames(repeatAlign), ]
  talsimForDendo <- as.matrix(reshape2::acast(talsimForDendo, TAL1 ~ TAL2, value.var = "Sim")) # melt then unmelt ...
  taldist <- 100 - talsimForDendo
  taldist <- taldist[rownames(repeatAlign), ]
  taldist <- taldist[, rownames(repeatAlign)]
  taleshclust <- stats::hclust(as.dist(taldist))
}

# 'maximum' tible columns:
# y (position of the array) | yLabel (array label) | x (position of domain in align) | textLabel (what will be displayed as text in domain labels) |
# textColor (?) | fillColor (repeat sim relative to ref, repeat cluster)
# Need a 'consensus' tible
# Add a symbol to designate the reference
if (exists("refTaleId")) {
  repeatAlignLong$arrayID[repeatAlignLong$arrayID == refTaleId] <-  paste0(
    repeatAlignLong$arrayID[repeatAlignLong$arrayID == refTaleId],
    "_#"
  )
}
if (exists("refTaleId") & exists("taleshclust")) {
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

# Add aesthetics as requested
if (countOfTales == 1) {
  if (!is.null(rvdAlign)) {
    finalPlot <- p
  } else {
    finalPlot <- p
  }
  print(finalPlot)
  return(finalPlot)
}

if (!is.null(repeatSim) & !is.null(rvdAlign)) {
  p <- bp + 
    ggplot2::scale_fill_gradient(name = "Similarity relative to reference",
                                 limits = c(70, 100),
                                 low = "red", high = "lightgrey"
    ) +
    ggplot2::scale_color_manual(name = "Match consensus RVD",
                                values = c(`TRUE` = "black", `FALSE` = "red")
    ) +
    ggplot2::geom_label(mapping = ggplot2::aes(fill = repeatSimVsRef, label = rvd,
                                               color = matchConsensusRvd),
                        label.size = NA,
                        family = "mono",
                        size = 2.5,
                        na.rm = TRUE
    )
  
} else if (!is.null(repeatSim) & is.null(rvdAlign)) {
  p <- bp + 
    ggplot2::scale_fill_discrete(name = "Repeats cluster",
                                 drop = TRUE,
                                 na.translate = FALSE) +
    ggplot2::geom_label(mapping = ggplot2::aes(fill = repeatClusterId,
                                               label = domCode),
                        color = "cyan",
                        label.size = NA,
                        family = "mono",
                        size = 2.5,
                        na.rm = TRUE
    )
  
  
} else if (is.null(repeatSim) & !is.null(rvdAlign)) {
  
  
  
} else if (is.null(repeatSim) & is.null(rvdAlign)) {
  
} else {
  
}



forCol <- colorRampPalette(c("#421727", "#6e2742", "#9a365c", "#b03e69", "#ffffff"))


# Merge tree and align
t <- ggtree::ggtree(ape::as.phylo(taleshclust))
  finalPlot <- (p) %>% aplot::insert_left(t, width = .1)

print(finalPlot)
return(finalPlot)
# Output : plot and return aplot plot


##############""












