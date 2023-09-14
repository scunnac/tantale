
load(file.path("/home/cunnac/TEMP/test_tantale/mining.RData"))
saveRDS(distalr_deci_output, file = testthat::test_path("data_for_tests", "sampleDistalrOutput.rds"))
saveRDS(grp, file = testthat::test_path("data_for_tests", "sampleDistalrGroups.rds"))



distalrOut <- readRDS(file = testthat::test_path("data_for_tests", "sampleDistalrOutput.rds"))
taleGroups <- readRDS(file = testthat::test_path("data_for_tests", "sampleDistalrGroups.rds"))

taleGroups$tellTaleId %in% distalrOut$tal.similarity

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

plot.type = 
talsim = 
repeatAlign
rvdAlign = NULL
repeatSim
repeat.clust.h.cut = 90
refgrep = NULL
consensusSeq = FALSE
noteColSet = NULL
save.path
...


# Getting repeat align

# If 0 sequences error
# if one sequence, no tree

# Getting rvd align if available and joinning


# joining repeat cluster if possible and requested


# joining repeat similarity relative to ref if possible and requested


# Tale rvd text color if possible and requested


# 'maximum' tible columns:
# y (position of the array) | yLabel (array label) | x (position of domain in align) | textLabel (what will be displayed as text in domain labels) |
# textColor (?) | fillColor (repeat sim relative to ref, repeat cluster)
# Need a 'consensus' tible


# Building TALE tree
###!!! THIS FAILS IF repeatAlign has a single sequence
talsim <- talsim[talsim$TAL1 %in% rownames(forCellNote), ]
talsim <- talsim[talsim$TAL2 %in% rownames(forCellNote), ]
talsim <- as.matrix(reshape2::acast(talsim, TAL1 ~ TAL2, value.var = "Sim")) # melt then unmelt ...


taldist <- 100 -talsim
taldist <- taldist[rownames(forCellNote), ]
taldist <- taldist[, rownames(forCellNote)]

clust <- hclust(as.dist(taldist))
forRowv <- as.dendrogram(clust)

# Create base plor


# Add aesthetics as requested

# Merge tree and align

# Output : plot and return aplot plot


##############""


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















