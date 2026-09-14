# Retired: msa_heatmap()
#
# Superseded by plot_tales_msa(), which returns a ggplot rather than drawing to
# a device. Its six plot_type values were combinations of features that
# plot_tales_msa() exposes as orthogonal arguments; the one thing it still did
# that plot_tales_msa() could not -- draw the consensus row -- was implemented
# there, which is what unblocked this.
#
# What is lost, and why it did not justify keeping 285 lines:
#   save_path     a returned ggplot is ggsave()-able
#   note_colors   a ggplot caller adds a scale instead
#   ...           passthrough to gplots::heatmap.2, backend-specific
#
# Kept rather than deleted because it is the only record of the heatmap.2
# layout, should anyone want to compare renderings.
#
# Moved out of R/ on 2026-09-14. See dev/restructuring-notes.md section 4.

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
#'  \code{tal_sim}, \code{forCellNote} - repeat/rvd alignment, and ref_pattern optionally.
#'  
#'  But if you don't have these alignments, you can provide \strong{repeat alignment}
#'  to \code{forMatrix} with \code{repeat_sim}, the repeat similarity data frame,
#'  the function will convert it into repeatClustID/Similarity alignment depending
#'  on the plot type. In case of \emph{repeat.clusters}, you may want to adjust
#'  the param \code{h_cut} to decrease/increase the number of repeat clusters.
#' 
#' 
#' 
#' 
#' @param tal_sim a \emph{three columns Tals similarity table} as obtained
#'  with \code{\link{tales_compare}} in the \code{tale_distances} element of the returned object.
#' @param repeat_align a multiple Tal repeat sequences alignment in the
#'  form of a matrix as returned by \code{\link{tales_align}}.
#' @param repeat_sim A long, three columns data frame with pairwise similarity
#' scores between repeats as available in the \code{domain_distances} element
#' of the object returned by the \code{\link{tales_compare}} function.
#' @param plot_type Either \code{"repeat.similarity"}, \code{"repeat.clusters"} ,
#'  \code{"repeat.clusters.with.rvd"}. Defines the type of plot that will be produced
#'   by the function. See below for details.
#' @param h_cut height for tree cutting when plot type
#'  in "repeat.clusters".
#' @param rvd_align (optional) when the rvds need to be labeled in the
#'  plot (plot_type = "repeat.similarity" or "repeat.clusters.with.rvd",
#'  a multiple Tal repeat sequences alignment in the form of a matrix as
#'  returned by \code{\link{tales_align}}.
#' @param ref_pattern regular expression pattern that will be used to search Tal names
#' to select the reference in the alignment.
#' @param consensus (logical) whether to display the consensus sequence when 
#' the plot type is "repeat.clusters.with.rvd".
#' @param note_colors In case rvdSim = NULL, vector of 2 colors for rvd alignment,
#' the first color is for matched rvds, and the second color is for mismatched ones.
#' In the other case, more colors should be supplied.
#' @param save_path file path to save the plot. If save_path is NULL, the heatmap 
#' will be printed. If save_path is specified, the image file will be created with
#' the format based on file extension.
#' @param ... any other arguments of \code{\link[gplots:heatmap.2]{heatmap.2}}
#' 
#' @return the return value of \code{\link[gplots:heatmap.2]{heatmap.2}}
#' 
#' @export
#' @family TALE plots
msa_heatmap <- function(tal_sim, repeat_align, rvd_align = NULL,
                        repeat_sim, h_cut = 10, ref_pattern = NULL,
                        consensus = FALSE, note_colors = NULL,
                        plot_type, save_path, ...) {
  
  
  if (startsWith(plot_type, "repeat.clusters")) {
    forMatrix <- .repeat_to_cluster_align(repeat_align = repeat_align, repeat_sim = repeat_sim, h_cut = h_cut)
  } else if (plot_type == "repeat.similarity") {
    forMatrix <- .repeat_to_sim_align(repeat_align = repeat_align, repeat_sim = repeat_sim, ref_tag = ref_pattern)
  } else if (!hasArg(plot_type) || is.null(plot_type) || is.na(plot_type)) {
    stop("Missing plot_type")
  } else{
    stop(glue::glue("\"{plot_type}\" plot is not available."))
  }
  
  
  if (plot_type == "repeat.clusters") {
    forCellNote <- repeat_align
  } else {
    if (is.null(rvd_align)) stop("\"{plot_type}\" plot requires rvd_align!")
    forCellNote <- rvd_align
  }
  
  forCellNote <- forCellNote[rownames(forMatrix),]
  
  
  
  # for 'Rowv' = dend
  
  ###!!! THIS FAILS IF repeat_align has a single sequence
  
  tal_sim <- tal_sim[tal_sim$TAL1 %in% rownames(forCellNote), ]
  tal_sim <- tal_sim[tal_sim$TAL2 %in% rownames(forCellNote), ]
  tal_sim <- as.matrix(reshape2::acast(tal_sim, TAL1 ~ TAL2, value.var = "Sim")) # melt then unmelt ...
  taldist <- 100 -tal_sim
  taldist <- taldist[rownames(forCellNote), ]
  taldist <- taldist[, rownames(forCellNote)]
  
  clust <- hclust(as.dist(taldist))
  forRowv <- as.dendrogram(clust)
  
  
  
  if (startsWith(plot_type, "repeat.clusters")) { # in case of repcode plotting
    # define 'col'
    # forCol <- function(x) scales::hue_pal(l = 55)(n=100)[0:x]
    forCol <- function(x) viridis::inferno(n=100, end = .9)[0:x]
    forBreaks <- 0:max(forMatrix, na.rm = T)
    
    # define 'key'
    forKeyxlab <- NA
    forKey <- FALSE
    forNoteCex <- 1
    
    # define 'notecol'
    if (is.null(note_colors)) {
      note_colors <- list(matched = "white", mismatched = "#01FFFF")
    } else {
      note_colors <- list(matched = note_colors[1], mismatched = note_colors[2])
    }
    if (endsWith(plot_type, "with.rvd")) {
      # consensus rvds
      # rvdsAlignedStrings <- apply(forCellNote, 1, function(x) paste(x, collapse = "-")) %>% BStringSet()
      # consensusRVD <- Biostrings::consensusString(rvdsAlignedStrings, ambiguityMap = "+")
      
      consensusRVD <- sapply(1:ncol(forCellNote), function(x) {
        allRVDs <- forCellNote[,x]
        freq <- sapply(unique(allRVDs), function(p) S4Vectors::countMatches(p, allRVDs))
        unique(allRVDs)[which.max(freq)]
      })
      
      # notecol
      rvdcol <- forCellNote
      # if (is.null(ref_pattern)) {
      #   reftalID <- which.max(apply(forCellNote, 1, function(x) length(x[!is.na(x)])))
      # } else {
      #   reftalID <- which(grepl(ref_pattern, rownames(forCellNote)))
      # }
      # rownames(forMatrix)[reftalID] <- paste0(rownames(forMatrix)[reftalID], "_")
      # for (k in 1:ncol(rvdcol)){
      #   rvd <- as.character(forCellNote[reftalID, k])
      #   if (is.na(rvd)) {
      #     rvdcol[,k] <- note_colors$mismatched
      #   } else {
      #     rvdcol[,k] <- ifelse(toupper(rvdcol[,k]) == toupper(rvd), note_colors$matched, note_colors$mismatched)
      #   }
      # }
      for (k in 1:ncol(rvdcol)){
        rvd <- consensusRVD[k]
        if (is.na(rvd)) {
          rvdcol[,k] <- note_colors$mismatched
        } else {
          rvdcol[,k] <- ifelse(toupper(rvdcol[,k]) == toupper(rvd), note_colors$matched, note_colors$mismatched)
        }
      }
      rvdcol <- rvdcol[order.dendrogram(forRowv), ]
      forNoteCol <- t(as.matrix(rvdcol))
      
      
      
    } else if (endsWith(plot_type, "repeat.clusters")) {
      forNoteCol <- "white"
    } else {
      stop(glue::glue("\"{plot_type}\" plot is not available."))
    }
  }
  else if (plot_type == "repeat.similarity") { # in case of rvd plotting
    # define 'col'
    # forCol <- colorRampPalette(c("dodgerblue4", "dodgerblue3", "dodgerblue", "deepskyblue", "white"))
    forCol <- colorRampPalette(c("#421727", "#6e2742", "#9a365c", "#b03e69", "#ffffff"))
    forBreaks <- 0:100
    
    rvdSim <- NULL
    # define 'notecol'
    if (!is.null(rvdSim) && is.matrix(rvdSim) && is.numeric(rvdSim)) { # in case rvd similarity matrix is provided
      if (is.null(note_colors)) {
        note_colors <- viridis::viridis(n=200, direction = -1)
      }
      col_range <- function(x) note_colors[as.integer(x)]
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
      if (is.null(note_colors)) {
        note_colors <- list(matched = "black", mismatched = "red")
      } else {
        note_colors <- list(matched = note_colors[1], mismatched = note_colors[2])
      }
      
      if (is.null(ref_pattern)) {
        sim_len <- apply(forMatrix, 1, function(x) {length(grep("100", x))})
        reftal <- names(which.max(sim_len))
        reftalID <- which(rownames(forCellNote) == reftal)
      } else {
        reftalID <- which(grepl(ref_pattern, rownames(forCellNote)))
      }
      rownames(forMatrix)[reftalID] <- paste0(rownames(forMatrix)[reftalID], "_#")
      rvdcol <- forCellNote
      for (k in 1:ncol(rvdcol)){
        rvd <- as.character(forCellNote[reftalID, k])
        if (is.na(rvd)) {
          rvdcol[,k] <- note_colors$mismatched
        } else {
          rvdcol[,k] <- ifelse(toupper(rvdcol[,k]) == toupper(rvd), note_colors$matched, note_colors$mismatched)
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
  extra.key <- function(x = NULL, check_plot_type = plot_type, hei = hei_top * 2.54, wid = wid_left * 2.54, col_set = note_colors) {
    if (check_plot_type == "repeat.similarity" || endsWith(check_plot_type, "with.rvd")) {
      if (!is.null(x) && is.matrix(x) && is.numeric(x)) {
        par(mai = c(hei*.4, 0, hei*.2, wid*.1), mgp = c(2, 1, 0))
        image(z = matrix(seq(-1, 1, by = .01), ncol = 1), col = col_set, yaxt = "n", xaxt = "n", xlab = "RVD similarity")
        axis(1, at = seq(0, 1, by = .25), labels = c("-1", NA, "0", NA, "1"))
      } else {
        par(mai = c(hei*.2, wid*.25, hei*.2,  wid*.25), mgp = c(2, 1, 0))
        image(z = matrix(c(0, 1), ncol = 2), col = "grey50", yaxt = "n", xaxt = "n")
        abline(h = 0.5, col = "grey", lwd = 1.5)
        text(0, 1, labels = ifelse(check_plot_type == "repeat.similarity", "reference", "consensus"), col = col_set$matched, font = 2)
        text(0, 0, labels = "other", col = col_set$mismatched, font = 2)
        mtext(side = 1, at = 0, text = "RVD alignment", cex = .75, col = "black", padj = 0.5)
        if (check_plot_type == "repeat.clusters.with.rvd" && consensus) {
          par(mar = c(0,0,0,0))
          image(z = matrix(1:length(consensusRVD), ncol = 1), col = "grey50", bg = "grey", yaxt = "n", xaxt = "n")
          for (i in 1:length(consensusRVD)) {
            abline(v = (i-.5)/(length(consensusRVD)-1), col = "grey")
            text((i-1)/(length(consensusRVD)-1), 0, labels = consensusRVD[i], font = 2, col = col_set$matched, cex = 1.2)
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
  
  # save plot to "save_path" if provided, the format of image depends on extension of "save_path"
  if (hasArg(save_path)) {
    img_format <- gsub(".*\\.", "", basename(save_path))
    img_size <-  list(save_path, width = (wid_left + wid_right + wid_left) * 2, height = (hei_top + hei_bottom) * 2)
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
                                    ifelse(consensus, hei_bottom/nrow(forMatrix), 0), hei_bottom), ##
                           lwid = c(wid_left, wid_right, wid_left), ##
                           
                           # trace
                           trace = "none")
  
  custom_arg_list <- as.list(substitute(list(...)))[-1L]
  
  default_arg_list <- default_arg_list[!names(default_arg_list) %in% names(custom_arg_list)]
  
  heatmap_plot <- do.call(gplots::heatmap.2, c(default_arg_list, custom_arg_list))
  
  if (hasArg(save_path)) {
    dev.off()
  }
  return(invisible(heatmap_plot))
}
