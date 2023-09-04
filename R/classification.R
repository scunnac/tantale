
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
    if (is.null(k)) message("WE SHOULD BE DOING SOMETHING")
    if (!is.numeric(numGroups) || length(numGroups) != 1) stop("'k' must be specified as a number!")
    # if (method == "euclidean") {
    #   taleTree <- hclust(d = dist(distMat, method = "euclidean"), method = "ward.D")
    # } else if (method == "distal") {
    #   taleTree <- hclust(as.dist(distMat), method = "ward.D")
    # }
    
    taleTree <- hclust(d = dist(distMat, method = "euclidean"), method = "ward.D")
    
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
      ggtree::geom_tiplab(ggtree::aes(label = label),
                          hjust = 1,
                          angle = 90,
                          align = FALSE,
                          color='black',
                          offset = -2,
                          
      ) +
      # ggplot2::scale_color_brewer("Groups", palette="BrBG") + # allowed maximum for palette BrBG is 11
      viridis::scale_color_viridis(discrete = T, option = "C", breaks = 1:numGroups) +
      ggplot2::geom_vline(xintercept = -(cutOff/2), linetype = 2) +
      ggtree::geom_text(x = (cutOff/2 - max(taleTree$height)/50),
                        y = 8, label = paste("half of 'cutOff' value: ", sprintf("%.2f", cutOff/2)),
                        color = "darkgrey", fontface = "plain") +
      ggplot2::xlab("Height/2") +
      ggtree::theme_dendrogram(plot.margin = ggplot2::margin(6,6,150,6))
  }
  
  
  
  
  
  if(plotTree == TRUE) {print(p)}
  return(taleGroups)
  
  
}






#' Heatmap plotting of rvd sequence variants
#' @description The function creates a graphical presentation from a tale annotation table. The output is like a heatmap that presents rvd sequence variants in Tal groups as column and respective strains as rows (or vice versa). It is different from a typical heatmap that it can display more than one value in a cell; for example, if one strain has 2 rvd sequence variants belong to 1 group, it will be displayed by 2 colors in 1 cell.
#' @param tale_annotation a data frame containing at least 3 columns for Tal groups, strain names, and rvd seqs, and 1 row is 1 Tal.
#' @param col "character", column name of \code{tale_annotation} to be displayed as rows in the heatmap (e.g. tal groups).
#' @param row "character", column name of \code{tale_annotation} to be displayed as columns in the heatmap (e.g. strain names).
#' @param value "character", column name for rvdseqs in the \code{tale_annotation}
#' @param truncTaleLab (optional, default = NULL) "character", column name of \code{tale_annotation} labeling the truncTales by TRUE/FALSE value. The truncTales are labeled by "T" in the heatmap cells, but if this argument is called.
#' @param extraCol (optional, default = NULL) "character", column name of \code{tale_annotation} containing other information (e.g. origin). It will be presented in a side bar on the right of the heatmap.
#' @param x.lab,y.lab,title character for x axix, y axis names and title 
#' @param mapcol character vector of colors for the cells.
#' @param mar.side margin of the heatmap for row dendrogram, col dendrogram, rownames, colnames, respectively. (by default, c(5, 5, 3, 3)).
#' @param sepwid numeric value for the width of separator between adjacent cells.
#' @param sepcol character of color for the separator between adjacent cells
#' @param inner_sepcol character of color for the separator between colors within 1 cell if there are more than 1.
#' @param save.path (optional) file path to save the plot, format of the image depends on the file extension. If save.path is NULL, the heatmap will be printed. If save.path is specified, the image file will be created.
#' @export
heatmap_talomes <- function(tale_annotation, col, row, value, truncTaleLab = NULL, extraCol = NULL,
                            x.lab = "TALE Group", y.lab = "Strain", title = "RVD sequences variants",
                            plot.type = "all",
                            mapcol = viridis::viridis(10), mar.side = c(5, 5, 3, 3),
                            sepwid = 5, sepcol = "white", inner_sepcol = "white", save.path = NULL) {
  
  ## rename colnames of tale annotation
  colnames(tale_annotation)[which(colnames(tale_annotation) == col)] <- "group"
  colnames(tale_annotation)[which(colnames(tale_annotation) == row)] <- "strain"
  colnames(tale_annotation)[which(colnames(tale_annotation) == value)] <- "rvdseq"
  if (!is.null(extraCol)) colnames(tale_annotation)[which(colnames(tale_annotation) == extraCol)] <- "extraCol"
  if (!is.null(truncTaleLab)) colnames(tale_annotation)[which(colnames(tale_annotation) == truncTaleLab)] <- "truncTale"
  
  tale_annotation %<>% dplyr::mutate(aberrantRepeat = ifelse(grepl("[a-z]", rvdseq), TRUE, FALSE))
  tale_annotation$rvdseq <- toupper(tale_annotation$rvdseq)
  
  if(is.numeric(tale_annotation$group)) tale_annotation$group <- paste0("G", tale_annotation$group)
  
  ## convert rvdseq to rvdfac, ranking of abundance
  tale_annotation %<>% dplyr::group_by(group) %>% dplyr::mutate(rvdfac = do.call((function(x) {
    levels(x) <- nlevels(x) + 1 - rank(table(x), ties.method = "first")
    return(as.integer(as.character(x)))
  }), list(as.factor(rvdseq))))
  
  ## get number of rvdseqs per group per strain
  ## it is used for the layout of heatmap cell
  numAlleles <- reshape2::dcast(tale_annotation, strain ~ group, value.var = "rvdseq", drop = F, fun.aggregate = function(x) length(x))
  
  rownames(numAlleles) <- numAlleles$strain
  numAlleles <- numAlleles[,-1]
  
  ## count number of alleles
  ## for column label
  variantsCount <- sapply(sort(unique(tale_annotation$group)), function(g) {
    length(unique(tale_annotation[tale_annotation$group == g,]$rvdseq))
  }, USE.NAMES = T)
  colnames(numAlleles) <- paste0(colnames(numAlleles), " #", variantsCount)
  
  tale_annotation %<>% dplyr::group_by(group, strain) %>% dplyr::mutate(reprsntRVDfac = ifelse(1 %in% rvdfac, 1, rvdfac[1]))
  
  ## representative alleles per group per strain
  ## in case a strain has more than 1 rvd seqs in 1 group
  ## find the most abundant rvdseqin that group
  ## if this strain has that rvdseq, take it as representative rvdseq
  ## if not, take randomly 1 rvdseq that strain has
  reprsntAlleles <- reshape2::dcast(tale_annotation, strain ~ group, value.var = "reprsntRVDfac", drop = F, fun.aggregate = function(x) as.integer(x[1]))
  reprsntAlleles <- as.data.frame(lapply(X = reprsntAlleles[,-1], FUN = as.integer))
  rownames(reprsntAlleles) <- rownames(numAlleles)
  colnames(reprsntAlleles) <- colnames(numAlleles)
  reprsntAlleles <- apply(reprsntAlleles, 2, function(x) { # allele '0' if missing
    x[is.na(x)] <- 0
    return(x)})
  
  ## compute dendrograms with the representative alleles
  Strain_HC <- hclust(ape::dist.gene(x = reprsntAlleles), method = "average")
  
  strainAlleles <- apply(reprsntAlleles, 1, as.integer)
  rownames(strainAlleles) <- colnames(numAlleles)
  TALE_HC <- hclust(ape::dist.gene(x = strainAlleles), method = "average")
  
  
  df <- par(no.readonly = T)
  ## plot sizes
  widleft <- mar.side[1]
  heitop <- mar.side[2]
  widright <- mar.side[3]
  heibot <- mar.side[4]
  
  if (plot.type == "single") { # plot representative alleles
    codedAlleles1 <- apply(reprsntAlleles, 2, function(x) ifelse(x == 0, NA, x))
    if (hasArg(save.path)) {
      img_format <- gsub(".*\\.", "", basename(save.path))
      img_size <-  list(save.path, width = (widleft + ncol(codedAlleles1) + widright)/2.54, height = (heitop + nrow(codedAlleles1) + heibot)/2.54)
      if (img_format %in% c("bmp", "jpeg", "png", "tiff")) {
        img_size <- c(img_size, units = "in", res = 1440)
      }
      do.call(img_format, img_size)
    }
    gplots::heatmap.2(as.matrix(codedAlleles1),
                      trace = "none",
                      col = mapcol[1:max(codedAlleles1, na.rm = T)],
                      breaks = 0:max(codedAlleles1, na.rm = T),
                      density.info = "none",
                      key = F,
                      xlab = x.lab,
                      ylab = y.lab,
                      margins = c(7, 7),
                      colsep = 0:(ncol(codedAlleles1)-0),
                      rowsep = 0:(nrow(codedAlleles1)-0),
                      sepcolor = sepcol,
                      sepwidth = rep(sepwid/100, 2),
                      Rowv = as.dendrogram(Strain_HC),
                      Colv = as.dendrogram(TALE_HC),
                      main = title,
                      na.color = "grey",
                      lmat = rbind(c(4, 3), c(2, 1)),
                      lhei = c(heitop, nrow(codedAlleles1) + heibot), ##
                      lwid = c(widleft, ncol(codedAlleles1) + widright)
    )
    
  } else if (plot.type == "all") { # plot all alleles
    uniqueRVD <- numAlleles
    rorder <- order.dendrogram(as.dendrogram(Strain_HC))
    corder <- order.dendrogram(as.dendrogram(TALE_HC))
    uniqueRVD <- uniqueRVD[rev(rorder), corder]
    
    
    colmat <- mapcol
    ## plot layout
    nplots <- nrow(uniqueRVD)*ncol(uniqueRVD)
    mainmat <- matrix(1:nplots, nrow = nrow(uniqueRVD))
    left.mat <- matrix(rep(nplots+1, nrow(uniqueRVD)), ncol = 1)
    top.mat <- matrix(c(0, rep(nplots+2, ncol(uniqueRVD))), nrow = 1, byrow = T)
    # right.mat <- matrix(rep(c(rep(0, heitop), (nplots+3):(nrow(uniqueRVD)+nplots+2)),2), ncol = 2, byrow = F)
    # bottom.mat <- matrix(c(rep(0, widleft), (nrow(uniqueRVD)+nplots+3):(nrow(uniqueRVD)+nplots+2+ncol(uniqueRVD)), 0, 0), nrow = 1)
    extcol <- matrix(c(0, rep(nplots+3, nrow(uniqueRVD))), ncol = 1, byrow = F)
    right.mat <- matrix(c(0, rep(nplots+4, nrow(uniqueRVD))), ncol = 1, byrow = F)
    legend.mat <- matrix(c(0, rep(nplots+5, nrow(uniqueRVD))), ncol = 1, byrow = F)
    bottom.mat <- matrix(c(0, rep(nplots+6, ncol(uniqueRVD)), 0, 0, 0), nrow = 1, byrow = T)
    laymat <- rbind(top.mat, cbind(left.mat, mainmat))
    laymat <- rbind(cbind(laymat, extcol, right.mat, legend.mat), bottom.mat)
    titmat <- matrix(c(0, rep(nplots+7, ncol(mainmat)), 0, 0, 0), nrow = 1, byrow = T)
    laymat <- rbind(titmat, laymat)
    on.exit(par(no.readonly = TRUE))
    
    if (hasArg(save.path)) {
      img_format <- gsub(".*\\.", "", basename(save.path))
      img_size <-  list(save.path, width = sum(widleft, rep(1, ncol(uniqueRVD)), widright)/2.54, height = sum(2, heitop, rep(1, nrow(uniqueRVD)), heibot)/2.54)
      if (img_format %in% c("bmp", "jpeg", "png", "tiff")) {
        img_size <- c(img_size, units = "in", res = 1440)
      }
      do.call(img_format, img_size)
    }
    
    layout(laymat, widths = c(widleft, rep(1, ncol(uniqueRVD)), ifelse(is.null(extraCol), 0.1, .7), widright, ifelse(is.null(extraCol), 0.1, 3)), heights = c(2, heitop, rep(1, nrow(uniqueRVD)), heibot))
    # layout.show(nplots+7)
    
    
    ## plot heatmap
    for (c in 1:ncol(uniqueRVD)) {
      gname <- gsub(" \\#\\d+", "", colnames(uniqueRVD)[c])
      # gname <- gsub("G", "", gname)
      g1 <- tale_annotation[tale_annotation$group == gname,]
      # g1$rvdseq <- as.integer(as.factor(g1$rvdseq))
      for (r in 1:nrow(uniqueRVD)) {
        par(mar = rep(0, 4))
        nelements <- uniqueRVD[r, c]
        if (nelements == 0) {
          image(z = matrix(0), col = "grey", axes = F)
          if (sepwid > 0) box(lwd = sepwid/2, col = sepcol)
        } else {
          g1s1 <- g1[g1$strain == rownames(uniqueRVD[r,]),]
          rvd.factor <- g1s1$rvdfac
          if (is.null(truncTaleLab)) {
            truncTale <- NA
          } else {
            truncTale <- sapply(g1s1$truncTale, function(l) ifelse(l, "T", NA))
          }
          # plot.index <- r + (c-1)*nrow(uniqueRVD)
          color.elements <- colmat[rvd.factor]
          par(mar = rep(0, 4))
          image(z = matrix(1:(nelements), ncol = 1), col = color.elements, axes = F)
          text(seq(0,1, length.out = nelements), 0, labels = truncTale, cex = 1, font = 2, col = "black")
          if (nelements > 1) {
            abline(v = seq(0.5/(nelements-1), 1-.5/(nelements-1), length.out = nelements -1), col = inner_sepcol, lty = 1)
          }
          if (sepwid > 0) {
            par(mar = rep(0, 4))
            box(lwd = sepwid/2, col = sepcol)}
        }
        
      }
    }
    
    ## dendrograms
    par(mai = rep(0, 4))
    plot(as.dendrogram(Strain_HC), horiz = TRUE, axes = FALSE, leaflab = "none", yaxs = "i")
    par(mai = rep(0, 4))
    plot(as.dendrogram(TALE_HC), axes = FALSE, leaflab = "none", xaxs = "i")
    
    if (!is.null(extraCol)) {
      ## extra column
      # extra.bar <- sample(c("Hanoi", "Hatay", "Namdinh", NA), nrow(uniqueRVD), replace = T)
      extra.bar <- sapply(rownames(uniqueRVD), function(s) {
        ext <- unique(tale_annotation$extraCol[tale_annotation$strain == s])
        return(ext)
      })
      rextra.bar <- unique(extra.bar[!is.na(extra.bar)])
      rextra.bar <- data.frame("lab" = sort(rextra.bar), "fac" = 1:length(rextra.bar))
      extra.col <- viridis::viridis(n = nrow(rextra.bar))[sapply(extra.bar, function(e) ifelse(is.na(e), NA, rextra.bar$fac[rextra.bar$lab == e]), simplify = T)]
      extra.col[is.na(extra.col)] <- "gray"
      par(mar = c(0, 0, 0, 0))
      image(z = matrix(1:nrow(uniqueRVD), nrow = 1), col = rev(extra.col), yaxt = "n", xaxt = "n", axes = F)
    } else {
      par(mar = c(0, 0, 0, 0))
      image(z = matrix(1:nrow(uniqueRVD), nrow = 1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    }
    
    
    ## rownames
    rownames(uniqueRVD) <- paste0(rownames(uniqueRVD),"  #", rowSums(uniqueRVD, na.rm = T))
    par(mar = c(0, 0.5, 0, 0))
    image(z = matrix(1:nrow(uniqueRVD), nrow = 1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    text(-1, seq(0, 1, length.out = nrow(uniqueRVD)), labels = rev(rownames(uniqueRVD)), font = 1, col = "black", cex = 1.2, adj = 0)
    mtext(side = 4, at = .5, text = "Strain", col = "black", padj = 0, line = -1)
    
    if (!is.null(extraCol)) {
      ## legend column
      par(mar = c(0,0.5,1,0))
      plot(rep(0, nrow(rextra.bar)), -seq(from = 0, by = .8, length.out = nrow(rextra.bar)), type = "p", pch = 15, col = viridis::viridis(n = nrow(rextra.bar)), axes = F, main = extraCol, xlab = NA, ylab = NA, cex = 4, ylim = c(-nrow(uniqueRVD), 0), xlim = c(0,2))
      text(rep(.3, nrow(rextra.bar)), -seq(from = 0, by = .8, length.out = nrow(rextra.bar)), labels = rextra.bar$lab, font = 1, col = "black", bg = "red", cex = 1.2, adj = 0)
    } else {
      par(mar = c(0,0,0,0))
      image(matrix(0), col = "white", axes = F)
    }
    
    
    ## colnames
    par(mar = c(0.5, 0, 0.5, 0))
    image(z = matrix(1:ncol(uniqueRVD), ncol = 1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    text(seq(0, 1, length.out = ncol(uniqueRVD)), 1, labels = colnames(uniqueRVD), font = 1, col = "black", bg = "red", cex = 1.2, srt = 90, adj = 1)
    mtext(side = 1, at = 0.5, text = "TALE Group", col = "black", padj = 0, line = -1)
    
    
    ## title
    par(xpd = T, mai = rep(0, 4))
    image(z = matrix(1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    text(0, 1/3, labels = "RVD sequences variants", font = 2, col = "black", cex = 2, pos = 1)
    
    
    par(df)
  }
  if (hasArg(save.path)) {
    dev.off()
  }
}



