# Draw the site's two maps of the package, from the recording of which
# object each exported function takes and returns
# (dev/function-graph-dataflow.tsv, written by
# dev/function-graph-dataflow.R):
#
# - pkgdown/assets/function_map.html: every exported function, interactive.
#   The reference index shows it in its first section (_pkgdown.yml);
#   dev/function-graph.qmd draws the same widget.
# - pkgdown/assets/function_map.svg: the main functions only, on the home
#   page (pkgdown/index.md).
#
# pkgdown copies pkgdown/assets/ to the site root. Rerun after rerunning the
# recorder, from the package root:
#
#   Rscript pkgdown/function_map.R
#
# Needs visNetwork, htmlwidgets (with pandoc, for a self-contained file) and
# ggplot2. dev/function-graph.qmd sources this file for its functions only.

# The recording, and what NAMESPACE exports now. Links between classes and
# functions; constructors share their class's name (`tales()` makes a
# `tales`), so node ids carry a prefix.
function_map_data <- function(root = ".") {
  tsv <- file.path(root, "dev", "function-graph-dataflow.tsv")
  flow <- read.delim(tsv, comment.char = "#", colClasses = "character",
                     na.strings = "NA")
  recorded_on <- sub("^# ", "", readLines(tsv, n = 1))
  ns_now <- readLines(file.path(root, "NAMESPACE"))
  exported_now <- sub("^export\\((.*)\\)$", "\\1",
                      grep("^export\\(", ns_now, value = TRUE))
  s3_now <- sub("^S3method\\((.*),(.*)\\)$", "\\1.\\2",
                grep("^S3method\\(", ns_now, value = TRUE))
  s3_now <- gsub("\"|^[a-z]+::", "", s3_now)

  package_classes <- c("tales", "tales_msa", "pairwise_distances",
                       "tale_distances", "domain_distances")
  bioc_classes <- c("BStringSet", "DNAStringSet", "AAStringSet")
  drawn_classes <- c(package_classes, bioc_classes)

  left_out <- "^is_|^validate_|_assert_|^print\\.|^format\\.|^\\[\\.|^dplyr_|^as_tales\\."
  shown <- flow[!grepl(left_out, flow$fn) & flow$fn %in% c(exported_now, s3_now), ]

  cls_id <- function(x) paste0("class:", x)
  fn_id <- function(x) paste0("fn:", x)
  split_pairs <- function(x) {
    x <- x[nzchar(x)]
    if (!length(x)) return(matrix(character(), 0, 2))
    do.call(rbind, strsplit(unlist(strsplit(x, ";")), "="))
  }
  links <- list()
  for (i in seq_len(nrow(shown))) {
    r <- shown[i, ]
    if (!is.na(r$arg_class) && r$arg_class %in% drawn_classes)
      links[[length(links) + 1]] <- c(cls_id(r$arg_class), fn_id(r$fn), "")
    for (p in seq_len(nrow(op <- split_pairs(r$other_args))))
      if (op[p, 2] %in% drawn_classes)
        links[[length(links) + 1]] <- c(cls_id(op[p, 2]), fn_id(r$fn), op[p, 1])
    if (r$return_class %in% drawn_classes)
      links[[length(links) + 1]] <- c(fn_id(r$fn), cls_id(r$return_class), "")
    for (p in seq_len(nrow(rp <- split_pairs(r$return_parts))))
      if (rp[p, 2] %in% drawn_classes)
        links[[length(links) + 1]] <- c(fn_id(r$fn), cls_id(rp[p, 2]), paste0("$", rp[p, 1]))
  }
  links <- unique(as.data.frame(do.call(rbind, links)))
  names(links) <- c("from", "to", "label")

  # The one link through files on disk.
  stopifnot(all(c("tell_tales", "tales_from_telltales") %in% exported_now))
  links <- rbind(links, data.frame(from = fn_id("tell_tales"),
                                   to = fn_id("tales_from_telltales"),
                                   label = "results folder"))

  ids <- unique(c(links$from, links$to))
  fns <- sort(sub("^fn:", "", grep("^fn:", ids, value = TRUE)))
  io_summary <- function(f) {
    r <- flow[flow$fn == f, ]
    paste(sprintf("%s &rarr; %s", ifelse(is.na(r$arg_class), "(no argument)", r$arg_class),
                  r$return_class), collapse = "<br>")
  }
  fn_nodes <- data.frame(
    id = fn_id(fns), label = paste0(fns, "()"),
    group = ifelse(fns %in% exported_now, "exported function", "S3 method"),
    shape = "dot", size = 12,
    title = vapply(fns, function(f) sprintf("<b>%s</b><br>%s", f, io_summary(f)), "")
  )
  class_nodes <- data.frame(
    id = cls_id(drawn_classes), label = drawn_classes,
    group = ifelse(drawn_classes %in% package_classes, "tantale class", "Biostrings class"),
    shape = "box", size = 20,
    title = sprintf("<b>%s</b>", drawn_classes)
  )
  class_nodes <- class_nodes[class_nodes$id %in% ids, ]

  list(flow = flow, recorded_on = recorded_on, exported_now = exported_now,
       s3_now = s3_now, nodes = rbind(class_nodes, fn_nodes),
       links = links)
}

# Every exported function, interactive, in the colours of R/palette.R
# (Tol muted), as in the figure below.
function_map_widget <- function(d, height = "80vh") {
  edges <- data.frame(
    from = d$links$from, to = d$links$to, label = d$links$label, arrows = "to",
    dashes = d$links$label == "results folder",
    color = "#7a7a7a", font.size = 11, font.color = "#555555"
  )
  visNetwork::visNetwork(d$nodes, edges, width = "100%", height = height) |>
    visNetwork::visGroups(groupname = "tantale class", color = "#332288",
                          font = list(color = "white", size = 16)) |>
    visNetwork::visGroups(groupname = "Biostrings class", color = "#88CCEE",
                          font = list(color = "black", size = 14)) |>
    visNetwork::visGroups(groupname = "exported function", color = "#117733") |>
    visNetwork::visGroups(groupname = "S3 method", color = "#DDCC77") |>
    visNetwork::visEdges(smooth = FALSE, arrows = list(to = list(scaleFactor = 0.5))) |>
    visNetwork::visPhysics(solver = "forceAtlas2Based",
                           forceAtlas2Based = list(gravitationalConstant = -120, springLength = 140),
                           stabilization = list(iterations = 600)) |>
    # Once laid out, the graph stops simulating, so an open page does no
    # work in the background.
    visNetwork::visEvents(stabilizationIterationsDone =
                            "function() { this.setOptions({physics: false}); }") |>
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1,
                                                   algorithm = "hierarchical"),
                           nodesIdSelection = list(enabled = TRUE, main = "Select a node")) |>
    visNetwork::visInteraction(navigationButtons = TRUE, tooltipDelay = 150) |>
    visNetwork::visLegend(useGroups = FALSE, position = "left", width = 0.12,
                          addNodes = data.frame(
                            label = c("tantale\nclass", "Biostrings\nclass", "exported", "S3\nmethod"),
                            color = c("#332288", "#88CCEE", "#117733", "#DDCC77"),
                            shape = c("box", "box", "dot", "dot"), font.size = 14,
                            font.color = c("white", "black", "black", "black")))
}

# The main functions only, top to bottom in the order of an analysis. Each
# node's place is set by hand (x in columns, y in rows), and so are the
# links, which keep only the main path; each must be in the recording. The
# grouping functions return the tales object they were given, with a group
# column; it is drawn a second time below them, so the figure reads
# downwards.
function_map_figure <- function(d) {
  place <- read.table(header = TRUE, text = "
    id                          x     y
    fn:tell_tales               0     0
    fn:tales_from_telltales     0     1
    class:tales                 0     2
    fn:plot.tales              -1.9   3
    fn:tales_anomalies         -1.3   3
    fn:tales_domain_distances  -0.47  3
    class:domain_distances     -0.47  4
    fn:tales_tale_distances    -0.47  5
    fn:tales_compare_distal     0.39  5
    fn:tales_compare_functal    1.27  5
    class:tale_distances        0     6
    fn:tales_group_hclust      -0.45  7
    fn:tales_group_kmedoids     0.45  7
    class:tales_grouped         0     8
    fn:tales_align             -0.66  9
    fn:talomes_heatmap          0     9
    fn:tales_predict_targets    0.81  9
    class:tales_msa            -0.66 10
    fn:plot.tales_msa          -0.66 11")
  e <- read.table(header = TRUE, text = "
    from                        to
    fn:tell_tales               fn:tales_from_telltales
    fn:tales_from_telltales     class:tales
    class:tales                 fn:plot.tales
    class:tales                 fn:tales_anomalies
    class:tales                 fn:tales_domain_distances
    class:tales                 fn:tales_compare_distal
    class:tales                 fn:tales_compare_functal
    fn:tales_domain_distances   class:domain_distances
    class:domain_distances      fn:tales_tale_distances
    fn:tales_tale_distances     class:tale_distances
    fn:tales_compare_distal     class:tale_distances
    fn:tales_compare_functal    class:tale_distances
    class:tale_distances        fn:tales_group_hclust
    class:tale_distances        fn:tales_group_kmedoids
    fn:tales_group_hclust       class:tales_grouped
    fn:tales_group_kmedoids     class:tales_grouped
    class:tales_grouped         fn:tales_align
    class:tales_grouped         fn:talomes_heatmap
    class:tales_grouped         fn:tales_predict_targets
    fn:tales_align              class:tales_msa
    class:tales_msa             fn:plot.tales_msa")
  recorded <- paste(d$links$from, d$links$to)
  as_recorded <- paste(sub("tales_grouped", "tales", e$from),
                       sub("tales_grouped", "tales", e$to))
  # Its tests give tales_group_kmedoids() the similarity table that
  # tales_compare_distal() returns; the articles give it a tale_distances.
  unrecorded <- "class:tale_distances fn:tales_group_kmedoids"
  stopifnot(all(as_recorded %in% c(recorded, unrecorded)),
            all(place$id %in% c(e$from, e$to)))
  e$label <- ifelse(e$from == "fn:tell_tales", "results folder", "")

  n <- merge(place, rbind(d$nodes, data.frame(
    id = "class:tales_grouped", label = "tales", group = "tantale class",
    shape = "box", size = 20, title = "")), by = "id")
  n$is_fn <- startsWith(n$id, "fn:")

  # Sizes in points: columns 150 pt apart, rows 34 pt; box widths from the
  # text width at the drawn font size.
  font_pt <- 10
  col_pt <- 150
  row_pt <- 34
  box_h <- 18
  n$w <- vapply(n$label, function(s) grid::convertWidth(
    grid::grobWidth(grid::textGrob(s, gp = grid::gpar(fontsize = font_pt))),
    "points", valueOnly = TRUE), numeric(1)) + 12
  n$px <- n$x * col_pt
  n$py <- -n$y * row_pt

  # Edges as vertical cubic curves from the bottom of one box to the top of
  # the next.
  bezier <- function(x0, y0, x1, y1, k = 0.5, m = 30) {
    t <- seq(0, 1, length.out = m)
    dy <- (y0 - y1) * k
    data.frame(
      x = (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * x0 + 3 * (1 - t) * t^2 * x1 + t^3 * x1,
      y = (1 - t)^3 * y0 + 3 * (1 - t)^2 * t * (y0 - dy) + 3 * (1 - t) * t^2 * (y1 + dy) + t^3 * y1)
  }
  paths <- do.call(rbind, lapply(seq_len(nrow(e)), function(i) {
    a <- n[n$id == e$from[i], ]
    b <- n[n$id == e$to[i], ]
    p <- bezier(a$px, a$py - box_h / 2, b$px, b$py + box_h / 2 + 1)
    p$edge <- i
    p$dashed <- e$label[i] == "results folder"
    p
  }))

  col <- c("tantale class" = "#332288", "Biostrings class" = "#88CCEE",
           "exported function" = "#117733", "S3 method" = "#DDCC77")
  n$fill <- ifelse(n$is_fn, "white", col[n$group])
  n$border <- col[n$group]
  n$text <- ifelse(n$group == "tantale class", "white", "black")

  ggplot2::ggplot() +
    ggplot2::geom_path(ggplot2::aes(x, y, group = edge,
                                    linetype = ifelse(dashed, "dashed", "solid")),
                       data = paths, colour = "grey55", linewidth = 0.35,
                       arrow = grid::arrow(length = grid::unit(4, "pt"), type = "closed")) +
    ggplot2::geom_rect(ggplot2::aes(xmin = px - w / 2, xmax = px + w / 2,
                                    ymin = py - box_h / 2, ymax = py + box_h / 2),
                       data = n, fill = n$fill, colour = n$border,
                       linewidth = ifelse(n$is_fn, 0.7, 0.3)) +
    ggplot2::geom_text(ggplot2::aes(px, py, label = label), data = n,
                       colour = n$text, size = font_pt, size.unit = "pt") +
    ggplot2::scale_linetype_identity() +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::theme_void()
}

if (sys.nframe() == 0L) {
  d <- function_map_data(".")
  htmlwidgets::saveWidget(function_map_widget(d, height = "95vh"),
                          "pkgdown/assets/function_map.html",
                          selfcontained = TRUE, title = "tantale function map")
  unlink("pkgdown/assets/function_map_files", recursive = TRUE)

  fig <- function_map_figure(d)
  b <- ggplot2::ggplot_build(fig)
  r <- b$layout$panel_params[[1]]
  w_pt <- diff(r$x.range) + 20
  h_pt <- diff(r$y.range) + 20
  ggplot2::ggsave("pkgdown/assets/function_map.svg", fig,
                  width = w_pt / 72, height = h_pt / 72, units = "in")
}
