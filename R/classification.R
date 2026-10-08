
##### Tale classification ####
# tales_group() (the joint student-code rework, ledger §11) split into two
# functions, one per clustering algorithm, because the two never shared much
# beyond "cluster distMat and attach the result": k_range/seed only ever
# meant anything for k-medoids, hclust's plot_tree only ever meant anything
# for hclust, and cramming both under one `method` switch is what made the
# original hard to read. Both still return `x` with `group` filled, via the
# shared correspondence check in .tales_attach_groups() below.


#' Validate x and build the pairwise distance matrix both grouping methods share
#' @return A square numeric matrix.
#' @noRd
.tales_group_distmat <- function(x, dists) {
  if (!is_tales(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  # Coercing accepts both a tale_distances and a plain table; as.matrix() then
  # replaces the hand-written acast() and asserts squareness on the way.
  # A parameter named `dists`, not `tale_distances` like the public-facing
  # argument this receives: R resolves tale_distances(tale_distances) fine
  # (a same-named local never masks a function in call position), but it
  # reads as a mistake, so this internal-only helper avoids it.
  as.matrix(tale_distances(dists))
}


#' Group TALEs by hierarchical clustering of their pairwise distance
#'
#' Clusters TALE arrays with \code{\link[stats:hclust]{hclust}} on their
#' pairwise distance, cuts the tree into exactly \code{k} groups, and returns
#' the \code{\link{tales}} object the distances were computed from with the
#' result attached as its \code{group} column.
#'
#' @details
#' The clustering is computed from `tale_distances`, but the result belongs on the
#' `tales` object the distances were computed from, so that is what comes
#' back. `group` is a recognised `tales` column, validated as constant within
#' an array -- it is an array-level property, like `seqnames`.
#'
#' Taking `x` rather than returning a bare lookup table is what makes the
#' correspondence checkable: the array names in `tale_distances` must be the array
#' names in `x`, and this is the only place that can be verified. A mismatch
#' is an error rather than a silent `NA` group, because a partly-grouped
#' object is the kind of thing that fails much later and confusingly.
#'
#' The bare mapping is still one line away if you want it:
#' `unique(out[c("array_id", "group")])`.
#'
#' The tree is built directly on the distances
#' (\code{stats::hclust(stats::as.dist(distMat))}), matching the original
#' DisTAL clustering this package reimplements. The tree is cut with
#' \code{stats::cutree(tree, k = k)}, which always succeeds, including when
#' a tie in merge heights would make a height-based cut ambiguous.
#'
#' @param x A [tales] object -- the one whose comparison produced `tale_distances`.
#' @param tale_distances A [tale_distances] object, as returned by [tales_compare_distal()].
#'   A plain data frame with the same columns is also accepted and coerced.
#' @param k Integer, the number of groups to cut the tree into.
#' @param plot_tree Logical, whether to draw the `ggtree` dendrogram: colored by
#'   group, with a dashed line at the cut height. `FALSE` by default.
#' @return `x` with an added (or replaced) `group` column.
#' @seealso [tales_group_kmedoids()], the alternative method;
#'   [tales_compare_distal()], which produces both inputs.
#' @export
#' @family pairwise distances
#' @examples
#' x <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
#'                                       package = "tantale"))
#' cmp <- tales_compare_distal(x)
#' grouped <- tales_group_hclust(cmp$tales, cmp$tale_distances, k = 2)
#' unique(grouped[c("array_id", "group")])
tales_group_hclust <- function(x, tale_distances, k = NULL, plot_tree = FALSE) {
  distMat <- .tales_group_distmat(x, tale_distances)

  if (!is.numeric(k) || length(k) != 1L || is.na(k)) {
    cli::cli_abort(
      "{.arg k} must be a single number, the number of groups to cut the tree into.",
      class = c("tantale_error_group_hclust_k", "tantale_error")
    )
  }

  taleTree <- stats::hclust(stats::as.dist(distMat), method = "ward.D")
  treeCuts <- stats::cutree(taleTree, k = k)
  taleGroups <- data.frame(name = names(treeCuts), group = treeCuts, row.names = NULL)

  if (isTRUE(plot_tree)) print(.tales_group_hclust_plot(taleTree, treeCuts, k))

  .tales_attach_groups(x, taleGroups)
}

#' The dendrogram plot for tales_group_hclust(), built only when asked for
#' @noRd
.tales_group_hclust_plot <- function(taleTree, treeCuts, k) {
  # The display cutoff: a height between the (n-k)-th and (n-k+1)-th merges,
  # i.e. exactly the height cutree(taleTree, h = cutOff) would need to
  # reproduce cutree(taleTree, k = k) -- derived directly from the sorted
  # merge heights rather than searched for.
  n <- length(taleTree$height) + 1L
  h <- sort(taleTree$height)
  cutOff <- mean(h[c(n - k, n - k + 1L)])

  g <- split(names(treeCuts), treeCuts)
  p <- ggtree::ggtree(taleTree)
  # tidytree::MRCA() emits a cli "Invalid edge matrix for <phylo>" message
  # for some subtree shapes, harmlessly -- it still returns the right MRCA
  # (a <tbl_df> internal fallback, not an error); see ledger section 6.
  clades <- vapply(g, function(nms) suppressMessages(tidytree::MRCA(p, nms)),
                   numeric(1))
  p <- tidytree::groupClade(p, clades, group_name = "subtree") +
    ggtree::aes(color = subtree)

  # k often runs to 25 groups or more, beyond what any palette keeps apart,
  # so colour only separates neighbouring clades (cycled in leaf order) and
  # the group number printed under each clade identifies it.
  leafOrder <- unique(treeCuts[taleTree$labels[taleTree$order]])
  cycle <- unname(.tol_muted[c("indigo", "rose", "teal", "wine", "olive",
                               "cyan", "purple", "green", "sand")])
  cladeColours <- stats::setNames(rep_len(cycle, k), leafOrder)
  groupLabels <- p$data %>%
    dplyr::filter(isTip) %>%
    dplyr::mutate(group = treeCuts[label]) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(y = mean(y), .groups = "drop")

  p + ggtree::layout_dendrogram() +
    ggtree::geom_tiplab(ggtree::aes(label = label),
                        hjust = 1, angle = 90, align = FALSE,
                        color = "black", offset = -2) +
    ggplot2::scale_color_manual(values = c(`0` = "grey60", cladeColours), guide = "none") +
    ggplot2::geom_vline(xintercept = -cutOff, linetype = 2, colour = "grey50") +
    ggplot2::geom_label(data = groupLabels,
                        mapping = ggplot2::aes(x = max(taleTree$height) * 0.04, y = y, label = group),
                        inherit.aes = FALSE, size = 3, fontface = "bold",
                        label.size = 0, fill = "white", colour = "grey20") +
    ggplot2::labs(subtitle = sprintf("Cut into %d groups at height %.2f (dashed line)", k, cutOff)) +
    ggplot2::xlab("Height") +
    ggtree::theme_dendrogram(plot.margin = ggplot2::margin(6, 6, 150, 6))
}


#' Group TALEs by k-medoids clustering of their pairwise distance
#'
#' Clusters TALE arrays with \code{\link[cluster:pam]{pam}} for every
#' candidate \code{k} in \code{k_range}, then returns the \code{\link{tales}}
#' object the distances were computed from with the chosen clustering's
#' result attached as its \code{group} column.
#'
#' @details
#' The clustering is computed from `tale_distances`, but the result belongs on the
#' `tales` object the distances were computed from, so that is what comes
#' back. `group` is a recognised `tales` column, validated as constant within
#' an array -- it is an array-level property, like `seqnames`.
#'
#' Taking `x` rather than returning a bare lookup table is what makes the
#' correspondence checkable: the array names in `tale_distances` must be the array
#' names in `x`, and this is the only place that can be verified. A mismatch
#' is an error rather than a silent `NA` group, because a partly-grouped
#' object is the kind of thing that fails much later and confusingly.
#'
#' The bare mapping is still one line away if you want it:
#' `unique(out[c("array_id", "group")])`.
#'
#' Unlike a hierarchical tree, PAM has no single structure that can be cut at
#' an arbitrary `k` after the fact: `k` is a parameter of the clustering
#' itself. So one clustering is computed per candidate in `k_range`, and `k`
#' picks which of those to keep:
#'
#' - `k` a single number uses that candidate directly.
#' - `k = "auto"` picks the elbow of the silhouette-vs-k curve, the point
#'   after which adding groups stops improving the fit much. It uses the
#'   first step of the Kneedle algorithm only, a practical heuristic; check
#'   the silhouette plot when the choice matters.
#' - `k = NULL` (the default) shows the silhouette plot and asks for a number
#'   at the console, but only when \code{\link{interactive}()} is `TRUE`.
#'   In a script, a test or a vignette render, `k = NULL` errors instead of
#'   blocking on input that will never arrive.
#'
#' @inheritParams tales_group_hclust
#' @param k_range Integer vector of candidate values of `k` to evaluate. Each
#'   must be between 1 and one less than the number of arrays being grouped
#'   (`cluster::pam()`'s own requirement).
#' @param k See Details.
#' @param seed Passed to \code{set.seed()} before every \code{cluster::pam()}
#'   call, so the same candidate always clusters the same way.
#' @param plot_silhouette Logical, whether to draw the silhouette-vs-k plot.
#' @return `x` with an added (or replaced) `group` column.
#' @seealso [tales_group_hclust()], the alternative method;
#'   [tales_compare_distal()], which produces both inputs.
#' @export
#' @family pairwise distances
#' @examples
#' x <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
#'                                       package = "tantale"))
#' cmp <- tales_compare_distal(x)
#' grouped <- tales_group_kmedoids(cmp$tales, cmp$tale_distances,
#'                                 k_range = 2:3, k = 2)
#' unique(grouped[c("array_id", "group")])
tales_group_kmedoids <- function(x, tale_distances, k_range = NULL, k = NULL,
                                 seed = 7, plot_silhouette = TRUE) {
  distMat <- .tales_group_distmat(x, tale_distances)

  if (is.null(k_range) || !is.numeric(k_range)) {
    cli::cli_abort(
      "{.arg k_range} must be a numeric vector of candidate {.arg k} values.",
      class = c("tantale_error_group_kmedoids_krange", "tantale_error")
    )
  }

  # cluster::pam() requires k in {1, ..., n - 1}; checked here, against the
  # actual number of arrays being clustered, so a bad k_range fails with a
  # message naming the arrays it came from -- not pam()'s own generic
  # "Number of clusters 'k' must be in {1,2, .., n-1}; hence n >= 2" surfacing
  # from inside lapply(), which is what a k_range = 2:4 example against a
  # 4-array fixture actually produced before this check existed.
  n <- nrow(distMat)
  if (any(k_range < 1) || any(k_range > n - 1)) {
    cli::cli_abort(
      c("{.arg k_range} must be between 1 and {n - 1} for {n} array{?s}.",
        "x" = "Got {.val {k_range}}."),
      class = c("tantale_error_group_kmedoids_krange", "tantale_error")
    )
  }

  allPam <- lapply(k_range, function(kpam) {
    set.seed(seed)
    as.list(cluster::pam(stats::as.dist(distMat), kpam))
  })
  silhVals <- vapply(allPam, function(a) a$silinfo$avg.width, numeric(1))

  if (is.null(k)) {
    if (!interactive()) {
      cli::cli_abort(
        c("{.arg k} is required outside an interactive session.",
          "i" = "Pass a number, or {.val auto} to pick one from the silhouette curve."),
        class = c("tantale_error_group_kmedoids_k", "tantale_error")
      )
    }
    if (isTRUE(plot_silhouette)) .tales_group_kmedoids_plot(k_range, silhVals)
    cat("Choose a number of groups:\t")
    numGroups <- as.numeric(readLines(con = stdin(), 1))
  } else if (identical(k, "auto")) {
    numGroups <- .tales_group_kmedoids_elbow(silhVals, k_range)
    if (isTRUE(plot_silhouette)) .tales_group_kmedoids_plot(k_range, silhVals, numGroups)
    cli::cli_inform("The number of groups is automatically decided based on the silhouette value: {numGroups}")
  } else if (is.numeric(k) && length(k) == 1L) {
    numGroups <- k
    if (isTRUE(plot_silhouette)) .tales_group_kmedoids_plot(k_range, silhVals, numGroups)
    cli::cli_inform("Number of groups is decided based on the provided value of k: {numGroups}")
  } else {
    cli::cli_abort(
      c("{.arg k} must be {.code NULL}, {.val auto}, or a single number.",
        "x" = "Got {.obj_type_friendly {k}}."),
      class = c("tantale_error_group_kmedoids_k", "tantale_error")
    )
  }

  group <- allPam[[which(k_range == numGroups)]]$clustering
  taleGroups <- data.frame(name = names(group), group = group, row.names = NULL)
  .tales_attach_groups(x, taleGroups)
}

#' The silhouette-vs-k plot for tales_group_kmedoids(), built only when asked for
#' @noRd
.tales_group_kmedoids_plot <- function(k_range, silhVals, highlight = NULL) {
  col <- if (is.null(highlight)) .tol_muted[["cyan"]] else
    ifelse(k_range != highlight, .tol_muted[["cyan"]], .tol_muted[["wine"]])
  plot(k_range, silhVals, pch = 19, col = col,
       xlab = "number of groups", ylab = "average silhouette values")
}

#' The elbow of the silhouette-vs-k curve
#'
#' A partial, first-step application of the Kneedle algorithm
#' (<https://raghavan.usc.edu/papers/kneedle-simplex11.pdf>): good enough on
#' the silhouette curves this has been tried on, not a full implementation.
#' @noRd
.tales_group_kmedoids_elbow <- function(v, k) {
  stopifnot(length(v) == length(k))
  n <- length(v)
  a <- (v[n] - v[1]) / (k[n] - k[1])
  b <- -1
  c <- (v[1] * k[n] - v[n] * k[1]) / (k[n] - k[1])
  d <- vapply(seq_len(n), function(i) abs(a * i + b * v[i] + c) / sqrt(a * a + b * b), numeric(1))
  k[which.max(d)]
}


#' Put a name->group mapping onto the tales object it was computed from
#'
#' Kept separate because it is the part with the invariant in it: every array
#' in `x` must be grouped, and every grouped name must be an array of `x`.
#' Either direction failing means the distances did not come from this object.
#' @noRd
.tales_attach_groups <- function(x, groups) {
  arrays <- unique(x$array_id)
  missing <- setdiff(arrays, groups$name)
  extra   <- setdiff(groups$name, arrays)
  if (length(missing) || length(extra)) {
    cli::cli_abort(
      c("{.arg tale_distances} does not describe the arrays in {.arg x}.",
        "x" = if (length(missing))
          "In {.arg x} but not grouped: {.val {utils::head(missing, 5)}}{if (length(missing) > 5) ' ...' else ''}",
        "x" = if (length(extra))
          "Grouped but not in {.arg x}: {.val {utils::head(extra, 5)}}{if (length(extra) > 5) ' ...' else ''}",
        "i" = "{.arg tale_distances} must come from comparing this same object."),
      class = c("tantale_error_group_mismatch", "tantale_error"))
  }
  x$group <- groups$group[match(x$array_id, groups$name)]
  x
}






#' Cross-tabulate row_var x col_var, applying fun_aggregate even to empty cells
#'
#' Replacement for `reshape2::dcast(df, row_var ~ col_var, value.var =,
#' drop = FALSE, fun.aggregate =)` (ledger section 14) -- verified
#' `identical()` to it, including the specific behaviour `talomes_heatmap()`
#' depends on: `drop = FALSE` means every row_var x col_var combination gets
#' a cell even when no rows match it, with `fun_aggregate` called on that
#' *empty* subset (giving `length(x) = 0` for the allele-count call,
#' `NA` for the representative-allele call) rather than a cell silently
#' skipped. `dplyr::summarise(.drop = FALSE)` is what reproduces that; a
#' plain `count()`/`pivot_wider(values_fill=)` would not, because it can
#' only fill a *constant*, not call the real aggregate function on nothing.
#' Column order follows `sort(levels(factor(col_var)))`, matching dcast's
#' own factor-based ordering -- load-bearing here, since the caller
#' zips `colnames(result)` positionally against a separately-sorted vector
#' right after.
#' @noRd
.dcast_count_matrix <- function(df, row_var, col_var, value_var, fun_aggregate) {
  df[[col_var]] <- factor(df[[col_var]])
  out <- df %>%
    dplyr::group_by(.data[[row_var]], .data[[col_var]], .drop = FALSE) %>%
    dplyr::summarise(..v = fun_aggregate(.data[[value_var]]), .groups = "drop") %>%
    tidyr::pivot_wider(id_cols = dplyr::all_of(row_var), names_from = dplyr::all_of(col_var),
                       values_from = "..v") %>%
    dplyr::arrange(.data[[row_var]]) %>%
    as.data.frame()
  out[, c(row_var, levels(df[[col_var]]))]
}

#' Heatmap of RVD sequence variants across strains and TALE groups
#'
#' @description
#' Draws a talome overview: one column per TALE group, one row per strain,
#' each cell showing which RVD sequence variant that strain carries in that
#' group. A strain's talome is its whole complement of TALEs, so the plot
#' shows at a glance which groups each strain has and where strains carry
#' different variants of the same TALE.
#'
#' Within a group, variants are ranked by how many strains carry them, and
#' the cell colour is that rank: the most common variant is the palest and
#' rarer ones are darker.
#' The \code{#} after each group label counts its distinct
#' variants. A white cell means the strain has no member in that group. A
#' cell can hold several colours side by side when a strain carries more
#' than one variant in the same group. Dendrograms order strains and groups
#' by the similarity of their variant profiles.
#'
#' @param tale_annotation A data frame with one row per TALE and at least a
#'   group, a strain and an RVD-sequence column. Or a \code{\link{tales}}
#'   object carrying group and strain columns, with one value per array (for
#'   instance a grouped \code{tales} to which a strain column was added):
#'   the RVD sequences are then computed with
#'   \code{\link{tales_rvd_strings}}, and arrays without repeats are left
#'   out.
#' @param group_col Name of the \code{tale_annotation} column holding TALE
#'   groups, drawn as columns (e.g. the \code{group} from
#'   \code{\link{tales_group_kmedoids}}).
#' @param strain_col Name of the column holding strain names, drawn as rows.
#' @param rvd_col Name of the column holding RVD sequences (e.g. from
#'   \code{\link{tales_rvd_strings}}). Not needed when
#'   \code{tale_annotation} is a \code{tales} object.
#' @param trunc_tales_col Optional name of a logical column marking
#'   truncated TALEs; those are labelled "T" in their cell (with
#'   \code{plot_type = "all"}).
#' @param extra_col Optional name of a column with further information about
#'   each strain (e.g. origin), drawn as a side bar on the right.
#' @param x_lab,y_lab,title Axis names and plot title.
#' @param colors Colours from the most common variant to the rarest. They
#'   are interpolated over the ranks present, so the most common variant always
#'   gets the first colour and the rarest the last. The default runs from pale
#'   to dark wine, so a rare variant stands out.
#' @param margins Space for the row dendrogram, column dendrogram, row
#'   names and column names, in that order, counted in heatmap cells.
#'   \code{NULL} (default) sizes the dendrograms to a quarter of the
#'   heatmap's width and height, between 1.5 and 5 cells. With
#'   \code{plot_type = "all"}, it gives the names and the title the room
#'   their text takes, whatever the size of the device, and the cells share
#'   the rest; with \code{plot_type = "single"}, the names get 3 cells, and
#'   the column dendrogram also holds the title and gets at least 3.5 cells.
#' @param sep_width Width of the separator between adjacent cells.
#' @param sep_color Colour of the separator between adjacent cells.
#' @param inner_sep_color Colour of the separator between variants within
#'   one cell.
#' @param save_path Optional file path; the image format follows the file
#'   extension. If \code{NULL} (default), the heatmap is drawn on the
#'   current device.
#' @param plot_type Either \code{"all"} to draw every allele, or
#'   \code{"single"} to draw one representative allele per group.
#' @return \code{NULL}, invisibly. Called for the side effect of drawing the
#'   heatmap (or, if \code{save_path} is given, writing it to a file).
#' @export
#' @family TALE plots
#' @examples
#' ann <- data.frame(
#'   group = c("G1", "G1", "G1", "G2", "G2"),
#'   strain = c("S1", "S2", "S3", "S1", "S2"),
#'   rvdseq = c("NI-HD-NG", "NI-HD-NG", "NN-HD-NG",
#'             "HD-NI-NG-NG", "HD-NI-NG-NG")
#' )
#' talomes_heatmap(ann, group_col = "group", strain_col = "strain",
#'                 rvd_col = "rvdseq")
#'
#' # From a tales object: one group and one strain per array
#' x <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
#'                                       package = "tantale"))
#' x$group <- c(ROI_00001 = 1, ROI_00002 = 1, ROI_00003 = 2, ROI_00004 = 2)[x$array_id]
#' x$strain <- c(ROI_00001 = "S1", ROI_00002 = "S2", ROI_00003 = "S1",
#'               ROI_00004 = "S2")[x$array_id]
#' talomes_heatmap(x, group_col = "group", strain_col = "strain")
talomes_heatmap <- function(tale_annotation, group_col, strain_col, rvd_col, trunc_tales_col = NULL, extra_col = NULL,
                            x_lab = "TALE Group", y_lab = "Strain", title = "RVD sequences variants",
                            plot_type = "all",
                            colors = .tantale_colours$variant_ranks, margins = NULL,
                            sep_width = 5, sep_color = "white", inner_sep_color = "white", save_path = NULL) {
  plot_type <- match.arg(plot_type, c("all", "single"))

  if (is_tales(tale_annotation)) {
    # one row per array: its RVD string and the per-array columns asked for
    per_array <- c(group_col, strain_col, trunc_tales_col, extra_col)
    missing <- setdiff(per_array, names(tale_annotation))
    if (length(missing) > 0L) {
      cli::cli_abort("{.arg tale_annotation} has no column{?s} {.field {missing}}.",
                     class = c("tantale_error_talome_column", "tantale_error"))
    }
    arrays <- unique(as.data.frame(tale_annotation)[c("array_id", per_array)])
    if (anyDuplicated(arrays$array_id)) {
      cli::cli_abort(c("{.field {per_array}} must hold one value per array.",
                       "x" = "{.val {arrays$array_id[duplicated(arrays$array_id)][1]}} has several."),
                     class = c("tantale_error_talome_column", "tantale_error"))
    }
    rvd <- tales_rvd_strings(tale_annotation)
    arrays$rvdseq <- as.character(rvd)[match(arrays$array_id, names(rvd))]
    tale_annotation <- arrays[!is.na(arrays$rvdseq), ]
    rvd_col <- "rvdseq"
  }

  ## rename colnames of tale annotation
  colnames(tale_annotation)[which(colnames(tale_annotation) == group_col)] <- "group"
  colnames(tale_annotation)[which(colnames(tale_annotation) == strain_col)] <- "strain"
  colnames(tale_annotation)[which(colnames(tale_annotation) == rvd_col)] <- "rvdseq"
  if (!is.null(extra_col)) colnames(tale_annotation)[which(colnames(tale_annotation) == extra_col)] <- "extra_col"
  if (!is.null(trunc_tales_col)) colnames(tale_annotation)[which(colnames(tale_annotation) == trunc_tales_col)] <- "truncTale"
  
  tale_annotation <- tale_annotation %>% dplyr::mutate(aberrantRepeat = grepl("[a-z]", rvdseq))
  tale_annotation$rvdseq <- toupper(tale_annotation$rvdseq)
  
  if(is.numeric(tale_annotation$group)) tale_annotation$group <- paste0("G", tale_annotation$group)
  
  ## convert rvdseq to rvdfac, ranking of abundance
  tale_annotation <- tale_annotation %>% dplyr::group_by(group) %>% dplyr::mutate(rvdfac = do.call((function(x) {
    levels(x) <- nlevels(x) + 1 - rank(table(x), ties.method = "first")
    return(as.integer(as.character(x)))
  }), list(as.factor(rvdseq))))
  
  ## get number of rvdseqs per group per strain
  ## it is used for the layout of heatmap cell
  numAlleles <- .dcast_count_matrix(tale_annotation, "strain", "group", "rvdseq", function(x) length(x))
  
  rownames(numAlleles) <- numAlleles$strain
  numAlleles <- numAlleles[,-1]
  
  ## count number of alleles
  ## for column label
  variantsCount <- vapply(sort(unique(tale_annotation$group)), function(g) {
    length(unique(tale_annotation[tale_annotation$group == g,]$rvdseq))
  }, integer(1), USE.NAMES = TRUE)
  colnames(numAlleles) <- paste0(colnames(numAlleles), " #", variantsCount)
  
  tale_annotation <- tale_annotation %>% dplyr::group_by(group, strain) %>% dplyr::mutate(reprsntRVDfac = ifelse(1 %in% rvdfac, 1, rvdfac[1]))
  
  ## representative alleles per group per strain
  ## in case a strain has more than 1 rvd seqs in 1 group
  ## find the most abundant rvdseqin that group
  ## if this strain has that rvdseq, take it as representative rvdseq
  ## if not, take randomly 1 rvdseq that strain has
  reprsntAlleles <- .dcast_count_matrix(tale_annotation, "strain", "group", "reprsntRVDfac", function(x) as.integer(x[1]))
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
  
  
  ## plot sizes
  auto_names <- is.null(margins)
  if (is.null(margins)) {
    # heatmap.2() ("single") draws the title inside the column-dendrogram
    # panel, which then needs room for it
    topMin <- if (plot_type == "single") 3.5 else 1.5
    margins <- c(min(5, max(1.5, ncol(numAlleles) / 4)),
                 min(5, max(topMin, nrow(numAlleles) / 4)), 3, 3)
  }
  widleft <- margins[1]
  heitop <- margins[2]
  widright <- margins[3]
  heibot <- margins[4]
  
  rankColours <- grDevices::colorRampPalette(colors)(max(tale_annotation$rvdfac))

  if (plot_type == "single") { # plot representative alleles
    codedAlleles1 <- apply(reprsntAlleles, 2, function(x) ifelse(x == 0, NA, x))
    if (!is.null(save_path)) {
      img_format <- gsub(".*\\.", "", basename(save_path))
      img_size <-  list(save_path, width = (widleft + ncol(codedAlleles1) + widright)/2.54, height = (heitop + nrow(codedAlleles1) + heibot)/2.54)
      if (img_format %in% c("bmp", "jpeg", "png", "tiff")) {
        img_size <- c(img_size, units = "in", res = 300)
      }
      do.call(img_format, img_size)
    }
    gplots::heatmap.2(as.matrix(codedAlleles1),
                      trace = "none",
                      col = rankColours[seq_len(max(codedAlleles1, na.rm = T))],
                      breaks = 0:max(codedAlleles1, na.rm = T),
                      density.info = "none",
                      key = F,
                      xlab = x_lab,
                      ylab = y_lab,
                      margins = c(7, 7),
                      colsep = 0:(ncol(codedAlleles1)-0),
                      rowsep = 0:(nrow(codedAlleles1)-0),
                      sepcolor = sep_color,
                      sepwidth = rep(sep_width/100, 2),
                      Rowv = as.dendrogram(Strain_HC),
                      Colv = as.dendrogram(TALE_HC),
                      main = title,
                      na.color = .tantale_colours$absent,
                      lmat = rbind(c(4, 3), c(2, 1)),
                      lhei = c(heitop, nrow(codedAlleles1) + heibot), ##
                      lwid = c(widleft, ncol(codedAlleles1) + widright)
    )
    
  } else if (plot_type == "all") { # plot all alleles
    uniqueRVD <- numAlleles
    rorder <- order.dendrogram(as.dendrogram(Strain_HC))
    corder <- order.dendrogram(as.dendrogram(TALE_HC))
    uniqueRVD <- uniqueRVD[rev(rorder), corder]
    
    
    colmat <- rankColours
    rowLabels <- paste0(rownames(uniqueRVD), "  #", rowSums(uniqueRVD, na.rm = TRUE))
    titleHeight <- 2
    if (auto_names) {
      # The name panels and the title get room in cm, measured from their
      # text, so that they fit whatever the size of the device; the cells
      # share what is left. Measured on the device drawn on, or, when the
      # file device is not open yet, on a throwaway one.
      if (!is.null(save_path)) grDevices::pdf(NULL)
      # layout() draws text at 0.66 of its size once the grid has three rows
      # or columns (see mfrow in ?par), and this one always has.
      line_cm <- graphics::par("csi") * 2.54 * 0.66
      text_cm <- function(s) max(graphics::strwidth(s, units = "inches", cex = 1.2 * 0.66)) * 2.54
      widright <- text_cm(rowLabels) + 2.5 * line_cm
      heibot <- text_cm(colnames(uniqueRVD)) + 2.5 * line_cm
      titleHeight <- 4 * line_cm
      if (!is.null(save_path)) grDevices::dev.off()
    }
    inCm <- function(x) if (auto_names) graphics::lcm(x) else x

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
    if (!is.null(save_path)) {
      img_format <- gsub(".*\\.", "", basename(save_path))
      img_size <-  list(save_path, width = sum(widleft, rep(1, ncol(uniqueRVD)), widright)/2.54, height = sum(titleHeight, heitop, rep(1, nrow(uniqueRVD)), heibot)/2.54)
      if (img_format %in% c("bmp", "jpeg", "png", "tiff")) {
        img_size <- c(img_size, units = "in", res = 300)
      }
      do.call(img_format, img_size)
    }
    # Read after the output device is chosen: par() on no device opens one.
    df <- par(no.readonly = TRUE)

    layout(laymat, widths = c(widleft, rep(1, ncol(uniqueRVD)), ifelse(is.null(extra_col), 0.1, 0.7), inCm(widright), ifelse(is.null(extra_col), 0.1, 3)),
           heights = c(inCm(titleHeight), heitop, rep(1, nrow(uniqueRVD)), inCm(heibot)))
    # layout.show(nplots+7)
    
    
    ## plot heatmap
    for (c in seq_len(ncol(uniqueRVD))) {
      gname <- gsub(" \\#\\d+", "", colnames(uniqueRVD)[c])
      # gname <- gsub("G", "", gname)
      g1 <- tale_annotation[tale_annotation$group == gname,]
      # g1$rvdseq <- as.integer(as.factor(g1$rvdseq))
      for (r in seq_len(nrow(uniqueRVD))) {
        par(mar = rep(0, 4))
        nelements <- uniqueRVD[r, c]
        if (nelements == 0) {
          image(z = matrix(0), col = .tantale_colours$absent, axes = F)
          if (sep_width > 0) box(lwd = sep_width/2, col = sep_color)
        } else {
          g1s1 <- g1[g1$strain == rownames(uniqueRVD[r,]),]
          rvd.factor <- g1s1$rvdfac
          if (is.null(trunc_tales_col)) {
            truncTale <- NA
          } else {
            truncTale <- ifelse(g1s1$truncTale, "T", NA)
          }
          # plot.index <- r + (c-1)*nrow(uniqueRVD)
          color.elements <- colmat[rvd.factor]
          par(mar = rep(0, 4))
          image(z = matrix(1:(nelements), ncol = 1), col = color.elements, axes = F)
          text(seq(0,1, length.out = nelements), 0, labels = truncTale, cex = 1, font = 2,
               col = .text_colour_on(color.elements))
          if (nelements > 1) {
            abline(v = seq(0.5/(nelements-1), 1-0.5/(nelements-1), length.out = nelements -1), col = inner_sep_color, lty = 1)
          }
          if (sep_width > 0) {
            par(mar = rep(0, 4))
            box(lwd = sep_width/2, col = sep_color)}
        }
        
      }
    }
    
    ## dendrograms
    par(mai = rep(0, 4))
    plot(as.dendrogram(Strain_HC), horiz = TRUE, axes = FALSE, leaflab = "none", yaxs = "i")
    par(mai = rep(0, 4))
    plot(as.dendrogram(TALE_HC), axes = FALSE, leaflab = "none", xaxs = "i")
    
    if (!is.null(extra_col)) {
      ## extra column
      # extra.bar <- sample(c("Hanoi", "Hatay", "Namdinh", NA), nrow(uniqueRVD), replace = T)
      extra.bar <- tale_annotation$extra_col[match(rownames(uniqueRVD), tale_annotation$strain)]
      rextra.bar <- unique(extra.bar[!is.na(extra.bar)])
      rextra.bar <- data.frame("lab" = sort(rextra.bar), "fac" = seq_along(rextra.bar))
      extra.col <- rep_len(.tol_light, nrow(rextra.bar))[match(extra.bar, rextra.bar$lab)]
      extra.col[is.na(extra.col)] <- .tantale_colours$no_value
      par(mar = c(0, 0, 0, 0))
      image(z = matrix(seq_len(nrow(uniqueRVD)), nrow = 1), col = rev(extra.col), yaxt = "n", xaxt = "n", axes = F)
    } else {
      par(mar = c(0, 0, 0, 0))
      image(z = matrix(seq_len(nrow(uniqueRVD)), nrow = 1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    }
    
    
    ## rownames
    par(mar = c(0, 0.5, 0, 0))
    image(z = matrix(seq_len(nrow(uniqueRVD)), nrow = 1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    text(-1, seq(0, 1, length.out = nrow(uniqueRVD)), labels = rev(rowLabels), font = 1, col = "black", cex = 1.2, adj = 0)
    mtext(side = 4, at = 0.5, text = y_lab, col = "black", padj = 0, line = -1)
    
    if (!is.null(extra_col)) {
      ## legend column
      par(mar = c(0,0.5,1,0))
      plot(rep(0, nrow(rextra.bar)), -seq(from = 0, by = 0.8, length.out = nrow(rextra.bar)), type = "p", pch = 15, col = rep_len(.tol_light, nrow(rextra.bar)), axes = F, main = extra_col, xlab = NA, ylab = NA, cex = 4, ylim = c(-nrow(uniqueRVD), 0), xlim = c(0,2))
      text(rep(0.3, nrow(rextra.bar)), -seq(from = 0, by = 0.8, length.out = nrow(rextra.bar)), labels = rextra.bar$lab, font = 1, col = "black", bg = "red", cex = 1.2, adj = 0)
    } else {
      par(mar = c(0,0,0,0))
      image(matrix(0), col = "white", axes = F)
    }
    
    
    ## colnames
    par(mar = c(0.5, 0, 0.5, 0))
    image(z = matrix(seq_len(ncol(uniqueRVD)), ncol = 1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    text(seq(0, 1, length.out = ncol(uniqueRVD)), 1, labels = colnames(uniqueRVD), font = 1, col = "black", bg = "red", cex = 1.2, srt = 90, adj = 1)
    mtext(side = 1, at = 0.5, text = x_lab, col = "black", padj = 0, line = -1)
    
    
    ## title
    par(xpd = T, mai = rep(0, 4))
    image(z = matrix(1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    text(0, 1/3, labels = title, font = 2, col = "black", cex = 2, pos = 1)
    
    
    par(df)
  }
  if (!is.null(save_path)) {
    grDevices::dev.off()
  }
  invisible(NULL)
}



