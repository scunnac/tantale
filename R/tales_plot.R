#### Plotting tales and tales_msa objects ####
#
# The two plot methods and everything only they use. They lived apart --
# plot.tales in distalr.R and plot.tales_msa in msa.R -- which are both
# legacy files named after the functions that used to dominate them rather
# than after what they now hold.
#
# Organised by method family rather than by class, matching tales_print.R
# and tales_summary.R. For those two the arrangement is forced, since
# format.tales and format.tales_msa share their layout helpers; here it is a
# choice, made so the sibling methods can be read against each other.
#
# The internals below the methods are all reached only from plot.tales_msa().
# plot.tales has none of its own: it calls .tales_require() from
# requirements.R, which the whole package shares.


#' Plot the domain composition of a set of TALE arrays
#'
#' @description
#' A compact, information-rich view of the arrays in a \code{tales} object: one
#' point per part, positioned by its place in the array, coloured by domain type
#' and filled by amino-acid length, with the RVD printed on each repeat.
#'
#' A \code{\link{tales_msa}} dispatches to \code{\link{plot.tales_msa}}
#' instead, being the more specific class.
#'
#' @param x A \code{\link{tales}} object, as returned by
#'   \code{\link{tales_from_telltale}} or in the \code{tales} element of
#'   \code{\link{tales_compare}}'s output. A legacy \code{tale_parts} data
#'   frame is accepted and converted.
#' @param position Which coordinate to lay the parts out on. \code{"array"}
#'   (default) uses \code{position_in_array}, so each array starts at 1 and runs
#'   contiguously. \code{"alignment"} uses \code{alignment_position}, which
#'   requires an aligned object (or one demoted from a \code{\link{tales_msa}},
#'   which keeps the column): gaps then appear as empty columns and shared
#'   features line up. Aberrant repeats, for instance, are visible as a column
#'   in the aligned layout and scattered in the unaligned one.
#' @param ... Unused, present for compatibility with the \code{plot} generic.
#' @return The ggplot object, invisibly printed as a side effect.
#' @method plot tales
#' @export
#' @family TALE plots
plot.tales <- function(x, position = c("array", "alignment"), ...) {
  position <- match.arg(position)
  if (!is_tales(x)) x <- tales(x)
  .tales_require(x, "plot.tales")
  if (identical(position, "alignment") && !"alignment_position" %in% names(x)) {
    cli::cli_abort(
      c("{.code position = \"alignment\"} needs the {.field alignment_position} column.",
        "i" = "Align first with {.fn tales_align}, or use {.code position = \"array\"}."),
      class = c("tantale_error_projection_column", "tantale_error")
    )
  }
  partsForPlots <- x %>%
    dplyr::mutate(label = dplyr::if_else(domain_type == "repeat", rvd, ""),
                  aa_length = factor(nchar(aa_seq)),
                  .x = if (identical(position, "alignment")) .data$alignment_position
                       else .data$position_in_array)

  p <- partsForPlots %>%
    ggplot2::ggplot(mapping = ggplot2::aes(fill = aa_length,
                                           color = domain_type,
                                           label = label,
                                           y = array_id,
                                           x = .x)) +
    ggplot2::scale_color_viridis_d(option = "rocket") +
    ggplot2::scale_fill_discrete() +
    ggplot2::scale_x_continuous(
      name = if (identical(position, "alignment")) "Position in alignment" else "Position in array",
      breaks = 1:100, minor_breaks = NULL) +
    ggplot2::geom_point(shape = 21, size = 5, stroke = 0.9) +
    ggnewscale::new_scale_color() +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_text(size = 2.1, color = "white") +
    ggplot2::labs(title = "Overview of TALE composition by genome") +
    ggplot2::theme_light()

  # seqnames groups arrays by source contig. It is optional in a tales, so the
  # facet is added only when it is there -- a fasta-derived object has none.
  if ("seqnames" %in% names(x)) {
    p <- p + ggplot2::facet_grid(seqnames ~ ., scales = "free_y", space = "free")
  }

  print(p)
  invisible(p)
}


#### The actual method for msa ploting ####

#' Plot a multiple alignment of TALEs
#'
#' @description Draws the alignment as a heatmap: one row per array, one
#'   column per alignment position, with cells coloured by one layer and
#'   optionally labelled with another.
#'
#' @details
#' A \code{tales_msa} carries every layer at once -- \code{rvd},
#' \code{dom_code} and whatever else the object holds -- so \code{fill} and
#' \code{label} name two of them rather than being passed as separate
#' matrices.
#'
#' Three things are decided independently, and it helps to read the figure
#' that way: what each cell *says*, what colour that text is, and what colour
#' the block behind it is.
#'
#' \strong{Cell text} is whatever \code{label} names, or nothing when
#' \code{label = NULL}. Termini are relabelled \code{N-} and \code{-C}; an
#' unidentified terminus keeps its \code{XXXXX} code, which is deliberately
#' not mistakable for an RVD. Repeat codes are padded to three characters so
#' columns line up.
#'
#' \strong{Text colour} always answers one question: does this element match
#' the consensus of its column? Cyan for yes, pink for no. The consensus is
#' the most frequent element in the column (\code{\link{tales_consensus}}),
#' taken over the labelled layer -- so the text colour and the text itself
#' always describe the same thing.
#'
#' \strong{Block fill} is what \code{fill_type} selects, and it is the only
#' part that can be unavailable:
#'
#' \tabular{lll}{
#'   \strong{fill_type} \tab \strong{shows} \tab \strong{needs} \cr
#'   \code{"repeat_clust"} \tab which cluster the repeat falls in, cut at \code{h_cut} \tab \code{domain_sim} \cr
#'   \code{"repeat_sim"} \tab protein-sequence similarity to the reference, 0-100 \tab \code{domain_sim} \cr
#'   \code{"rvd_sim"} \tab how alike the RVD's DNA-binding preference is to the reference's, -1 to 1 \tab a \code{label} layer \cr
#' }
#'
#' With no \code{domain_sim} and no \code{label}, every block is flat grey:
#' the text still carries the consensus comparison, but there is nothing to
#' colour blocks by.
#'
#' A cell with no value for the chosen layer keeps its text and loses its
#' colour. In \code{"rvd_sim"} that is the termini, which have no DNA-binding
#' preference and so no position on a specificity scale.
#'
#' \strong{The reference} matters for both similarity fills.
#' \code{ref_pattern} is matched against the array names and must identify
#' exactly one, otherwise the default is used with a warning; by default it is
#' the array with the most non-gap parts, ties broken alphabetically. The
#' reference row is marked with a trailing \code{_#}.
#'
#' \strong{Two panels may be attached.} Supplying \code{tal_sim} with more
#' than one array adds a dendrogram panel on the left; \code{consensus = TRUE}
#' adds a consensus panel on top. When either is present the return value is an
#' \code{aplot} composition rather than a single ggplot, so modify the
#' alignment through its \code{plotlist} element rather than adding layers to
#' the result directly.
#'
#' @param x A \code{\link{tales_msa}} object.
#' @param fill Layer whose values colour the cells. Defaults to
#'   \code{"dom_code"} when present, otherwise the first available residue
#'   layer.
#' @param label Layer whose values are written in the cells. Left unset it
#'   defaults to \code{"rvd"} when that is not already the \code{fill}; pass
#'   \code{NULL} explicitly for an unlabelled heatmap.
#' @param tal_sim Pairwise distances between whole TALEs, as the
#'   \code{tale_distances} element of a \code{\link{tales_compare}} result.
#'   Used to build the tree panel that orders the alignment rows.
#' @param domain_sim Pairwise distances between repeat units, as the
#'   \code{domain_distances} element of a \code{\link{tales_compare}}
#'   result. Used to group repeats into clusters, and to score each repeat
#'   against the reference TALE\'s repeat at the same alignment column.
#' @param h_cut Height at which the repeat tree is cut to define clusters.
#'   Interpreted on a distance scale, so 0 means identical.
#' @param ref_pattern Regular expression matched against the array names to
#'   choose the reference TALE. Must identify exactly one.
#' @param consensus Whether to add a consensus panel above the alignment. The
#'   consensus is the most frequent element in each column of the labelled
#'   layer, so it always matches what the cells say.
#' @param fill_type One of \code{"repeat_clust"}, \code{"repeat_sim"} or
#'   \code{"rvd_sim"}. The first two colour cells by repeat cluster or by
#'   protein-sequence similarity to the reference. \code{"rvd_sim"} colours
#'   them instead by how alike each RVD's *DNA-binding preference* is to the
#'   reference TALE's RVD at that position, on a diverging scale over
#'   \code{[-1, 1]}. The repeat- and RVD-level views genuinely differ:
#'   \code{HD} and \code{ND} are distinct repeats with identical
#'   specificity, while repeats differing only at positions 12-13 are
#'   near-identical proteins targeting different bases.
#' @param ... Unused, present for compatibility with the \code{plot} generic.
#'
#' @return An \code{\link[aplot:insert_left]{aplot}} object.
#' @method plot tales_msa
#' @export
#' @family TALE plots
plot.tales_msa <- function(x, fill = NULL, label = NULL,
                           tal_sim = NULL, domain_sim = NULL,
                           h_cut = 10,
                           ref_pattern = NULL,
                           consensus = FALSE,
                           fill_type = "repeat_clust",
                           ...) {

  # Resolve the two layers against what this alignment actually carries.
  available <- intersect(TALES_RESIDUE_COLS, names(x))
  if (is.null(fill)) fill <- if ("dom_code" %in% available) "dom_code" else available[1]
  # missing(), not is.null(): an explicit label = NULL asks for no labels at
  # all, which is different from not having said which layer to use.
  if (missing(label) && "rvd" %in% available && !identical(fill, "rvd")) label <- "rvd"
  for (layer in c(fill, label)) {
    if (!is.null(layer) && !layer %in% names(x)) {
      cli::cli_abort(
        c("This alignment has no {.field {layer}} layer.",
          "i" = "Available: {.field {available}}"),
        class = c("tantale_error_msa_layer", "tantale_error")
      )
    }
  }

  arrayNames <- unique(x$array_id)
  countOfTales <- length(arrayNames)
  if (countOfTales < 1L) {
    cli::cli_abort(
      "Cannot plot an alignment with no arrays in it.",
      class = c("tantale_error_msa_empty", "tantale_error")
    )
  }

  # The drawing code below works on one matrix per layer. as.matrix() always
  # returns a matrix, including for a single array, so the layers cannot
  # arrive malformed the way a hand-assembled matrix could.
  repeat_align <- as.matrix(x, value = fill)
  rvd_align <- if (!is.null(label)) as.matrix(x, value = label) else NULL
  repeat_sim <- domain_sim
  arrayNames <- rownames(repeat_align)
  countOfTales <- nrow(repeat_align)

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
                                                          na.value = "grey35",
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





#' Build the one-row consensus panel used by plot.tales_msa()
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





#### Turning an alignment of something into an alignment of something else ####
#
# Each of these takes the matrix plot.tales_msa() is drawing and returns a
# matrix of the same shape holding a different quantity: a cluster id, a
# similarity to the reference, an RVD specificity score. They are the fill
# layers, and plot.tales_msa() is their only caller, so they live beside it.

.repeat_to_sim_align <- function(repeat_align, repeat_sim, ref_tag = NULL) {
  # A function that substitute the repeatIDs with the aa similarity relative to a
  # reference repeat for each column. The ref repeat is the one from a TALE that
  # is defined as a reference in the alignment. This function takes as input, the
  # repeat alignment and the df output by `.format_repeat_dist_mat()` This
  # function outputs the modified alignment matrix
  
  refRowIdx <- match(.pick_ref_name(repeat_align, ref_tag = ref_tag), rownames(repeat_align))
  simAlign <- apply(repeat_align, 2,
                    function(column) {
                      refState <- column[refRowIdx]
                      relevantSims <- subset(repeat_sim, subset = id1 == refState)
                      sim <- 100 - relevantSims$dissim[match(column, relevantSims$id2, nomatch = NA)]
                      if (is.na(refState)) sim[!is.na(column)] <- 0 # if reference repeat is NA, set the aligned repeat sim = 0
                      return(sim)
                    }
  )
  simAlign <- matrix(simAlign, nrow = nrow(repeat_align)) # in case of 1-row matrix
  rownames(simAlign) <- rownames(repeat_align)
  colnames(simAlign) <- colnames(repeat_align)
  return(simAlign)
}

#' Convert repeat alignment to clusterID alignment
#'
#' @param repeat_sim A long, three columns data frame with pairwise similarity scores between repeats as available in the \code{domain_distances} element of the object returned by the \code{\link{tales_compare}} function.
#' @param repeat_align a multiple Tal repeat sequences alignment in the form of a matrix as returned by \code{\link{tales_align}}.
#' @param h_cut a numeric value indicating the height at which to cut the hclust tree of repeats. Interpreted on a distance scale (0 = identical).
#' @return a matrix with exactly the same dimension as the input \code{repeat_sim} but containing clusterID instead of
#' repeatID.
#' @noRd
.repeat_to_cluster_align <- function(repeat_sim, repeat_align, h_cut = 10) {
  # as.dist() expects a DISTANCE, which is what the class stores.
  repeat_dissim <- as.matrix(reshape2::acast(repeat_sim, id1 ~ id2, value.var = "dissim"))
  dist_clust <- hclust(as.dist(repeat_dissim))
  dist_cut <- as.data.frame(cbind(RepID = dist_clust$labels, Rep_clust = cutree(dist_clust, h = h_cut)))
  clustIDAlign <- apply(repeat_align, 2,
                        function(column){
                          as.numeric(dist_cut$Rep_clust[match(column, dist_cut$RepID)])
                        })
  clustIDAlign <- matrix(clustIDAlign, nrow = nrow(repeat_align)) # in case of 1-row matrix
  rownames(clustIDAlign) <- rownames(repeat_align)
  colnames(clustIDAlign) <- colnames(repeat_align)
  return(clustIDAlign)
}

#' Recode an RVD alignment as similarity to a reference row
#'
#' Substitutes each RVD with a score expressing how similar its DNA-binding
#' preference is to the RVD of a reference TALE, column by column. This is the
#' RVD-level counterpart of \code{.repeat_to_sim_align()}, which works on
#' protein sequence similarity instead: the two come apart, since repeats can be
#' sequence-divergent yet share an RVD, or near-identical yet differ at
#' positions 12-13.
#'
#' Currently unwired: no \code{fill_type} in either plotting function requests
#' an RVD-level layer. It is the only consumer of the internal
#' \code{rvdSimDf} dataset.
#'
#' @param rvd_align A character matrix of aligned RVDs.
#' @param rvd_sims A data frame of pairwise RVD similarity with columns
#'   \code{rvd1}, \code{rvd2} and \code{Cor}. Defaults to the package's
#'   internal \code{rvdSimDf}.
#' @param ref_tag Pattern selecting the reference row; see
#'   \code{.pick_ref_name()}.
#' @return A numeric matrix with the dimensions and dimnames of
#'   \code{rvd_align}.
#' @keywords internal
.rvd_to_match_align <- function(rvd_align, rvd_sims = rvdSimDf, ref_tag = NULL) {
  refRowIdx <- match(.pick_ref_name(rvd_align, ref_tag = ref_tag), rownames(rvd_align))
  simAlign <- apply(rvd_align, 2,
                    function(column) {
                      refState <- column[refRowIdx]
                      relevantSims <- subset(rvd_sims, subset = rvd1 == refState)
                      relevantSims$Cor[match(column, relevantSims$rvd2, nomatch = NA)]
                    }
  )
  simAlign <- matrix(simAlign, nrow = nrow(rvd_align)) # in case of 1-row matrix
  rownames(simAlign) <- rownames(rvd_align)
  colnames(simAlign) <- colnames(rvd_align)
  return(simAlign)
}


#### Choosing the reference row ####

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
