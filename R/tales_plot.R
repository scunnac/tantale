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
#' point per part, positioned by its place in the array, outlined by domain
#' type and filled by amino-acid length. Each repeat carries its RVD, each
#' terminus the short name of its code (\code{N-}, \code{-C}, or \code{??}
#' for a terminus that does not match its TALE domain profile; see
#' \code{\link{tales_anchor_codes}}). Arrays are listed from the top in
#' alphabetical order of \code{array_id}, compared byte by byte as in every
#' projection of a \code{tales} object, so the order does not depend on the
#' locale.
#'
#' Each length gets one colour, whatever the part. The canonical 34-aa repeat
#' and the 20-aa half-repeat that ends every array always get calm colours
#' (sand and pale blue), so a repeat of any other length (an aberrant repeat,
#' for instance) stands out. The other lengths take strong colours in
#' increasing order of length, from Paul Tol's schemes, which stay distinct
#' for colour-blind readers; past 17 distinct lengths the colours repeat.
#'
#' A \code{\link{tales_msa}} dispatches to \code{\link{plot.tales_msa}}
#' instead, being the more specific class.
#'
#' @param x A \code{\link{tales}} object, as returned by
#'   \code{\link{tales_from_telltales}} or in the \code{tales} element of
#'   \code{\link{tales_compare_distal}}'s output. A data frame of parts is
#'   accepted and converted with \code{\link{tales}}.
#' @param position Which coordinate to lay the parts out on. \code{"array"}
#'   (default) uses \code{position_in_array}, so each array starts at 1 and runs
#'   contiguously. \code{"alignment"} uses \code{alignment_position}, which
#'   requires an aligned object (or one demoted from a \code{\link{tales_msa}},
#'   which keeps the column): gaps then appear as empty columns and shared
#'   features line up. Aberrant repeats, for instance, are visible as a column
#'   in the aligned layout and scattered in the unaligned one.
#' @param facet_by Names of the columns whose values split the plot into
#'   panels, stacked in rows whose height follows the number of arrays.
#'   Each column must hold one value per array: \code{"seqnames"} (the
#'   default, one panel per source sequence) or \code{"strain"} for a set of
#'   genomes, for instance, and \code{c("strain", "seqnames")} for both.
#'   \code{NULL} draws a single panel, as does the default when \code{x}
#'   has no \code{seqnames} column.
#' @param ... Unused, present for compatibility with the \code{plot} generic.
#' @return A ggplot object. Like any ggplot, it is drawn when printed, which
#'   happens automatically at the console; inside a loop or a function, call
#'   \code{print()} on it.
#' @method plot tales
#' @export
#' @family TALE plots
#' @examples
#' x <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
#'                                       package = "tantale"))
#' plot(x)
#' plot(x, facet_by = NULL)
plot.tales <- function(x, position = c("array", "alignment"), facet_by = "seqnames", ...) {
  position <- match.arg(position)
  if (!is_tales(x)) x <- tales(x)
  .tales_require(x, "plot.tales")
  if (missing(facet_by) && !"seqnames" %in% names(x)) facet_by <- NULL
  if (!is.null(facet_by)) {
    if (!is.character(facet_by) || !length(facet_by) || anyNA(facet_by)) {
      cli::cli_abort("{.arg facet_by} must be a character vector of column names, or {.code NULL}.",
                     class = c("tantale_error_plot_facet", "tantale_error"))
    }
    absent <- setdiff(facet_by, names(x))
    if (length(absent)) {
      cli::cli_abort("{.arg facet_by} names {?a column/columns} absent from {.arg x}: {.field {absent}}.",
                     class = c("tantale_error_plot_facet", "tantale_error"))
    }
    # a column varying within an array would split that array across panels
    varying <- facet_by[vapply(facet_by, function(nm) {
      any(tapply(x[[nm]], x$array_id, function(z) length(unique(z))) > 1L)
    }, logical(1))]
    if (length(varying)) {
      cli::cli_abort("{.arg facet_by} {?column/columns} {.field {varying}} must hold one value per array.",
                     class = c("tantale_error_plot_facet", "tantale_error"))
    }
  }
  if (identical(position, "alignment") && !"alignment_position" %in% names(x)) {
    cli::cli_abort(
      c("{.code position = \"alignment\"} needs the {.field alignment_position} column.",
        "i" = "Align first with {.fn tales_align}, or use {.code position = \"array\"}."),
      class = c("tantale_error_projection_column", "tantale_error")
    )
  }
  anchors <- tales_anchor_codes()
  partsForPlots <- x %>%
    dplyr::mutate(
      # repeats show their RVD, termini their code's short name (N-, -C, ??)
      label = dplyr::coalesce(names(anchors)[match(rvd, anchors)],
                              dplyr::if_else(domain_type == "repeat", rvd, "")),
      aa_length = nchar(aa_seq),
      .x = if (identical(position, "alignment")) .data$alignment_position
           else .data$position_in_array)

  # One colour per length, whatever the part (ledger §49): the canonical
  # 34-aa repeat and the final 20-aa half-repeat keep calm colours, the other
  # lengths take strong ones in increasing order, recycled past 17.
  lengths <- sort(unique(partsForPlots$aa_length))
  others <- lengths[!lengths %in% c(20L, 34L)]
  pool <- c(unname(.tol_muted[c("rose", "indigo", "purple", "green", "olive", "wine", "teal")]),
            .tol_light)
  lengthColours <- stats::setNames(character(length(lengths)), lengths)
  lengthColours[as.character(others)] <- pool[(seq_along(others) - 1L) %% length(pool) + 1L]
  lengthColours[names(lengthColours) == "34"] <- .tol_muted[["sand"]]
  lengthColours[names(lengthColours) == "20"] <- .tol_muted[["cyan"]]
  partsForPlots <- partsForPlots %>%
    dplyr::mutate(length_aa = factor(aa_length, lengths),
                  colour = unname(lengthColours[as.character(aa_length)]),
                  label_colour = .text_colour_on(colour),
                  # byte order (ledger §40), first array at the top: a
                  # discrete y axis puts its first level at the bottom
                  array_id = factor(array_id,
                                    levels = rev(levels(.array_factor(array_id)))))

  p <- partsForPlots %>%
    ggplot2::ggplot(mapping = ggplot2::aes(fill = length_aa,
                                           color = domain_type,
                                           label = label,
                                           y = array_id,
                                           x = .x)) +
    ggplot2::scale_color_manual(name = "Domain type",
                                breaks = TALES_DOMAIN_TYPES,
                                values = c("N-terminus" = .tol_muted[["wine"]],
                                           "repeat" = "#BBBBBB",
                                           "C-terminus" = .tol_muted[["indigo"]])) +
    ggplot2::scale_fill_manual(name = "Length (aa)", values = lengthColours) +
    ggplot2::scale_x_continuous(
      name = if (identical(position, "alignment")) "Position in alignment" else "Position in array",
      breaks = 1:100, minor_breaks = NULL) +
    ggplot2::geom_point(shape = 21, size = 5, stroke = 0.9) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_text(mapping = ggplot2::aes(color = label_colour), size = 2.1) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(title = "Overview of TALE composition by genome") +
    ggplot2::theme_light() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = .tantale_colours$strip, colour = NA),
                   strip.text = ggplot2::element_text(colour = "grey20", face = "bold"))
  
  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(!!!rlang::syms(facet_by)),
                                 scales = "free_y", space = "free")
  }
  
  p
}


#### The actual method for msa ploting ####

#' Plot a multiple alignment of TALEs
#'
#' @description Draws the alignment as a heatmap: one row per array, one
#'   column per alignment position, with cells coloured by one layer and
#'   optionally labelled with another.
#'
#' @details
#' A \code{tales_msa} carries every layer at once (\code{rvd},
#' \code{dom_code} and whatever else the object holds), so \code{fill} and
#' \code{label} simply name two of them.
#'
#' Three things are decided independently, and it helps to read the figure
#' that way: what each cell *says*, what colour that text is, and what colour
#' the block behind it is.
#'
#' \strong{Cell text} is whatever \code{label} names, or nothing when
#' \code{label = NULL}. Termini are relabelled \code{N-} and \code{-C}; an
#' unidentified terminus keeps its \code{XXXXX} code.
#' Domain codes are padded to three characters so columns line up.
#'
#' \strong{Text colour} always answers one question: does this element match
#' the consensus of its column? Black for yes, red for no, grey where the
#' column has no consensus (a gap, or a tie). The block fills are all pale,
#' so the text reads on every one of them. The consensus is the most
#' frequent element in the column (\code{\link{tales_consensus}}), taken over
#' the labelled layer, so the text colour and the text itself always describe
#' the same thing.
#'
#' \strong{Block fill} is what \code{fill_type} selects, and it is the only
#' part that can be unavailable:
#'
#' \tabular{lll}{
#'   \strong{fill_type} \tab \strong{shows} \tab \strong{needs} \cr
#'   \code{"domain_clust"} \tab which cluster the domain falls in, cut at \code{h_cut} \tab \code{domain_distances} \cr
#'   \code{"domain_sim"} \tab protein-sequence similarity to the reference, 0-100 \tab \code{domain_distances} \cr
#'   \code{"rvd_sim"} \tab how alike the RVD's DNA-binding preference is to the reference's, -1 to 1 \tab a \code{label} layer \cr
#' }
#'
#' Every value here scores across the whole \code{tales_msa}, termini
#' included: \code{domain_distances} covers every distinct part sequence,
#' and \code{"domain_clust"}/\code{"domain_sim"} colour a terminus cell
#' exactly like a repeat cell.
#'
#' With no \code{domain_distances} and no \code{label}, every block is flat
#' grey: the text still carries the consensus comparison, but there is
#' nothing to colour blocks by.
#'
#' A cell with no value for the chosen layer keeps its text and is filled
#' grey. In \code{"rvd_sim"} that is the termini, which have no DNA-binding
#' preference and so no position on a specificity scale. The RVD similarity
#' table is TALVEZ's and covers 17 RVDs; a rarer RVD (\code{NV}, say) scores
#' 1 where it is identical to the reference's RVD and is grey elsewhere.
#'
#' \strong{The reference} matters for both similarity fills.
#' \code{ref_pattern} is matched against the array names and must identify
#' exactly one, otherwise the default is used with a warning; by default it is
#' the array with the most non-gap parts, ties broken alphabetically. The
#' reference row is marked with a trailing \code{_#}.
#'
#' \strong{Two panels may be attached.} Supplying \code{tale_distances} with
#' more than one array adds a dendrogram panel on the left; \code{consensus =
#' TRUE} adds a consensus panel on top. When either is present the return
#' value is an \code{aplot} composition: modify the alignment through its
#' \code{plotlist} element, since layers added to the composition itself do
#' not reach it.
#'
#' @param x A \code{\link{tales_msa}} object.
#' @param fill Layer whose values colour the cells. Defaults to
#'   \code{"dom_code"} when present, otherwise the first available residue
#'   layer.
#' @param label Layer whose values are written in the cells. Left unset it
#'   defaults to \code{"rvd"} when that is not already the \code{fill}; pass
#'   \code{NULL} explicitly for an unlabelled heatmap.
#' @param tale_distances Pairwise distances between whole TALEs, as the
#'   \code{tale_distances} element of a \code{\link{tales_compare_distal}} result.
#'   Used to build the tree panel that orders the alignment rows.
#' @param domain_distances Pairwise distances between distinct domains, as
#'   the \code{domain_distances} element of a \code{\link{tales_compare_distal}}
#'   result. Used to group domains into clusters, and to score each domain
#'   against the reference TALE\'s domain at the same alignment column.
#' @param h_cut Height at which the domain tree is cut to define clusters.
#'   Interpreted on a distance scale, so 0 means identical.
#' @param ref_pattern Regular expression matched against the array names to
#'   choose the reference TALE. Must identify exactly one.
#' @param consensus Whether to add a consensus panel above the alignment. The
#'   consensus is the most frequent element in each column of the labelled
#'   layer, so it always matches what the cells say.
#' @param fill_type One of \code{"domain_clust"}, \code{"domain_sim"} or
#'   \code{"rvd_sim"}. The first two colour cells by domain cluster or by
#'   protein-sequence similarity to the reference, for every distinct part,
#'   termini included (see Details). \code{"rvd_sim"} colours them instead
#'   by how alike each RVD's *DNA-binding preference* is to the reference
#'   TALE's RVD at that position, on a diverging scale over \code{[-1, 1]}.
#'   The domain- and RVD-level views genuinely differ: repeats carrying
#'   \code{HD} and \code{ND} differ in sequence yet both favour cytosine,
#'   while repeats differing only at positions 12-13 are near-identical
#'   proteins targeting different bases.
#' @param ... Unused, present for compatibility with the \code{plot} generic.
#'
#' @return The alignment plot, returned invisibly after being printed as a
#'   side effect: a \code{ggplot} normally, or, when \code{tale_distances}
#'   and/or \code{consensus} add extra panels, an
#'   \code{\link[aplot:insert_left]{aplot}} composition -- see Details.
#' @method plot tales_msa
#' @export
#' @family TALE plots
#' @examples
#' aligned <- data.frame(
#'   array_id = c("A1", "A1", "A1", "A2", "A2"),
#'   position_in_array = c(1L, 2L, 3L, 1L, 2L),
#'   alignment_position = c(1L, 2L, 3L, 1L, 3L),
#'   rvd = c("NTERM", "HD", "CTERM", "NTERM", "CTERM")
#' )
#' msa <- tales_msa(aligned)
#' plot(msa)
plot.tales_msa <- function(x, fill = NULL, label = NULL,
                           tale_distances = NULL, domain_distances = NULL,
                           h_cut = 10,
                           ref_pattern = NULL,
                           consensus = FALSE,
                           fill_type = "domain_clust",
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
  
  # The drawing code below still needs one matrix per layer for the
  # domain-cluster/similarity helpers below, which genuinely need a
  # domain-by-position layout to substitute values column by column.
  # as.matrix() always returns a matrix, including for a single array, so
  # the layers cannot arrive malformed the way a hand-assembled matrix could.
  domain_align <- as.matrix(x, value = fill)
  rvd_align <- if (!is.null(label)) as.matrix(x, value = label) else NULL
  arrayNames <- rownames(domain_align)
  countOfTales <- nrow(domain_align)
  
  # Both tables are addressed as id1/id2/dissim below; pairwise_distances()
  # validates them and folds a sim column into dissim.
  if (!is.null(tale_distances))   tale_distances   <- tibble::as_tibble(pairwise_distances(tale_distances))
  if (!is.null(domain_distances)) domain_distances <- tibble::as_tibble(pairwise_distances(domain_distances))
  
  
  # domainAlignLong/rvdAlignLong's second column used to be named
  # "position_in_array" -- a pre-existing mislabel: as.matrix.tales_msa()
  # indexes its columns by x$alignment_position, not each array's own
  # position_in_array, so that was always what these matrices' columns (and
  # everything melted from them below) actually carried. Renamed throughout
  # this function to match (ledger §21's closing note; class-design.md
  # §4.6). The x-axis title says "Position in alignment" to match (§22).
  
  # Getting domain align
  if (!is.null(domain_align)) {
    domainAlignLong <- .matrix_to_long(domain_align) %>%
      dplyr::as_tibble()
    colnames(domainAlignLong) <- c("array_id", "alignment_position", "dom_code")
    domainAlignLong %<>% dplyr::mutate(array_id = as.character(array_id),
                                       dom_code = stringr::str_pad(dom_code, 3, "left"))
    # Consensus match computed straight off x, the tales_msa itself: no need to
    # wait for domain_align/domainAlignLong above, since
    # .tales_consensus_match_long() takes a tales_msa directly and does its own
    # sparse-to-complete accounting for gaps -- this is the specific piece of
    # plot.tales_msa() that no longer touches the array-by-position matrix at all.
    domainMatchConsensusLong <- .tales_consensus_match_long(x, value_col = fill)
    colnames(domainMatchConsensusLong) <- c("array_id", "alignment_position", "matchConsensusDomain")
    domainAlignLong %<>% dplyr::left_join(domainMatchConsensusLong,
                                          by = dplyr::join_by(array_id, alignment_position))
  }
  
  
  # Getting rvd align if available
  if (!is.null(rvd_align)) {
    rvdAlignLong <- .matrix_to_long(rvd_align) %>%
      dplyr::as_tibble()
    colnames(rvdAlignLong) <- c("array_id", "alignment_position", "rvd")
    rvdAlignLong %<>% dplyr::mutate(array_id = as.character(array_id),
                                    rvd = gsub("NTERM", "N-", rvd),
                                    rvd = gsub("CTERM", "-C", rvd)
    )
    
    # Tale rvd text color if possible: coloring of RVDs in alignment depending
    # on whether they match the consensus at the position. Computed off x
    # directly, same as the domain case above -- see the comment there.
    rvdMatchConsensusLong <- .tales_consensus_match_long(x, value_col = label)
    colnames(rvdMatchConsensusLong) <- c("array_id", "alignment_position", "matchConsensusRvd")
    
    # RVD similarity relative to the reference, the RVD-level counterpart of
    # domainSimVsRef: that one scores protein sequence similarity, this one
    # how alike two RVDs' DNA-binding preferences are. They come apart -- HD
    # and ND are different domains with identical specificity, while domains
    # differing only at positions 12-13 are nearly identical proteins
    # targeting different bases. The domain block below overrides refTaleId
    # when domain distances are given.
    refTaleId <- .pick_ref_name(align = rvd_align, ref_tag = ref_pattern)
    rvdSimAlignLong <- .rvd_to_match_align(rvd_align = rvd_align,
                                           ref_tag = ref_pattern) %>%
      .matrix_to_long() %>%
      dplyr::as_tibble()
    colnames(rvdSimAlignLong) <- c("array_id", "alignment_position", "rvdSimVsRef")
    rvdSimAlignLong %<>% dplyr::mutate(array_id = as.character(array_id))
    
    rvdAlignLong %<>%
      dplyr::left_join(rvdMatchConsensusLong,
                       by = dplyr::join_by(array_id, alignment_position)) %>%
      dplyr::left_join(rvdSimAlignLong,
                       by = dplyr::join_by(array_id, alignment_position))
  }
  
  # Assign main alignment object in long format
  if (!is.null(domain_align) & !is.null(rvd_align)) {
    domainAlignLong %<>% dplyr::inner_join(rvdAlignLong,
                                           by = dplyr::join_by(array_id, alignment_position),
                                           unmatched = "error",
                                           relationship = "one-to-one")
  } else if (!is.null(domain_align) & is.null(rvd_align)) {
    domainAlignLong <- domainAlignLong
  } else if (is.null(domain_align)) {
    domainAlignLong <- rvdAlignLong
  } else {
    cli::cli_abort(
      c("Cannot build the requested plot from these arguments.",
        "i" = "Check {.arg fill}, {.arg label} and {.arg fill_type} against the object's layers."),
      class = c("tantale_error_bad_argument", "tantale_error"))
  }
  
  # joining domain cluster if possible
  # joining domain similarity relative to ref
  if (!is.null(domain_distances) & !is.null(domain_align)) {
    # .domain_to_cluster_align()/.domain_to_sim_align() still take a
    # parameter literally named domain_sim -- untouched here, §21 items 1-4
    # reassigned to the maintainer (ledger §21 closing note, §19). Only the
    # outer variable being passed in is renamed.
    domainClusterAlignLong <- .domain_to_cluster_align(domain_align = domain_align,
                                                       domain_sim = domain_distances,
                                                       h_cut = h_cut) %>%
      .matrix_to_long() %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(value = as.character(value))
    colnames(domainClusterAlignLong) <- c("array_id", "alignment_position", "domainClusterId")
    
    refTaleId <- .pick_ref_name(align = domain_align, ref_tag = ref_pattern)
    domainSimAlignLong <- .domain_to_sim_align(domain_align = domain_align,
                                               domain_sim = domain_distances,
                                               ref_tag = ref_pattern) %>%
      .matrix_to_long() %>%
      dplyr::as_tibble()
    colnames(domainSimAlignLong) <- c("array_id", "alignment_position", "domainSimVsRef")
    # Join with main tible
    domainAlignLong %<>%
      dplyr::left_join(domainClusterAlignLong,
                       by = dplyr::join_by(array_id, alignment_position)) %>%
      dplyr::left_join(domainSimAlignLong,
                       by = dplyr::join_by(array_id, alignment_position))
  }
  
  
  # Building TALE tree if possible
  if (!is.null(tale_distances) & countOfTales > 1) {
    taleDistForDendo <- tale_distances[tale_distances$id1 %in% arrayNames, ]
    taleDistForDendo <- taleDistForDendo[taleDistForDendo$id2 %in% arrayNames, ]
    taleDist <- .pairwise_long_to_matrix(taleDistForDendo, "dissim")
    taleDist <- taleDist[arrayNames, ]
    taleDist <- taleDist[, arrayNames]
    taleshclust <- stats::hclust(as.dist(taleDist))
  }
  
  
  # Add a symbol to designate the reference if necessary
  if (exists("refTaleId")) { # in the tibble
    domainAlignLong$array_id[domainAlignLong$array_id == refTaleId] <-  paste0(
      domainAlignLong$array_id[domainAlignLong$array_id == refTaleId],
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
  bp <- domainAlignLong %>% ggplot2::ggplot(mapping = ggplot2::aes(
    x = alignment_position, y = array_id)
  ) +
    ggplot2::scale_x_discrete(
      name = "Position in alignment",
      limits = factor(1:max(domainAlignLong$alignment_position))
    ) +
    ggplot2::scale_y_discrete(name = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  # Colours (R/palette.R). Every fill is pale, because the text colour
  # carries the consensus match and has to read on all of them.
  # scale_fill_manual() takes `values`, not `palette`: the name collided with
  # the `palette` discrete_scale() supplies internally, so every call to this
  # function aborted with "formal argument 'palette' matched by multiple actual
  # arguments" regardless of fill_type. discrete_scale() is the scale that
  # actually accepts a palette *function*, which is what is wanted here.
  domainClusterFillScale <- ggplot2::discrete_scale(aesthetics = "fill",
                                                    name = "Repeats cluster",
                                                    palette = grDevices::colorRampPalette(.tantale_colours$clusters),
                                                    drop = TRUE,
                                                    na.translate = FALSE,
                                                    guide = NULL)
  # Low similarity is the strong colour, so a divergent domain stands out.
  domainSimFillScale <- ggplot2::scale_fill_gradientn(
    name = "Similarity relative to reference",
    colours = rev(.tantale_colours$sequential),
    na.value = .tantale_colours$no_value,
    guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(10, "lines")))
  # The RVD score is a correlation on [-1, 1], so it wants a diverging scale
  # centred on zero rather than the sequential one used for domain similarity.
  rvdSimFillScale <- ggplot2::scale_fill_gradientn(
    name = "RVD specificity vs reference\n(grey: no score)",
    colours = .tantale_colours$diverging, limits = c(-1, 1),
    na.value = .tantale_colours$no_value,
    guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(10, "lines")))
  labelConsensusColorScale <- ggplot2::scale_color_manual(name = "Match consensus?",
                                                          na.value = .tantale_colours$no_consensus,
                                                          values = c(`FALSE` = .tantale_colours$mismatch,
                                                                     `TRUE` = .tantale_colours$match)
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
                          linewidth = NA,
                          family = "mono",
                          size = 3, fontface = "bold",
                          na.rm = TRUE
      )
  } else if (!is.null(domain_distances) & !is.null(rvd_align)) {
    if (fill_type == "domain_sim") {
      p <- bp +
        domainSimFillScale +
        labelConsensusColorScale +
        ggplot2::geom_label(mapping = ggplot2::aes(fill = domainSimVsRef,
                                                   label = rvd,
                                                   color = matchConsensusRvd),
                            linewidth = NA,
                            family = "mono",
                            size = 3, fontface = "bold",
                            na.rm = TRUE
        )
    } else if (fill_type == "domain_clust") {
      p <- bp +
        domainClusterFillScale +
        labelConsensusColorScale +
        ggplot2::geom_label(mapping = ggplot2::aes(fill = domainClusterId,
                                                   label = rvd,
                                                   color = matchConsensusRvd),
                            linewidth = NA,
                            family = "mono",
                            size = 3, fontface = "bold",
                            na.rm = TRUE
        )
    } else {
      cli::cli_abort("{.arg fill_type} must be {.val domain_clust}, {.val domain_sim} or {.val rvd_sim}.", class = c("tantale_error_msa_layer", "tantale_error"))
    }
  } else if (!is.null(domain_distances) & is.null(rvd_align)) {
    if (fill_type == "domain_sim") {
      p <- bp +
        domainSimFillScale +
        labelConsensusColorScale +
        ggplot2::geom_label(mapping = ggplot2::aes(fill = domainSimVsRef,
                                                   label = dom_code,
                                                   color = matchConsensusDomain),
                            linewidth = NA,
                            family = "mono",
                            size = 3, fontface = "bold",
                            na.rm = TRUE
        )
    } else if (fill_type == "domain_clust") {
      p <- bp +
        domainClusterFillScale +
        labelConsensusColorScale +
        ggplot2::geom_label(mapping = ggplot2::aes(fill = domainClusterId,
                                                   label = dom_code,
                                                   color = matchConsensusDomain),
                            linewidth = NA,
                            family = "mono",
                            size = 3, fontface = "bold",
                            na.rm = TRUE
        )
    } else {
      cli::cli_abort("{.arg fill_type} must be {.val domain_clust}, {.val domain_sim} or {.val rvd_sim}.", class = c("tantale_error_msa_layer", "tantale_error"))
    }
  } else if (is.null(domain_distances) & !is.null(rvd_align)) {
    p <- bp +
      labelConsensusColorScale +
      ggplot2::geom_label(mapping = ggplot2::aes(label = rvd,
                                                 color = matchConsensusRvd),
                          fill = .tantale_colours$no_value,
                          linewidth = NA,
                          family = "mono",
                          size = 3, fontface = "bold",
                          na.rm = TRUE
      )
  } else if (is.null(domain_distances) & is.null(rvd_align)) {
    p <- bp +
      labelConsensusColorScale +
      ggplot2::geom_label(mapping = ggplot2::aes(label = dom_code,
                                                 color = matchConsensusDomain),
                          fill = .tantale_colours$no_value,
                          linewidth = NA,
                          family = "mono",
                          size = 3, fontface = "bold",
                          na.rm = TRUE
      )
  } else {
    cli::cli_abort(
      c("Cannot build the requested plot from these arguments.",
        "i" = "Check {.arg fill}, {.arg label} and {.arg fill_type} against the object's layers."),
      class = c("tantale_error_bad_argument", "tantale_error"))
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
    consensusAlign <- if (!is.null(rvd_align)) rvd_align else domain_align
    # aplot's height is a *ratio* of the main plot, so a fixed value would grow
    # with the number of arrays -- several rows tall for a large group. Scale it
    # so the consensus stays about one alignment row high whatever the count.
    consensusHeight <- max(0.08, min(0.45, 1.0 / countOfTales))
    finalPlot <- aplot::insert_top(
      finalPlot,
      .consensus_panel(consensusAlign,
                       n_positions = max(domainAlignLong$alignment_position),
                       pad = is.null(rvd_align)),
      height = consensusHeight
    )
  }
  
  # aplot composes panels with patchwork's guide collection, which places the
  # collected legend by the *combined* object's own theme, not by any
  # individual panel's -- so the "bottom" position set on the main panel above
  # is silently ignored once a tree and/or consensus panel is attached, and
  # the legend falls back to patchwork's default of "right". Printing a
  # patchwork-converted copy is the only way to make that take effect; the
  # value returned to the caller stays the aplot, so its $plotlist/[i,j]
  # panel-indexing API is unchanged for anyone composing further.
  if (inherits(finalPlot, "aplot")) {
    print(aplot::as.patchwork(finalPlot) &
            ggplot2::theme(legend.position = "bottom"))
  } else {
    print(finalPlot)
  }
  invisible(finalPlot)
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
#' Still matrix-based, unlike the rest of \code{plot.tales_msa()}'s consensus
#' handling: \code{domain_align}/\code{rvd_align} already exist regardless
#' (the domain-cluster/similarity helpers still need them), so there is
#' nothing to save by also rebuilding this one call site around
#' \code{\link{.tales_consensus_long}}.
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
                        linewidth = NA, family = "mono",
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
#
# "domain", not "repeat": dom_code names any distinct part sequence, and the
# two termini are scored by these helpers exactly like a repeat is (ledger
# section 19) -- the old "repeat_*" names for these internals claimed a
# specificity the code never had.

.domain_to_sim_align <- function(domain_align, domain_sim, ref_tag = NULL) {
  # A function that substitute the domain codes with the aa similarity relative to a
  # reference domain for each column. The ref domain is the one from a TALE that
  # is defined as a reference in the alignment. This function takes as input, the
  # domain alignment and the df output by `.format_domain_distances_mat()` This
  # function outputs the modified alignment matrix
  
  refRowIdx <- match(.pick_ref_name(domain_align, ref_tag = ref_tag), rownames(domain_align))
  simAlign <- apply(domain_align, 2,
                    function(column) {
                      refState <- column[refRowIdx]
                      relevantSims <- subset(domain_sim, subset = id1 == refState)
                      sim <- 100 - relevantSims$dissim[match(column, relevantSims$id2, nomatch = NA)]
                      if (is.na(refState)) sim[!is.na(column)] <- 0 # if reference domain is NA, set the aligned domain sim = 0
                      return(sim)
                    }
  )
  simAlign <- matrix(simAlign, nrow = nrow(domain_align)) # in case of 1-row matrix
  rownames(simAlign) <- rownames(domain_align)
  colnames(simAlign) <- colnames(domain_align)
  return(simAlign)
}

#' Convert domain alignment to clusterID alignment
#'
#' @param domain_sim A long, three columns data frame with pairwise similarity scores between domains as available in the \code{domain_distances} element of the object returned by the \code{\link{tales_compare_distal}} function.
#' @param domain_align a multiple Tal domain-code sequences alignment in the form of a matrix as returned by \code{\link{tales_align}}.
#' @param h_cut a numeric value indicating the height at which to cut the hclust tree of domains. Interpreted on a distance scale (0 = identical).
#' @return a matrix with exactly the same dimension as the input \code{domain_align} but containing clusterID instead of
#' dom_code.
#' @noRd
.domain_to_cluster_align <- function(domain_sim, domain_align, h_cut = 10) {
  # as.dist() expects a DISTANCE, which is what the class stores.
  domain_dissim <- .pairwise_long_to_matrix(domain_sim, "dissim")
  dist_clust <- hclust(as.dist(domain_dissim))
  dist_cut <- as.data.frame(cbind(DomID = dist_clust$labels, Dom_clust = cutree(dist_clust, h = h_cut)))
  clustIDAlign <- apply(domain_align, 2,
                        function(column){
                          as.numeric(dist_cut$Dom_clust[match(column, dist_cut$DomID)])
                        })
  clustIDAlign <- matrix(clustIDAlign, nrow = nrow(domain_align)) # in case of 1-row matrix
  rownames(clustIDAlign) <- rownames(domain_align)
  colnames(clustIDAlign) <- colnames(domain_align)
  return(clustIDAlign)
}

#' Recode an RVD alignment as similarity to a reference row
#'
#' Substitutes each RVD with a score expressing how similar its DNA-binding
#' preference is to the RVD of a reference TALE, column by column. This is the
#' RVD-level counterpart of \code{.domain_to_sim_align()}, which works on
#' protein sequence similarity instead: the two come apart, since domains can be
#' sequence-divergent yet share an RVD, or near-identical yet differ at
#' positions 12-13.
#'
#' Wired through \code{fill_type = "rvd_sim"} in \code{\link{plot.tales_msa}}
#' -- the RVD-level counterpart of \code{"domain_sim"}. The internal
#' \code{rvdSimDf} dataset it reads also feeds \code{.rvd_score_table()}.
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
                      sims <- relevantSims$Cor[match(column, relevantSims$rvd2, nomatch = NA)]
                      # An RVD outside rvd_sims (rare RVDs such as NV) still
                      # scores 1 against an identical reference RVD, as in
                      # .rvd_score_table(). Gaps and termini stay NA.
                      same <- is.na(sims) & !is.na(refState) & !is.na(column) &
                        column == refState & !refState %in% tales_anchor_codes()
                      sims[same] <- 1
                      sims
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
      cli::cli_warn(
        c("{.arg ref_tag} matches {sum(match)} sequences, not one.",
          "i" = "Falling back to the default reference, the longest array."),
        class = "tantale_warning_ambiguous_ref")
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
