#### Summarising tales and tales_msa objects ####
#
# print() has to be cheap: it fires every time an object is echoed. summary()
# does not, so the things worth computing but too slow or too detailed for a
# print method go here -- chiefly the anomaly report, which otherwise nobody
# sees unless they think to call tales_anomalies().
#
# Following the R convention, summary() returns an object and a print method
# renders it. That way the numbers are available, not merely visible:
# summary(x)$n_distinct_repeats beats re-deriving it.


#' Emit label/value pairs with a label column wide enough for all of them
#' @noRd
.summary_lines <- function(pairs) {
  pairs <- pairs[!vapply(pairs, function(p) is.null(p[[2]]), logical(1))]
  if (!length(pairs)) return(invisible(NULL))
  w <- max(nchar(vapply(pairs, `[[`, character(1), 1)))
  for (p in pairs) cat(sprintf("  %-*s  %s\n", w, p[[1]], p[[2]]))
  invisible(NULL)
}


#' Summarise a tales object
#'
#' @description
#' The measures worth knowing about a set of TALE arrays before doing anything
#' with them, and the ones too expensive for \code{\link{print.tales}}.
#'
#' @details
#' \strong{Distinct domains} counts different \code{dom_code}s against the
#' number of parts, broken down by domain type. \emph{Domains}, not repeats:
#' a \code{dom_code} identifies any distinct part sequence, and the two
#' termini are parts like the repeats are. On the reference fixture 251
#' distinct codes cover 180 repeats and 71 termini, so calling the total a
#' repeat count would overstate it by nearly a third.
#'
#' The ratio to parts is a property of the biology: TALEs reuse repeats
#' heavily, within an array and between arrays. It is also what decides the
#' cost of \code{\link{tales_compare_distal}}, whose pairwise comparison runs over
#' distinct domains rather than over parts.
#'
#' \strong{Repeats per array} counts repeats only, excluding the termini, so
#' it is the number that determines how long a target box each TALE
#' recognises. Note that \code{tell_tales()}'s log reports array length as the
#' number of domain hits, which includes the two termini and is therefore
#' larger by up to two.
#'
#' \strong{Complete arrays} are those where both an N- and a C-terminus were
#' identified. An incomplete one is not necessarily wrong -- the array may sit
#' at the end of a contig, or a terminus may simply not have been detected --
#' but it is the precondition several downstream functions depend on.
#'
#' \strong{Anomalies} come from \code{\link{tales_anomalies}}, and are the
#' reason this is a \code{summary()} rather than part of printing: they are
#' worth computing, and too slow to compute every time an object is echoed.
#'
#' @param object A \code{\link{tales}} object.
#' @param ... Unused.
#' @return An object of class \code{summary.tales}: a list of the measures,
#'   so they can be used as well as read.
#' @method summary tales
#' @export
#' @family tales objects
summary.tales <- function(object, ...) {
  arrays <- unique(object$array_id)
  has <- function(col) col %in% names(object)

  repeats <- if (has("domain_type")) {
    tapply(object$domain_type == "repeat", object$array_id, sum)
  } else NULL

  complete <- if (has("domain_type")) {
    tapply(object$domain_type, object$array_id,
           function(v) all(c("N-terminus", "C-terminus") %in% v))
  } else NULL

  structure(
    list(
      class_name = class(object)[1],
      n_arrays = length(arrays),
      n_parts = nrow(object),
      n_distinct_domains = if (has("dom_code")) length(unique(object$dom_code)) else NULL,
      n_distinct_by_type = if (has("dom_code") && has("domain_type")) {
        vapply(split(object$dom_code, object$domain_type),
               function(v) length(unique(v)), integer(1))
      } else NULL,
      # repeats only, and NA is a missing RVD rather than a kind of one
      n_distinct_rvds = if (has("rvd")) {
        v <- if (has("domain_type")) object$rvd[object$domain_type == "repeat"] else object$rvd
        length(unique(v[!is.na(v)]))
      } else NULL,
      repeats_per_array = if (!is.null(repeats) && length(repeats)) {
        c(min = min(repeats), median = stats::median(repeats), max = max(repeats))
      } else NULL,
      n_complete = if (!is.null(complete)) sum(complete) else NULL,
      n_seqnames = if (has("seqnames")) length(unique(object$seqnames)) else NULL,
      n_groups = if (has("group")) length(unique(object$group)) else NULL,
      namespace = tales_namespace(object),
      anomalies = tales_anomalies(object)
    ),
    class = "summary.tales")
}


#' @param x A \code{summary.tales} object.
#' @rdname summary.tales
#' @method print summary.tales
#' @export
print.summary.tales <- function(x, ...) {
  cat(cli::format_inline("{.cls {x$class_name}} summary"), "\n", sep = "")
  .summary_lines(list(
    list("arrays / parts", paste(x$n_arrays, "/", x$n_parts)),
    list("distinct domains", if (!is.null(x$n_distinct_domains))
      paste0(x$n_distinct_domains, " of ", x$n_parts, " parts",
             if (!is.null(x$n_distinct_by_type))
               paste0("  (", paste(names(x$n_distinct_by_type),
                                   x$n_distinct_by_type, sep = " ",
                                   collapse = ", "), ")"))),
    list("distinct RVDs", x$n_distinct_rvds),
    list("repeats per array", if (!is.null(x$repeats_per_array))
      sprintf("min %g   median %g   max %g", x$repeats_per_array[["min"]],
              x$repeats_per_array[["median"]], x$repeats_per_array[["max"]])),
    list("arrays with both termini", if (!is.null(x$n_complete))
      paste0(x$n_complete, " of ", x$n_arrays)),
    list("source sequences", x$n_seqnames),
    list("groups", x$n_groups),
    list("dom_code namespace", if (!is.null(x$namespace) && nzchar(x$namespace)) x$namespace),
    list("anomalies", if (nrow(x$anomalies) == 0L) "none" else
      paste0(nrow(x$anomalies), " in ",
             length(unique(x$anomalies$array_id)), " array(s)"))
  ))
  if (nrow(x$anomalies) > 0L) {
    byCheck <- table(x$anomalies$check)
    .summary_lines(lapply(names(byCheck), function(k) list(paste0("  ", k), byCheck[[k]])))
    cat("  ", cli::format_inline("Inspect with {.fn tales_anomalies}."), "\n", sep = "")
  }
  invisible(x)
}


#' Summarise a tales_msa object
#'
#' @description
#' What an alignment looks like, in numbers: how gappy it is, and how much the
#' arrays actually agree.
#'
#' @details
#' \strong{Columns with no consensus} is the measure worth having. A column
#' has no consensus when no element is strictly more common than the rest --
#' see \code{\link{tales_consensus}} -- so the count says how much of the
#' alignment is agreement and how much is three different answers.
#'
#' It is reported per layer, and the two usually differ in a way that means
#' something. Repeats that are distinct proteins can carry the same RVD, so a
#' column can lack a \code{dom_code} consensus while having a clear
#' \code{rvd} one: the arrays disagree about the repeat and agree about the
#' base it binds. The reverse is rarer and more interesting.
#'
#' @param object A \code{\link{tales_msa}} object.
#' @param ... Unused.
#' @return An object of class \code{summary.tales_msa}.
#' @method summary tales_msa
#' @export
#' @family tales objects
summary.tales_msa <- function(object, ...) {
  # dom_code first, matching the order print() previews them in
  layers <- intersect(c("dom_code", "rvd"), names(object))
  width <- tales_width(object)
  arrays <- unique(object$array_id)

  perLayer <- lapply(layers, function(lyr) {
    m <- as.matrix(object, value = lyr)
    consensus <- tales_consensus(m)
    matches <- tales_consensus_match(m, long = FALSE)
    list(layer = lyr,
         n_no_consensus = sum(is.na(consensus)),
         n_matching = sum(matches, na.rm = TRUE),
         n_scorable = sum(!is.na(matches)))
  })
  names(perLayer) <- layers

  m1 <- if (length(layers)) as.matrix(object, value = layers[1]) else NULL

  structure(
    list(
      class_name = class(object)[1],
      n_arrays = length(arrays),
      width = width,
      n_gaps = if (!is.null(m1)) sum(is.na(m1)) else NULL,
      n_cells = if (!is.null(m1)) length(m1) else NULL,
      n_ungapped_columns = if (!is.null(m1)) sum(colSums(is.na(m1)) == 0L) else NULL,
      layers = perLayer,
      namespace = tales_namespace(object),
      anomalies = tales_anomalies(object)
    ),
    class = "summary.tales_msa")
}


#' @param x A \code{summary.tales_msa} object.
#' @rdname summary.tales_msa
#' @method print summary.tales_msa
#' @export
print.summary.tales_msa <- function(x, ...) {
  cat(cli::format_inline("{.cls {x$class_name}} summary"), "\n", sep = "")
  pairs <- list(
    list("arrays / width", paste(x$n_arrays, "/", x$width)),
    list("gaps", if (!is.null(x$n_gaps))
      sprintf("%d of %d cells (%.0f%%)", x$n_gaps, x$n_cells,
              100 * x$n_gaps / x$n_cells)),
    list("columns with no gap", if (!is.null(x$n_ungapped_columns))
      paste0(x$n_ungapped_columns, " of ", x$width))
  )
  for (l in x$layers) {
    pairs <- c(pairs, list(
      list(paste0("no consensus (", l$layer, ")"),
           paste0(l$n_no_consensus, " of ", x$width, " columns")),
      list(paste0("matching consensus (", l$layer, ")"),
           paste0(l$n_matching, " of ", l$n_scorable, " scorable cells"))))
  }
  pairs <- c(pairs, list(
    list("dom_code namespace", if (!is.null(x$namespace) && nzchar(x$namespace)) x$namespace),
    list("anomalies", if (nrow(x$anomalies) == 0L) "none" else nrow(x$anomalies))))
  .summary_lines(pairs)
  invisible(x)
}
