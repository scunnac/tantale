##### TALE relatedness by predicted DNA-binding specificity (FuncTAL) #####
# A from-scratch reimplementation of QueTAL's FuncTAL, not a port (ledger
# §12b). The original (inst/tools/QueTAL_v1.1/FuncTAL/FuncTAL_v.1.1.pl) could
# not be kept working: it needs Bio::Perl, which BioPerl dropped in its 1.7
# reorganisation. Rather than patch around that, the one thing FuncTAL does
# that nothing else in the package does -- compare TALEs by the DNA sequence
# their repeats are predicted to bind, not by repeat sequence identity -- is
# rebuilt on `universalmotif`, an actively maintained Bioconductor package
# for exactly this kind of comparison.
#
# This is a deliberate divergence, not a bug: `universalmotif::compare_motifs()`
# computes Pearson correlation per aligned column and combines the column
# scores (`score.strat`); FuncTAL flattened the whole padded, overlapping
# region -- both positions and the four bases at each -- into one vector and
# took a single correlation over that. These are different statistics, not
# two settings of the same one, verified empirically on real fixture data
# before writing this file, not assumed. Scores and any downstream tree from
# `tales_compare_functal()` will not match a FuncTAL run on the same input.


#' Build one array's position weight matrix from its RVD sequence
#'
#' One column per repeat, in order; each column is that repeat's RVD's row
#' from \code{\link{rvd_dna_specificity}} (or \code{"XX"}'s flat row, for an
#' RVD not in it -- deliberately uninformative, so an unrecognised RVD does
#' not pull a comparison toward false similarity).
#' @noRd
.functal_pwm <- function(rvd_seq) {
  spec <- rvd_dna_specificity
  rows <- match(rvd_seq, spec$rvd, nomatch = match("XX", spec$rvd))
  m <- t(as.matrix(spec[rows, c("A", "C", "G", "T")]))
  rownames(m) <- c("A", "C", "G", "T")
  m
}

# Metrics compare_motifs() reports as similarities (higher = more similar);
# everything else it reports as a distance already (closer to zero = more
# similar) -- see ?compare_motifs. Needed to build a `dissim` column that
# means the same thing regardless of which `method` was used.
.FUNCTAL_SIMILARITY_METHODS <- c("PCC", "WPCC", "SW", "ALLR", "ALLR_LL", "BHAT")


#' Compare TALEs by predicted DNA-binding specificity (FuncTAL)
#'
#' @description
#' Quantifies how TALE arrays relate by the DNA sequence their repeats are
#' predicted to bind, rather than by repeat sequence identity
#' (\code{\link{tales_compare_distal}}). Each array's repeats are turned into
#' a position weight matrix (PWM) over the RVD-to-base specificity code, and
#' PWMs are compared pairwise with \code{\link[universalmotif]{compare_motifs}}.
#'
#' @details
#' This is a reimplementation, not a port, of QueTAL's FuncTAL: the original
#' Perl tool could not be kept working (it needs \code{Bio::Perl}, dropped by
#' BioPerl's 1.7 reorganisation; see \code{dev/restructuring-notes.md} §12b),
#' so the comparison itself was rebuilt on \code{universalmotif} rather than
#' patched. \strong{Results diverge from the original FuncTAL tool, and this
#' is by design, not an approximation to be improved away.} FuncTAL scored
#' two RVD arrays by flattening their entire padded, overlapping alignment
#' (positions and bases together) into one vector and taking a single
#' Pearson correlation. \code{compare_motifs()} instead correlates matched
#' columns individually and combines the column scores (\code{score.strat}).
#' Verified empirically to disagree on real data before writing this
#' function, not assumed to differ only in magnitude.
#'
#' Only \code{rvd}, in repeat order, drives the comparison (via
#' \code{\link{tales_rvd_strings}}, which drops the two termini by default --
#' DNA-binding specificity is a property of the repeat region, and FuncTAL
#' itself never scored termini either, since its RVD extraction only ever
#' found repeats). An RVD absent from \code{\link{rvd_dna_specificity}} is
#' scored with a flat, uninformative row rather than dropped, so an unusual
#' RVD costs a comparison specificity rather than an error.
#'
#' Only a handful of \code{\link[universalmotif]{compare_motifs}}'s many
#' options are exposed here, chosen for what actually varies across TALE
#' arrays (their differing repeat counts) rather than mirrored wholesale; see
#' \code{?compare_motifs} for the rest, several of which are worth revisiting
#' -- ledger §12b lists them as follow-ups.
#'
#' @param x A \code{\link{tales}} object carrying an \code{rvd} column.
#' @param method One of \code{compare_motifs()}'s comparison metrics
#'   (\code{"PCC"} default). \code{"PCC"}, \code{"WPCC"}, \code{"SW"},
#'   \code{"ALLR"}, \code{"ALLR_LL"} and \code{"BHAT"} are similarities
#'   (higher is more similar); the rest are already distances. Either kind
#'   is turned into a proper \code{dissim} in the result.
#' @param tryRC Also score each pair's reverse-complement alignment and keep
#'   the better one. \code{FALSE} by default: an RVD array's specificity
#'   code has a fixed reading direction (5' to 3' target, N- to C-terminal
#'   repeat order), so comparing against a reverse complement asks a
#'   different, narrower biological question -- do these two TALEs target
#'   opposite strands of overlapping sites -- worth asking on purpose, not
#'   folded silently into every comparison.
#' @param min.overlap Minimum aligned width to accept, as `compare_motifs()`
#'   defines it. Defaults to \code{1} (any overlap at all), not
#'   \code{compare_motifs()}'s own default of \code{6} -- generic TF motifs
#'   are usually longer than 6 positions, but a TALE array legitimately has
#'   as few as half a dozen repeats, so the upstream default would silently
#'   refuse to compare some real, short arrays.
#' @param normalise.scores Penalise alignments that leave much of either
#'   motif unaligned. \code{TRUE} by default: TALE array lengths vary widely
#'   (see \code{min.overlap}), so an unnormalised score would favour matches
#'   that are mostly one array hanging off the end of a much longer one.
#' @param score.strat How \code{compare_motifs()} combines per-column scores
#'   into one alignment score. \code{"a.mean"} (its own default): a sum
#'   would scale with array length, favouring long-vs-long comparisons
#'   for no biological reason.
#' @param nthreads Passed to \code{compare_motifs()}.
#' @return A \code{\link{tale_distances}} object -- interchangeable with
#'   \code{\link{tales_compare_distal}}'s, so it can be handed directly to
#'   \code{\link{tales_group_hclust}}/\code{\link{tales_group_kmedoids}}.
#'   No tree is built here, matching \code{tales_compare_distal()}'s own
#'   division of labour: compare here, cluster/tree there.
#' @seealso [tales_compare_distal()], comparing by repeat sequence instead.
#' @export
#' @family pairwise distances
#' @examples
#' x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
#'                                      package = "tantale"))
#' fn <- tales_compare_functal(x)
#' fn
tales_compare_functal <- function(x, method = "PCC", tryRC = FALSE,
                                  min.overlap = 1, normalise.scores = TRUE,
                                  score.strat = "a.mean", nthreads = 1) {
  if (!is_tales(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }

  rvd_strings <- tales_rvd_strings(x)
  # tales_rvd_strings() silently drops an array with zero repeats (all
  # termini) rather than erroring -- it has nothing left to render for that
  # one array, but plenty for the object as a whole. A comparison must not
  # silently proceed over fewer arrays than it was given.
  missing_arrays <- setdiff(unique(x$array_id), names(rvd_strings))
  if (length(missing_arrays) > 0L) {
    cli::cli_abort(
      c("Every array needs at least one repeat to build a specificity PWM.",
        "x" = "No repeats found for {.val {utils::head(missing_arrays, 5)}}."),
      class = c("tantale_error_functal_empty", "tantale_error")
    )
  }
  rvd_list <- strsplit(as.character(rvd_strings), "-", fixed = TRUE)
  names(rvd_list) <- names(rvd_strings)

  motifs <- lapply(names(rvd_list), function(id) {
    universalmotif::create_motif(.functal_pwm(rvd_list[[id]]),
                                 alphabet = "DNA", type = "PCM", name = id)
  })

  sim <- universalmotif::compare_motifs(
    motifs, method = method, tryRC = tryRC, min.overlap = min.overlap,
    normalise.scores = normalise.scores, score.strat = score.strat,
    nthreads = nthreads
  )

  ids <- rownames(sim)
  long <- data.frame(
    id1 = rep(ids, times = length(ids)),
    id2 = rep(ids, each = length(ids)),
    score = as.vector(sim)
  )
  long$dissim <- if (method %in% .FUNCTAL_SIMILARITY_METHODS) {
    1 - long$score
  } else {
    long$score
  }

  tale_distances(long[c("id1", "id2", "dissim")])
}
