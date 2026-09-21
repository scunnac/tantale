#' RVD-to-DNA-binding-specificity weights
#'
#' One row per RVD, giving its relative binding weight for each of the four
#' bases. This is the table \code{\link{tales_compare_functal}} stacks, one
#' row per repeat, to build an array's position weight matrix.
#'
#' @format A tibble with 404 rows and 5 columns:
#' \describe{
#'   \item{rvd}{The two-letter RVD code (or the anchor codes \code{"XX"},
#'     an unrecognised RVD's fallback, and \code{"H*"}/\code{"N*"}).}
#'   \item{A, C, G, T}{Relative binding weight for that base. Not
#'     normalised to sum to one -- \code{\link{tales_compare_functal}} hands
#'     them to \code{\link[universalmotif:create_motif]{create_motif}} as
#'     raw counts (\code{type = "PCM"}), which normalises internally.}
#' }
#'
#' @details
#' Ported verbatim from QueTAL FuncTAL's own table
#' (\code{inst/tools/QueTAL_v1.1/FuncTAL/Info/2014mat18}) -- the values are
#' unchanged, only a header and column names added. Conceptually related to
#' but distinct from the internal \code{rvdSimDf} used by
#' \code{\link{tales_align}}'s RVD scoring: that one is a *derived*
#' RVD-vs-RVD similarity (a correlation between two RVDs' base-preference
#' profiles, itself built from TALVEZ's much smaller 17-RVD \code{mat1}),
#' used to score repeat *substitutions* during sequence alignment. This one
#' is the *raw*, per-RVD base preference itself, one level upstream, used to
#' build a whole array's binding-specificity model for comparison against
#' another array's -- a different consumer, and not interchangeable with
#' \code{rvdSimDf}.
#'
#' @references
#' Pérez-Quintero A.L. et al. (2015). QueTAL: a suite of tools to classify
#' and compare TAL effectors functionally and phylogenetically.
#' \emph{Frontiers in Plant Science} \strong{6}, 545.
#' \doi{10.3389/fpls.2015.00545}
#'
#' @seealso \code{\link{tales_compare_functal}}, its only consumer.
#' @family pairwise distances
#' @examples
#' rvd_dna_specificity[rvd_dna_specificity$rvd %in% c("HD", "NI", "NG", "NN"), ]
"rvd_dna_specificity"
