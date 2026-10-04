#' RVD-to-DNA-binding-specificity weights
#'
#' One row per RVD, giving its relative binding weight for each of the four
#' bases. This is the table \code{\link{tales_compare_functal}} stacks, one
#' row per repeat, to build an array's position weight matrix.
#'
#' @format A tibble with 404 rows and 5 columns:
#' \describe{
#'   \item{rvd}{The two-letter RVD code. \code{"N*"} and \code{"H*"} are
#'     RVDs whose residue 13 is missing; \code{"OO"} is position 0, the base
#'     just before the first repeat's target, where TALEs prefer a T;
#'     \code{"XX"} is the flat fallback row used for an unrecognised RVD.}
#'   \item{A, C, G, T}{Relative binding weight for that base. The weights
#'     do not sum to one: \code{\link{tales_compare_functal}} hands
#'     them to \code{\link[universalmotif:create_motif]{create_motif}} as
#'     raw counts (\code{type = "PCM"}), which normalises internally.}
#' }
#'
#' @details
#' Taken verbatim from QueTAL FuncTAL's own table (\code{Info/2014mat18} in
#' FuncTAL 1.1, by A. L. Pérez-Quintero, redistributed with permission): the
#' values are unchanged, only a header and column names added. Conceptually related to
#' but distinct from the internal \code{rvdSimDf} used by
#' \code{\link{tales_align}}'s RVD scoring: that one is a *derived*
#' RVD-vs-RVD similarity (a correlation between two RVDs' base-preference
#' profiles, itself built from TALVEZ's much smaller 17-RVD \code{mat1}),
#' used to score repeat *substitutions* during sequence alignment. This one
#' is the *raw*, per-RVD base preference itself, one level upstream, used to
#' build a whole array's binding-specificity model for comparison against
#' another array's. The two tables serve different consumers and cannot
#' replace each other.
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


#' Curated TALE annotations for ten published *Xanthomonas oryzae* genomes
#'
#' The TALEs of ten complete, published genomes, named as the literature
#' names them and annotated by hand. It is the reference this package's own
#' discovery can be checked against: three of the ten genomes, MAI1, BAI3
#' and PXO86, are the ones [tantale_genome()] installs, so a
#' [tell_tales()] run on them can be compared array by array with what is
#' recorded here.
#'
#' @format A tibble with 128 rows and 10 columns. One row is one TALE gene.
#' \describe{
#'   \item{strain}{The strain, e.g. \code{"PXO99A"}. Ten of them.}
#'   \item{label, tal_name}{The TALE's name in the literature, capitalised
#'     (\code{"Tal4a"}) and as written in the source publication
#'     (\code{"tal4a"}). \code{strain} and \code{label} together identify
#'     a row: a label repeats across strains, never within one.
#'     \code{tal_name} is missing where the source gives no name.}
#'   \item{annotale_group}{The class AnnoTALE assigns the TALE
#'     (\code{"TalAH30"}), where it was recorded.}
#'   \item{replicon_id, genome_id}{The RefSeq replicon and the GenBank
#'     assembly the TALE comes from, e.g. \code{"NC_010717"} and
#'     \code{"GCA_000019585.2"}. Present for every row: this is what lets a
#'     reader fetch the sequence.}
#'   \item{pubmed}{PubMed identifier(s) of the publication(s) reporting the
#'     genome or the TALE, comma-separated where there is more than one.}
#'   \item{trunc_tale}{\code{TRUE} for the nine truncTALEs, TALEs whose
#'     C-terminal region is naturally shortened so that they no longer
#'     activate transcription. See the truncTALE article for what this does
#'     to discovery and to frameshift correction.}
#'   \item{rvd_seq}{The repeat-variable diresidues in array order,
#'     dash-separated (\code{"NI-HD-NG-..."}), the form
#'     \code{\link{tales_predict_targets}} takes. A \emph{lowercase} RVD
#'     (\code{"ng"}, \code{"n*"}) marks a repeat of non-standard length,
#'     the usual convention; 16 of the 128 arrays carry at least one, so
#'     do not upper-case this column before comparing it with anything.}
#'   \item{unusual_feature}{A note, for the nine arrays that carry one, on
#'     what makes the gene atypical: repeats of non-standard length,
#'     deletions, duplications, a premature stop.}
#' }
#'
#' @details
#' The table is sparse where the sources are: \code{tal_name} is absent for
#' 19 rows, \code{annotale_group} for 54, and \code{pubmed} for 19. The
#' gaps are the state of the curation, not placeholders to be filled by a
#' guess.
#'
#' Identifiers that name an array within a single run are deliberately not
#' included: the \code{\link{tell_tales}} array id and AnnoTALE's
#' \code{tempTALE} name both depend on the run that produced them, so
#' neither would survive a rerun. Matching a row to a fresh
#' \code{\link{tell_tales}} result means matching on what the TALE is,
#' its \code{rvd_seq}, rather than on what some run called it.
#'
#' No column of DisTAL groups is included. The working file carried one, but
#' how its values had been computed was not recorded, so they could not be
#' vouched for. \code{\link{tales_compare_distal}} followed by
#' \code{\link{tales_group_hclust}} computes such groups from
#' \code{rvd_seq}.
#'
#' @source
#' Compiled by Bao Tram Vi during the work that became this package, from
#' the genome records and publications cited in each row; published here for
#' the first time. The genome sequences themselves are the GenBank
#' assemblies named in \code{genome_id}, at
#' <https://www.ncbi.nlm.nih.gov/datasets/genome/>.
#'
#' @family TALE discovery
#' @examples
#' # the TALEs of one strain
#' subset(tale_annotations, strain == "PXO86",
#'        select = c("label", "annotale_group", "trunc_tale"))
#'
#' # the naturally truncated ones, across all ten genomes
#' table(tale_annotations$strain, tale_annotations$trunc_tale)
#'
#' # what makes the atypical arrays atypical
#' subset(tale_annotations, !is.na(unusual_feature),
#'        select = c("strain", "label", "unusual_feature"))
"tale_annotations"
