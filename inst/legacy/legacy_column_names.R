# Retired: the legacy column-name converters of tales() and
# pairwise_distances()
#
# They renamed camelCase columns (arrayID, positionInArray, ...) and the old
# distance-table spellings (TAL1/TAL2, RepU1/RepU2, Sim, Dissim, arlemScore,
# maxLength, normArlemScore) to the current schema. No code in R/ produced
# those spellings any more and no shipped data carried them; they served only
# tables saved by old versions. Retired before 1.0.0 at the maintainer's
# request (ledger §37). To read such a table, rename its columns first; the
# vectors below map each old name to the current one.
#
# .as_mafft_score_table() lost its RepU1/RepU2/Sim branch at the same time.
#
# Moved out of R/ on 2026-10-01.

# Legacy camelCase -> target snake_case. Kept so the class can be used against
# output of the not-yet-renamed pipeline; see class-design.md §1.1.
# `seqnames` is deliberately absent: it keeps its Bioconductor spelling.
TALES_LEGACY_NAMES <- c(
  arrayID         = "array_id",
  positionInArray = "position_in_array",
  positionInCrd   = "position_in_crd",
  domainType      = "domain_type",
  aaSeq           = "aa_seq",
  dnaSeq          = "dna_seq",
  domCode         = "dom_code",
  sourceDirectory = "source_directory"
)

#' Rename legacy camelCase columns to the target schema
#' @param x A data frame.
#' @return The same data frame with recognised legacy names replaced.
#' @noRd
.tales_rename_legacy <- function(x) {
  hit <- intersect(names(x), names(TALES_LEGACY_NAMES))
  if (length(hit) == 0L) return(x)
  clash <- intersect(unname(TALES_LEGACY_NAMES[hit]), names(x))
  if (length(clash) > 0L) {
    cli::cli_abort(
      c("Cannot rename legacy columns: the target name{?s} {.field {clash}} {?is/are} already present.",
        "i" = "Drop or rename {cli::qty(clash)}the duplicate{?s} before calling {.fn tales}."),
      class = c("tantale_error_tales_name_clash", "tantale_error")
    )
  }
  names(x)[match(hit, names(x))] <- unname(TALES_LEGACY_NAMES[hit])
  x
}

# Legacy spellings -> canonical. Both entity vocabularies collapse onto the
# same id columns; renaming by name also fixes repeat.similarity listing its
# ids in the order RepU2, RepU1.
PAIRWISE_DISTANCES_LEGACY_NAMES <- c(
  TAL1           = "id1",
  TAL2           = "id2",
  RepU1          = "id1",
  RepU2          = "id2",
  Sim            = "sim",
  Dissim         = "dissim",
  arlemScore     = "arlem_score",
  maxLength      = "max_length",
  normArlemScore = "norm_arlem_score"
)

#' @noRd
.pairwise_distances_rename_legacy <- function(x) {
  hit <- intersect(names(x), names(PAIRWISE_DISTANCES_LEGACY_NAMES))
  if (length(hit) == 0L) return(x)
  target <- unname(PAIRWISE_DISTANCES_LEGACY_NAMES[hit])
  clash <- intersect(target, names(x))
  if (length(clash) > 0L) {
    cli::cli_abort(
      c("Cannot rename legacy columns: the target name{?s} {.field {clash}} {?is/are} already present.",
        "i" = "Drop or rename {cli::qty(clash)}the duplicate{?s} first."),
      class = c("tantale_error_distances_name_clash", "tantale_error")
    )
  }
  names(x)[match(hit, names(x))] <- target
  x
}
