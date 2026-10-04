# Builds the exported `tale_annotations` dataset from Bao Tram Vi's curated
# table of the TALEs of ten published Xanthomonas oryzae genomes
# (data-raw/ref_tale_annotation.tsv, her working file kept here unchanged).
# See ?tale_annotations and dev/restructuring-notes.md §57.
#
# Three columns are dropped and one renamed:
#  - telltale2_id and annotale_id name an array within one particular run
#    (`ROI_00001`, `tempTALE1`): neither is reproducible, so neither belongs
#    in a reference table (maintainer, 2026-10-04). strain + label is the
#    key, and it is unique over all 128 rows. annotale_group is kept: those
#    are AnnoTALE's own class names (`TalAH30`), which do carry across runs.
#  - distal_group is dropped: it is empty in every row of this file, and how
#    the values in the sibling file were computed is not recorded, so they
#    cannot be vouched for (ledger §57, Q128). tales_compare_distal() and
#    tales_group_hclust() recompute such groups from `rvd_seq`.
#  - truncTALE -> trunc_tale: columns are snake_case (dev/CLAUDE.md).

path <- file.path("data-raw", "ref_tale_annotation.tsv")

tale_annotations <- readr::read_tsv(
  path,
  col_types = readr::cols(.default = readr::col_character(),
                          truncTALE = readr::col_logical())
) |>
  dplyr::select(-"distal_group", -"telltale2_id", -"annotale_id") |>
  dplyr::rename(trunc_tale = "truncTALE") |>
  dplyr::relocate("strain", "label", "tal_name") |>
  dplyr::arrange(.data$strain, .data$label)

stopifnot(
  nrow(tale_annotations) == 128L,
  dplyr::n_distinct(tale_annotations$strain) == 10L,
  # strain + label is the key: labels repeat across strains, not within one
  !anyDuplicated(tale_annotations[c("strain", "label")]),
  !anyNA(tale_annotations$label),
  # every row says which genome it came from, which is what lets a reader
  # fetch the sequence (rOpenSci asks a dataset to name its source)
  !anyNA(tale_annotations$genome_id),
  !anyNA(tale_annotations$replicon_id)
)

usethis::use_data(tale_annotations, overwrite = TRUE)
