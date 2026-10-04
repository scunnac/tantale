# Builds the exported `tale_annotations` dataset from Bao Tram Vi's curated
# table of the TALEs of ten published Xanthomonas oryzae genomes
# (data-raw/ref_tale_annotation.tsv, her working file kept here unchanged).
# See ?tale_annotations and dev/restructuring-notes.md §57.
#
# Four columns of the working file are dropped:
#  - telltale2_id and annotale_id name an array within one particular run
#    (`ROI_00001`, `tempTALE1`): neither is reproducible, so neither belongs
#    in a reference table (maintainer, 2026-10-04). strain + label is the
#    key, and it is unique over all 128 rows.
#  - distal_group is empty in every row of this file, and how the values in
#    the sibling file were computed is not recorded, so they cannot be
#    vouched for (Q128). tales_compare_distal() and tales_group_hclust()
#    recompute such groups from `rvd_seq`.
#  - annotale_group is replaced by annotale_class, derived below. The
#    working file held values like "TalAH30", but the trailing number is
#    the member's index within the class on the day it was assigned: the
#    catalogue grows and members renumber, so all 73 of Tram's numbers had
#    already drifted (TalAH30 is TalAH25 today) while the class letters
#    agreed in every one. Only the class is kept (Q140).
#
# `truncTALE` keeps its spelling: it is the field's term for these TALEs,
# as `iTALE` is, and the snake_case rule is for R code, not for a name the
# literature writes that way (maintainer, 2026-10-04).

# AnnoTALE's class catalogue, from `AnnoTALEcli-1.5.jar loadAndView`,
# downloaded 2026-10-04 (17 minutes, 440 MB of output). Its
# List_of_classes.txt is kept here gzipped, so this script needs no
# network and the dataset can be rebuilt as it was. Format: per class, a
# header line, then one line per member -- an aligned RVD row ("--" for a
# gap), a tab, then "<TALE id> <species> <strain>[ (Pseudo)]".
.annotale_catalogue <- function(
    path = file.path("data-raw", "annotale_List_of_classes.txt.gz")) {
  rows <- grep("\t", readLines(path, warn = FALSE), value = TRUE, fixed = TRUE)
  parts <- strsplit(rows, "\t", fixed = TRUE)
  tag <- trimws(vapply(parts, `[`, "", 2))
  strain <- sub("\\s*\\(.*$", "", sub("^\\S+\\s+", "", tag))
  rvd <- vapply(parts, function(x) {
    v <- strsplit(trimws(x[[1]]), " +")[[1]]
    paste(v[v != "--"], collapse = "-")
  }, "")
  data.frame(
    # the id is "<class><member index>"; only the class survives a rerun
    annotale_class = sub("[0-9]+$", "", sub(" .*$", "", tag)),
    strain = sub("^(Xoo|Xoc|Xo|Xac|X)\\s+", "", strain),
    rvd = toupper(rvd),
    stringsAsFactors = FALSE
  )
}

tale_annotations <- readr::read_tsv(
  file.path("data-raw", "ref_tale_annotation.tsv"),
  col_types = readr::cols(.default = readr::col_character(),
                          truncTALE = readr::col_logical())
) |>
  dplyr::select(-"distal_group", -"telltale2_id", -"annotale_id",
                -"annotale_group") |>
  dplyr::relocate("strain", "label", "tal_name") |>
  dplyr::arrange(.data$strain, .data$label)

# One correction to the curated values (maintainer, Q148). MAI1 TalH was
# recorded with 12 repeats; the genome it came from has 16 at that locus.
# A tell_tales() run on the MAI1 genome tantale_genome() installs finds it
# as one array of 16 (ROI_00007, 2026-10-05), and AnnoTALE's catalogue
# records the same 16 for MAI1 (TalDN23). The 12 are an exact prefix: the
# last four repeats, NG-NI-NG-NI, had been lost. The source TSV is kept
# as Tram left it; the correction is made here, where it can be read.
.talh <- with(tale_annotations, strain == "MAI1" & label == "TalH")
stopifnot(sum(.talh) == 1L,
          tale_annotations$rvd_seq[.talh] == "NN-HD-NV-HD-NI-NG-NI-NN-NS-HD-HD-NI")
tale_annotations$rvd_seq[.talh] <-
  "NN-HD-NV-HD-NI-NG-NI-NN-NS-HD-HD-NI-NG-NI-NG-NI"

# Matched on strain plus the RVD string, case-insensitively: lowercase
# marks an aberrant repeat here and the catalogue upper-cases. The strain
# has to be part of the key -- a TALE shared between strains would
# otherwise take whichever copy the catalogue happened to list first, and
# with RVD alone 127 rows match but some to the wrong strain's copy.
catalogue <- .annotale_catalogue()
tale_annotations$annotale_class <- catalogue$annotale_class[
  match(paste(tale_annotations$strain, toupper(tale_annotations$rvd_seq)),
        paste(catalogue$strain, catalogue$rvd))
]
tale_annotations <- dplyr::relocate(tale_annotations, "annotale_class",
                                    .after = "tal_name")

stopifnot(
  nrow(tale_annotations) == 128L,
  dplyr::n_distinct(tale_annotations$strain) == 10L,
  # strain + label is the key: labels repeat across strains, not within one
  !anyDuplicated(tale_annotations[c("strain", "label")]),
  !anyNA(tale_annotations$label),
  # every row says which genome it came from, which is what lets a reader
  # fetch the sequence (rOpenSci asks a dataset to name its source)
  !anyNA(tale_annotations$genome_id),
  !anyNA(tale_annotations$replicon_id),
  # 127 of 128 carry a class; only PXO99A Tal7b does not, being the
  # 5-repeat allele its own unusual_feature describes (§57).
  sum(!is.na(tale_annotations$annotale_class)) == 127L,
  # MAI1 TalH, once corrected, matches MAI1's own catalogue entry
  tale_annotations$annotale_class[.talh] == "TalDN",
  # a class is letters only, with no member index
  all(grepl("^Tal[A-Z]+$", stats::na.omit(tale_annotations$annotale_class)))
)

usethis::use_data(tale_annotations, overwrite = TRUE)
