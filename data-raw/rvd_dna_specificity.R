# Builds the exported `rvd_dna_specificity` dataset from QueTAL FuncTAL's own
# table (inst/legacy/QueTAL_v1.1/FuncTAL/Info/2014mat18), values unchanged --
# only a header and column names added. See ?rvd_dna_specificity and
# dev/restructuring-notes.md §12b for what it feeds (tales_compare_functal()).

path <- system.file("legacy", "QueTAL_v1.1", "FuncTAL", "Info", "2014mat18",
                    package = "tantale", mustWork = TRUE)
rvd_dna_specificity <- readr::read_tsv(
  path, col_names = c("rvd", "A", "C", "G", "T"), show_col_types = FALSE,
  # "NA" is an RVD (Asn-Ala); readr's default would read it as missing.
  na = character()
)

usethis::use_data(rvd_dna_specificity, overwrite = TRUE)
