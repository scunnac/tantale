# Build the tell_tales() output directories used by the tests and examples.
#
#   inst/extdata/tellTaleExampleOutput          full tell_tales() run on the
#   tests/testthat/data_for_tests/example_output  shipped BAI3 sample (4 TALEs)
#
# and three copies of the latter, each with one inconsistency planted in its
# ROI_00004, for the tests of .tale_parts():
#
#   err_missing_nterm  N-terminus removed from both AnnoTALE part files: the
#                      array has no N-terminus part, positions start at the
#                      first repeat
#   err_missing_dna    N-terminus removed from the protein parts only: the two
#                      files disagree, the array is left out
#   err_array_count    last RVD removed from TALE_RVDs.fasta: RVDs and repeat
#                      parts no longer correspond, an error
#
# The err_* copies keep only what .tale_parts() reads.
#
# Rerun after any change to what tell_tales() writes.
#
# Run with:  Rscript data-raw/make_telltale_test_fixtures.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

TEST_DIR <- "tests/testthat/data_for_tests"
ROI      <- "ROI_00004"

subject <- system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                       package = "tantale", mustWork = TRUE)
run <- file.path(tempdir(), "tell_tales_example")
unlink(run, recursive = TRUE)
tell_tales(subject_file = subject, output_dir = run)

replace_dir <- function(from, to) {
  unlink(to, recursive = TRUE)
  dir.create(to, recursive = TRUE)
  file.copy(list.files(from, full.names = TRUE), to, recursive = TRUE)
}
replace_dir(run, "inst/extdata/tellTaleExampleOutput")
replace_dir(run, file.path(TEST_DIR, "example_output"))

# a copy holding only what .tale_parts() reads
minimal_copy <- function(name) {
  to <- file.path(TEST_DIR, name)
  unlink(to, recursive = TRUE)
  dir.create(file.path(to, "annotale"), recursive = TRUE)
  file.copy(file.path(run, c("array_report.tsv", "hits_report.tsv")), to)
  for (roi in list.files(file.path(run, "annotale"))) {
    dir.create(file.path(to, "annotale", roi))
    parts <- list.files(file.path(run, "annotale", roi), "^TALE_.*\\.fasta$", full.names = TRUE)
    file.copy(parts, file.path(to, "annotale", roi))
  }
  file.path(to, "annotale", ROI)
}

drop_nterm <- function(file) {
  seqs <- Biostrings::readBStringSet(file)
  Biostrings::writeXStringSet(seqs[!grepl("N-terminus", names(seqs), fixed = TRUE)], file)
}

roi <- minimal_copy("err_missing_nterm")
drop_nterm(file.path(roi, "TALE_Protein_parts.fasta"))
drop_nterm(file.path(roi, "TALE_DNA_parts.fasta"))

roi <- minimal_copy("err_missing_dna")
drop_nterm(file.path(roi, "TALE_Protein_parts.fasta"))

roi <- minimal_copy("err_array_count")
rvdFile <- file.path(roi, "TALE_RVDs.fasta")
rvds <- Biostrings::readBStringSet(rvdFile)
rvds[[1]] <- Biostrings::BString(sub("-[^-]+$", "", as.character(rvds[[1]])))
Biostrings::writeXStringSet(rvds, rvdFile)
