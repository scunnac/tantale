# tell_tales() gives up at three points, and each must say why rather than
# letting a later stage fail on an empty table. The third of these used to be
# missing: the run carried on and died inside Bioconductor with "Rle of type
# 'NULL' is not supported", naming nothing the caller could act on.

subject <- function() {
  system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
              package = "tantale", mustWork = TRUE)
}

test_that("no hmmer hit at all is reported, not crashed on", {
  fasta <- tempfile(fileext = ".fa")
  set.seed(1)
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(paste(sample(Biostrings::DNA_BASES, 10000, replace = TRUE),
                                   collapse = "")),
    filepath = fasta)
  expect_warning(tell_tales(subject_file = fasta, output_dir = tempfile()),
                 regexp = "found no TALE cds hit")
})

test_that("a score threshold that rejects everything is reported", {
  out <- tempfile()
  expect_warning(
    suppressMessages(tell_tales(subject_file = subject(), output_dir = out,
                                nterm_min_score = 1e6,
                                repeat_min_score = 1e6,
                                cterm_min_score = 1e6)),
    regexp = "No record remains after filtering")
})

test_that("min_domain_hits rejecting every subject sequence is reported", {
  # Regression, restructuring-notes.md 8.0: this filter had no guard.
  out <- tempfile()
  expect_warning(
    suppressMessages(tell_tales(subject_file = subject(), output_dir = out,
                                min_domain_hits = 1000)),
    regexp = "No subject sequence carries more than")
})

test_that("giving up still returns the output directory, invisibly", {
  out <- tempfile()
  expect_invisible(
    res <- suppressWarnings(suppressMessages(
      tell_tales(subject_file = subject(), output_dir = out,
                 min_domain_hits = 1000))))
  expect_identical(res, out)
})
