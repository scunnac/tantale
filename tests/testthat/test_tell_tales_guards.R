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
    regexp = "No subject sequence carries at least")
})

test_that("giving up still returns the output directory, invisibly", {
  out <- tempfile()
  expect_invisible(
    res <- suppressWarnings(suppressMessages(
      tell_tales(subject_file = subject(), output_dir = out,
                 min_domain_hits = 1000))))
  expect_identical(res, out)
})


#### min_array_length ####

# The fixture's four arrays carry 26, 14, 26 and 22 repeats.

test_that("min_array_length defaults to keeping everything", {
  out <- tempfile()
  suppressWarnings(suppressMessages(tell_tales(subject_file = subject(), output_dir = out)))
  expect_length(list.dirs(file.path(out, "annotale"), recursive = FALSE), 4L)
})

test_that("min_array_length drops only the arrays below it", {
  # 20 is between the 14-repeat array and the next smallest, at 22
  out <- tempfile()
  suppressWarnings(suppressMessages(
    tell_tales(subject_file = subject(), output_dir = out, min_array_length = 20)))
  expect_length(list.dirs(file.path(out, "annotale"), recursive = FALSE), 3L)

  report <- readr::read_tsv(file.path(out, "arrayReport.tsv"),
                            show_col_types = FALSE, progress = FALSE)
  expect_equal(nrow(report), 3L)
})

test_that("min_array_length counts repeats, not all hits", {
  # The 14-repeat array also carries two termini, so it has 16 hits. A
  # threshold of 15 must drop it: were all hits counted, 16 would survive.
  out <- tempfile()
  suppressWarnings(suppressMessages(
    tell_tales(subject_file = subject(), output_dir = out, min_array_length = 15)))
  expect_length(list.dirs(file.path(out, "annotale"), recursive = FALSE), 3L)
})

test_that("min_array_length rejecting every array is reported, not crashed on", {
  out <- tempfile()
  expect_warning(
    suppressMessages(tell_tales(subject_file = subject(), output_dir = out,
                                min_array_length = 500)),
    regexp = "No TALE array has at least")
})

test_that("min_domain_hits is inclusive, as documented", {
  # talRegion6 carries 24 hits. At 24 it must be kept; the filter used to be
  # strictly greater-than, which dropped it.
  out <- tempfile()
  suppressWarnings(suppressMessages(
    tell_tales(subject_file = subject(), output_dir = out, min_domain_hits = 24)))
  report <- readr::read_tsv(file.path(out, "hitsReport.tsv"),
                            show_col_types = FALSE, progress = FALSE)
  expect_true("talRegion6" %in% report$seqnames)
})
