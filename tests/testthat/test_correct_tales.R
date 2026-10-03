# The genome itself no longer ships with the package (ledger §50); these
# tests run on its TALE loci, cut out with 3 kb on either side by
# data-raw/make_bai311_tale_loci.R. On the whole genome TALEcorrection makes
# 70 corrections; 69 of them fall in these windows, at the same positions
# with the same edits, and the 70th (an insertion at 585387) lies outside
# every TALE array tell_tales() finds.
excerpt <- function() test_path("data_for_tests", "bai311_tale_loci.fa")

test_that("correct_tales output a tibble of expected shape when return_corrections is TRUE", {
  t <- correct_tales(uncorrected_path = excerpt(),
                corrected_path = tempfile(), return_corrections = TRUE,
                conda_bin = "auto")
  # On the whole genome the count was 70, not the pre-fix 63: correct_tales()
  # was feeding the repeat and C-terminus nHMMER outputs to
  # TALEcorrection.jar's r=/c= flags swapped relative to what the tool itself
  # documents those flags as wanting (confirmed by running it and reading its
  # own printed usage) -- fixed 2026-09-21, and the correction count changed
  # as a real consequence, not a regression.
  expect_setequal(dim(t), c(69, 4))
})

test_that("correct_tales output a the path of an existing file when return_corrections is FALSE", {
  f <- correct_tales(uncorrected_path = excerpt(),
                corrected_path = tempfile(fileext = ".fa"), return_corrections = FALSE,
                conda_bin = "auto")
  expect_true(fs::file_exists(f))
  unlink(f)
})
