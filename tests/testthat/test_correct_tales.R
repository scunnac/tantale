

test_that("correct_tales output a tibble of expected shape when return_corrections is TRUE", {
  t <- correct_tales(uncorrected_path = "/home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/BAI3-1-1.fa",
                corrected_path = tempfile(), return_corrections = TRUE,
                conda_bin = "auto")
  # 70, not the pre-fix 63: correct_tales() was feeding the repeat and
  # C-terminus nHMMER outputs to TALEcorrection.jar's r=/c= flags swapped
  # relative to what the tool itself documents those flags as wanting
  # (confirmed by running it and reading its own printed usage) -- fixed
  # 2026-09-21, and the correction count changed as a real consequence,
  # not a regression.
  expect_setequal(dim(t), c(70, 4))
})

test_that("correct_tales output a the path of an existing file when return_corrections is FALSE", {
  f <- correct_tales(uncorrected_path = "/home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/BAI3-1-1.fa",
                corrected_path = tempfile(tmpdir = "~"), return_corrections = FALSE,
                conda_bin = "auto")
  expect_true(fs::file_exists(f))
  unlink(f)
})
