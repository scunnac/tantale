

test_that("correct_tales output a tibble of expected shape when return_corrections is TRUE", {
  t <- correct_tales(uncorrected_path = "/home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/BAI3-1-1.fa",
                corrected_path = tempfile(), return_corrections = TRUE,
                conda_bin = "auto")
  expect_setequal(dim(t), c(63, 4))
})

test_that("correct_tales output a the path of an existing file when return_corrections is FALSE", {
  f <- correct_tales(uncorrected_path = "/home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/BAI3-1-1.fa",
                corrected_path = tempfile(tmpdir = "~"), return_corrections = FALSE,
                conda_bin = "auto")
  expect_true(fs::file_exists(f))
  unlink(f)
})
