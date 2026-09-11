

test_that("correcTales output a tibble of expected shape when returnCorrectionsTble is TRUE", {
  t <- correcTales(uncorrectedAssemblyPath = "/home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/BAI3-1-1.fa",
                correctedAssemblyPath = tempfile(), returnCorrectionsTble = TRUE,
                condaBinPath = "auto")
  expect_setequal(dim(t), c(63, 4))
})

test_that("correcTales output a the path of an existing file when returnCorrectionsTble is FALSE", {
  f <- correcTales(uncorrectedAssemblyPath = "/home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/BAI3-1-1.fa",
                correctedAssemblyPath = tempfile(tmpdir = "~"), returnCorrectionsTble = FALSE,
                condaBinPath = "auto")
  expect_true(fs::file_exists(f))
  unlink(f)
})
