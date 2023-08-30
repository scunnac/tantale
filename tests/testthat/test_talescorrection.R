

test_that("corTales output a tibble of expected shape when returnCorrectionsTble is TRUE", {
  t <- corTales(uncorrectedAssemblyPath = "/home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/BAI3-1-1.fa",
                correctedAssemblyPath = tempfile(), returnCorrectionsTble = TRUE,
                condaBinPath = "/home/cunnac/bin/miniconda3/condabin/conda")
  expect_setequal(dim(t), c(63, 4))
})

test_that("corTales output a the path of an existing file when returnCorrectionsTble is FALSE", {
  f <- corTales(uncorrectedAssemblyPath = "/home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/BAI3-1-1.fa",
                correctedAssemblyPath = tempfile(tmpdir = "~"), returnCorrectionsTble = FALSE,
                condaBinPath = "/home/cunnac/bin/miniconda3/condabin/conda")
  expect_true(fs::file_exists(f))
  unlink(f)
})
