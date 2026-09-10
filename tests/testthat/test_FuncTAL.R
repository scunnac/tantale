
test_that("FuncTAL() defaults outputDir to the current working directory", {
  expect_identical(formals(FuncTAL)$outputDir, quote(getwd()))
})

test_that("FuncTAL() errors when the bundled FuncTAL perl script cannot run", {
  # FuncTAL_v.1.1.pl currently requires List::MoreUtils and Bio::Perl, which
  # are not part of the 'tantale' conda environment (see project notes on
  # inst/tools/QueTAL_v1.1/FuncTAL). Until that's resolved, FuncTAL() is
  # expected to fail; this test documents that known, current failure so a
  # silent regression (e.g. the exit code check being dropped again) would
  # be caught, and it will start failing loudly on its own once the conda
  # environment gains the missing Perl dependencies.
  skip_on_cran()
  outDir <- tempfile("FuncTAL")
  dir.create(outDir)
  expect_error(
    FuncTAL(TALfile = test_path("data_for_tests", "sample_funcTAL_input.txt"),
           outputDir = outDir),
    "FuncTAL failed"
  )
  expect_length(list.files(outDir), 0)
})
