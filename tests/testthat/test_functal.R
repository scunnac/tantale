
test_that("functal() defaults output_dir to the current working directory", {
  expect_identical(formals(functal)$output_dir, quote(getwd()))
})

test_that("functal() errors when the bundled functal perl script cannot run", {
  # FuncTAL_v.1.1.pl currently requires List::MoreUtils and Bio::Perl, which
  # are not part of the 'tantale' conda environment (see project notes on
  # inst/tools/QueTAL_v1.1/functal). Until that's resolved, functal() is
  # expected to fail; this test documents that known, current failure so a
  # silent regression (e.g. the exit code check being dropped again) would
  # be caught, and it will start failing loudly on its own once the conda
  # environment gains the missing Perl dependencies.
  skip_on_cran()
  output_dir <- tempfile("functal")
  dir.create(output_dir)
  expect_error(
    functal(tal_file = test_path("data_for_tests", "sample_funcTAL_input.txt"),
           output_dir = output_dir),
    "functal failed"
  )
  expect_length(list.files(output_dir), 0)
})
