# plot() had no test and was broken: it calls mutate() and
# ggplot() unqualified, and neither dplyr nor ggplot2 was imported into the
# package namespace, so it only worked if the *user* had attached them.
# R CMD check's "no visible global function definition" NOTE was pointing at
# real breakage, not the usual NSE false positive.

test_that("plot() on a tales returns a ggplot", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  p <- plot(tales(out$tale_parts))
  expect_s3_class(p, "ggplot")
})

test_that("plot() works with only the package attached", {
  # The regression this guards: unqualified calls resolving through the user's
  # search path rather than the package namespace.
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  expect_no_error(plot(tantale::tales(out$tale_parts)))
})

test_that("seqnames is optional: the facet is added only when present", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(out$tale_parts)
  withFacet <- plot(x)
  withoutFacet <- plot(x[, setdiff(names(x), "seqnames")])
  expect_s3_class(withFacet, "ggplot")
  expect_s3_class(withoutFacet, "ggplot")
  expect_s3_class(withFacet$facet, "FacetGrid")
  expect_s3_class(withoutFacet$facet, "FacetNull")
})

test_that("it reports the columns it needs rather than failing obscurely", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(out$tale_parts)
  expect_error(plot(x[, setdiff(names(x), "aa_seq")]),
               class = "tantale_error_projection_column")
})

test_that("a legacy tale_parts data frame is still accepted", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  expect_s3_class(suppressWarnings(plot(tales(out$tale_parts))), "ggplot")
})

test_that("plot() dispatches to the composition plot for a tales", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(out$tale_parts)
  p <- plot(x)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "Overview of TALE composition by genome")
})

test_that("a tales_msa still dispatches to plot.tales_msa, not plot.tales", {
  # tales_msa inherits from tales, so the more specific method must win. If
  # plot.tales ever shadowed it, an alignment would silently render as a
  # composition plot.
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(out$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:3], ]
  msa <- suppressWarnings(suppressMessages(tales_align(sub, residue_col = "rvd")))
  p <- suppressWarnings(suppressMessages(plot(msa)))
  expect_false(identical(p$labels$title, "Overview of TALE composition by genome"))
})

test_that("position = 'alignment' lays parts out on the alignment coordinate", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(out$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:6], ]
  msa <- suppressWarnings(suppressMessages(tales_align(sub, residue_col = "rvd")))
  back <- suppressWarnings(suppressMessages(as_tales(msa)))

  arrayLayout <- plot(back, position = "array")
  alignLayout <- plot(back, position = "alignment")
  # the aligned layout spans the alignment width; the array one only the longest array
  expect_equal(max(alignLayout$data$.x), tales_width(msa))
  expect_equal(max(arrayLayout$data$.x), max(back$position_in_array))
  expect_gt(max(alignLayout$data$.x), max(arrayLayout$data$.x))
})

test_that("position = 'alignment' needs an aligned object", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(out$tale_parts)
  expect_error(plot(x, position = "alignment"),
               class = "tantale_error_projection_column")
})

test_that("the default layout is unchanged", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(out$tale_parts)
  expect_identical(plot(x)$data$.x, x$position_in_array)
})
