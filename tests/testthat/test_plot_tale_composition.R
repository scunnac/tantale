# plot_tale_composition() had no test and was broken: it calls mutate() and
# ggplot() unqualified, and neither dplyr nor ggplot2 was imported into the
# package namespace, so it only worked if the *user* had attached them.
# R CMD check's "no visible global function definition" NOTE was pointing at
# real breakage, not the usual NSE false positive.

test_that("plot_tale_composition() returns a ggplot", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  p <- plot_tale_composition(out$tale_parts)
  expect_s3_class(p, "ggplot")
})

test_that("plot_tale_composition() works with only the package attached", {
  # The regression this guards: unqualified calls resolving through the user's
  # search path rather than the package namespace.
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  expect_no_error(tantale::plot_tale_composition(out$tale_parts))
})
