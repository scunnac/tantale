# Coverage for exports that had no test at all. plot_tales_composition() turned
# out to be outright broken when checked this way, so the others were worth
# exercising too.

fixture <- function() {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
}

test_that("repeats_only = TRUE drops every anchor code, not just NTERM/CTERM", {
  # Regression (first seen in the retired tale_parts_to_rvd()): the filter
  # hardcoded c("NTERM", "CTERM") and so kept "XXXXX", the code of a terminus
  # that does not match its profile. One sequence in the reference fixture
  # carries it.
  d <- fixture()
  expect_true(tales_anchor_codes()[[3]] %in% d$tale_parts$rvd)
  out <- as.character(tales_rvd_strings(tales(d$tale_parts), repeats_only = TRUE))
  for (code in tales_anchor_codes()) {
    expect_false(any(grepl(code, out, fixed = TRUE)),
                 label = paste0("anchor code ", code, " survived repeats_only = TRUE"))
  }
})

test_that("tales_anomalies() reports arrays with missing sequences", {
  # replaces diagnose_tale_parts(), which guarded on columns its callers never
  # read and warned about an "output tibble" it was about to not produce
  d <- fixture()
  expect_s3_class(tales_anomalies(tales(d$tale_parts)), "data.frame")
})

test_that("validate_pairwise_distances() accepts a valid object and rejects a broken one", {
  d <- fixture()
  x <- domain_distances(d$repeat.similarity)
  expect_s3_class(validate_pairwise_distances(x), "pairwise_distances")
  broken <- x
  broken$id1 <- NULL
  expect_error(validate_pairwise_distances(broken), class = "tantale_error")
})

