# Coverage for exports that had no test at all. plot_tales_composition() turned
# out to be outright broken when checked this way, so the others were worth
# exercising too.

fixture <- function() {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
}

test_that("tale_parts_to_rvd() renders one RVD string per array", {
  d <- fixture()
  out <- tale_parts_to_rvd(d$tale_parts)
  expect_s4_class(out, "BStringSet")
  expect_length(out, length(unique(d$tale_parts$arrayID)))
})

test_that("rvd_only = TRUE drops every anchor code, not just NTERM/CTERM", {
  # Regression: the filter hardcoded c("NTERM", "CTERM") and so kept "XXXXX",
  # the sentinel for a terminus detected but not identified. One sequence in
  # the reference fixture carries it.
  d <- fixture()
  out <- as.character(tale_parts_to_rvd(d$tale_parts, rvd_only = TRUE))
  for (code in tales_anchor_codes()) {
    expect_false(any(grepl(code, out, fixed = TRUE)),
                 label = paste0("anchor code ", code, " survived rvd_only = TRUE"))
  }
})

test_that("rvd_only = FALSE keeps the termini", {
  d <- fixture()
  out <- as.character(tale_parts_to_rvd(d$tale_parts, rvd_only = FALSE))
  expect_true(any(grepl("NTERM", out, fixed = TRUE)))
})

test_that("diagnose_tale_parts() reports arrays with missing sequences", {
  d <- fixture()
  expect_s3_class(diagnose_tale_parts(d$tale_parts), "data.frame")
})

test_that("repeat_to_rvd_map_distalr() maps each repeat code to exactly one RVD", {
  d <- fixture()
  m <- repeat_to_rvd_map_distalr(d$tale_parts)
  expect_named(m, c("repeatID", "RVD"))
  # the mapping must be a function: one RVD per repeat code, never two
  expect_identical(anyDuplicated(m$repeatID), 0L)
  expect_identical(nrow(m), length(unique(d$tale_parts$domCode)))
})

test_that("validate_pairwise_distances() accepts a valid object and rejects a broken one", {
  d <- fixture()
  x <- domain_distances(d$repeat.similarity)
  expect_s3_class(validate_pairwise_distances(x), "pairwise_distances")
  broken <- x
  broken$id1 <- NULL
  expect_error(validate_pairwise_distances(broken), class = "tantale_error")
})

test_that("repeat_to_rvd_map() builds the code -> RVD mapping", {
  repeat_vecs <- list(a = c("r1", "r2"), b = c("r2", "r3"))
  rvd_vecs    <- list(a = c("NI", "HD"), b = c("HD", "NG"))
  m <- repeat_to_rvd_map(repeat_vecs, rvd_vecs)
  expect_named(m, c("repeatID", "RVD"))
  expect_identical(nrow(m), 3L)          # r1, r2, r3 -- r2 shared, counted once
  expect_identical(anyDuplicated(m$repeatID), 0L)
})

test_that("repeat_to_rvd_map() rejects a repeat code with two different RVDs", {
  # Ledger 2 keeps this function on the retirement list but insists its
  # assertion be migrated first: it is the only place the package enforces
  # that a repeat code maps to exactly one RVD.
  repeat_vecs <- list(a = c("r1", "r2"), b = c("r2", "r3"))
  ambiguous   <- list(a = c("NI", "HD"), b = c("NN", "NG"))  # r2 -> HD and NN
  expect_error(repeat_to_rvd_map(repeat_vecs, ambiguous))
})
