# Errors raised by the package carry classed conditions, so they can be
# asserted by class rather than by matching message text. See ledger 9.5:
# before the cli conversion many of these were bare stop() calls whose
# condition message was the empty string.

fixture_tales <- function() {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  tales(out$tale_parts)
}

test_that("every error the package raises carries the tantale_error parent", {
  # a representative sample across the different subsystems
  expect_error(tales(data.frame(nope = 1)), class = "tantale_error")
  expect_error(pairwise_distances(data.frame(id1 = "a")), class = "tantale_error")
  expect_error(tales_compare(data.frame(a = 1)), class = "tantale_error")
})


#### tales_compare() guards ####

test_that("an unknown aln_method is rejected by name", {
  x <- fixture_tales()
  # This guard was unreachable before: logger::log_errors() && stop(...) short
  # circuited, so the user saw a complaint about calling handlers instead.
  expect_error(tales_compare(x, aln_method = "not_a_method"),
               class = "tantale_error_aln_method")
  expect_error(tales_compare(x, aln_method = "not_a_method"), "not_a_method")
})

test_that("tales_compare() refuses a tales with no aa_seq", {
  x <- fixture_tales()
  x$aa_seq <- NULL
  expect_error(tales_compare(x), class = "tantale_error_compare_no_aa")
})


#### the RVD conversion internals repaired in ledger 2 ####

rvd_align_fixture <- function() {
  matrix(c("NI", "HD", NA, "NI", NA, "NG"), nrow = 2, byrow = TRUE,
         dimnames = list(c("a", "b"), NULL))
}

test_that(".rvd_to_repeat_align() rejects a row with no repeat vector", {
  m <- rvd_align_fixture()
  expect_error(tantale:::.rvd_to_repeat_align(m, list(a = c("r1", "r2"))),
               class = "tantale_error_rvd_repeat_missing")
})

test_that(".rvd_to_repeat_align() rejects a non-gap/repeat count mismatch", {
  m <- rvd_align_fixture()
  bad <- list(a = c("r1", "r2", "r3"), b = c("r3", "r4"))
  expect_error(tantale:::.rvd_to_repeat_align(m, bad),
               class = "tantale_error_rvd_repeat_length")
})

test_that(".rvd_to_repeat_align() back-maps positionally when counts agree", {
  m <- rvd_align_fixture()
  ok <- list(a = c("r1", "r2"), b = c("r3", "r4"))
  out <- tantale:::.rvd_to_repeat_align(m, ok)
  expect_equal(dim(out), dim(m))
  expect_identical(rownames(out), rownames(m))
  expect_identical(as.vector(t(out)), c("r1", "r2", NA, "r3", NA, "r4"))
})

test_that(".rvd_to_match_align() runs against the internal rvdSimDf", {
  # Its default argument used to be tantale::rvdSimDf, which errors because
  # rvdSimDf lives in sysdata.rda -- the function could never run (ledger 2).
  m <- rvd_align_fixture()
  m[is.na(m)] <- "NG"
  out <- tantale:::.rvd_to_match_align(m)
  expect_true(is.matrix(out))
  expect_equal(dim(out), dim(m))
  expect_type(out, "double")
})


#### name clashes between legacy and canonical vocabularies ####

test_that("a table carrying both spellings of a column is refused", {
  df <- data.frame(id1 = c("a", "b"), id2 = c("b", "a"),
                   dissim = c(0, 5), Dissim = c(0, 5))
  expect_error(pairwise_distances(df), class = "tantale_error_distances_name_clash")
})

test_that("tales() also refuses a legacy/canonical name clash", {
  df <- data.frame(arrayID = "a1", array_id = "a1", positionInArray = 1L,
                   rvd = "NI", stringsAsFactors = FALSE)
  expect_error(tales(df), class = "tantale_error_tales_name_clash")
})

test_that("the name-clash messages are formattable", {
  # Regression: both bullets carried a {?s} plural marker with no quantity to
  # count, so cli failed with "Cannot pluralize without a quantity" and the
  # classed condition was replaced by a bare simpleError.
  d1 <- data.frame(id1 = c("a", "b"), id2 = c("b", "a"),
                   dissim = c(0, 5), Dissim = c(0, 5))
  expect_error(pairwise_distances(d1), "already present")
  d2 <- data.frame(arrayID = "a1", array_id = "a1", positionInArray = 1L,
                   rvd = "NI", stringsAsFactors = FALSE)
  expect_error(tales(d2), "already present")
})


#### Empty-input guards ####

test_that("projecting a tales with nothing left to render errors", {
  x <- fixture_tales()
  # keep only the terminus sentinels, then ask for repeats only
  anchors <- x[x$rvd %in% tales_anchor_codes(), ]
  expect_error(tales_rvd_strings(anchors, rvd_only = TRUE),
               class = "tantale_error_projection_empty")
})

test_that("aligning an empty tales errors", {
  x <- fixture_tales()
  expect_error(tales_align(x[0, ], residue_col = "rvd"),
               class = "tantale_error_msa_empty")
})

test_that("a part with no amino acid sequence is reported with its array", {
  x <- fixture_tales()
  bad <- unique(x$array_id)[1]
  x$aa_seq[x$array_id == bad][1] <- NA_character_
  expect_error(suppressWarnings(tales_compare(x)),
               class = "tantale_error_parts_no_aa")
  # the offending array must be named, not just counted
  expect_error(suppressWarnings(tales_compare(x)), bad, fixed = TRUE)
})

test_that("an empty-string aa_seq is reported too, not just NA", {
  # The guard covers is.na() | == "", but the list of affected arrays used to
  # collect only the NA ones, so an empty string produced a message naming
  # nothing at all.
  x <- fixture_tales()
  bad <- unique(x$array_id)[1]
  x$aa_seq[x$array_id == bad][1] <- ""
  expect_error(suppressWarnings(tales_compare(x)), bad, fixed = TRUE)
})
