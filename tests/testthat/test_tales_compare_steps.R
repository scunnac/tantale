# The three steps tales_compare() composes (ledger 8.5b).
#
# The golden baseline already pins that the composition produces what the
# monolith produced. These test the properties that only matter now that the
# steps are separately callable -- chiefly that they refuse the ways they can
# be misused, since exporting step 1 makes the run-dependence of dom_code
# part of the public API.

# A tales as it looks *before* comparison. The stored fixture already
# carries a dom_code from the run that produced it, so it is dropped here --
# otherwise the "needs codes assigned first" tests below would be handed
# codes and pass for the wrong reason.
cmp_fixture <- function(n = 3) {
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- suppressWarnings(tales(d$tale_parts))
  x <- x[x$array_id %in% unique(x$array_id)[seq_len(n)], ]
  x[setdiff(names(x), "dom_code")]
}


#### step 1: assigning codes ####

test_that("tales_assign_domain_codes() gives one code per distinct sequence", {
  x <- cmp_fixture()
  out <- tales_assign_domain_codes(x)

  expect_s3_class(out, "tales")
  expect_true("dom_code" %in% names(out))
  expect_false(anyNA(out$dom_code))
  # the defining property: code <-> sequence is one to one
  expect_equal(dplyr::n_distinct(out$dom_code), dplyr::n_distinct(out$aa_seq))
  bySeq <- tapply(out$dom_code, out$aa_seq, function(z) length(unique(z)))
  expect_true(all(bySeq == 1L))
})

test_that("tales_assign_domain_codes() codes termini as well as repeats", {
  # A dom_code names a distinct *domain*, not a repeat -- the termini are
  # parts like the repeats are, and get codes too.
  out <- tales_assign_domain_codes(cmp_fixture())
  coded <- out$dom_code[out$domain_type != "repeat"]
  expect_true(length(coded) > 0L)
  expect_false(anyNA(coded))
})

test_that("tales_assign_domain_codes() stamps a namespace", {
  out <- tales_assign_domain_codes(cmp_fixture())
  expect_true(nzchar(tales_namespace(out)))
})

test_that("codes are run-dependent, and the namespace says so", {
  # This is the hazard the export makes public, so it is pinned rather than
  # only documented: a different set of arrays is a different vocabulary.
  three <- tales_assign_domain_codes(cmp_fixture(3))
  two   <- tales_assign_domain_codes(cmp_fixture(2))
  expect_false(identical(tales_namespace(three), tales_namespace(two)))
})

test_that("tales_assign_domain_codes() refuses parts with no aa_seq", {
  x <- cmp_fixture()
  x$aa_seq[1] <- NA_character_
  expect_error(tales_assign_domain_codes(x),
               class = "tantale_error_parts_no_aa")
})

test_that("tales_assign_domain_codes() needs a tales with aa_seq", {
  x <- cmp_fixture()
  expect_error(tales_assign_domain_codes(tibble::as_tibble(x)),
               class = "tantale_error_tales_type")
  expect_error(tales_assign_domain_codes(x[setdiff(names(x), "aa_seq")]),
               class = "tantale_error_compare_no_aa")
})


#### step 2: domain distances ####

test_that("tales_domain_distances() compares distinct domains, not parts", {
  x <- tales_assign_domain_codes(cmp_fixture())
  dd <- suppressMessages(tales_domain_distances(x))

  expect_s3_class(dd, "domain_distances")
  n <- dplyr::n_distinct(x$dom_code)
  expect_equal(nrow(dd), n^2)
  expect_setequal(unique(dd$id1), unique(x$dom_code))
})

test_that("tales_domain_distances() carries the namespace through", {
  x <- tales_assign_domain_codes(cmp_fixture())
  dd <- suppressMessages(tales_domain_distances(x))
  expect_identical(tales_namespace(dd), tales_namespace(x))
})

test_that("tales_domain_distances() needs codes assigned first", {
  x <- cmp_fixture()
  expect_error(suppressMessages(tales_domain_distances(x)),
               class = "tantale_error_projection_column")
})

test_that("tales_domain_distances() rejects an unknown aln_method", {
  x <- tales_assign_domain_codes(cmp_fixture())
  expect_error(suppressMessages(tales_domain_distances(x, aln_method = "nope")),
               class = "tantale_error_aln_method")
})


#### step 3: TALE distances ####

test_that("tales_tale_distances() returns one score per ordered pair", {
  x <- tales_assign_domain_codes(cmp_fixture())
  dd <- suppressMessages(tales_domain_distances(x))
  td <- suppressMessages(tales_tale_distances(x, dd))

  expect_s3_class(td, "tale_distances")
  n <- dplyr::n_distinct(x$array_id)
  expect_equal(nrow(td), n^2)
  expect_setequal(unique(td$id1), unique(x$array_id))
})

test_that("tales_tale_distances() refuses distances from another run", {
  # The reason both objects are stamped: step 3 indexes domains by code, so
  # distances keyed by another call's codes would silently compare the wrong
  # domains rather than fail.
  x <- tales_assign_domain_codes(cmp_fixture(3))
  other <- tales_assign_domain_codes(cmp_fixture(2))
  ddOther <- suppressMessages(tales_domain_distances(other))

  expect_error(suppressMessages(tales_tale_distances(x, ddOther)),
               class = "tantale_error_namespace_mismatch")
})

test_that("tales_tale_distances() needs codes assigned first", {
  x <- tales_assign_domain_codes(cmp_fixture())
  dd <- suppressMessages(tales_domain_distances(x))
  expect_error(suppressMessages(tales_tale_distances(cmp_fixture(), dd)),
               class = "tantale_error_projection_column")
})


#### the composition ####

test_that("the three steps reproduce tales_compare()", {
  # The claim the decomposition rests on. The golden baseline pins the
  # values; this pins that the wrapper is exactly its parts.
  x <- cmp_fixture()

  whole <- suppressWarnings(suppressMessages(tales_compare(x)))

  coded <- tales_assign_domain_codes(x)
  dd <- suppressMessages(tales_domain_distances(coded))
  td <- suppressMessages(tales_tale_distances(coded, dd))

  expect_equal(as.data.frame(whole$tales), as.data.frame(coded))
  expect_equal(as.data.frame(whole$domain_distances), as.data.frame(dd))
  expect_equal(as.data.frame(whole$tale_distances), as.data.frame(td))
})
