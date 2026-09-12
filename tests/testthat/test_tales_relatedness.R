# Tests for tales_relatedness(). See dev/class-design.md §1.2 and
# dev/restructuring-notes.md §1 (what the returned list shrank to, and why).
#
# These run the real pipeline (pairwise protein alignment + ARLEM), so they are
# slower than the rest of the suite.

example_tales <- function() {
  suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
}

relatedness_once <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- suppressWarnings(suppressMessages(
        tales_relatedness(example_tales())
      ))
    }
    cached
  }
})


test_that("tales_relatedness() returns exactly three typed objects", {
  res <- relatedness_once()
  expect_named(res, c("tales", "repeat_sim", "tale_sim"))
  expect_s3_class(res$tales, "tales")
  expect_s3_class(res$repeat_sim, "repeat_sim")
  expect_s3_class(res$tale_sim, "tale_sim")
})

test_that("the three dropped slots are gone", {
  res <- relatedness_once()
  expect_false(any(c("repeats.code", "coded.repeats.str", "repeats.cluster")
                   %in% names(res)))
})

test_that("dom_code is minted onto the returned tales", {
  res <- relatedness_once()
  x <- example_tales()
  expect_true("dom_code" %in% names(res$tales))
  expect_false("dom_code" %in% names(x))
  expect_identical(nrow(res$tales), nrow(x))
})

test_that("dom_code is in one-to-one correspondence with aa_seq", {
  # the invariant validate_tales() enforces, here on freshly minted codes
  res <- relatedness_once()
  tp <- res$tales
  expect_equal(dplyr::n_distinct(tp$dom_code), dplyr::n_distinct(tp$aa_seq))
  expect_silent(validate_tales(tp))
})


#### Namespace stamping ####

test_that("the tales and repeat_sim carry the same namespace", {
  res <- relatedness_once()
  expect_type(tales_namespace(res$tales), "character")
  expect_identical(tales_namespace(res$repeat_sim), tales_namespace(res$tales))
})

test_that("tale_sim is deliberately unstamped, being keyed by array_id", {
  res <- relatedness_once()
  expect_null(tales_namespace(res$tale_sim))
})

test_that("the namespace is the hash of the part set, so a rerun matches", {
  res <- relatedness_once()
  expect_identical(
    tales_namespace(res$tales),
    tantale:::.tales_dom_code_namespace(example_tales()$aa_seq)
  )
})


#### The similarity tables ####

test_that("repeat_sim is keyed by dom_code and is square", {
  res <- relatedness_once()
  expect_setequal(unique(res$repeat_sim$id1), unique(res$tales$dom_code))
  expect_silent(sim_assert_square(res$repeat_sim))
})

test_that("tale_sim is keyed by array_id and is square", {
  res <- relatedness_once()
  expect_setequal(unique(res$tale_sim$id1), unique(res$tales$array_id))
  expect_silent(sim_assert_square(res$tale_sim))
})

test_that("both tables put their id columns first, in order", {
  res <- relatedness_once()
  expect_identical(names(res$repeat_sim)[1:3], c("id1", "id2", "sim"))
  expect_identical(names(res$tale_sim)[1:3], c("id1", "id2", "sim"))
})

test_that("self-similarity is 100", {
  res <- relatedness_once()
  selfRep <- res$repeat_sim$sim[res$repeat_sim$id1 == res$repeat_sim$id2]
  expect_true(all(selfRep == 100))
})


#### Preconditions and guards ####

test_that("a plain data frame is refused", {
  expect_error(tales_relatedness(data.frame(a = 1)),
               class = "tantale_error_tales_type")
})

test_that("a tales without aa_seq is refused", {
  x <- as_tales(test_path("data_for_tests", "tellTaleExampleOutput",
                          "rvdSequences.fas"), sep = "-")
  expect_error(tales_relatedness(x),
               class = "tantale_error_relatedness_no_aa")
})

test_that("re-minting over an existing dom_code warns", {
  res <- relatedness_once()
  expect_warning(
    suppressMessages(tales_relatedness(res$tales)),
    class = "tantale_warning_relatedness_remint"
  )
})


#### The deprecated distalr() still behaves as before ####

test_that("distalr() is deprecated but returns its original six slots", {
  w <- NULL
  out <- suppressWarnings(suppressMessages(withCallingHandlers(
    distalr(suppressWarnings(
      tantale:::.tale_parts(test_path("data_for_tests", "tellTaleExampleOutput"))
    )),
    deprecatedWarning = function(c) w <<- conditionMessage(c)
  )))
  expect_match(w, "tales_relatedness")
  expect_named(out, c("tale_parts", "repeats.code", "coded.repeats.str",
                      "repeat.similarity", "tal.similarity", "repeats.cluster"))
  # legacy column vocabulary preserved for existing callers
  expect_true("domCode" %in% names(out$tale_parts))
  expect_setequal(names(out$repeat.similarity), c("RepU1", "RepU2", "Dissim", "Sim"))
})
