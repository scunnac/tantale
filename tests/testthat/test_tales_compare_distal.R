# Tests for tales_compare_distal(). See dev/class-design.md §1.2 and
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
        tales_compare_distal(example_tales())
      ))
    }
    cached
  }
})


test_that("tales_compare_distal() returns exactly three typed objects", {
  res <- relatedness_once()
  expect_named(res, c("tales", "domain_distances", "tale_distances"))
  expect_s3_class(res$tales, "tales")
  expect_s3_class(res$domain_distances, "domain_distances")
  expect_s3_class(res$tale_distances, "tale_distances")
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

test_that("the tales and domain_distances carry the same namespace", {
  res <- relatedness_once()
  expect_type(tales_namespace(res$tales), "character")
  expect_identical(tales_namespace(res$domain_distances), tales_namespace(res$tales))
})

test_that("tale_distances is deliberately unstamped, being keyed by array_id", {
  res <- relatedness_once()
  expect_null(tales_namespace(res$tale_distances))
})

test_that("the namespace is the hash of the part set, so a rerun matches", {
  res <- relatedness_once()
  expect_identical(
    tales_namespace(res$tales),
    tantale:::.tales_dom_code_namespace(example_tales()$aa_seq)
  )
})


#### The similarity tables ####

test_that("domain_distances is keyed by dom_code and is square", {
  res <- relatedness_once()
  expect_setequal(unique(res$domain_distances$id1), unique(res$tales$dom_code))
  expect_silent(distances_assert_square(res$domain_distances))
})

test_that("tale_distances is keyed by array_id and is square", {
  res <- relatedness_once()
  expect_setequal(unique(res$tale_distances$id1), unique(res$tales$array_id))
  expect_silent(distances_assert_square(res$tale_distances))
})

test_that("both tables put their id columns first, in order", {
  res <- relatedness_once()
  expect_identical(names(res$domain_distances)[1:3], c("id1", "id2", "dissim"))
  expect_identical(names(res$tale_distances)[1:3], c("id1", "id2", "dissim"))
})

test_that("self-comparison is distance 0", {
  res <- relatedness_once()
  selfRep <- res$domain_distances$dissim[res$domain_distances$id1 == res$domain_distances$id2]
  expect_true(all(selfRep == 0))
})


#### Preconditions and guards ####

test_that("a plain data frame is refused", {
  expect_error(tales_compare_distal(data.frame(a = 1)),
               class = "tantale_error_tales_type")
})

test_that("a tales without aa_seq is refused", {
  x <- as_tales(test_path("data_for_tests", "tellTaleExampleOutput",
                          "rvd_sequences.fas"), sep = "-")
  expect_error(tales_compare_distal(x),
               class = "tantale_error_compare_no_aa")
})

test_that("re-minting over an existing dom_code warns", {
  res <- relatedness_once()
  expect_warning(
    suppressMessages(tales_compare_distal(res$tales)),
    class = "tantale_warning_relatedness_remint"
  )
})


#### aa_seq or dna_seq ####

test_that(".translate_parts() reproduces the stored aa_seq exactly", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  tp <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))$tale_parts
  expect_identical(tantale:::.translate_parts(tp[["dna_seq"]]), tp[["aa_seq"]])
})

test_that("translation needs no.init.codon: repeats start on CTG/TTG", {
  # The default forces the first codon to M. TALE repeats begin on alternative
  # start codons, so without the flag every repeat gains a spurious leading M.
  ctg <- "CTGACCCCGGAACAGGTG"          # L T P E Q V
  expect_identical(tantale:::.translate_parts(ctg), "LTPEQV")
  expect_false(identical(
    as.character(Biostrings::translate(Biostrings::DNAStringSet(ctg))), "LTPEQV"))
})

test_that("a trailing stop codon is stripped", {
  expect_identical(tantale:::.translate_parts("ATGGATTAA"), "MD")
})

test_that("an out-of-frame sequence is refused rather than silently truncated", {
  expect_error(tantale:::.translate_parts("ATGGA"),
               class = "tantale_error_translate_frame")
})

test_that("tales_compare_distal() falls back to dna_seq and keeps the derived column", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  noAa <- x[, setdiff(names(x), "aa_seq")]
  expect_warning(res <- suppressMessages(tales_compare_distal(noAa)),
                 class = "tantale_warning_translated_aa")
  expect_true("aa_seq" %in% names(res$tales))
})

test_that("the dna_seq path gives the same distances as the aa_seq path", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  viaAa  <- suppressWarnings(suppressMessages(tales_compare_distal(x)))
  viaDna <- suppressWarnings(suppressMessages(
    tales_compare_distal(x[, setdiff(names(x), "aa_seq")])))
  expect_equal(as.data.frame(viaDna$domain_distances),
               as.data.frame(viaAa$domain_distances))
  expect_equal(as.data.frame(viaDna$tale_distances),
               as.data.frame(viaAa$tale_distances))
})

test_that("neither aa_seq nor dna_seq is still an error", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  expect_error(tales_compare_distal(x[, setdiff(names(x), c("aa_seq", "dna_seq"))]),
               class = "tantale_error_compare_no_aa")
})
