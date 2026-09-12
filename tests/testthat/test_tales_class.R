# Tests for the 'tales' class. See dev/class-design.md §2.

# A minimal valid object: two arrays, N-term + 2 repeats + C-term each.
minimal_tales_df <- function() {
  tibble::tibble(
    array_id = rep(c("a1", "a2"), each = 4),
    position_in_array = rep(1:4, times = 2),
    domain_type = rep(c("N-terminus", "repeat", "repeat", "C-terminus"), 2),
    position_in_crd = rep(c(NA, 1L, 2L, NA), 2),
    rvd = rep(c("NTERM", "NI", "HD", "CTERM"), 2),
    aa_seq = rep(c("MDP", "LTPA", "LTPD", "SIVA"), 2),
    dom_code = rep(c("1", "2", "3", "4"), 2),
    seqnames = rep(c("c1", "c2"), each = 4)
  )
}

test_that("tales() accepts a valid table and returns a tibble subclass", {
  x <- tales(minimal_tales_df())
  expect_s3_class(x, "tales")
  expect_s3_class(x, "tbl_df")
  expect_true(is_tales(x))
})

test_that("tales() renames legacy camelCase columns", {
  df <- minimal_tales_df()
  names(df)[names(df) == "array_id"] <- "arrayID"
  names(df)[names(df) == "position_in_array"] <- "positionInArray"
  names(df)[names(df) == "domain_type"] <- "domainType"
  names(df)[names(df) == "position_in_crd"] <- "positionInCrd"
  names(df)[names(df) == "aa_seq"] <- "aaSeq"
  names(df)[names(df) == "dom_code"] <- "domCode"

  x <- tales(df)
  expect_true(all(c("array_id", "position_in_array", "domain_type",
                    "position_in_crd", "aa_seq", "dom_code") %in% names(x)))
  expect_false(any(c("arrayID", "positionInArray") %in% names(x)))
})

test_that("tales() coerces a double position_in_array to integer", {
  df <- minimal_tales_df()
  df$position_in_array <- as.numeric(df$position_in_array)
  expect_type(tales(df)$position_in_array, "integer")
})

test_that("seqnames keeps its Bioconductor spelling", {
  x <- tales(minimal_tales_df())
  expect_true("seqnames" %in% names(x))
})


#### Column contract ####

test_that("a missing key column is an error", {
  df <- minimal_tales_df()
  df$array_id <- NULL
  expect_error(tales(df), class = "tantale_error_tales_missing_column")
})

test_that("at least one residue column is required", {
  df <- minimal_tales_df()
  df$rvd <- NULL
  df$dom_code <- NULL
  expect_error(tales(df), class = "tantale_error_tales_missing_column")
})

test_that("either residue column alone is enough", {
  df <- minimal_tales_df()
  expect_s3_class(tales(df[setdiff(names(df), "dom_code")]), "tales")
  expect_s3_class(tales(df[setdiff(names(df), c("rvd", "aa_seq"))]), "tales")
})

test_that("extra columns are preserved without complaint", {
  df <- minimal_tales_df()
  df$strain <- "X"
  expect_silent(x <- tales(df))
  expect_true("strain" %in% names(x))
})

test_that("a wrong column type is an error", {
  df <- minimal_tales_df()
  df$array_id <- seq_len(nrow(df))
  expect_error(tales(df), class = "tantale_error_tales_type")
})


#### Hard invariants ####

test_that("an empty tales is valid", {
  x <- tales(minimal_tales_df()[0, ])
  expect_s3_class(x, "tales")
  expect_equal(nrow(x), 0L)
})

test_that("a duplicated key is an error", {
  df <- minimal_tales_df()
  df$position_in_array[2] <- 1L
  expect_error(tales(df), class = "tantale_error_tales_duplicate_key")
})

test_that("NA or non-positive position_in_array is an error", {
  df <- minimal_tales_df()
  df$position_in_array[1] <- NA_integer_
  expect_error(tales(df), class = "tantale_error_tales_position")

  df2 <- minimal_tales_df()
  df2$position_in_array[1] <- 0L
  expect_error(tales(df2), class = "tantale_error_tales_position")
})

test_that("NA in a residue column is an error", {
  df <- minimal_tales_df()
  df$rvd[1] <- NA_character_
  expect_error(tales(df), class = "tantale_error_tales_na")
})


#### Conditional invariants ####

test_that("an unknown domain_type is an error", {
  df <- minimal_tales_df()
  df$domain_type[2] <- "middle"
  expect_error(tales(df), class = "tantale_error_tales_domain_type")
})

test_that("at most one terminus of each kind per array", {
  df <- minimal_tales_df()
  df$domain_type[3] <- "C-terminus"
  df$position_in_crd[3] <- NA_integer_
  expect_error(tales(df), class = "tantale_error_tales_terminus")
})

test_that("zero termini is allowed - a repeats-only subset stays valid", {
  x <- tales(minimal_tales_df())
  repeats_only <- dplyr::filter(x, domain_type == "repeat")
  expect_s3_class(repeats_only, "tales")
  expect_silent(validate_tales(repeats_only))
  # position_in_array no longer starts at 1, and that is fine
  expect_false(min(repeats_only$position_in_array) == 1L)
})

test_that("position_in_crd must be NA exactly on non-repeat parts", {
  df <- minimal_tales_df()
  df$position_in_crd[1] <- 1L
  expect_error(tales(df), class = "tantale_error_tales_crd")
})

test_that("position_in_crd must agree with position_in_array on order", {
  df <- minimal_tales_df()
  # swap the two repeats' crd positions: same set, wrong order
  df$position_in_crd[df$domain_type == "repeat"] <- c(2L, 1L, 2L, 1L)
  expect_error(tales(df), class = "tantale_error_tales_crd")
})

test_that("a shifted but correctly ordered position_in_crd is valid", {
  # This is what a repeats-only subset looks like, so it must not error:
  # the exact offset relation is a precondition, not an invariant.
  df <- minimal_tales_df()
  df$position_in_crd[df$domain_type == "repeat"] <- c(3L, 4L, 3L, 4L)
  expect_s3_class(tales(df), "tales")
})

test_that("tales_assert_complete() catches what the validator deliberately allows", {
  x <- tales(minimal_tales_df())
  expect_silent(tales_assert_complete(x))

  repeats_only <- dplyr::filter(x, domain_type == "repeat")
  expect_s3_class(repeats_only, "tales")                       # still valid
  expect_error(tales_assert_complete(repeats_only),            # but not complete
               class = "tantale_error_tales_incomplete")
})

test_that("the aa_seq <-> dom_code correspondence must be bijective", {
  df <- minimal_tales_df()
  df$dom_code[2] <- "99"          # same aa_seq in array 2 keeps code "2"
  expect_error(tales(df), class = "tantale_error_tales_dom_code")
})

test_that("seqnames must be constant within an array", {
  df <- minimal_tales_df()
  df$seqnames[1] <- "other"
  expect_error(tales(df), class = "tantale_error_tales_inconsistent")
})

test_that("an empty dna_seq warns but does not fail", {
  df <- minimal_tales_df()
  df$dna_seq <- c("ATG", rep("", 7))
  expect_warning(x <- tales(df), class = "tantale_warning_tales_missing_dna")
  expect_s3_class(x, "tales")
})


#### dplyr policy ####

test_that("row subsetting preserves the class and the invariants", {
  x <- tales(minimal_tales_df())
  expect_s3_class(dplyr::filter(x, array_id == "a1"), "tales")
  expect_s3_class(x[1:4, ], "tales")
  expect_s3_class(dplyr::arrange(x, dplyr::desc(position_in_array)), "tales")
  expect_s3_class(dplyr::slice(x, 1), "tales")
})

test_that("filtering to no rows is allowed", {
  x <- tales(minimal_tales_df())
  out <- dplyr::filter(x, array_id == "nope")
  expect_s3_class(out, "tales")
  expect_equal(nrow(out), 0L)
})

test_that("dropping a required column degrades silently to a tibble", {
  x <- tales(minimal_tales_df())
  out <- dplyr::select(x, -array_id)
  expect_false(is_tales(out))
  expect_s3_class(out, "tbl_df")
})

test_that("dropping one residue column keeps the class while the other remains", {
  x <- tales(minimal_tales_df())
  expect_s3_class(dplyr::select(x, -dom_code), "tales")
})

test_that("`[` degrades on column subsetting too", {
  # dplyr's dplyr_col_select() only calls dplyr_reconstruct() for plain
  # data.frame/data.table, so `[.tales` is what makes select() degrade.
  # This test pins that behaviour in case the dplyr internals change.
  x <- tales(minimal_tales_df())
  expect_s3_class(x[, c("array_id", "position_in_array", "rvd")], "tales")
  expect_false(is_tales(x[, c("position_in_array", "rvd")]))
  expect_false(is_tales(x["rvd"]))
})

test_that("mutate that collides the key is an error", {
  x <- tales(minimal_tales_df())
  expect_error(
    dplyr::mutate(x, position_in_array = 1L),
    class = "tantale_error_tales_duplicate_key"
  )
})

test_that("mutate that does not touch the key is fine", {
  x <- tales(minimal_tales_df())
  expect_s3_class(dplyr::mutate(x, note = "hello"), "tales")
})


#### Namespace tag ####

test_that("the dom_code namespace is carried through subsetting", {
  x <- tales(minimal_tales_df(), dom_code_namespace = "run-abc")
  expect_identical(tales_namespace(x), "run-abc")
  expect_identical(tales_namespace(dplyr::filter(x, array_id == "a1")), "run-abc")
  expect_identical(tales_namespace(dplyr::mutate(x, note = 1)), "run-abc")
})


#### Anchor codes ####

test_that("tales_anchor_codes() covers all three terminus sentinels", {
  expect_setequal(tales_anchor_codes(), c("NTERM", "CTERM", "XXXXX"))
})


#### Against real pipeline output ####

test_that("a real distalr tale_parts table validates as tales", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  tp <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))$tale_parts
  x <- tales(tp)
  expect_s3_class(x, "tales")
  expect_equal(nrow(x), nrow(tp))
  expect_true(all(c("array_id", "position_in_array", "dom_code") %in% names(x)))
})
