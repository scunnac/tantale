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

test_that("NA in a residue column is an anomaly, not an error", {
  df <- minimal_tales_df()
  df$rvd[1] <- NA_character_
  expect_warning(x <- tales(df), class = "tantale_warning_tales_anomalous")
  expect_true("missing_rvd" %in% tales_anomalies(x)$check)
})


#### Conditional invariants ####

test_that("an unknown domain_type is an anomaly, not an error", {
  df <- minimal_tales_df()
  df$domain_type[2] <- "middle"
  expect_warning(x <- tales(df), class = "tantale_warning_tales_anomalous")
  expect_true("domain_type_unknown" %in% tales_anomalies(x)$check)
})

test_that("more than one terminus of a kind is an anomaly", {
  df <- minimal_tales_df()
  df$domain_type[3] <- "C-terminus"
  df$position_in_crd[3] <- NA_integer_
  expect_warning(x <- tales(df), class = "tantale_warning_tales_anomalous")
  expect_true(any(grepl("terminus", tales_anomalies(x)$check)))
})

test_that("zero termini is allowed - a repeats-only subset stays valid", {
  x <- tales(minimal_tales_df())
  repeats_only <- dplyr::filter(x, domain_type == "repeat")
  expect_s3_class(repeats_only, "tales")
  expect_silent(validate_tales(repeats_only))
  # position_in_array no longer starts at 1, and that is fine
  expect_false(min(repeats_only$position_in_array) == 1L)
})

test_that("position_in_crd misplacement is an anomaly", {
  df <- minimal_tales_df()
  # a value on a terminus, chosen not to collide with any repeat's coordinate
  # (a collision would trip the structural uniqueness check first, correctly)
  df$position_in_crd[1] <- 99L
  expect_warning(x <- tales(df), class = "tantale_warning_tales_anomalous")
  expect_true("crd_placement" %in% tales_anomalies(x)$check)
})

test_that("position_in_crd disagreeing on order is an anomaly", {
  df <- minimal_tales_df()
  # swap the two repeats' crd positions: same set, wrong order
  df$position_in_crd[df$domain_type == "repeat"] <- c(2L, 1L, 2L, 1L)
  expect_warning(x <- tales(df), class = "tantale_warning_tales_anomalous")
  expect_true("crd_order" %in% tales_anomalies(x)$check)
})

test_that("a duplicated position_in_crd is still a structural error", {
  # uniqueness is the part of the CRD contract that stays hard: a repeated
  # coordinate makes the array unindexable, not merely odd
  df <- minimal_tales_df()
  df$position_in_crd[df$domain_type == "repeat"] <- c(1L, 1L, 1L, 2L)
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

test_that("seqnames varying within an array is an anomaly", {
  df <- minimal_tales_df()
  df$seqnames[1] <- "other"
  expect_warning(x <- tales(df), class = "tantale_warning_tales_anomalous")
  expect_true("seqnames_inconsistent" %in% tales_anomalies(x)$check)
  # and sanitize removes the offending array
  expect_false(unique(df$array_id)[1] %in% suppressWarnings(tales(df, sanitize = TRUE))$array_id)
})

test_that("an empty dna_seq warns but does not fail", {
  df <- minimal_tales_df()
  df$dna_seq <- c("ATG", rep("", 7))
  expect_warning(x <- tales(df), class = "tantale_warning_tales_anomalous")
  expect_s3_class(x, "tales")
  expect_true("missing_dna_seq" %in% tales_anomalies(x)$check)
})

test_that("the same aa_seq paired with two different rvd values is an anomaly", {
  df <- minimal_tales_df()
  # a1 and a2 share the "LTPA" aa_seq at position 2; give a2's copy a
  # different rvd without touching dom_code, so only this check fires.
  df$rvd[df$array_id == "a2" & df$aa_seq == "LTPA"] <- "HD"
  expect_warning(x <- tales(df), class = "tantale_warning_tales_anomalous")
  an <- tales_anomalies(x)
  expect_true("aa_seq_rvd_inconsistent" %in% an$check)
  expect_setequal(an$array_id[an$check == "aa_seq_rvd_inconsistent"], c("a1", "a2"))
  # both arrays share the offending aa_seq, so sanitize drops both
  clean <- suppressWarnings(tales(df, sanitize = TRUE))
  expect_false(any(c("a1", "a2") %in% clean$array_id))
})

test_that("different aa_seq sharing one rvd is not an anomaly", {
  df <- minimal_tales_df()
  # a2's repeat at position 2 becomes a different protein with a fresh
  # dom_code, but keeps a1's rvd at that position -- legitimate: distinct
  # repeats routinely share a binding specificity.
  df$aa_seq[df$array_id == "a2" & df$position_in_array == 2] <- "LTPZ"
  df$dom_code[df$array_id == "a2" & df$position_in_array == 2] <- "5"
  x <- tales(df)
  expect_equal(nrow(tales_anomalies(x)), 0L)
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


#### as_tales() ####

test_that("as_tales() builds a tales from an RVD fasta", {
  x <- as_tales(test_path("data_for_tests", "tellTaleExampleOutput", "rvd_sequences.fas"),
                sep = "-")
  expect_s3_class(x, "tales")
  expect_setequal(names(x), c("array_id", "position_in_array", "rvd"))
  # one row per element of each sequence
  expect_setequal(as.integer(table(x$array_id)), c(28L, 16L, 28L, 24L))
})

test_that("as_tales() numbers position_in_array from sequence order", {
  x <- as_tales(test_path("data_for_tests", "tellTaleExampleOutput", "rvd_sequences.fas"),
                sep = "-")
  first <- dplyr::filter(x, array_id == x$array_id[1])
  expect_identical(first$position_in_array, seq_len(nrow(first)))
})

test_that("as_tales() accepts a BStringSet and a list, matching the file path result", {
  p <- test_path("data_for_tests", "tellTaleExampleOutput", "rvd_sequences.fas")
  from_path <- as_tales(p, sep = "-")
  from_set <- as_tales(Biostrings::readBStringSet(p), sep = "-")
  from_list <- as_tales(as.list(as.character(Biostrings::readBStringSet(p))), sep = "-")
  expect_equal(from_set, from_path)
  expect_equal(from_list, from_path)
})

test_that("as_tales() puts repeat codes in dom_code when asked", {
  x <- suppressWarnings(
    as_tales(test_path("data_for_tests", "Out_CodedRepeats.fa"),
             sep = " ", residue_col = "dom_code")
  )
  expect_s3_class(x, "tales")
  expect_true("dom_code" %in% names(x))
  expect_false("rvd" %in% names(x))
})

test_that("as_tales() errors on unnamed sequences", {
  expect_error(
    as_tales(list("NI-HD-NG"), sep = "-"),
    class = "tantale_error_tales_unnamed"
  )
})

test_that("as_tales() on a data frame is tales()", {
  df <- minimal_tales_df()
  expect_equal(as_tales(df), tales(df))
})


#### tales_from_telltale() ####

test_that("tales_from_telltale() returns a validated tales", {
  x <- suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
  expect_s3_class(x, "tales")
  expect_true(all(c("array_id", "position_in_array", "rvd", "domain_type",
                    "position_in_crd", "aa_seq", "seqnames") %in% names(x)))
  # dom_code is minted later, by the relatedness computation
  expect_false("dom_code" %in% names(x))
})

test_that("tales_from_telltale() output holds complete arrays", {
  x <- suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
  expect_silent(tales_assert_complete(x))
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





test_that("a missing residue is now a biological anomaly, not a structural error", {
  # .tales_check_residue_na() was absorbed into .tales_anomalies() when the
  # class was relaxed: a part with no RVD is odd biology, not a corrupt table,
  # so it loads with a warning rather than aborting.
  df <- tibble::tibble(array_id = c("a", "a"), position_in_array = 1:2,
                       rvd = c("NI", NA))
  expect_warning(x <- tales(df), class = "tantale_warning_tales_anomalous")
  expect_s3_class(x, "tales")
  expect_identical(tales_anomalies(x)$array_id, "a")
  expect_identical(tales_anomalies(x)$check, "missing_rvd")
  expect_equal(nrow(suppressWarnings(tales(df, sanitize = TRUE))), 0L)
})
