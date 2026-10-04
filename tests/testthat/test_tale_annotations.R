# The curated reference table (ledger §57). These guard the promises the
# documentation makes about it, since a reader joins against them.

test_that("tale_annotations has the documented shape and key", {
  expect_s3_class(tale_annotations, "tbl_df")
  expect_identical(dim(tale_annotations), c(128L, 10L))
  expect_named(tale_annotations,
               c("strain", "label", "tal_name", "annotale_group",
                 "replicon_id", "genome_id", "pubmed", "truncTALE",
                 "rvd_seq", "unusual_feature"))
  # strain + label identifies a row, and neither is ever missing
  expect_false(anyNA(tale_annotations$strain))
  expect_false(anyNA(tale_annotations$label))
  expect_identical(anyDuplicated(tale_annotations[c("strain", "label")]), 0L)
  expect_identical(dplyr::n_distinct(tale_annotations$strain), 10L)
})

test_that("every row names the genome it came from", {
  # this is what lets a reader fetch the sequence, so it is not optional
  expect_false(anyNA(tale_annotations$genome_id))
  expect_false(anyNA(tale_annotations$replicon_id))
  expect_true(all(grepl("^GC[AF]_[0-9]+\\.[0-9]+$", tale_annotations$genome_id)))
})

test_that("the RVD strings parse as RVDs", {
  expect_false(anyNA(tale_annotations$rvd_seq))
  rvds <- unlist(strsplit(tale_annotations$rvd_seq, "-", fixed = TRUE))
  # two characters each: a residue pair, or a residue plus "*" when the
  # 13th is missing. Case carries meaning, so it is not normalised away.
  expect_true(all(nchar(rvds) == 2L))
  expect_true(all(grepl("^[A-Za-z][A-Za-z*]$", rvds)))
  # the four common RVDs are the bulk of any real talome
  expect_true(all(c("HD", "NG", "NI", "NN") %in% rvds))
})

test_that("lowercase RVDs mark aberrant repeats and are kept as such", {
  # A lowercase RVD is the convention for a repeat of non-standard length;
  # tantale reads it that way elsewhere (grepl("[a-z]", ...) in
  # classification.R), so upper-casing this column would destroy
  # information. 16 of the 128 arrays carry at least one.
  aberrant <- grepl("[a-z]", tale_annotations$rvd_seq)
  expect_identical(sum(aberrant), 16L)
  rvds <- unlist(strsplit(tale_annotations$rvd_seq, "-", fixed = TRUE))
  expect_true(any(grepl("[a-z]", rvds)))
})

test_that("truncTALE is logical and flags the documented nine", {
  # the column keeps the field's spelling, not snake_case (ledger §57)
  expect_type(tale_annotations$truncTALE, "logical")
  expect_false(anyNA(tale_annotations$truncTALE))
  expect_identical(sum(tale_annotations$truncTALE), 9L)
})

test_that("the three genomes tantale installs are all present", {
  expect_true(all(c("MAI1", "BAI3", "PXO86") %in% tale_annotations$strain))
})
