# tales_consensus() and tales_consensus_match(), on hand-built alignments
# whose answers are known, plus their tales_msa-native counterparts
# (ledger §21), which must agree with them on a real alignment.

aln <- matrix(c("HD", "HD", "NI",     # 1: HD is a strict majority
                "NG", "NI", "HD",     # 2: three-way tie
                "NN", NA,   NA,       # 3: gaps outnumber the one residue
                "NI", "NI", NA,       # 4: a majority despite a gap
                "HD", "hd", "HD"),    # 5: an aberrant (lower-case) HD
              nrow = 3, dimnames = list(c("A1", "A2", "A3"), NULL))

test_that("the consensus needs a strict majority, and a gap majority gives none", {
  expect_identical(tales_consensus(aln), c("HD", NA, NA, "NI", "HD"))
})

test_that("the consensus does not depend on the order of the arrays", {
  for (i in 1:5) {
    shuffled <- aln[sample(nrow(aln)), , drop = FALSE]
    expect_identical(tales_consensus(shuffled), tales_consensus(aln))
  }
})

test_that("matches are TRUE/FALSE, NA where there is no consensus, gaps never match", {
  m <- tales_consensus_match(aln, long = FALSE)
  expect_type(m, "logical")
  expect_identical(dimnames(m), dimnames(aln))
  expect_identical(unname(m[, 1]), c(TRUE, TRUE, FALSE))
  expect_true(all(is.na(m[, 2])))
  expect_true(all(is.na(m[, 3])))
  expect_identical(unname(m[, 4]), c(TRUE, TRUE, FALSE))  # the gap is not a match
  expect_identical(unname(m[, 5]), c(TRUE, TRUE, TRUE))   # case is ignored
})

test_that("the long form holds the same values, keyed by alignment position", {
  wide <- tales_consensus_match(aln, long = FALSE)
  long <- tales_consensus_match(aln)
  expect_named(long, c("array_id", "alignment_position", "tales_consensus_match"))
  expect_identical(nrow(long), length(wide))
  idx <- cbind(match(long$array_id, rownames(wide)), long$alignment_position)
  expect_identical(long$tales_consensus_match, wide[idx])
})

test_that("the tales_msa-native versions agree with the matrix versions", {
  msa <- readRDS(test_path("data_for_tests", "sampleTalesMsa.rds"))
  for (layer in c("dom_code", "rvd")) {
    mat <- as.matrix(msa, value = layer)
    native <- .tales_consensus_long(msa, layer)
    expect_identical(native[[layer]][order(native$alignment_position)],
                     unname(tales_consensus(mat)), info = layer)

    wide <- tales_consensus_match(mat, long = FALSE)
    native_match <- .tales_consensus_match_long(msa, layer)
    idx <- cbind(match(native_match$array_id, rownames(wide)),
                 native_match$alignment_position)
    expect_identical(native_match$match, unname(wide[idx]), info = layer)
  }
})
