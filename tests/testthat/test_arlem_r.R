# The R implementation of ARLEM's alignment (R/arlem.R, ledger 33), which
# replaced the ARLEM 1.0 executable tantale used to bundle (ledger 28).
#
# Expected scores here all come from that executable: the toy cases were
# read off it by hand, and arlem_reference_scores.rds holds its answers for
# 36 cases, 1126 pairs (data-raw/make_arlem_reference_scores.R).

toy_cost <- function(values) {
  # upper triangle, row by row, of a symmetric matrix over types 1..k
  k <- (1 + sqrt(1 + 8 * length(values))) / 2
  m <- matrix(0, k, k, dimnames = list(seq_len(k), seq_len(k)))
  m[lower.tri(m)] <- values
  m <- t(m)
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

toy_coded <- function(...) Biostrings::BStringSet(c(...))

arlem_quiet <- function(...) suppressMessages(.arlem_scores_r(...))

scores_of <- function(coded, cost, ...) arlem_quiet(coded, cost, ...)$arlem_score


#### the model, pinned to the binary's own answers ####

test_that("matching, and a leading unit explained as an insertion", {
  cost <- toy_cost(c(5, 8, 5))
  coded <- toy_coded("1 1 2 3", "1 2 3 3", "2 2 2")
  # 0-1: pure substitutions (1/2, 2/3); 0-2 and 1-2: one unit inserted
  # before the first match, at the insertion cost
  expect_equal(scores_of(coded, cost, dup = 10, ins = 10), c(10, 20, 20))
})

test_that("duplication and insertion costs each enter where they should", {
  cost <- toy_cost(c(50, 80, 50))
  coded <- toy_coded("1 2 3", "1 3", "2", "1 1 1 3")
  expect_equal(scores_of(coded, cost, dup = 10, ins = 10),
               c(10, 20, 30, 30, 20, 50))
  # expensive insertions: runs must be grown by duplication instead
  expect_equal(scores_of(coded, cost, dup = 10, ins = 1000),
               c(60, 120, 60, 140, 20, 160))
  # cheap duplications: "1 1 1" costs two copies of 1, not two insertions
  expect_equal(scores_of(coded, cost, dup = 1, ins = 10),
               c(10, 20, 12, 30, 2, 32))
})

test_that("scores are reported like the binary's: one row per pair, 0-based", {
  out <- arlem_quiet(toy_coded("1", "2", "1 2"), toy_cost(5))
  expect_named(out, c("id1", "id2", "arlem_score"))
  expect_equal(out$id1, c(0, 0, 1))
  expect_equal(out$id2, c(1, 2, 2))
})

test_that("a single array has no pairs to score", {
  out <- arlem_quiet(toy_coded("1 2"), toy_cost(5))
  expect_equal(nrow(out), 0L)
})

test_that("a domain code missing from the cost matrix is refused", {
  expect_error(arlem_quiet(toy_coded("1 2", "1 7"), toy_cost(5)),
               class = "tantale_error_arlem_codes")
})


#### agreement with the ARLEM executable's recorded answers ####

reference <- readRDS(test_path("data_for_tests", "arlem_reference_scores.rds"))

test_that("scores match the executable's on random maps", {
  expect_length(reference$random, 35)
  for (i in seq_along(reference$random)) {
    case <- reference$random[[i]]
    out <- arlem_quiet(Biostrings::BStringSet(case$maps), case$cost,
                           dup = case$dup, ins = case$ins, insert = case$insert)
    expect_identical(out$arlem_score, case$scores$arlem_score,
                     info = paste("case", i, "insert", case$insert))
    expect_equal(out$id1, case$scores$id1)
    expect_equal(out$id2, case$scores$id2)
  }
})

test_that("scores match the executable's on 44 real TALE arrays", {
  case <- reference$real
  expect_equal(nrow(case$scores), choose(44, 2))
  out <- arlem_quiet(Biostrings::BStringSet(case$maps), case$cost)
  expect_identical(out$arlem_score, case$scores$arlem_score)
})

test_that("tales_tale_distances() reproduces the executable's distances", {
  # the whole step 3, from the fixture's own codes and domain distances
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  parts <- tibble::as_tibble(d$tale_parts)
  coded <- .coded_seq_set(parts)
  cost <- suppressMessages(.arlem_cost_matrix(domain_distances(d$repeat.similarity)))
  expect_identical(cost, reference$real$cost)

  x <- suppressWarnings(tales_quietly(d$tale_parts))
  dd <- domain_distances(d$repeat.similarity,
                         dom_code_namespace = tales_namespace(x))
  td <- suppressMessages(tales_tale_distances(x, dd))
  expected <- tale_distances(.normalise_arlem_scores(
    .arlem_scores_long(reference$real$scores, length(coded)), parts, coded))
  expect_identical(tibble::as_tibble(td), tibble::as_tibble(expected))
})
