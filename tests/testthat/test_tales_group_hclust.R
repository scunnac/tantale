# Tests for tales_group_hclust(). See ledger §11's 2026-09-19 joint session.

distalrOut <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
tal_sim <- distalrOut$tal.similarity
nTales <- length(unique(tal_sim$TAL1))
tls <- suppressWarnings(tales(distalrOut$tale_parts))

expect_grouped <- function(out, k = NULL) {
  testthat::expect_s3_class(out, "tales")
  testthat::expect_true("group" %in% names(out))
  testthat::expect_false(anyNA(out$group))
  # group is an array-level property: constant within each array
  perArray <- tapply(out$group, out$array_id, function(z) length(unique(z)))
  testthat::expect_true(all(perArray == 1L))
  testthat::expect_equal(length(unique(out$array_id)), nTales)
  if (!is.null(k)) testthat::expect_equal(length(unique(out$group)), k)
  invisible(out)
}

test_that("an explicit k returns the expected number of groups", {
  grp <- tales_group_hclust(tls, tal_sim, k = 3)
  expect_grouped(grp, k = 3)
})

test_that("k must be a single number", {
  expect_error(tales_group_hclust(tls, tal_sim),
               class = "tantale_error_group_hclust_k")
  expect_error(tales_group_hclust(tls, tal_sim, k = c(2, 3)),
               class = "tantale_error_group_hclust_k")
  expect_error(tales_group_hclust(tls, tal_sim, k = "nope"),
               class = "tantale_error_group_hclust_k")
})

test_that("cutree(k=) never fails on a tied tree, unlike the old bisection search", {
  # Reproduces the tie case worked through with the maintainer: two merges at
  # the exact same height used to make the old height-search hit lo >= hi and
  # throw "Cannot determine k groups!". cutree(k=) has no such failure mode.
  d <- matrix(c(
    0, 2, 5, 5,
    2, 0, 5, 5,
    5, 5, 0, 2,
    5, 5, 2, 0
  ), 4, 4, dimnames = list(c("T1", "T2", "T3", "T4"), c("T1", "T2", "T3", "T4")))
  tree <- stats::hclust(stats::as.dist(d), method = "ward.D")
  expect_equal(tree$height[1], tree$height[2]) # the tie this test relies on

  cuts <- stats::cutree(tree, k = 3)
  expect_length(unique(cuts), 3L)
})

test_that("clusters the distance matrix directly, not the Euclidean distance between rows", {
  # A, B, C mutually distance 1 (obviously one group); D is a distant
  # outlier. If tales_group_hclust() were still using
  # dist(distMat, "euclidean") on the rows, this would not necessarily hold
  # -- as.dist(distMat) guarantees it, since A/B/C's raw mutual distance is
  # the smallest possible.
  labs <- c("A", "B", "C", "D")
  d <- matrix(0, 4, 4, dimnames = list(labs, labs))
  d["A", "B"] <- 1; d["A", "C"] <- 1; d["B", "C"] <- 1
  d["A", "D"] <- 9; d["B", "D"] <- 9; d["C", "D"] <- 9
  d[lower.tri(d)] <- t(d)[lower.tri(d)]

  # id1/id2/dissim directly -- not the legacy TAL1/TAL2/Sim spelling, which
  # is a *similarity* and gets inverted to dissim = 100 - sim on the way in.
  fake_sim <- data.frame(
    id1 = rep(labs, each = 4), id2 = rep(labs, 4),
    dissim = as.vector(d)
  )
  fake_tls <- suppressWarnings(tales(data.frame(
    array_id = labs, position_in_array = 1L, rvd = "HD"
  )))
  grp <- tales_group_hclust(fake_tls, fake_sim, k = 2)
  out <- unique(as.data.frame(grp)[c("array_id", "group")])
  abc <- out$group[out$array_id %in% c("A", "B", "C")]
  expect_true(length(unique(abc)) == 1L)
  expect_false(out$group[out$array_id == "D"] %in% abc)
})

test_that("plot_tree defaults to FALSE and skips the ggtree machinery", {
  # No warning/message from ggtree/tidytree should fire when not asked for.
  expect_no_condition(tales_group_hclust(tls, tal_sim, k = 3))
})

test_that("plot_tree = TRUE draws without erroring", {
  grDevices::pdf(NULL) # no display needed in a test run
  on.exit(grDevices::dev.off())
  expect_no_error(suppressWarnings(tales_group_hclust(tls, tal_sim, k = 3, plot_tree = TRUE)))
})


#### the tales <-> distances correspondence (ledger 5.2, shared helper) ####

test_that("rejects distances that describe other arrays", {
  other <- tal_sim
  other$TAL1 <- paste0("not_", other$TAL1)
  other$TAL2 <- paste0("not_", other$TAL2)
  expect_error(tales_group_hclust(tls, other, k = 3),
               class = "tantale_error_group_mismatch")
})

test_that("rejects a tales that is missing some compared arrays", {
  sub <- tls[tls$array_id %in% unique(tls$array_id)[1:2], ]
  expect_error(tales_group_hclust(sub, tal_sim, k = 3),
               class = "tantale_error_group_mismatch")
})

test_that("needs a tales object", {
  expect_error(tales_group_hclust(distalrOut$tale_parts, tal_sim, k = 3),
               class = "tantale_error_tales_type")
})

test_that("the bare mapping is recoverable from the returned object", {
  out <- tales_group_hclust(tls, tal_sim, k = 3)
  mapping <- unique(as.data.frame(out)[c("array_id", "group")])
  expect_equal(nrow(mapping), nTales)
})
