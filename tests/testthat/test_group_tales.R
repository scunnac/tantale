
distalrOut <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
tal_sim <- distalrOut$tal.similarity
nTales <- length(unique(tal_sim$TAL1))
tls <- suppressWarnings(tales(distalrOut$tale_parts))

# tales_group() now returns the tales object with `group` filled, rather than
# a bare name/group table (ledger 5.2). These assert the shape once so the
# per-method tests below can concentrate on the clustering.
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

test_that("tales_group() with k-medoids and an explicit k returns the expected number of groups", {
  grp <- tales_group(tls, tal_sim = tal_sim, k = 3, k_range = 2:6, method = "k-medoids")
  expect_grouped(grp, k = 3)
})

test_that("tales_group() with k-medoids and k = 'auto' picks a number of groups automatically", {
  expect_message(
    grp <- tales_group(tls, tal_sim = tal_sim, k = "auto", k_range = 2:6, method = "k-medoids"),
    "automatically decided"
  )
  expect_grouped(grp)
})

test_that("tales_group() with k-medoids forces plot_tree off and warns via message", {
  expect_message(
    tales_group(tls, tal_sim = tal_sim, k = 3, k_range = 2:6, method = "k-medoids", plot_tree = TRUE),
    "tale tree will not be plotted"
  )
})

test_that("tales_group() with k-medoids errors on invalid k_range", {
  expect_error(
    tales_group(tls, tal_sim = tal_sim, k = 3, k_range = "not numeric", method = "k-medoids"),
    "invalid k values"
  )
  expect_error(
    tales_group(tls, tal_sim = tal_sim, k = 3, k_range = NULL, method = "k-medoids"),
    "invalid k values"
  )
})

test_that("tales_group() with k-medoids errors on invalid k", {
  # A length > 1 k trips `if (k == "auto")` with R's own "condition has length
  # > 1" error before ever reaching the intended stop("invalid k value!") -
  # still errors, just not with the intended message.
  expect_error(
    tales_group(tls, tal_sim = tal_sim, k = c(2, 3), k_range = 2:6, method = "k-medoids")
  )
})

test_that("tales_group() with hclust and an explicit k returns the expected number of groups", {
  grp <- tales_group(tls, tal_sim = tal_sim, k = 3, method = "hclust")
  expect_grouped(grp, k = 3)
})

test_that("tales_group() with hclust errors when k is not a single number", {
  expect_error(tales_group(tls, tal_sim = tal_sim, method = "hclust"))
  expect_error(tales_group(tls, tal_sim = tal_sim, k = c(2, 3), method = "hclust"))
})


#### the tales <-> distances correspondence (ledger 5.2) ####

test_that("tales_group() rejects distances that describe other arrays", {
  # The reason the function takes x at all: this is the only place the
  # correspondence can be checked.
  other <- tal_sim
  other$TAL1 <- paste0("not_", other$TAL1)
  other$TAL2 <- paste0("not_", other$TAL2)
  expect_error(tales_group(tls, tal_sim = other, k = 3, method = "hclust"),
               class = "tantale_error_group_mismatch")
})

test_that("tales_group() rejects a tales that is missing some compared arrays", {
  sub <- tls[tls$array_id %in% unique(tls$array_id)[1:2], ]
  expect_error(tales_group(sub, tal_sim = tal_sim, k = 3, method = "hclust"),
               class = "tantale_error_group_mismatch")
})

test_that("tales_group() needs a tales object", {
  expect_error(tales_group(distalrOut$tale_parts, tal_sim = tal_sim, k = 3,
                           method = "hclust"),
               class = "tantale_error_tales_type")
})

test_that("the bare mapping is recoverable from the returned object", {
  out <- tales_group(tls, tal_sim = tal_sim, k = 3, method = "hclust")
  mapping <- unique(as.data.frame(out)[c("array_id", "group")])
  expect_equal(nrow(mapping), nTales)
})
