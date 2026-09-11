
distalrOut <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
tal_sim <- distalrOut$tal.similarity
nTales <- length(unique(tal_sim$TAL1))

test_that("group_tales() with k-medoids and an explicit k returns the expected number of groups", {
  grp <- group_tales(tal_sim = tal_sim, k = 3, k_range = 2:6, method = "k-medoids")
  expect_s3_class(grp, "data.frame")
  expect_named(grp, c("name", "group"))
  expect_equal(nrow(grp), nTales)
  expect_equal(length(unique(grp$group)), 3)
})

test_that("group_tales() with k-medoids and k = 'auto' picks a number of groups automatically", {
  expect_message(
    grp <- group_tales(tal_sim = tal_sim, k = "auto", k_range = 2:6, method = "k-medoids"),
    "automatically decided"
  )
  expect_s3_class(grp, "data.frame")
  expect_named(grp, c("name", "group"))
  expect_equal(nrow(grp), nTales)
})

test_that("group_tales() with k-medoids forces plot_tree off and warns via message", {
  expect_message(
    group_tales(tal_sim = tal_sim, k = 3, k_range = 2:6, method = "k-medoids", plot_tree = TRUE),
    "tale tree will not be plotted"
  )
})

test_that("group_tales() with k-medoids errors on invalid k_range", {
  expect_error(
    group_tales(tal_sim = tal_sim, k = 3, k_range = "not numeric", method = "k-medoids"),
    "invalid k values"
  )
  expect_error(
    group_tales(tal_sim = tal_sim, k = 3, k_range = NULL, method = "k-medoids"),
    "invalid k values"
  )
})

test_that("group_tales() with k-medoids errors on invalid k", {
  # A length > 1 k trips `if (k == "auto")` with R's own "condition has length
  # > 1" error before ever reaching the intended stop("invalid k value!") -
  # still errors, just not with the intended message.
  expect_error(
    group_tales(tal_sim = tal_sim, k = c(2, 3), k_range = 2:6, method = "k-medoids")
  )
})

test_that("group_tales() with hclust and an explicit k returns the expected number of groups", {
  grp <- group_tales(tal_sim = tal_sim, k = 3, method = "hclust")
  expect_s3_class(grp, "data.frame")
  expect_named(grp, c("name", "group"))
  expect_equal(nrow(grp), nTales)
  expect_equal(length(unique(grp$group)), 3)
})

test_that("group_tales() with hclust errors when k is not a single number", {
  expect_error(group_tales(tal_sim = tal_sim, method = "hclust"))
  expect_error(group_tales(tal_sim = tal_sim, k = c(2, 3), method = "hclust"))
})
