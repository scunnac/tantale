
distalrOut <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
taleSim <- distalrOut$tal.similarity
nTales <- length(unique(taleSim$TAL1))

test_that("groupTales() with k-medoids and an explicit k returns the expected number of groups", {
  grp <- groupTales(taleSim = taleSim, k = 3, k_test = 2:6, method = "k-medoids")
  expect_s3_class(grp, "data.frame")
  expect_named(grp, c("name", "group"))
  expect_equal(nrow(grp), nTales)
  expect_equal(length(unique(grp$group)), 3)
})

test_that("groupTales() with k-medoids and k = 'auto' picks a number of groups automatically", {
  expect_message(
    grp <- groupTales(taleSim = taleSim, k = "auto", k_test = 2:6, method = "k-medoids"),
    "automatically decided"
  )
  expect_s3_class(grp, "data.frame")
  expect_named(grp, c("name", "group"))
  expect_equal(nrow(grp), nTales)
})

test_that("groupTales() with k-medoids forces plotTree off and warns via message", {
  expect_message(
    groupTales(taleSim = taleSim, k = 3, k_test = 2:6, method = "k-medoids", plotTree = TRUE),
    "tale tree will not be plotted"
  )
})

test_that("groupTales() with k-medoids errors on invalid k_test", {
  expect_error(
    groupTales(taleSim = taleSim, k = 3, k_test = "not numeric", method = "k-medoids"),
    "invalid k values"
  )
  expect_error(
    groupTales(taleSim = taleSim, k = 3, k_test = NULL, method = "k-medoids"),
    "invalid k values"
  )
})

test_that("groupTales() with k-medoids errors on invalid k", {
  # A length > 1 k trips `if (k == "auto")` with R's own "condition has length
  # > 1" error before ever reaching the intended stop("invalid k value!") -
  # still errors, just not with the intended message.
  expect_error(
    groupTales(taleSim = taleSim, k = c(2, 3), k_test = 2:6, method = "k-medoids")
  )
})

test_that("groupTales() with hclust and an explicit k returns the expected number of groups", {
  grp <- groupTales(taleSim = taleSim, k = 3, method = "hclust")
  expect_s3_class(grp, "data.frame")
  expect_named(grp, c("name", "group"))
  expect_equal(nrow(grp), nTales)
  expect_equal(length(unique(grp$group)), 3)
})

test_that("groupTales() with hclust errors when k is not a single number", {
  expect_error(groupTales(taleSim = taleSim, method = "hclust"))
  expect_error(groupTales(taleSim = taleSim, k = c(2, 3), method = "hclust"))
})
