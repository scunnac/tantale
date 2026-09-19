# Tests for tales_group_kmedoids(). See ledger §11's 2026-09-19 joint session.

distalrOut <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
tal_sim <- distalrOut$tal.similarity
nTales <- length(unique(tal_sim$TAL1))
tls <- suppressWarnings(tales(distalrOut$tale_parts))

expect_grouped <- function(out, k = NULL) {
  testthat::expect_s3_class(out, "tales")
  testthat::expect_true("group" %in% names(out))
  testthat::expect_false(anyNA(out$group))
  perArray <- tapply(out$group, out$array_id, function(z) length(unique(z)))
  testthat::expect_true(all(perArray == 1L))
  testthat::expect_equal(length(unique(out$array_id)), nTales)
  if (!is.null(k)) testthat::expect_equal(length(unique(out$group)), k)
  invisible(out)
}

# No display needed for base plot() calls in a test run.
local_no_display <- function(envir = parent.frame()) {
  withr::local_pdf(NULL, .local_envir = envir)
}


#### grouping itself ####

test_that("an explicit k returns the expected number of groups", {
  local_no_display()
  grp <- tales_group_kmedoids(tls, tal_sim, k_range = 2:6, k = 3)
  expect_grouped(grp, k = 3)
})

test_that("k = 'auto' picks a number of groups automatically", {
  local_no_display()
  expect_message(
    grp <- tales_group_kmedoids(tls, tal_sim, k_range = 2:6, k = "auto"),
    "automatically decided"
  )
  expect_grouped(grp)
})

test_that("the k-medoids clustering returned is a real PAM assignment, not a placeholder", {
  # Pinning the point raised in the joint session: k-medoids returns actual
  # cluster::pam() cluster membership for the chosen k, the same kind of
  # real result tales_group_hclust() returns via cutree() -- not merely a
  # suggested k value.
  local_no_display()
  grp <- tales_group_kmedoids(tls, tal_sim, k_range = 2:6, k = 3)
  distMat <- as.matrix(tale_distances(tal_sim))
  set.seed(7)
  expected <- cluster::pam(stats::as.dist(distMat), 3)$clustering
  out <- unique(as.data.frame(grp)[c("array_id", "group")])
  expect_equal(out$group[match(names(expected), out$array_id)], unname(expected))
})


#### k_range ####

test_that("errors on invalid k_range", {
  expect_error(tales_group_kmedoids(tls, tal_sim, k_range = "not numeric", k = 3),
               class = "tantale_error_group_kmedoids_krange")
  expect_error(tales_group_kmedoids(tls, tal_sim, k_range = NULL, k = 3),
               class = "tantale_error_group_kmedoids_krange")
})


#### k: the three-type contract ####

test_that("k = NULL outside an interactive session errors instead of blocking on stdin", {
  expect_false(interactive()) # the premise this test relies on, under testthat
  expect_error(tales_group_kmedoids(tls, tal_sim, k_range = 2:6),
               class = "tantale_error_group_kmedoids_k")
})

test_that("an invalid k value or type is a classed error, including length > 1", {
  # Previously `k == "auto"` on a length-2 k tripped R's own "condition has
  # length > 1" before ever reaching the intended message (old §11
  # observation) -- identical(k, "auto") fixes that; both now hit the same
  # classed error.
  expect_error(tales_group_kmedoids(tls, tal_sim, k_range = 2:6, k = c(2, 3)),
               class = "tantale_error_group_kmedoids_k")
  expect_error(tales_group_kmedoids(tls, tal_sim, k_range = 2:6, k = TRUE),
               class = "tantale_error_group_kmedoids_k")
})


#### seed ####

test_that("seed is documented and overridable, defaulting to the old hardcoded 7", {
  local_no_display()
  distMat <- as.matrix(tale_distances(tal_sim))

  set.seed(7)
  baseline <- cluster::pam(stats::as.dist(distMat), 3)$clustering
  grp_default <- tales_group_kmedoids(tls, tal_sim, k_range = 2:6, k = 3)
  out_default <- unique(as.data.frame(grp_default)[c("array_id", "group")])
  expect_equal(out_default$group[match(names(baseline), out_default$array_id)],
               unname(baseline))

  expect_equal(formals(tales_group_kmedoids)$seed, 7)
})


#### plotting ####

test_that("plot_silhouette = FALSE draws nothing", {
  # With no display device open, any attempt to plot would error.
  expect_no_error(
    tales_group_kmedoids(tls, tal_sim, k_range = 2:6, k = 3, plot_silhouette = FALSE)
  )
})


#### the tales <-> distances correspondence (ledger 5.2, shared helper) ####

test_that("rejects distances that describe other arrays", {
  local_no_display()
  other <- tal_sim
  other$TAL1 <- paste0("not_", other$TAL1)
  other$TAL2 <- paste0("not_", other$TAL2)
  expect_error(tales_group_kmedoids(tls, other, k_range = 2:6, k = 3),
               class = "tantale_error_group_mismatch")
})

test_that("rejects a tales that is missing some compared arrays", {
  local_no_display()
  sub <- tls[tls$array_id %in% unique(tls$array_id)[1:2], ]
  expect_error(tales_group_kmedoids(sub, tal_sim, k_range = 2:6, k = 3),
               class = "tantale_error_group_mismatch")
})

test_that("needs a tales object", {
  expect_error(tales_group_kmedoids(distalrOut$tale_parts, tal_sim, k_range = 2:6, k = 3),
               class = "tantale_error_tales_type")
})

test_that("the bare mapping is recoverable from the returned object", {
  local_no_display()
  out <- tales_group_kmedoids(tls, tal_sim, k_range = 2:6, k = 3)
  mapping <- unique(as.data.frame(out)[c("array_id", "group")])
  expect_equal(nrow(mapping), nTales)
})
