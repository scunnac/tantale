# Tests for the 'pairwise_distances' class. See dev/class-design.md §3.

# A complete 3x3 distance, long form. Self-comparison is 0.
minimal_distances_df <- function() {
  ids <- c("a", "b", "c")
  g <- expand.grid(id1 = ids, id2 = ids, stringsAsFactors = FALSE)
  tibble::tibble(
    id1 = g$id1,
    id2 = g$id2,
    dissim = ifelse(g$id1 == g$id2, 0, 50)
  )
}

# The same thing in the legacy similarity vocabulary.
legacy_sim_df <- function() {
  d <- minimal_distances_df()
  tibble::tibble(id1 = d$id1, id2 = d$id2, sim = 100 - d$dissim)
}

test_that("pairwise_distances() accepts a complete square table", {
  x <- pairwise_distances(minimal_distances_df())
  expect_s3_class(x, "pairwise_distances")
  expect_s3_class(x, "tbl_df")
  expect_true(is_pairwise_distances(x))
})

test_that("the subclasses add semantics but no structure", {
  df <- minimal_distances_df()
  expect_s3_class(tale_distances(df), "tale_distances")
  expect_s3_class(tale_distances(df), "pairwise_distances")
  expect_s3_class(domain_distances(df), "domain_distances")
  expect_s3_class(domain_distances(df), "pairwise_distances")
  # identical payload, differing only in class
  expect_equal(as.data.frame(tale_distances(df)), as.data.frame(domain_distances(df)))
})


#### Column contract ####

test_that("the required columns are enforced", {
  df <- minimal_distances_df()
  expect_error(pairwise_distances(df[c("id1", "id2")]),
               class = "tantale_error_distances_missing_column")
  expect_error(pairwise_distances(df[c("id1", "dissim")]),
               class = "tantale_error_distances_missing_column")
})

test_that("dissim must be numeric and ids must not be NA", {
  df <- minimal_distances_df()
  df$dissim <- as.character(df$dissim)
  expect_error(pairwise_distances(df), class = "tantale_error_distances_type")

  df2 <- minimal_distances_df()
  df2$id1[1] <- NA_character_
  expect_error(pairwise_distances(df2), class = "tantale_error_distances_na")
})

test_that("extra columns are preserved", {
  df <- minimal_distances_df()
  df$arlem_score <- 1
  df$note <- "x"
  x <- pairwise_distances(df)
  expect_true(all(c("arlem_score", "note") %in% names(x)))
})

test_that("a non-square table is still valid - squareness is a precondition", {
  df <- minimal_distances_df()
  asymmetric <- df[df$id1 == "a", ]
  expect_s3_class(pairwise_distances(asymmetric), "pairwise_distances")
})


#### Legacy renaming ####

test_that("TAL1/TAL2 and Sim are renamed", {
  df <- legacy_sim_df()
  names(df) <- c("TAL1", "TAL2", "Sim")
  x <- tale_distances(df)
  expect_setequal(names(x), c("id1", "id2", "dissim"))
  # a legacy Sim is converted to the distance, and not kept alongside it
  expect_false("sim" %in% names(x))
  expect_equal(x$dissim, 100 - df$Sim)
})

test_that("RepU1/RepU2 are renamed despite their reversed column order", {
  df <- legacy_sim_df()
  # repeat.similarity really does list RepU2 first
  df <- df[c("id2", "id1", "sim")]
  names(df) <- c("RepU2", "RepU1", "Sim")
  x <- domain_distances(df)
  expect_setequal(names(x), c("id1", "id2", "dissim"))
  # renaming is by name, so the reversed order does not swap the ids
  expect_identical(x$id1, as.character(df$RepU1))
})

test_that("the tal.similarity extras are renamed to snake_case", {
  df <- minimal_distances_df()
  names(df) <- c("TAL1", "TAL2", "Sim")
  df$arlemScore <- 1
  df$maxLength <- 2
  df$normArlemScore <- 3
  x <- tale_distances(df)
  expect_true(all(c("arlem_score", "max_length", "norm_arlem_score") %in% names(x)))
})


#### Preconditions ####

test_that("distances_assert_square() accepts a complete table", {
  expect_silent(distances_assert_square(pairwise_distances(minimal_distances_df())))
})

test_that("distances_assert_square() rejects an asymmetric filter", {
  x <- pairwise_distances(minimal_distances_df())
  expect_error(distances_assert_square(x[x$id1 == "a", ]),
               class = "tantale_error_distances_not_square")
})

test_that("filtering both id columns keeps a table square", {
  x <- pairwise_distances(minimal_distances_df())
  both <- x[x$id1 %in% c("a", "b") & x$id2 %in% c("a", "b"), ]
  expect_silent(distances_assert_square(both))
})


#### Views ####

test_that("as.matrix() builds the square matrix with sorted dimnames", {
  x <- pairwise_distances(minimal_distances_df())
  m <- as.matrix(x)
  expect_equal(dim(m), c(3L, 3L))
  expect_identical(rownames(m), c("a", "b", "c"))
  expect_identical(colnames(m), c("a", "b", "c"))
  expect_identical(diag(m), c(a = 0, b = 0, c = 0))
})

test_that("as.matrix() matches the acast() call it replaces", {
  x <- pairwise_distances(minimal_distances_df())
  expected <- reshape2::acast(as.data.frame(x), id1 ~ id2, value.var = "dissim")
  expect_equal(as.matrix(x), expected)
})

test_that("as.matrix() errors on a non-square table and on a missing column", {
  x <- pairwise_distances(minimal_distances_df())
  expect_error(as.matrix(x[x$id1 == "a", ]),
               class = "tantale_error_distances_not_square")
  expect_error(as.matrix(x, value = "sim"),
               class = "tantale_error_distances_layer")
})

test_that("distances_restrict() keeps the table square", {
  x <- pairwise_distances(minimal_distances_df())
  out <- distances_restrict(x, c("a", "b"))
  expect_s3_class(out, "pairwise_distances")
  expect_identical(nrow(out), 4L)
  expect_silent(distances_assert_square(out))
})

test_that("distances_restrict() errors on an unknown id", {
  x <- pairwise_distances(minimal_distances_df())
  expect_error(distances_restrict(x, c("a", "zzz")),
               class = "tantale_error_distances_unknown_id")
})


#### dplyr policy ####

test_that("row subsetting preserves the class", {
  x <- pairwise_distances(minimal_distances_df())
  expect_s3_class(dplyr::filter(x, dissim == 0), "pairwise_distances")
  expect_s3_class(x[1:3, ], "pairwise_distances")
})

test_that("dropping a required column degrades silently to a tibble", {
  x <- tale_distances(minimal_distances_df())
  out <- dplyr::select(x, -dissim)
  expect_false(is_pairwise_distances(out))
  expect_s3_class(out, "tbl_df")
})

test_that("the namespace tag is carried", {
  x <- domain_distances(minimal_distances_df(), dom_code_namespace = "run-abc")
  expect_identical(tales_namespace(x), "run-abc")
  expect_identical(tales_namespace(dplyr::filter(x, dissim == 0)), "run-abc")
})


#### Against real pipeline output ####

test_that("the real distalr similarity tables validate", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))

  rs <- domain_distances(out$repeat.similarity)
  ts <- tale_distances(out$tal.similarity)

  expect_s3_class(rs, "domain_distances")
  expect_s3_class(ts, "tale_distances")
  expect_setequal(names(rs), c("id1", "id2", "dissim"))
  expect_true(all(c("arlem_score", "max_length", "norm_arlem_score") %in% names(ts)))

  # both are complete squares, as recorded in the design doc
  expect_silent(distances_assert_square(rs))
  expect_silent(distances_assert_square(ts))
  expect_equal(dim(as.matrix(rs)), c(251L, 251L))
  expect_equal(dim(as.matrix(ts)), c(44L, 44L))
})

test_that("as.matrix() reproduces what the existing call sites compute", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))

  # classification.R:24 does 100 - acast(tal_sim, TAL1 ~ TAL2, value.var = "Sim")
  legacy <- 100 - reshape2::acast(out$tal.similarity, TAL1 ~ TAL2, value.var = "Sim")
  # as.matrix() now yields the distance directly, so no inversion is needed
  viaClass <- as.matrix(tale_distances(out$tal.similarity))
  expect_equal(viaClass, legacy)
})
