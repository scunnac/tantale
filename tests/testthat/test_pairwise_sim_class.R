# Tests for the 'pairwise_sim' class. See dev/class-design.md §3.

# A complete 3x3 similarity, long form.
minimal_sim_df <- function() {
  ids <- c("a", "b", "c")
  g <- expand.grid(id1 = ids, id2 = ids, stringsAsFactors = FALSE)
  tibble::tibble(
    id1 = g$id1,
    id2 = g$id2,
    sim = ifelse(g$id1 == g$id2, 100, 50)
  )
}

test_that("pairwise_sim() accepts a complete square table", {
  x <- pairwise_sim(minimal_sim_df())
  expect_s3_class(x, "pairwise_sim")
  expect_s3_class(x, "tbl_df")
  expect_true(is_pairwise_sim(x))
})

test_that("the subclasses add semantics but no structure", {
  df <- minimal_sim_df()
  expect_s3_class(tale_sim(df), "tale_sim")
  expect_s3_class(tale_sim(df), "pairwise_sim")
  expect_s3_class(repeat_sim(df), "repeat_sim")
  expect_s3_class(repeat_sim(df), "pairwise_sim")
  # identical payload, differing only in class
  expect_equal(as.data.frame(tale_sim(df)), as.data.frame(repeat_sim(df)))
})


#### Column contract ####

test_that("the required columns are enforced", {
  df <- minimal_sim_df()
  expect_error(pairwise_sim(df[c("id1", "id2")]),
               class = "tantale_error_sim_missing_column")
  expect_error(pairwise_sim(df[c("id1", "sim")]),
               class = "tantale_error_sim_missing_column")
})

test_that("sim must be numeric and ids must not be NA", {
  df <- minimal_sim_df()
  df$sim <- as.character(df$sim)
  expect_error(pairwise_sim(df), class = "tantale_error_sim_type")

  df2 <- minimal_sim_df()
  df2$id1[1] <- NA_character_
  expect_error(pairwise_sim(df2), class = "tantale_error_sim_na")
})

test_that("extra columns are preserved", {
  df <- minimal_sim_df()
  df$dissim <- 100 - df$sim
  df$note <- "x"
  x <- pairwise_sim(df)
  expect_true(all(c("dissim", "note") %in% names(x)))
})

test_that("a non-square table is still valid - squareness is a precondition", {
  df <- minimal_sim_df()
  asymmetric <- df[df$id1 == "a", ]
  expect_s3_class(pairwise_sim(asymmetric), "pairwise_sim")
})


#### Legacy renaming ####

test_that("TAL1/TAL2 and Sim are renamed", {
  df <- minimal_sim_df()
  names(df) <- c("TAL1", "TAL2", "Sim")
  x <- tale_sim(df)
  expect_setequal(names(x), c("id1", "id2", "sim"))
})

test_that("RepU1/RepU2 are renamed despite their reversed column order", {
  df <- minimal_sim_df()
  # repeat.similarity really does list RepU2 first
  df <- df[c("id2", "id1", "sim")]
  names(df) <- c("RepU2", "RepU1", "Sim")
  x <- repeat_sim(df)
  expect_setequal(names(x), c("id1", "id2", "sim"))
  # renaming is by name, so the reversed order does not swap the ids
  expect_identical(x$id1, as.character(df$RepU1))
})

test_that("the tal.similarity extras are renamed to snake_case", {
  df <- minimal_sim_df()
  names(df) <- c("TAL1", "TAL2", "Sim")
  df$arlemScore <- 1
  df$maxLength <- 2
  df$normArlemScore <- 3
  x <- tale_sim(df)
  expect_true(all(c("arlem_score", "max_length", "norm_arlem_score") %in% names(x)))
})


#### Preconditions ####

test_that("sim_assert_square() accepts a complete table", {
  expect_silent(sim_assert_square(pairwise_sim(minimal_sim_df())))
})

test_that("sim_assert_square() rejects an asymmetric filter", {
  x <- pairwise_sim(minimal_sim_df())
  expect_error(sim_assert_square(x[x$id1 == "a", ]),
               class = "tantale_error_sim_not_square")
})

test_that("filtering both id columns keeps a table square", {
  x <- pairwise_sim(minimal_sim_df())
  both <- x[x$id1 %in% c("a", "b") & x$id2 %in% c("a", "b"), ]
  expect_silent(sim_assert_square(both))
})


#### Views ####

test_that("as.matrix() builds the square matrix with sorted dimnames", {
  x <- pairwise_sim(minimal_sim_df())
  m <- as.matrix(x)
  expect_equal(dim(m), c(3L, 3L))
  expect_identical(rownames(m), c("a", "b", "c"))
  expect_identical(colnames(m), c("a", "b", "c"))
  expect_identical(diag(m), c(a = 100, b = 100, c = 100))
})

test_that("as.matrix() matches the acast() call it replaces", {
  x <- pairwise_sim(minimal_sim_df())
  expected <- reshape2::acast(as.data.frame(x), id1 ~ id2, value.var = "sim")
  expect_equal(as.matrix(x), expected)
})

test_that("as.matrix() errors on a non-square table and on a missing column", {
  x <- pairwise_sim(minimal_sim_df())
  expect_error(as.matrix(x[x$id1 == "a", ]),
               class = "tantale_error_sim_not_square")
  expect_error(as.matrix(x, value = "dissim"),
               class = "tantale_error_sim_layer")
})

test_that("sim_restrict() keeps the table square", {
  x <- pairwise_sim(minimal_sim_df())
  out <- sim_restrict(x, c("a", "b"))
  expect_s3_class(out, "pairwise_sim")
  expect_identical(nrow(out), 4L)
  expect_silent(sim_assert_square(out))
})

test_that("sim_restrict() errors on an unknown id", {
  x <- pairwise_sim(minimal_sim_df())
  expect_error(sim_restrict(x, c("a", "zzz")),
               class = "tantale_error_sim_unknown_id")
})


#### dplyr policy ####

test_that("row subsetting preserves the class", {
  x <- pairwise_sim(minimal_sim_df())
  expect_s3_class(dplyr::filter(x, sim == 100), "pairwise_sim")
  expect_s3_class(x[1:3, ], "pairwise_sim")
})

test_that("dropping a required column degrades silently to a tibble", {
  x <- tale_sim(minimal_sim_df())
  out <- dplyr::select(x, -sim)
  expect_false(is_pairwise_sim(out))
  expect_s3_class(out, "tbl_df")
})

test_that("the namespace tag is carried", {
  x <- repeat_sim(minimal_sim_df(), dom_code_namespace = "run-abc")
  expect_identical(tales_namespace(x), "run-abc")
  expect_identical(tales_namespace(dplyr::filter(x, sim == 100)), "run-abc")
})


#### Against real pipeline output ####

test_that("the real distalr similarity tables validate", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))

  rs <- repeat_sim(out$repeat.similarity)
  ts <- tale_sim(out$tal.similarity)

  expect_s3_class(rs, "repeat_sim")
  expect_s3_class(ts, "tale_sim")
  expect_setequal(names(rs), c("id1", "id2", "dissim", "sim"))
  expect_true(all(c("arlem_score", "max_length", "norm_arlem_score") %in% names(ts)))

  # both are complete squares, as recorded in the design doc
  expect_silent(sim_assert_square(rs))
  expect_silent(sim_assert_square(ts))
  expect_equal(dim(as.matrix(rs)), c(251L, 251L))
  expect_equal(dim(as.matrix(ts)), c(44L, 44L))
})

test_that("as.matrix() reproduces what the existing call sites compute", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))

  # classification.R:24 does 100 - acast(tal_sim, TAL1 ~ TAL2, value.var = "Sim")
  legacy <- 100 - reshape2::acast(out$tal.similarity, TAL1 ~ TAL2, value.var = "Sim")
  viaClass <- 100 - as.matrix(tale_sim(out$tal.similarity))
  expect_equal(viaClass, legacy)
})
