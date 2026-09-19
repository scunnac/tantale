# Tests for tales_compare_functal(). See ledger §12b: a reimplementation of
# QueTAL's FuncTAL on universalmotif, not a port -- results are expected to
# differ from the original Perl tool, verified empirically before writing
# the function, not assumed.

example_tales <- function() {
  suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
}


#### shape ####

test_that("returns a tale_distances object, square, over all arrays", {
  x <- example_tales()
  out <- tales_compare_functal(x)
  expect_s3_class(out, "tale_distances")
  n <- dplyr::n_distinct(x$array_id)
  expect_equal(nrow(out), n^2)
  expect_setequal(unique(out$id1), unique(x$array_id))
  expect_silent(distances_assert_square(out))
})

test_that("self-comparison is distance 0", {
  out <- tales_compare_functal(example_tales())
  self <- out$dissim[out$id1 == out$id2]
  expect_true(all(abs(self) < 1e-8))
})

test_that("the result plugs directly into tales_group_hclust()", {
  # The point of returning a tale_distances object rather than a bare
  # matrix: it is interchangeable with tales_compare_distal()'s output.
  x <- example_tales()
  out <- tales_compare_functal(x)
  grp <- tales_group_hclust(x, out, k = 2)
  expect_s3_class(grp, "tales")
  expect_equal(length(unique(grp$group)), 2L)
})

test_that("needs a tales object", {
  expect_error(tales_compare_functal(data.frame(a = 1)),
               class = "tantale_error_tales_type")
})


#### RVD handling ####

test_that("only rvd drives the comparison, in repeat order, termini dropped", {
  # tales_rvd_strings()'s own default (rvd_only = TRUE) does the dropping;
  # this pins that tales_compare_functal() relies on that default rather
  # than passing rvd_only = FALSE itself.
  x <- example_tales()
  rvd <- tales_rvd_strings(x)
  expect_false(any(grepl("NTERM|CTERM", rvd)))
})

test_that("an RVD absent from the specificity table is scored, not dropped", {
  # A made-up RVD ("ZZ") falls back to the flat "XX" row rather than
  # erroring or silently vanishing from the array.
  df <- data.frame(
    array_id = c("A1", "A1", "A2", "A2"),
    position_in_array = c(1L, 2L, 1L, 2L),
    domain_type = c("repeat", "repeat", "repeat", "repeat"),
    rvd = c("NI", "ZZ", "NI", "HD")
  )
  x <- suppressWarnings(tales(df))
  expect_no_error(out <- tales_compare_functal(x))
  expect_equal(nrow(out), 4L)
})

test_that("an all-terminus array (no repeats) errors rather than being silently dropped", {
  # tales_rvd_strings() itself just omits an array with nothing to render --
  # fine for a projection, not fine for a comparison that must cover every
  # array it was given.
  df <- data.frame(
    array_id = c("A1", "A1", "A2"), position_in_array = c(1L, 2L, 1L),
    domain_type = c("N-terminus", "repeat", "N-terminus"),
    rvd = c("NTERM", "NI", "NTERM")
  )
  x <- suppressWarnings(tales(df))
  err <- expect_error(tales_compare_functal(x), class = "tantale_error_functal_empty")
  expect_match(conditionMessage(err), "A2")
})


#### method direction: similarity vs distance metrics ####

test_that("a similarity metric (PCC) is inverted to a proper dissim", {
  x <- example_tales()
  out <- tales_compare_functal(x, method = "PCC")
  # Distinct arrays: dissim should not be 0 (they are not identical PWMs
  # in this fixture), and self-comparisons are exactly 0 (tested above).
  cross <- out$dissim[out$id1 != out$id2]
  expect_true(all(cross > 0))
})

test_that("a distance metric (EUCL) is used as-is, not inverted", {
  x <- example_tales()
  out <- tales_compare_functal(x, method = "EUCL")
  self <- out$dissim[out$id1 == out$id2]
  expect_true(all(abs(self) < 1e-8))
  cross <- out$dissim[out$id1 != out$id2]
  expect_true(all(cross > 0))
})


#### divergence from the original FuncTAL is real, not a magnitude difference ####

test_that("compare_motifs() PCC and FuncTAL's own padded-flatten PCC disagree", {
  # Pins the empirical finding the roxygen docs assert: these are different
  # statistics, not two parameterisations of the same one. Reproduced small
  # here rather than trusted from memory.
  m1 <- matrix(c(50,10,5,5, 5,10,1,50, 1,1,1,1), nrow = 4,
              dimnames = list(c("A","C","G","T"), NULL))
  m2 <- matrix(c(50,10,5,5, 15,50,5,5), nrow = 4,
              dimnames = list(c("A","C","G","T"), NULL))
  mot1 <- universalmotif::create_motif(m1, alphabet = "DNA", type = "PCM", name = "m1")
  mot2 <- universalmotif::create_motif(m2, alphabet = "DNA", type = "PCM", name = "m2")
  um_pcc <- universalmotif::compare_motifs(list(mot1, mot2), method = "PCC",
                                           min.overlap = 1, tryRC = FALSE)[1, 2]

  functal_style_pcc <- function(m1, m2) {
    len1 <- ncol(m1); len2 <- ncol(m2)
    best <- -Inf
    for (offset2 in seq(len1 - 1, -(len2 - 1))) {
      mm1 <- cbind(matrix(0, 4, max(0, -offset2)), m1, matrix(0, 4, max(0, offset2 + len2 - len1)))
      mm2 <- cbind(matrix(0, 4, max(0, offset2)), m2)
      w <- max(ncol(mm1), ncol(mm2))
      if (ncol(mm1) < w) mm1 <- cbind(mm1, matrix(0, 4, w - ncol(mm1)))
      if (ncol(mm2) < w) mm2 <- cbind(mm2, matrix(0, 4, w - ncol(mm2)))
      r <- suppressWarnings(cor(as.vector(mm1), as.vector(mm2)))
      if (!is.na(r) && r > best) best <- r
    }
    best
  }
  ft_pcc <- functal_style_pcc(m1, m2)
  expect_false(isTRUE(all.equal(um_pcc, ft_pcc)))
})
