distalrOut <- readRDS(file = testthat::test_path("data_for_tests", "sampleDistalrOutput.rds"))
repeatMsaByGroup <- readRDS(file = testthat::test_path("data_for_tests", "sampleRepeatMsaByGroup.rds"))

# The alignment as the class holds it: one object carrying every layer, which
# is what plot() takes. Built by tales_align() from the same three arrays the
# matrix fixture covers.
msa <- readRDS(testthat::test_path("data_for_tests", "sampleTalesMsa.rds"))

# Matrices are still the input to the internal helpers (.consensus_panel(),
# .rvd_to_match_align()), so the matrix fixture stays for those.
repeat_align <- repeatMsaByGroup[[6]]
rvd_align <- repeat_to_rvd_align(repeat_align = repeat_align,
                                 rvd_map = repeat_to_rvd_map_distalr(distalrOut$tale_parts))


#### the arguments are named layers of the object ####

test_that("plot() refuses a layer the alignment does not carry", {
  expect_error(plot(msa, fill = "not_a_layer"), class = "tantale_error_msa_layer")
  expect_error(plot(msa, label = "not_a_layer"), class = "tantale_error_msa_layer")
})

test_that("plot() refuses an empty alignment in terms of the alignment", {
  # Previously this reported a problem with `repeat_align`, an argument a
  # plot() caller never supplies and cannot inspect.
  expect_error(plot(msa[0, ]), class = "tantale_error_msa_empty")
})

test_that("a single-array alignment plots", {
  # The old matrix interface could be handed a vector here, by a caller who
  # subset without drop = FALSE. as.matrix() on the object cannot produce one.
  one <- msa[msa$array_id == unique(msa$array_id)[1], ]
  expect_s3_class(suppressMessages(plot(one)), "ggplot")
})


# --- Assertions ------------------------------------------------------------
# Everything above is exploratory script kept from development. Until the
# ggplot2 4.x `palette` fix, the plotting code aborted on every call, so none of
# it could have asserted anything; these are the first real checks.

test_that("plot() returns a ggplot for both fill types", {
  for (ft in c("repeat_clust", "repeat_sim")) {
    p <- suppressWarnings(suppressMessages(plot(
      msa, domain_sim = distalrOut$repeat.similarity, fill_type = ft
    )))
    expect_s3_class(p, "ggplot")
  }
})

test_that("the returned plot actually renders", {
  p <- suppressWarnings(suppressMessages(plot(
    msa, domain_sim = distalrOut$repeat.similarity
  )))
  f <- withr::local_tempfile(fileext = ".png")
  suppressWarnings(suppressMessages(
    ggplot2::ggsave(f, p, width = 8, height = 3, dpi = 72)
  ))
  expect_true(file.exists(f))
  expect_gt(file.size(f), 1000)
})


#### consensus row ####

test_that("consensus = TRUE adds a panel and consensus = FALSE does not", {
  without <- suppressMessages(plot(msa, consensus = FALSE))
  with    <- suppressMessages(plot(msa, consensus = TRUE))
  # with a consensus the result is an aplot composition carrying an extra panel
  expect_gt(length(with), length(without))
})

test_that("the consensus panel labels the same layer as the cells", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleRepeatMsaByGroup.rds")))
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  m <- readRDS(test_path("data_for_tests", "sampleRepeatMsaByGroup.rds"))[[4]]
  rvd <- repeat_to_rvd_align(repeat_align = m,
                             rvd_map = repeat_to_rvd_map_distalr(d$tale_parts))
  # rvd_align supplied -> consensus must be of the RVDs, not the repeat codes
  panel <- tantale:::.consensus_panel(rvd, n_positions = ncol(rvd))
  expect_s3_class(panel, "ggplot")
  expect_identical(nrow(panel$data), ncol(rvd))
  expect_identical(unique(panel$data$array_id), "Consensus")
})

test_that(".consensus_panel() reproduces tales_consensus(), with terminus relabelling", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleRepeatMsaByGroup.rds")))
  m <- readRDS(test_path("data_for_tests", "sampleRepeatMsaByGroup.rds"))[[4]]
  expected <- tales_consensus(m)
  expected <- gsub("NTERM", "N-", expected)
  expected <- gsub("CTERM", "-C", expected)
  panel <- tantale:::.consensus_panel(m, n_positions = ncol(m))
  expect_identical(panel$data$label, expected)
})

test_that("the consensus panel pads repeat codes exactly as the cells do", {
  m <- matrix(c("1", "22", "333", "1", "22", "333"), nrow = 2, byrow = TRUE,
              dimnames = list(c("a", "b"), NULL))
  padded <- tantale:::.consensus_panel(m, n_positions = 3, pad = TRUE)
  expect_identical(padded$data$label, stringr::str_pad(c("1", "22", "333"), 3, "left"))
  bare <- tantale:::.consensus_panel(m, n_positions = 3, pad = FALSE)
  expect_identical(bare$data$label, c("1", "22", "333"))
})


#### fill_type = "rvd_sim" ####

fixture_rvd <- function() {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  skip_if_not(file.exists(test_path("data_for_tests", "sampleRepeatMsaByGroup.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  m <- readRDS(test_path("data_for_tests", "sampleRepeatMsaByGroup.rds"))[[4]]
  list(d = d, m = m,
       rvd = repeat_to_rvd_align(repeat_align = m,
                                 rvd_map = repeat_to_rvd_map_distalr(d$tale_parts)))
}

test_that("rvd_sim fills by RVD specificity relative to the reference", {
  f <- fixture_rvd()
  p <- suppressMessages(plot(msa, fill = "dom_code", label = "rvd",
                             fill_type = "rvd_sim"))
  layer <- if (is.null(p$plotlist)) p$data else p$plotlist[[1]]$data
  expect_true("rvdSimVsRef" %in% names(layer))
  # a correlation, so bounded and signed -- unlike the 0-100 repeat similarity
  expect_lte(max(layer$rvdSimVsRef, na.rm = TRUE), 1)
  expect_gte(min(layer$rvdSimVsRef, na.rm = TRUE), -1)
  expect_lt(min(layer$rvdSimVsRef, na.rm = TRUE), 0)   # some pair is anti-correlated
})

test_that("the reference row scores 1 against itself throughout", {
  f <- fixture_rvd()
  ref <- tantale:::.pick_ref_name(f$rvd, ref_tag = NULL)
  sc <- tantale:::.rvd_to_match_align(f$rvd)
  expect_true(all(sc[ref, ] == 1, na.rm = TRUE))
})

test_that("opposite specificities score strongly negative", {
  # NG binds T (5/10/1/50), NN binds A and G (30/10/30/1): the RVD view must
  # show these as opposed, where a repeat-level view shows only "different"
  f <- fixture_rvd()
  sc <- tantale:::.rvd_to_match_align(f$rvd)
  expect_lt(min(sc, na.rm = TRUE), -0.9)
})

test_that("rvd_sim needs a labelled RVD layer", {
  expect_error(plot(msa, label = NULL, fill_type = "rvd_sim"),
               class = "tantale_error_msa_layer")
})

test_that("an unknown fill_type is refused and names the valid ones", {
  expect_error(
    suppressMessages(plot(msa, domain_sim = distalrOut$repeat.similarity,
                          fill_type = "nonsense")),
    class = "tantale_error_msa_layer")
})

test_that("the tree panel is built from either distance vocabulary", {
  # plot() reads id1/id2/dissim; pairwise_distances() lets the older
  # TAL1/TAL2/Sim spelling in at the door. Both must give the same figure.
  d <- distalrOut
  legacy <- suppressMessages(plot(msa, tal_sim = d$tal.similarity))
  canonical <- suppressMessages(plot(msa, tal_sim = tale_distances(d$tal.similarity)))
  expect_s3_class(legacy, "aplot")
  # a tree panel was actually added, not silently skipped
  expect_true(any(vapply(legacy$plotlist, function(p) inherits(p, "ggtree"), logical(1))))

  tips <- function(p) {
    tr <- p$plotlist[[which(vapply(p$plotlist, function(q) inherits(q, "ggtree"), logical(1)))]]
    tr$data$label[tr$data$isTip][order(tr$data$y[tr$data$isTip])]
  }
  expect_identical(tips(legacy), tips(canonical))
  # the tree must order on distance, not on its inverse: the two nearest TALEs
  # are neighbouring leaves
  # the reference row is marked with a trailing _#
  expect_setequal(sub("_#$", "", tips(legacy)), unique(msa$array_id))
})

#### consensus is a property of the alignment, not of its row order ####

test_that("tales_consensus() does not depend on the order of the rows", {
  m <- readRDS(test_path("data_for_tests", "sampleRepeatMsaByGroup.rds"))[[4]]
  reference <- tales_consensus(m)
  set.seed(20260915)
  for (i in 1:8) {
    permuted <- m[sample(nrow(m)), , drop = FALSE]
    expect_identical(tales_consensus(permuted), reference)
  }
})

test_that("a column where every array differs is broken deterministically", {
  # Three arrays, three distinct repeats: no majority exists. Whatever is
  # reported must at least be the same for the same alignment.
  m <- matrix(c("b", "a", "c",
                "x", "x", "y"), nrow = 3,
              dimnames = list(c("t1", "t2", "t3"), NULL))
  expect_identical(tales_consensus(m)[1], "a")   # smallest in sort order
  expect_identical(tales_consensus(m)[2], "x")   # a genuine majority is unaffected
  expect_identical(tales_consensus(m[c(3, 1, 2), ]), tales_consensus(m))
})

#### tales_consensus_match() ####

test_that("tales_consensus_match() returns logicals, not the strings TRUE/FALSE", {
  # It used to assign TRUE into the character matrix it was handed, which
  # stores "TRUE"; sum(), which() and ! then all did the wrong thing, despite
  # the documented return being a logical matrix.
  m <- matrix(c("HD", "NI", "HD",
                "NG", "NG", "NG"), nrow = 3,
              dimnames = list(c("t1", "t2", "t3"), NULL))
  w <- tales_consensus_match(m, long = FALSE)
  expect_type(w, "logical")
  expect_identical(dim(w), dim(m))
  expect_identical(dimnames(w), dimnames(m))
  expect_identical(sum(w), 5L)          # only t2's NI differs from consensus HD
  expect_type(tales_consensus_match(m, long = TRUE)$tales_consensus_match, "logical")
})

test_that("a gap never counts as matching the consensus", {
  m <- matrix(c("HD", NA, "HD",
                NA, "NI", NA), nrow = 3,
              dimnames = list(c("t1", "t2", "t3"), NULL))
  w <- tales_consensus_match(m, long = FALSE)
  # column 1: consensus HD, so the gap in t2 is FALSE and never NA
  expect_identical(unname(w[, 1]), c(TRUE, FALSE, TRUE))
  expect_false(anyNA(w))
  # column 2 is mostly gap, so its consensus is a gap and nothing matches it
  expect_identical(unname(w[, 2]), c(FALSE, FALSE, FALSE))
})
