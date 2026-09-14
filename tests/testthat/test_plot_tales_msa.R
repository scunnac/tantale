
# load(file.path("/home/cunnac/TEMP/test_tantale/mining.RData"))
# saveRDS(distalr_deci_output, file = testthat::test_path("data_for_tests", "sampleDistalrOutput.rds"))
# saveRDS(grp, file = testthat::test_path("data_for_tests", "sampleDistalrGroups.rds"))
# saveRDS(repeatMsaByGroup_withSim, file = testthat::test_path("data_for_tests", "repeatMsaByGroup.rds"))

distalrOut <- readRDS(file = testthat::test_path("data_for_tests", "sampleDistalrOutput.rds"))
#taleGroups <- readRDS(file = testthat::test_path("data_for_tests", "sampleDistalrGroups.rds"))
repeatMsaByGroup <- readRDS(file = testthat::test_path("data_for_tests", "sampleRepeatMsaByGroup.rds"))

repeat_align <- repeatMsaByGroup[[6]]
rvd_align <- repeat_to_rvd_align(repeat_align = repeat_align,
                                   rvd_map = repeat_to_rvd_map_distalr(distalrOut$tale_parts))


tales_consensus(repeat_align)
tales_consensus(rvd_align)

tales_consensus_match(repeat_align)


try(plot_tales_msa(repeat_align = NULL,
                   tal_sim = NULL,
                   repeat_sim = NULL,
                   rvd_align = NULL,
                   h_cut = 90,
                   ref_pattern = NULL,
                   consensus = FALSE,
                   fill_type = "repeat_sim" # "repeatClust"
))

try(plot_tales_msa(repeat_align = repeat_align[3,],
                   tal_sim = NULL,
                   repeat_sim = distalrOut$repeat.similarity,
                   rvd_align = rvd_align[3, , drop = FALSE],
                   h_cut = 90,
                   ref_pattern = NULL,
                   consensus = FALSE,
                   fill_type = "repeat_sim" # "repeatClust"
))

try(plot_tales_msa(repeat_align = repeat_align[3,, drop = FALSE],
                   tal_sim = NULL,
                   repeat_sim = distalrOut$repeat.similarity,
                   rvd_align = rvd_align[3, ],
                   h_cut = 90,
                   ref_pattern = NULL,
                   consensus = FALSE,
                   fill_type = "repeat_sim" # "repeatClust"
))


##############""














# --- Assertions ------------------------------------------------------------
# Everything above is exploratory script kept from development. Until the
# ggplot2 4.x `palette` fix, plot_tales_msa() aborted on every call, so none of
# it could have asserted anything; these are the first real checks.

test_that("plot_tales_msa() returns a ggplot for both fill types", {
  m <- repeatMsaByGroup[[which(sapply(repeatMsaByGroup,
                                      function(z) is.matrix(z) && nrow(z) > 2))[1]]]
  for (ft in c("repeat_clust", "repeat_sim")) {
    p <- suppressWarnings(suppressMessages(plot_tales_msa(
      repeat_align = m,
      repeat_sim = distalrOut$repeat.similarity,
      fill_type = ft
    )))
    expect_s3_class(p, "ggplot")
  }
})

test_that("the returned plot actually renders", {
  m <- repeatMsaByGroup[[which(sapply(repeatMsaByGroup,
                                      function(z) is.matrix(z) && nrow(z) > 2))[1]]]
  p <- suppressWarnings(suppressMessages(plot_tales_msa(
    repeat_align = m, repeat_sim = distalrOut$repeat.similarity
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
  skip_if_not(file.exists(test_path("data_for_tests", "sampleRepeatMsaByGroup.rds")))
  m <- readRDS(test_path("data_for_tests", "sampleRepeatMsaByGroup.rds"))[[4]]
  without <- suppressMessages(plot_tales_msa(repeat_align = m, consensus = FALSE))
  with    <- suppressMessages(plot_tales_msa(repeat_align = m, consensus = TRUE))
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
  expect_identical(unique(panel$data$arrayID), "Consensus")
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
  p <- suppressMessages(plot_tales_msa(repeat_align = f$m, rvd_align = f$rvd,
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

test_that("rvd_sim needs an rvd_align", {
  f <- fixture_rvd()
  expect_error(plot_tales_msa(repeat_align = f$m, fill_type = "rvd_sim"),
               class = "tantale_error_msa_layer")
})

test_that("an unknown fill_type is refused and names the valid ones", {
  f <- fixture_rvd()
  expect_error(
    suppressMessages(plot_tales_msa(repeat_align = f$m, rvd_align = f$rvd,
                                    repeat_sim = f$d$repeat.similarity,
                                    fill_type = "nonsense")),
    class = "tantale_error_msa_layer")
})
