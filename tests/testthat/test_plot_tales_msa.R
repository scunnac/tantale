
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
