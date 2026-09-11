
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












