
# load(file.path("/home/cunnac/TEMP/test_tantale/mining.RData"))
# saveRDS(distalr_deci_output, file = testthat::test_path("data_for_tests", "sampleDistalrOutput.rds"))
# saveRDS(grp, file = testthat::test_path("data_for_tests", "sampleDistalrGroups.rds"))
# saveRDS(repeatMsaByGroup_withSim, file = testthat::test_path("data_for_tests", "repeatMsaByGroup.rds"))

load_all()
distalrOut <- readRDS(file = testthat::test_path("data_for_tests", "sampleDistalrOutput.rds"))
#taleGroups <- readRDS(file = testthat::test_path("data_for_tests", "sampleDistalrGroups.rds"))
repeatMsaByGroup <- readRDS(file = testthat::test_path("data_for_tests", "repeatMsaByGroup.rds"))

repeatAlign <- repeatMsaByGroup[[6]]
rvdAlign <- convertRepeat2RvdAlign(repeatAlign = repeatAlign,
                                   repeat2RvdMapping = getRepeat2RvdMappingFromDistalr(distalrOut$taleParts))


taleAlignConsensus(repeatAlign)
taleAlignConsensus(rvdAlign)

matchConsensus(repeatAlign, taleAlignConsensus(repeatAlign))


try(ggplotTalesMsa(repeatAlign = NULL,
                   talsim = NULL,
                   repeatSim = NULL,
                   rvdAlign = NULL,
                   repeat.clust.h.cut = 90,
                   refgrep = NULL,
                   consensusSeq = FALSE,
                   fillType = "repeatSim" # "repeatClust"
))

try(ggplotTalesMsa(repeatAlign = repeatAlign[3,],
                   talsim = NULL,
                   repeatSim = distalrOut$repeat.similarity,
                   rvdAlign = rvdAlign[3, , drop = FALSE],
                   repeat.clust.h.cut = 90,
                   refgrep = NULL,
                   consensusSeq = FALSE,
                   fillType = "repeatSim" # "repeatClust"
))

try(ggplotTalesMsa(repeatAlign = repeatAlign[3,, drop = FALSE],
                   talsim = NULL,
                   repeatSim = distalrOut$repeat.similarity,
                   rvdAlign = rvdAlign[3, ],
                   repeat.clust.h.cut = 90,
                   refgrep = NULL,
                   consensusSeq = FALSE,
                   fillType = "repeatSim" # "repeatClust"
))


##############""












