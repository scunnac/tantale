

test_that("talvez output a tibble with the expected dims", {
  talvezPreds <- talvez(rvdSeqs = system.file("extdata", "Sample_TALEs_RVDSeqs_AnnoTALE.fasta", package = "tantale", mustWork = T),
                        subjDnaSeqFile = system.file("extdata", "cladeIII_sweet_promoters.fasta", package = "tantale", mustWork = T),
                        optParam = "-t 0 -l 19",
                        condaBinPath = "/home/cunnac/bin/miniconda3/condabin/conda")
  expect_true(identical(dim(talvezPreds), c(90L,9L)))
})

