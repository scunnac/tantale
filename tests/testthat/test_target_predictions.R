

test_that("talvez output a tibble with the expected dims", {
  talvezPreds <- talvez(rvd_seqs = system.file("extdata", "Sample_TALEs_RVDSeqs_AnnoTALE.fasta",
                                              package = "tantale", mustWork = T),
                        subj_file = system.file("extdata", "cladeIII_sweet_promoters.fasta",
                                                     package = "tantale", mustWork = T),
                        opt_param = "-t 0 -l 19",
                        conda_bin = "auto")
  expect_true(identical(dim(talvezPreds), c(90L,9L)))
})

