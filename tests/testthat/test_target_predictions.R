

test_that("talvez output a tibble with the expected dims", {
  talvezPreds <- talvez(rvd_seqs = system.file("extdata", "Sample_TALEs_RVDSeqs_AnnoTALE.fasta",
                                              package = "tantale", mustWork = T),
                        subj_file = system.file("extdata", "cladeIII_sweet_promoters.fasta",
                                                     package = "tantale", mustWork = T),
                        opt_param = "-t 0 -l 19",
                        conda_bin = "auto")
  expect_true(identical(dim(talvezPreds), c(90L,9L)))
})



test_that("tales_predict_targets() accepts a tales object and records the method", {
  x <- suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
  preds <- suppressWarnings(suppressMessages(tales_predict_targets(
    x,
    subj_file = system.file("extdata", "cladeIII_sweet_promoters.fasta",
                            package = "tantale", mustWork = TRUE),
    method = "talvez"
  )))
  expect_s3_class(preds, "tbl_df")
  expect_true("method" %in% names(preds))
  expect_true(all(preds$method == "talvez"))
  # predictions are for the arrays we supplied
  expect_true(all(preds$taleId %in% unique(x$array_id)))
})

test_that("tales_predict_targets() matches calling the backend directly", {
  rvds <- system.file("extdata", "Sample_TALEs_RVDSeqs_AnnoTALE.fasta",
                      package = "tantale", mustWork = TRUE)
  subj <- system.file("extdata", "cladeIII_sweet_promoters.fasta",
                      package = "tantale", mustWork = TRUE)
  direct <- suppressMessages(talvez(rvd_seqs = rvds, subj_file = subj,
                                    opt_param = "-t 0 -l 19"))
  viaGeneric <- suppressMessages(tales_predict_targets(
    rvds, subj_file = subj, method = "talvez", opt_param = "-t 0 -l 19"
  ))
  expect_equal(viaGeneric[names(direct)], direct)
})
