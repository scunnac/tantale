# run_annotale_predict() and run_annotale_build() on the toy regions: two
# copies of one TALE, one of them carrying a frameshift
# (data-raw/make_toy_tale_regions.R). Needs a Java runtime; fails when it
# is missing.

toy_fasta <- test_path("data_for_tests", "toy_tal_regions.fasta")
predict_out <- file.path(tempdir(), "test_annotale_predict")
unlink(predict_out, recursive = TRUE)
predict_status <- suppressMessages(
  run_annotale_predict(toy_fasta, output_dir = predict_out, prefix = "toy")
)

test_that("run_annotale_predict() returns 0 and writes both stages", {
  expect_identical(predict_status, 0L)
  expect_true(file.exists(file.path(predict_out, "Predict", "protocol_predict.txt")))
  expect_true(all(file.exists(file.path(predict_out, "Analyze",
    c("TALE_RVDs.fasta", "TALE_Protein_parts.fasta", "TALE_DNA_parts.fasta")))))
})

test_that("predict finds both TALEs and flags the frameshifted copy", {
  dna <- Biostrings::readDNAStringSet(list.files(
    file.path(predict_out, "Predict"), "^TALE_DNA_sequences_", full.names = TRUE))
  expect_length(dna, 2L)
  expect_true(all(startsWith(names(dna), "toy-")))  # the prefix given
  expect_identical(sum(grepl("(Pseudo)", names(dna), fixed = TRUE)), 1L)
})

test_that("run_annotale_build() puts the two copies in one class", {
  predicted <- list.files(file.path(predict_out, "Predict"),
                          pattern = "^TALE_DNA_sequences_", full.names = TRUE)
  build_out <- withr::local_tempdir()
  expect_identical(suppressMessages(run_annotale_build(predicted, output_dir = build_out)),
                   0L)
  classes <- list.dirs(build_out, recursive = FALSE, full.names = FALSE)
  expect_identical(classes, "Class_1")
  expect_true(file.exists(file.path(build_out, "Class_builder.xml")))
})
