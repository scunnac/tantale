# Every external program the package runs must stop the call when it exits
# with a non-zero status, rather than leave empty or partial output for
# later stages to trip over (ledger START HERE, 2026-09-23). Each test makes
# the program fail on purpose: a jar path that does not exist, or an HMM
# directory with no profiles in it. They need Java and the tantale conda
# environment, like the functions they test, and fail when those are
# missing.

toy_fasta <- test_path("data_for_tests", "toy_tal_regions.fasta")
no_jar <- file.path(tempdir(), "no_such_tool.jar")

test_that(".annotale_exec() stops with AnnoTALE's class and shows its stderr", {
  err <- expect_error(
    .annotale_exec("sh -c 'echo \"Exception {oops}\" >&2; exit 3'", "analyze",
                   quiet = TRUE),
    class = "tantale_error_annotale_failed")
  msg <- conditionMessage(err)
  expect_match(msg, "AnnoTALE analyze failed with exit status 3", fixed = TRUE)
  # braces in the program's own text are shown, not interpolated
  expect_match(msg, "Exception {oops}", fixed = TRUE)
})

test_that(".annotale_exec() returns 0 invisibly on success", {
  expect_invisible(res <- .annotale_exec("true", "build"))
  expect_identical(res, 0L)
})

test_that("run_annotale_predict() stops when AnnoTALE fails", {
  expect_error(
    suppressMessages(run_annotale_predict(toy_fasta, output_dir = tempfile(),
                                          annotale_jar = no_jar)),
    class = "tantale_error_annotale_failed")
})

test_that("run_annotale_build() stops when AnnoTALE fails", {
  expect_error(
    suppressMessages(run_annotale_build(toy_fasta, output_dir = tempfile(),
                                        annotale_jar = no_jar)),
    class = "tantale_error_annotale_failed")
})

test_that("a java_args Java refuses stops each jar wrapper", {
  bad <- "-XX:+NoSuchOption"
  expect_error(
    suppressMessages(run_annotale_predict(toy_fasta, output_dir = tempfile(),
                                          java_args = bad)),
    class = "tantale_error_annotale_failed")
  expect_error(
    suppressMessages(run_annotale_build(toy_fasta, output_dir = tempfile(),
                                        java_args = bad)),
    class = "tantale_error_annotale_failed")
  expect_error(
    suppressMessages(preditale(Biostrings::BStringSet(c(a = "NI-HD-NG-NN-NI-HD")),
                               subj_file = toy_fasta, output_dir = tempfile(),
                               java_args = bad)),
    class = "tantale_error_exec_failed")
  expect_error(
    suppressMessages(correct_tales(toy_fasta, corrected_path = tempfile(),
                                   java_args = bad)),
    class = "tantale_error_exec_failed")
})

test_that("tell_tales()'s AnnoTALE step fails with a class its caller catches", {
  expect_error(
    .run_annotale_analyze(toy_fasta, output_dir = tempfile(),
                          annotale_jar = no_jar),
    class = "tantale_error_annotale_failed")
})

test_that("the nHMMER search stops when nhmmer fails", {
  out <- withr::local_tempdir()
  expect_error(
    suppressMessages(.run_nhmmer_search(
      subject_file = toy_fasta,
      hmm_file = file.path(out, "missing.hmm"),
      search_out_file = file.path(out, "search.txt"),
      readable_out_file = file.path(out, "readable.txt"))),
    class = "tantale_error_exec_failed")
})

test_that("preditale() stops when PrediTALE fails", {
  rvds <- Biostrings::BStringSet(c(a = "NI-HD-NG-NN-NI-HD"))
  expect_error(
    suppressMessages(preditale(rvds, subj_file = toy_fasta,
                               output_dir = tempfile(),
                               predictor_path = no_jar)),
    class = "tantale_error_exec_failed")
})

test_that("correct_tales() stops when any of its nHMMER searches fails", {
  expect_error(
    suppressMessages(correct_tales(toy_fasta, corrected_path = tempfile(),
                                   hmm_path = withr::local_tempdir())),
    class = "tantale_error_exec_failed")
})
