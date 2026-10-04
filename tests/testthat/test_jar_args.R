# opt_param and java_args of the jar wrappers (ledger §51). The checks run
# before any program starts, so these tests need neither Java nor the tools.

toy_fasta <- test_path("data_for_tests", "toy_tal_regions.fasta")

test_that("opt_param may not set a key the wrapper sets itself", {
  expect_error(run_annotale_predict(toy_fasta, output_dir = tempfile(),
                                    opt_param = "Sensitive=true s=other"),
               class = "tantale_error_jar_args")
  expect_error(run_annotale_build(toy_fasta, output_dir = tempfile(),
                                  opt_param = "c=3 outdir=/tmp/x"),
               class = "tantale_error_jar_args")
  err <- expect_error(preditale(toy_fasta, subj_file = toy_fasta,
                                opt_param = "s=other.fa"),
                      class = "tantale_error_jar_args")
  # the message names the argument to use instead
  expect_match(conditionMessage(err), "subj_file", fixed = TRUE)
})

test_that("a key's letters inside another option's value are not a key", {
  # build's own s= (significance level) is not one of its reserved keys
  expect_null(.check_jar_args("c=3 s=0.05", "", reserved = c(t = "fasta_file"),
                              fn = "run_annotale_build"))
  expect_null(.check_jar_args('Strand="forward strand" sl=1e-5', "",
                              reserved = c(s = "subj_file"), fn = "preditale"))
  expect_error(.check_jar_args('Strand="forward strand" s=x', "",
                               reserved = c(s = "subj_file"), fn = "preditale"),
               class = "tantale_error_jar_args")
})

test_that("opt_param and java_args must be single strings", {
  expect_error(run_annotale_build(toy_fasta, output_dir = tempfile(),
                                  opt_param = c("c=3", "s=0.05")),
               class = "tantale_error_jar_args")
  expect_error(correct_tales(toy_fasta, java_args = NA_character_),
               class = "tantale_error_jar_args")
  expect_error(preditale(toy_fasta, subj_file = toy_fasta, java_args = 2),
               class = "tantale_error_jar_args")
})
