

test_that("talvez output a tibble with the expected dims", {
  talvezPreds <- talvez(rvd_seqs = system.file("extdata", "Sample_TALEs_RVDSeqs_AnnoTALE.fasta",
                                              package = "tantale", mustWork = T),
                        subj_file = system.file("extdata", "cladeIII_sweet_promoters.fasta",
                                                     package = "tantale", mustWork = T),
                        opt_param = "-t 0 -l 19",
                        conda_bin = "auto")
  expect_identical(dim(talvezPreds), c(90L,9L))
  expect_named(talvezPreds, c("tale_id", "rvds", "subj_seq_id", "start", "end",
                              "strand", "ebe_seq", "score", "rank"))
})



test_that("tales_predict_targets() accepts a tales object and records the method", {
  x <- suppressWarnings(
    tales_from_telltales(test_path("data_for_tests", "example_output"))
  )
  preds <- suppressWarnings(suppressMessages(tales_predict_targets(
    x,
    subj_file = system.file("extdata", "cladeIII_sweet_promoters.fasta",
                            package = "tantale", mustWork = TRUE),
    method = "talvez"
  )))
  expect_s3_class(preds, "tbl_df")
  expect_named(preds, c("tale_id", "rvds", "subj_seq_id", "start", "end",
                        "strand", "ebe_seq", "score", "method", "rank"))
  expect_true(all(preds$method == "talvez"))
  # predictions are for the arrays we supplied
  expect_true(all(preds$tale_id %in% unique(x$array_id)))
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
  expect_identical(viaGeneric[names(direct)], direct)
})


# preditale() and plot_target_preds() -----------------------------------
# PrediTALE needs a Java runtime; these fail when it is missing.

sweet <- system.file("extdata", "cladeIII_sweet_promoters.fasta",
                     package = "tantale", mustWork = TRUE)
rvd_file <- system.file("extdata", "Sample_TALEs_RVDSeqs_AnnoTALE.fasta",
                        package = "tantale", mustWork = TRUE)
pt_preds <- suppressMessages(preditale(rvd_file, subj_file = sweet,
                                       output_dir = tempfile()))

test_that("preditale() returns the documented columns, one row per site", {
  expect_s3_class(pt_preds, "tbl_df")
  expect_named(pt_preds, c("tale_id", "rvds", "subj_seq_id", "start", "end",
                           "strand", "ebe_seq", "score", "pval"))
  expect_gt(nrow(pt_preds), 0)
  expect_true(all(pt_preds$tale_id %in% names(Biostrings::readBStringSet(rvd_file))))
})

test_that("each predicted site covers position 0 plus one base per RVD", {
  n_rvds <- lengths(strsplit(pt_preds$rvds, "-", fixed = TRUE))
  expect_identical(as.integer(pt_preds$end - pt_preds$start + 1L),
                   as.integer(n_rvds + 1L))
  expect_identical(nchar(pt_preds$ebe_seq), as.integer(n_rvds + 1L))
})

test_that("preditale()'s EBE sequences are the subject's bases at those coordinates", {
  gr <- GenomicRanges::makeGRangesFromDataFrame(pt_preds, seqnames.field = "subj_seq_id")
  expect_identical(as.character(BSgenome::getSeq(Biostrings::readDNAStringSet(sweet), gr)),
                   pt_preds$ebe_seq)
})

test_that("preditale() gives the same sites from a file or from sequences", {
  from_set <- suppressMessages(preditale(Biostrings::readBStringSet(rvd_file),
                                         subj_file = sweet, output_dir = tempfile()))
  key <- function(p) p[order(p$tale_id, p$subj_seq_id, p$start, p$strand), ]
  expect_equal(key(from_set), key(pt_preds), ignore_attr = TRUE)
})

test_that("preditale()'s opt_param reaches PrediTALE: forward strand only", {
  fwd <- suppressMessages(preditale(rvd_file, subj_file = sweet,
                                    opt_param = 'Strand="forward strand"',
                                    output_dir = tempfile()))
  expect_gt(nrow(fwd), 0)
  expect_true(all(fwd$strand == "+"))
  expect_true(any(pt_preds$strand == "-"))  # the default searches both
})

test_that("preditale() refuses a directory that already holds its results", {
  out <- tempfile()
  suppressMessages(preditale(rvd_file, subj_file = sweet, output_dir = out))
  expect_error(suppressMessages(preditale(rvd_file, subj_file = sweet, output_dir = out)),
               class = "tantale_error_output_not_empty")
})

test_that("preditale() refuses an rvd_seqs of the wrong type", {
  expect_error(preditale(42, subj_file = sweet),
               class = "tantale_error_rvd_seqs")
})

best <- pt_preds[order(-pt_preds$score), ][1, ]
window <- paste0(best$subj_seq_id, ":", best$start - 30, "-", best$end + 30)

test_that("plot_target_preds() draws one label per RVD and the DNA in the window", {
  p <- plot_target_preds(best, subj_file = sweet, filter_range = window)
  expect_s3_class(p, "ggplot")
  expect_no_warning(ggplot2::ggplot_build(p))
  # the RVD layer: one cell per RVD, plus position 0
  n_rvds <- length(strsplit(best$rvds, "-", fixed = TRUE)[[1]])
  expect_identical(nrow(p$data), n_rvds + 1L)
  expect_equal(range(p$data$xPos), c(best$start, best$end))
  # the sequence layer: both strands over the whole window
  dna <- p$layers[[2]]$data
  expect_equal(nrow(dna), 2 * (best$end - best$start + 61))
})

test_that("plot_target_preds() refuses predictions made on other sequences", {
  wrong <- best
  wrong$ebe_seq <- paste(rev(strsplit(wrong$ebe_seq, "", fixed = TRUE)[[1]]), collapse = "")
  expect_error(plot_target_preds(wrong, subj_file = sweet, filter_range = window),
               class = "tantale_error_ebe_mismatch")
})

test_that("plot_target_preds() wants exactly one range", {
  expect_error(plot_target_preds(best, subj_file = sweet,
                                 filter_range = c(window, window)),
               class = "tantale_error_bad_argument")
})
