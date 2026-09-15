# Regression baseline for the whole pipeline.
#
# These tests assert nothing about what the values *should* be. They record
# what the pipeline currently produces, so that a refactor which is supposed to
# change nothing can be shown to have changed nothing -- and so that one which
# does change something says exactly where.
#
# This replaces an ad-hoc baseline kept outside the repository. It earned its
# place: it caught tales_consensus() breaking ties by row order, a bug that no
# other test would have noticed, because the consensus was still *a* plausible
# value.
#
# When a change here is intended, inspect the diff, satisfy yourself that every
# line of it is something you meant, then accept it with
# testthat::snapshot_accept("golden").
#
# MAFFT and arlem are required. If they are missing these tests fail rather
# than skip: a baseline that quietly does not run is worse than none.

fx <- function() {
  path <- test_path("data_for_tests", "sampleDistalrOutput.rds")
  expect_true(file.exists(path))
  readRDS(path)
}


#### the tales object itself ####

test_that("golden: the tales column contract", {
  x <- tales(fx()$tale_parts)
  # Small and worth reading in full: this is the package's central contract.
  expect_golden(names(x))
  expect_golden(vapply(x, typeof, character(1)))
  expect_golden(fingerprint(x))
})

test_that("golden: anomalies reported for the reference fixture", {
  x <- tales(fx()$tale_parts)
  expect_golden(as.data.frame(tales_anomalies(x)))
})

test_that("golden: the consumer requirements table", {
  expect_golden(as.data.frame(tales_requirements()))
})


#### projections ####

test_that("golden: projections of a tales onto strings and maps", {
  d <- fx()
  x <- tales(d$tale_parts)
  expect_golden(fingerprint(tales_rvd_strings(x)))
  expect_golden(fingerprint(tales_coded_strings(x)))
  expect_golden(fingerprint(tales_domain_codes(x)))
  expect_golden(fingerprint(repeat_to_rvd_map_distalr(d$tale_parts)))
  expect_golden(fingerprint(tale_parts_to_rvd(d$tale_parts)))
})


#### the expensive paths ####

test_that("golden: tales_compare() on four arrays", {
  d <- fx()
  x <- tales(d$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:4], ]
  out <- suppressWarnings(suppressMessages(tales_compare(sub)))

  expect_named(out, c("tales", "domain_distances", "tale_distances"))
  expect_golden(fingerprint(out$tales))
  expect_golden(fingerprint(out$domain_distances))
  expect_golden(fingerprint(out$tale_distances))
})

test_that("golden: tales_align() on both residue layers", {
  d <- fx()
  x <- tales(d$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:4], ]
  byRvd  <- suppressWarnings(suppressMessages(tales_align(sub, residue_col = "rvd")))
  byCode <- suppressWarnings(suppressMessages(tales_align(sub, residue_col = "dom_code")))

  expect_golden(fingerprint(byRvd))
  expect_golden(fingerprint(byCode))
  # width is what the back-mapping depends on, so pin it separately
  expect_golden(c(rvd = tales_width(byRvd), dom_code = tales_width(byCode)))
})

test_that("golden: tales_group() partitions the arrays the same way", {
  expect_golden(fingerprint(suppressWarnings(suppressMessages(
      tales_group(fx()$tal.similarity, k = 4, method = "hclust")))))
})


#### what the figure is built from ####

test_that("golden: the data behind plot() on a tales_msa", {
  # Every fill and label layer at once, so a change to any of them shows up.
  d <- fx()
  msa <- readRDS(test_path("data_for_tests", "sampleTalesMsa.rds"))
  p <- suppressMessages(plot(msa, fill = "dom_code", label = "rvd",
                             domain_sim = d$repeat.similarity,
                             fill_type = "repeat_clust"))
  layer <- if (is.null(p$plotlist)) p$data else p$plotlist[[1]]$data
  expect_golden(fingerprint(layer))
})

test_that("golden: the consensus of the reference alignment", {
  # Small enough to read, and the place a silent change does the most damage:
  # a wrong consensus still looks like a consensus.
  msa <- readRDS(test_path("data_for_tests", "sampleTalesMsa.rds"))
  expect_golden(tales_consensus(as.matrix(msa, value = "rvd")))
  expect_golden(tales_consensus(as.matrix(msa, value = "dom_code")))
})


#### tell_tales(), the pipeline entry point ####

# The refactor of restructuring-notes.md 5.3 is meant to change nothing about
# what tell_tales() writes. This is what says so. Before it existed the
# function had two tests, both of which only checked that it did not error --
# no cover at all for 745 lines.
#
# Correction is left off here: it is 29x the rest of the pipeline and scales
# with the 1057-sequence reference set, not with the subject (ledger 8.1). The
# correction branch gets its own baseline once the toy fixture exists.

test_that("golden: tell_tales() writes the same files with the same contents", {
  out <- file.path(tempdir(), "golden_telltale")
  unlink(out, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  res <- suppressWarnings(suppressMessages(tell_tales(
    subject_file = system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                               package = "tantale", mustWork = TRUE),
    output_dir = out)))
  expect_identical(res, out)

  expect_golden(telltale_fingerprint(out))
})

test_that("golden: the tables tell_tales() writes, column by column", {
  # The digests above say "something changed"; these say which column, which
  # is what saves the time when it does.
  out <- file.path(tempdir(), "golden_telltale_tables")
  unlink(out, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  suppressWarnings(suppressMessages(tell_tales(
    subject_file = system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                               package = "tantale", mustWork = TRUE),
    output_dir = out)))

  for (f in c("hitsReport.tsv", "domainsReport.tsv", "arrayReport.tsv")) {
    tbl <- readr::read_tsv(file.path(out, f), show_col_types = FALSE,
                           progress = FALSE)
    expect_golden(fingerprint(as.data.frame(tbl)))
  }
})

test_that("golden: a tell_tales() run loads back as a tales object", {
  # The end-to-end contract: what the entry point writes is what the class
  # reads. A refactor that kept every file byte-identical but broke this would
  # still have broken the pipeline.
  out <- file.path(tempdir(), "golden_telltale_roundtrip")
  unlink(out, recursive = TRUE)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  suppressWarnings(suppressMessages(tell_tales(
    subject_file = system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                               package = "tantale", mustWork = TRUE),
    output_dir = out)))

  x <- suppressWarnings(tales_from_telltale(out))
  expect_s3_class(x, "tales")
  # source_directory records where the run happened, so it holds this
  # session's tempdir and changes every time. Keep the part that carries
  # signal -- which ROI each part came from -- and drop the prefix.
  x$source_directory <- basename(x$source_directory)
  expect_golden(fingerprint(x))
})
