# Does frameshift correction actually correct anything? (ledger 8.1)
#
# Until this file existed, nothing asserted that it did. The correction
# branch was covered only by a golden digest, which pins the code path but
# would happily keep passing if the correction silently stopped working --
# the digest would just record the new, wrong answer once accepted.
#
# The fixture is built by data-raw/make_toy_tale_regions.R and carries its
# own answer key. Three regions cut from the shipped BAI3 sequences:
#
#   toy_intact      one complete TALE, untouched
#   toy_frameshift  the same TALE with one nucleotide inserted mid-array
#   toy_no_tale     a stretch of the same genome with no TALE in it
#
# Two copies of the *same* TALE is the point: the intact one is the control,
# so any difference between them is attributable to the inserted base rather
# than to the two arrays being different TALEs.
#
# MAFFT, HMMER and arlem are required. If they are missing these fail rather
# than skip.

toy_subject <- function() test_path("data_for_tests", "toy_tal_regions.fasta")
toy_truth   <- function() {
  readr::read_tsv(test_path("data_for_tests", "toy_tal_regions_truth.tsv"),
                  show_col_types = FALSE, progress = FALSE)
}

# Both runs are needed by several expectations and cost a few seconds each,
# so each is computed once per file.
toy_run <- local({
  cache <- list()
  function(correct) {
    key <- if (correct) "corrected" else "plain"
    if (is.null(cache[[key]])) {
      out <- file.path(tempdir(), paste0("toy_", key))
      unlink(out, recursive = TRUE)
      args <- list(subject_file = toy_subject(), output_dir = out)
      if (correct) {
        args$correct_array <- TRUE
        args$correction_ref <- test_path("data_for_tests", "correction_ref_20.fa.gz")
      }
      suppressWarnings(suppressMessages(do.call(tell_tales, args)))
      r <- readr::read_tsv(file.path(out, "arrayReport.tsv"),
                           show_col_types = FALSE, progress = FALSE)
      cache[[key]] <<- list(dir = out, report = as.data.frame(r))
    }
    cache[[key]]
  }
})

# the real arrays, dropping the single-hit fragments each region also yields
toy_arrays <- function(report) {
  a <- report[report$n_domain_hits > 4, ]
  a[match(c("toy_intact", "toy_frameshift"), a$seqnames), ]
}


test_that("the fixture is what its answer key says it is", {
  truth <- toy_truth()
  seqs <- Biostrings::readDNAStringSet(toy_subject())
  expect_setequal(names(seqs), truth$seqname)
  # the frameshifted copy is its twin plus exactly one base
  w <- Biostrings::width(seqs)
  names(w) <- names(seqs)
  expect_equal(unname(w[["toy_frameshift"]]), unname(w[["toy_intact"]]) + 1L)
  # and identical to it up to the insertion point
  at <- truth$insertion_at[truth$seqname == "toy_frameshift"]
  expect_identical(
    as.character(Biostrings::subseq(seqs[["toy_frameshift"]], 1L, at)),
    as.character(Biostrings::subseq(seqs[["toy_intact"]], 1L, at)))
})


test_that("a region with no TALE yields no array", {
  # Pinned on real genomic sequence rather than the random DNA the older
  # test generated: random DNA is an easier negative than the real thing.
  report <- toy_run(FALSE)$report
  expect_false("toy_no_tale" %in% report$seqnames)
})


test_that("the inserted base truncates the ORF when correction is off", {
  # If this stops being true the fixture has lost its point, so it is
  # asserted rather than assumed.
  a <- toy_arrays(toy_run(FALSE)$report)
  intact <- a[a$seqnames == "toy_intact", ]
  shifted <- a[a$seqnames == "toy_frameshift", ]

  expect_gt(intact$longest_orf_length, shifted$longest_orf_length)
  expect_gt(intact$orf_coverage, shifted$orf_coverage)
  # both still found as arrays -- the frameshift breaks the ORF, not the
  # HMMER-level detection of the repeats
  expect_equal(intact$n_domain_hits, shifted$n_domain_hits)
})


test_that("correction recovers the intact TALE from the frameshifted copy", {
  # The assertion the correction branch never had: not that the call
  # returned, but that it produced the right answer.
  a <- toy_arrays(toy_run(TRUE)$report)
  intact <- a[a$seqnames == "toy_intact", ]
  shifted <- a[a$seqnames == "toy_frameshift", ]

  expect_identical(shifted$rvd_string, intact$rvd_string)
  expect_equal(shifted$longest_orf_length, intact$longest_orf_length)
  expect_equal(shifted$orf_coverage, intact$orf_coverage)
})


test_that("correction charges the extra base to the frameshifted copy", {
  # The intact control is not corrected to zero indels -- the 20-sequence
  # test reference is too small for that (ledger 8.1b), and this fixture is
  # not the place to relitigate it. What must hold is the *difference*: one
  # inserted base, one extra insertion called.
  a <- toy_arrays(toy_run(TRUE)$report)
  intact <- a[a$seqnames == "toy_intact", ]
  shifted <- a[a$seqnames == "toy_frameshift", ]

  expect_equal(shifted$predicted_ins_count, intact$predicted_ins_count + 1L)
  expect_equal(shifted$predicted_dels_count, intact$predicted_dels_count)
})
