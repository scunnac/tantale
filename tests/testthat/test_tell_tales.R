

test_that("send message if no hmmer hit", {
  # a random dna sequence file
  fasta <- tempfile()
  Biostrings::DNAStringSet(x = paste(sample(Biostrings::DNA_BASES, size = 10000, replace = TRUE), collapse = "")) %>%
  Biostrings::writeXStringSet(filepath = fasta)
  expect_warning(tell_tales(subject_file = fasta, output_dir = tempfile()),
                 class = "tantale_warning_no_hits")
})

test_that("telltale no correction runs without error", {
  expect_invisible(tell_tales(subject_file = system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                                                      package = "tantale", mustWork = T),
                            output_dir = tempfile()
                            )
                   )
})

# PXO86's ROI_00019 is a genuine truncTALE: its stop codon falls inside the
# C-terminal part, so AnnoTALE's record ends in "*". The report's lengths
# count residues, like the tales object, so the two must agree (ledger §32.4).
test_that("array_report.tsv terminus lengths do not count the stop codon", {
  out <- tempfile()
  suppressWarnings(suppressMessages(tell_tales(
    subject_file = test_path("data_for_tests", "pxo86_roi18_19_excerpt.fa"),
    output_dir = out)))
  report <- readr::read_tsv(file.path(out, "array_report.tsv"),
                            show_col_types = FALSE, progress = FALSE)
  expect_true(all(c("cterm_aa_length", "nterm_aa_length") %in% names(report)))
  x <- suppressWarnings(tales_from_telltale(out))
  cterm <- x[x$domain_type == "C-terminus", ]
  expect_gt(nrow(cterm), 0)
  expect_equal(report$cterm_aa_length[match(cterm$array_id, report$array_id)],
               nchar(cterm$aa_seq))
  expect_true(42 %in% report$cterm_aa_length)   # the truncated C-terminus
})

test_that(".aa_residue_count() ignores stop codons", {
  s <- Biostrings::AAStringSet(c(a = "RRKRS*", b = "SVGGTI", c = ""))
  expect_identical(.aa_residue_count(s), c(5L, 6L, 0L))
})
