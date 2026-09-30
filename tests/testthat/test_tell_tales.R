

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


#### Terminus check with the protein profiles ####

termini_cases <- function() {
  s <- Biostrings::readAAStringSet(test_path("data_for_tests", "termini_profile_cases.fa"))
  part <- sub("^\\S+ (\\S+) .*$", "\\1", names(s))
  names(s) <- sub(" .*$", "", names(s))
  list(`N-terminus` = s[part == "N-terminus"], `C-terminus` = s[part == "C-terminus"])
}

test_that("genuine termini match their profile, complete or truncated; unrelated segments do not", {
  hits <- .tale_termini_hmmsearch(termini_cases(), max_evalue = 1e-5,
                                  hmm_dir = system.file("extdata", "hmmProfile", package = "tantale"))
  hit <- function(id, col) hits[[col]][hits$array_id == id]
  expect_true(hit("MAI1_ROI_00006", "nterm_aa_hit"))
  expect_true(hit("PXO86_truncTALE", "cterm_aa_hit"))            # 42 aa
  expect_lt(hit("PXO86_truncTALE", "cterm_aa_evalue"), 1e-10)
  expect_false(hit("BAI3-1-1_ROI_00001", "cterm_aa_hit"))
  expect_false(hit("BAI3-1-1_ROI_00006", "nterm_aa_hit"))
  # no segment on the other side: NA, not FALSE
  expect_true(is.na(hit("MAI1_ROI_00006", "cterm_aa_hit")))
  expect_true(is.na(hit("PXO86_truncTALE", "nterm_aa_evalue")))
})

test_that("the terminus check needs the protein profiles", {
  expect_error(.tale_termini_hmmsearch(termini_cases(), max_evalue = 1e-5, hmm_dir = tempdir()),
               class = "tantale_error_hmm_missing")
})

test_that("tell_tales() codes termini from the protein profiles and drops DNA-only AnnoTALE output", {
  out <- tempfile()
  suppressWarnings(suppressMessages(tell_tales(
    subject_file = test_path("data_for_tests", "pxo86_roi18_19_excerpt.fa"),
    output_dir = out)))
  report <- readr::read_tsv(file.path(out, "array_report.tsv"),
                            show_col_types = FALSE, progress = FALSE)
  expect_true(all(c("n_dna_hits", "nterm_dna_hit", "cterm_dna_hit", "nterm_aa_evalue",
                    "cterm_aa_evalue", "nterm_aa_hit", "cterm_aa_hit") %in% names(report)))
  truncated <- report[report$cterm_aa_length == 42 & !is.na(report$cterm_aa_length), ]
  expect_equal(nrow(truncated), 1)
  # the nhmmer DNA search misses this truncated C-terminus; the protein profile does not
  expect_false(truncated$cterm_dna_hit)
  expect_true(truncated$cterm_aa_hit)
  expect_match(truncated$rvd_string, "-CTERM$")

  # AnnoTALE splits ROI_00001's DNA but cannot translate it: no parts file is left
  expect_false(file.exists(file.path(out, "annotale", "ROI_00001", "TALE_DNA_parts.fasta")))
  expect_false(file.exists(file.path(out, "annotale", "ROI_00001", "TALE_Protein_parts.fasta")))
  x <- tales_from_telltale(out)
  expect_false("ROI_00001" %in% x$array_id)
  expect_identical(x$rvd[x$array_id == truncated$array_id & x$domain_type == "C-terminus"], "CTERM")
})
