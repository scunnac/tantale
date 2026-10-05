

test_that("send message if no hmmer hit", {
  # a random dna sequence file
  fasta <- tempfile()
  Biostrings::DNAStringSet(x = c(random = paste(sample(Biostrings::DNA_BASES, size = 10000, replace = TRUE), collapse = ""))) %>%
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

test_that("the closing log is printed line by line, braces in paths as text", {
  # cli_inform() rewrapped the log into one paragraph and evaluated {...}
  out <- file.path(tempfile(), "run{1}")
  msgs <- testthat::capture_messages(suppressWarnings(tell_tales(
    subject_file = system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                               package = "tantale", mustWork = TRUE),
    output_dir = out)))
  printed <- unlist(strsplit(paste(msgs, collapse = ""), "\n", fixed = TRUE))
  log <- readLines(file.path(out, "tell_tales.log"))
  expect_true(grep("^nterm_min_score:", log, value = TRUE) %in% printed)
  expect_true(paste0("Output directory:\t", out) %in% printed)
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
  x <- suppressWarnings(tales_from_telltales(out))
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

test_that("canonical termini match their profile; truncated and unrelated segments do not", {
  hits <- .tale_termini_hmmsearch(termini_cases(), max_evalue = 1e-5,
                                  hmm_dir = system.file("extdata", "hmmProfile", package = "tantale"))
  hit <- function(id, col) hits[[col]][hits$array_id == id]
  expect_true(hit("MAI1_ROI_00006", "nterm_aa_hit"))
  # a truncTALE C-terminus (42 aa) matches strongly but covers too little
  # of the profile to be canonical (ledger §59)
  expect_false(hit("PXO86_truncTALE", "cterm_aa_hit"))
  expect_lt(hit("PXO86_truncTALE", "cterm_aa_evalue"), 1e-10)
  expect_false(hit("BAI3-1-1_ROI_00001", "cterm_aa_hit"))
  expect_false(hit("BAI3-1-1_ROI_00006", "nterm_aa_hit"))
  # both reach the profile end next to the repeats
  expect_identical(hit("MAI1_ROI_00006", "nterm_aa_profile_gap"), 0L)
  expect_identical(hit("PXO86_truncTALE", "cterm_aa_profile_gap"), 0L)
  # the far end tells them apart: the complete N-terminus covers its whole
  # profile in one stretch, the truncated C-terminus stops after 42 residues
  expect_identical(hit("MAI1_ROI_00006", "nterm_aa_far_gap"), 0L)
  expect_identical(hit("MAI1_ROI_00006", "nterm_aa_domains"), 1L)
  expect_gt(hit("MAI1_ROI_00006", "nterm_aa_cover"), 0.95)
  expect_gt(hit("PXO86_truncTALE", "cterm_aa_far_gap"), 200L)
  expect_lt(hit("PXO86_truncTALE", "cterm_aa_cover"), 0.2)
  # with a looser cover threshold it would count
  loose <- .tale_termini_hmmsearch(termini_cases(), max_evalue = 1e-5, min_cover = 0.1,
                                   hmm_dir = system.file("extdata", "hmmProfile", package = "tantale"))
  expect_true(loose$cterm_aa_hit[loose$array_id == "PXO86_truncTALE"])
  expect_gt(hit("PXO86_truncTALE", "cterm_aa_score"), 0)
  # no segment on the other side: NA, not FALSE
  expect_true(is.na(hit("MAI1_ROI_00006", "cterm_aa_hit")))
  expect_true(is.na(hit("PXO86_truncTALE", "nterm_aa_evalue")))
})

test_that("a terminus that a frameshift cuts short of the repeats is no match", {
  # BAI3-1-1 raw assembly: the 287-aa N-terminus matches the profile strongly
  # up to position 150 of 288, then reads in another frame (ledger §42)
  hits <- .tale_termini_hmmsearch(termini_cases(), max_evalue = 1e-5,
                                  hmm_dir = system.file("extdata", "hmmProfile", package = "tantale"))
  row <- hits[hits$array_id == "BAI3-1-1_ROI_00002", ]
  expect_lt(row$nterm_aa_evalue, 1e-80)
  expect_identical(row$nterm_aa_profile_gap, 138L)
  expect_false(row$nterm_aa_hit)
  # the tolerances decide it: allowing 138 positions, and a match over half
  # the profile, makes it a match
  loose <- .tale_termini_hmmsearch(termini_cases(), max_evalue = 1e-5, max_profile_gap = 138L,
                                   min_cover = 0.5,
                                   hmm_dir = system.file("extdata", "hmmProfile", package = "tantale"))
  expect_true(loose$nterm_aa_hit[loose$array_id == "BAI3-1-1_ROI_00002"])
})

test_that("the terminus check needs the protein profiles", {
  expect_error(.tale_termini_hmmsearch(termini_cases(), max_evalue = 1e-5, hmm_dir = tempdir()),
               class = "tantale_error_hmm_missing")
})

test_that("tell_tales() codes termini from the protein profiles and drops DNA-only AnnoTALE output", {
  out <- tempfile()
  # the terminus/repeat overlaps are normal, and the duplicate repeat hits
  # merged; the array on roi18_region is extended past its end, then clipped (§36)
  suppressWarnings(expect_no_warning(expect_no_warning(suppressMessages(tell_tales(
    subject_file = test_path("data_for_tests", "pxo86_roi18_19_excerpt.fa"),
    output_dir = out)), class = "tantale_warning_overlapping_hits"),
    message = "out-of-bound"))
  report <- readr::read_tsv(file.path(out, "array_report.tsv"),
                            show_col_types = FALSE, progress = FALSE)
  expect_true(all(c("n_dna_hits", "nterm_dna_hit", "cterm_dna_hit", "nterm_aa_evalue",
                    "cterm_aa_evalue", "nterm_aa_hit", "cterm_aa_hit") %in% names(report)))
  truncated <- report[report$cterm_aa_length == 42 & !is.na(report$cterm_aa_length), ]
  expect_identical(nrow(truncated), 1L)
  # the nhmmer DNA search misses this truncated C-terminus; the protein profile does not
  expect_false(truncated$cterm_dna_hit)
  expect_false(truncated$cterm_aa_hit)   # not canonical (ledger §59)
  # below the score threshold, the DNA still carries the two ends of a
  # C-terminus, split by a deletion (ledger §59)
  expect_true(all(c("nterm_dna_score", "cterm_dna_score", "nterm_dna_cover", "cterm_dna_cover",
                    "nterm_dna_pieces", "cterm_dna_pieces") %in% names(report)))
  expect_identical(truncated$cterm_dna_pieces, 2)
  expect_lt(truncated$cterm_dna_score, 200)
  expect_gt(truncated$cterm_dna_cover, 0.15)
  expect_lt(truncated$cterm_dna_cover, 0.3)
  hits <- readr::read_tsv(file.path(out, "hits_report.tsv"), show_col_types = FALSE, progress = FALSE)
  expect_true(all(c("score", "evalue", "hmm_from", "hmm_to") %in% names(hits)))
  expect_false(anyNA(hits$score))
  expect_match(truncated$rvd_string, "-XXXXX$")

  # AnnoTALE splits ROI_00001's DNA but cannot translate it: no parts file is left
  expect_false(file.exists(file.path(out, "annotale", "ROI_00001", "TALE_DNA_parts.fasta")))
  expect_false(file.exists(file.path(out, "annotale", "ROI_00001", "TALE_Protein_parts.fasta")))
  # and reading the directory says so (§49)
  expect_warning(x <- tales_from_telltales(out), "ROI_00001",
                 class = "tantale_warning_annotale_unparsed")
  expect_false("ROI_00001" %in% x$array_id)
  expect_identical(x$rvd[x$array_id == truncated$array_id & x$domain_type == "C-terminus"], "XXXXX")
})


#### Subject preparation and hit ranges (§36) ####

test_that("subject preparation keeps full headers and writes nothing next to the input", {
  dir <- withr::local_tempdir()
  fa <- file.path(dir, "subject.fa")
  writeLines(c(">ctg1 plasmid pXO1", "ACGTACGTAC", ">ctg2 chromosome", "ACGTAC"), fa)
  subject <- suppressMessages(.telltale_prepare_subject(fa))
  expect_identical(list.files(dir), "subject.fa")
  expect_identical(GenomeInfoDb::seqnames(subject$seqinfo), c("ctg1 plasmid pXO1", "ctg2 chromosome"))
  expect_identical(unname(GenomeInfoDb::seqlengths(subject$seqinfo)), c(10L, 6L))
})

test_that("subject preparation refuses duplicated or empty sequence names", {
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">ctg1", "ACGT", ">ctg1", "ACGT"), fa)
  expect_error(suppressMessages(.telltale_prepare_subject(fa)),
               class = "tantale_error_seqnames")
  writeLines(c(">ctg1", "ACGT", ">", "ACGT"), fa)
  expect_error(suppressMessages(.telltale_prepare_subject(fa)),
               class = "tantale_error_seqnames")
})

test_that("hit ranges carry the original names and lengths, quietly, when a sequence has no hit", {
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">ctg1 plasmid pXO1", strrep("ACGT", 50), ">ctg2 chromosome", strrep("ACGT", 25)), fa)
  subject <- suppressMessages(.telltale_prepare_subject(fa))
  # hits on seq1 only
  hits <- data.frame(target_name = "seq1", start = c(11, 51), end = c(40, 90),
                     strand = "+", query_name = "repeat", hit_id = c("DOM_1", "DOM_2"))
  expect_no_warning(gr <- .telltale_hits_to_ranges(hits, NULL, subject$seqlevels, subject$seqinfo))
  expect_identical(GenomeInfoDb::seqlevels(gr), "ctg1 plasmid pXO1")
  expect_identical(unname(GenomeInfoDb::seqlengths(gr)), 200L)
})


#### Overlapping hits within an array (§36) ####

overlap_case <- function(ranges, types) {
  gr <- GenomicRanges::GRanges("ctg1", IRanges::IRanges(ranges[, 1], ranges[, 2]), strand = "+",
                               query_name = types, hit_id = paste0("DOM_", seq_along(types)))
  list(gr = gr, seqs = Biostrings::DNAStringSet(c(ctg1 = strrep("ACGT", 100))),
       hmm = list(nterm = "nterm", repeats = "repeat", cterm = "cterm"))
}

test_that("a terminus hit overlapping the adjacent repeat hit is not reported", {
  x <- overlap_case(rbind(c(1, 100), c(97, 198), c(199, 300), c(281, 380)),
                    c("nterm", "repeat", "repeat", "cterm"))
  expect_no_warning(out <- .telltale_group_arrays(x$gr, min_gap = 50, subject_seqs = x$seqs, hmm = x$hmm))
  expect_identical(unname(S4Vectors::mcols(out$by_array)$n_dna_hits), 4L)
})

test_that("overlapping hits of the same domain type are reported", {
  x <- overlap_case(rbind(c(1, 100), c(101, 202), c(150, 251), c(252, 350)),
                    c("nterm", "repeat", "repeat", "cterm"))
  expect_warning(.telltale_group_arrays(x$gr, min_gap = 50, subject_seqs = x$seqs, hmm = x$hmm),
                 class = "tantale_warning_overlapping_hits")
})

