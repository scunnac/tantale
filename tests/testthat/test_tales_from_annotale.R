# tales_from_annotale() reads run_annotale_predict()'s output. The fixture is
# AnnoTALE predict + analyze on the BAI3 sample that tellTaleExampleOutput was
# built from (data-raw/make_annotale_example_output.R).

annotale_example <- function() {
  system.file("extdata", "annotaleExampleOutput", package = "tantale", mustWork = TRUE)
}

# a copy of the fixture's analyze files, optionally without predict's GFF3
annotale_copy <- function(with_gff = TRUE, env = parent.frame()) {
  to <- withr::local_tempdir(.local_envir = env)
  file.copy(file.path(annotale_example(), "Analyze"), to, recursive = TRUE)
  if (with_gff) file.copy(file.path(annotale_example(), "Predict"), to, recursive = TRUE)
  to
}

test_that("tales_from_annotale() reads AnnoTALE's TALEs with their contigs and termini", {
  x <- tales_from_annotale(annotale_example())
  expect_true(is_tales(x))
  expect_identical(sort(unique(x$array_id)),
                   paste0("bai3_sample_tal_genomic_regions-tempTALE", 1:4))
  expect_true("seqnames" %in% names(x))
  expect_setequal(x$seqnames, c("talRegion5", "talRegion6"))
  termini <- x[x$domain_type != "repeat", ]
  expect_identical(nrow(termini), 8L)
  expect_setequal(termini$rvd, c("NTERM", "CTERM"))
})

test_that("tales_from_annotale() and tales_from_telltales() agree on the same TALEs", {
  a <- tales_from_annotale(annotale_example())
  t <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput", package = "tantale"))
  expect_setequal(unname(as.character(tales_rvd_strings(a, repeats_only = FALSE))),
                  unname(as.character(tales_rvd_strings(t, repeats_only = FALSE))))
})

test_that("tales_from_annotale() works without predict's GFF3, leaving out seqnames", {
  x <- tales_from_annotale(annotale_copy(with_gff = FALSE))
  expect_false("seqnames" %in% names(x))
  expect_identical(length(unique(x$array_id)), 4L)
})

test_that("a terminus unlike the TALE profile is coded XXXXX", {
  dir <- annotale_copy()
  f <- file.path(dir, "Analyze", "TALE_Protein_parts.fasta")
  parts <- Biostrings::readAAStringSet(f)
  hit <- grep("tempTALE2 .*: C-terminus$", names(parts))
  expect_length(hit, 1L)
  parts[hit] <- Biostrings::AAStringSet("MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRG")
  Biostrings::writeXStringSet(parts, f)
  expect_warning(x <- tales_from_annotale(dir), class = "tantale_warning_tales_anomalous")
  expect_identical(tales_anomalies(x)$check, "terminus_unmatched")
  cterm <- x[x$domain_type == "C-terminus", ]
  expect_identical(cterm$rvd[cterm$array_id == "bai3_sample_tal_genomic_regions-tempTALE2"], "XXXXX")
  expect_identical(sum(cterm$rvd == "CTERM"), 3L)
})

test_that("tales_from_annotale() refuses a directory without AnnoTALE output", {
  expect_error(tales_from_annotale(withr::local_tempdir()),
               class = "tantale_error_annotale_missing")
})

test_that("a TALE in which AnnoTALE found no repeat is read with its two termini", {
  # AnnoTALE predict on BAI3-1-1 reports such a pseudogene: N- and C-terminus
  # parts, and an empty record in TALE_RVDs.fasta
  dir <- annotale_copy(with_gff = FALSE)
  analyze <- file.path(dir, "Analyze")
  id <- "BAI3-1-1-tempTALE8 (Pseudo) [207114-207183:-1]"
  prot <- Biostrings::readAAStringSet(file.path(analyze, "TALE_Protein_parts.fasta"))
  dna <- Biostrings::readDNAStringSet(file.path(analyze, "TALE_DNA_parts.fasta"))
  nterm <- grep("tempTALE1 .*: N-terminus$", names(prot))
  cterm <- grep("tempTALE1 .*: C-terminus$", names(prot))
  extra <- function(x) stats::setNames(x[c(nterm, cterm)],
                                       paste0(id, c(": N-terminus", ": C-terminus")))
  Biostrings::writeXStringSet(c(prot, extra(prot)), file.path(analyze, "TALE_Protein_parts.fasta"))
  Biostrings::writeXStringSet(c(dna, extra(dna)), file.path(analyze, "TALE_DNA_parts.fasta"))
  cat(">", id, "\n\n", sep = "", file = file.path(analyze, "TALE_RVDs.fasta"), append = TRUE)

  expect_warning(x <- tales_from_annotale(dir), class = "tantale_warning_tales_anomalous")
  pseudo <- x[x$array_id == "BAI3-1-1-tempTALE8", ]
  expect_identical(pseudo$domain_type, c("N-terminus", "C-terminus"))
  an <- tales_anomalies(x)
  expect_identical(an$check[an$array_id == "BAI3-1-1-tempTALE8"], "no_repeat")
  expect_identical(length(unique(x$array_id)), 5L)
})
