

# a random dna sequence file
fasta <- tempfile()
Biostrings::DNAStringSet(x = paste(sample(Biostrings::DNA_BASES, size = 10000, replace = TRUE), collapse = "")) %>%
  Biostrings::writeXStringSet(filepath = fasta)

test_that("send message if no hmmer hit", {
  expect_warning(tellTale(subjectFile = fasta, outputDir = tempfile()),
                 regexp = "NhmmerSearch found no TALE cds hit")
})

test_that("telltale no correction runs without error", {
  expect_invisible(tellTale(subjectFile = system.file("extdata", "bai3_sample_tal_regions.fasta", package = "tantale", mustWork = T),
                          outputDir = tempfile()
                          )
                 )
})
