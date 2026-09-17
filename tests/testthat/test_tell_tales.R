

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
