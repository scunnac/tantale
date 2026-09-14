

test_that("build_repeat_msa result with RVDs is of expected dims",
          {aln <- .build_repeat_msa(input_seqs = system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",package = "tantale", mustWork = TRUE),
                                 sep = "-", repeat_sims = NULL,
                                 mafft_opts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
                                 mafft_path = system.file("tools", "mafft-linux64",package = "tantale", mustWork = TRUE),
                                 gap_symbol = "-")
          expect_true(identical(dim(aln), c(11L, 26L)))}
)


test_that("build_repeat_msa result with coded parts is of expected dims",
          {aln <- .build_repeat_msa(input_seqs = system.file("extdata", "small_Out_CodedRepeats.fa",package = "tantale", mustWork = TRUE),
                                 sep = " ", repeat_sims = NULL,
                                 mafft_opts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
                                 mafft_path = system.file("tools", "mafft-linux64",package = "tantale", mustWork = TRUE),
                                 gap_symbol = "-")
          expect_true(identical(dim(aln), c(3L, 28L)))}
)


test_that("build_repeat_msa throughts a warning and return a NA matrix if provided with empty sequences",
          {seqs <- Biostrings::BStringSet()
          expect_warning(aln <- .build_repeat_msa(input_seqs = seqs,
                                               sep = " ", repeat_sims = NULL,
                                               mafft_opts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
                                               mafft_path = system.file("tools", "mafft-linux64",package = "tantale", mustWork = TRUE),
                                               gap_symbol = "-"))
          expect_equal(aln, matrix())}
)

test_that("build_repeat_msa deals properly with single sequence inputs",
          {seqs <- Biostrings::readBStringSet(system.file("extdata", "small_Out_CodedRepeats.fa",package = "tantale", mustWork = TRUE))
          aln <- .build_repeat_msa(input_seqs = seqs[1],
                                sep = " ", repeat_sims = NULL,
                                mafft_opts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
                                mafft_path = system.file("tools", "mafft-linux64",package = "tantale", mustWork = TRUE),
                                gap_symbol = "-")
          expect_true(identical(dim(aln), c(1L, 28L)))}
)



#### repeat_sims vocabulary ####

test_that(".as_mafft_score_table() accepts the canonical and legacy vocabularies", {
  # MAFFT scores matches, so it wants a similarity. A domain_distances stores
  # the distance, and must be inverted on the way in.
  canonical <- data.frame(id1 = c("a","b"), id2 = c("b","a"), dissim = c(0, 40))
  out <- tantale:::.as_mafft_score_table(canonical)
  expect_named(out, c("RepU1", "RepU2", "Sim"))
  expect_equal(out$Sim, c(100, 60))

  legacy <- data.frame(RepU1 = c("a","b"), RepU2 = c("b","a"), Sim = c(100, 60))
  expect_equal(tantale:::.as_mafft_score_table(legacy)$Sim, c(100, 60))
})

test_that(".as_mafft_score_table() refuses a table it cannot read", {
  expect_error(tantale:::.as_mafft_score_table(data.frame(a = 1, b = 2)),
               class = "tantale_error_msa_sim_table")
  expect_error(tantale:::.as_mafft_score_table(data.frame(id1 = "a", id2 = "b")),
               class = "tantale_error_msa_sim_table")
})

test_that("tales_align() accepts a domain_distances for repeat_sims", {
  # Regression: .build_repeat_msa() hardcoded RepU1/RepU2/Sim, so the
  # documented path -- handing it the typed object tales_compare() returns --
  # failed with "Can't subset columns that don't exist".
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:3], ]
  viaClass <- suppressWarnings(suppressMessages(
    tales_align(sub, residue_col = "dom_code",
                repeat_sims = domain_distances(d$repeat.similarity))))
  viaLegacy <- suppressWarnings(suppressMessages(
    tales_align(sub, residue_col = "dom_code",
                repeat_sims = d$repeat.similarity)))
  expect_s3_class(viaClass, "tales_msa")
  expect_equal(as.data.frame(viaClass), as.data.frame(viaLegacy))
})
