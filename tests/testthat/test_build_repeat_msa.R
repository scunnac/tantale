

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

