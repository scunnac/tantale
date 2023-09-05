

test_that("buildRepeatMsa result with RVDs is of expected dims",
          {aln <- buildRepeatMsa(inputSeqs = system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",package = "tantale", mustWork = TRUE),
                                 sep = "-", distalRepeatSims = NULL,
                                 mafftOpts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
                                 mafftPath = system.file("tools", "mafft-linux64",package = "tantale", mustWork = TRUE),
                                 gapSymbol = "-")
          expect_true(identical(dim(aln), c(11L, 26L)))}
)


test_that("buildRepeatMsa result with coded parts is of expected dims",
          {aln <- buildRepeatMsa(inputSeqs = system.file("extdata", "small_Out_CodedRepeats.fa",package = "tantale", mustWork = TRUE),
                                 sep = " ", distalRepeatSims = NULL,
                                 mafftOpts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
                                 mafftPath = system.file("tools", "mafft-linux64",package = "tantale", mustWork = TRUE),
                                 gapSymbol = "-")
          expect_true(identical(dim(aln), c(3L, 28L)))}
)


test_that("buildRepeatMsa throughts a warning and return a NA matrix if provided with empty sequences",
          {seqs <- Biostrings::BStringSet()
          expect_warning(aln <- buildRepeatMsa(inputSeqs = seqs,
                                               sep = " ", distalRepeatSims = NULL,
                                               mafftOpts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
                                               mafftPath = system.file("tools", "mafft-linux64",package = "tantale", mustWork = TRUE),
                                               gapSymbol = "-"))
          expect_equal(aln, matrix())}
)

test_that("buildRepeatMsa deals properly with single sequence inputs",
          {seqs <- Biostrings::readBStringSet(system.file("extdata", "small_Out_CodedRepeats.fa",package = "tantale", mustWork = TRUE))
          aln <- buildRepeatMsa(inputSeqs = seqs[1],
                                sep = " ", distalRepeatSims = NULL,
                                mafftOpts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
                                mafftPath = system.file("tools", "mafft-linux64",package = "tantale", mustWork = TRUE),
                                gapSymbol = "-")
          expect_true(identical(dim(aln), c(1L, 28L)))}
)

