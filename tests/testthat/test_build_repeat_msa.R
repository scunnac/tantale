

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
  expect_named(out, c("id1", "id2", "sim"))
  expect_equal(out$sim, c(100, 60))

  legacy <- data.frame(RepU1 = c("a","b"), RepU2 = c("b","a"), Sim = c(100, 60))
  expect_equal(tantale:::.as_mafft_score_table(legacy)$sim, c(100, 60))
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


#### RVD alignments get a scoring matrix ####

test_that(".rvd_score_table() covers every pair and fills XX neutrally", {
  t <- tantale:::.rvd_score_table(c("NI", "NN", "HD", "XX"))
  expect_identical(nrow(t), 16L)
  expect_false(anyNA(t$sim))
  g <- function(a, b) t$sim[t$id1 == a & t$id2 == b]
  # XX means "terminus detected, identity unknown": neutral against everything
  expect_identical(g("XX", "NI"), 0)
  expect_identical(g("XX", "HD"), 0)
  # ...but maximal against itself. Not a claim about knowledge: MAFFT produces
  # unusable output when the diagonal is not high, and a low value would assert
  # that an XX must *not* align with an XX.
  expect_identical(g("XX", "XX"), 1)
  expect_identical(g("NI", "NI"), 1)
  # and real RVD pairs keep their correlation
  expect_equal(g("NI", "NN"), 0.63, tolerance = 0.01)
})

test_that("RVD alignment opts in to the built-in matrix with repeat_sims = \"rvd\"", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:4], ]
  withMat <- suppressWarnings(suppressMessages(
    tales_align(sub, residue_col = "rvd", repeat_sims = "rvd")))
  without <- suppressWarnings(suppressMessages(tales_align(sub, residue_col = "rvd")))
  # the matrix must actually change the alignment, or it is not being used
  expect_false(identical(as.data.frame(withMat), as.data.frame(without)))
})

test_that("a repeat-keyed table is refused for an RVD alignment", {
  # it is keyed by dom_code, which has no meaning for RVDs; silently ignoring
  # it was the old behaviour and hid the mistake
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:4], ]
  expect_error(
    suppressMessages(tales_align(sub, residue_col = "rvd",
                                 repeat_sims = domain_distances(d$repeat.similarity))),
    class = "tantale_error_msa_sim_table")
})

test_that("NULL and FALSE both mean no matrix", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:4], ]
  a <- suppressWarnings(suppressMessages(tales_align(sub, residue_col = "dom_code")))
  b <- suppressWarnings(suppressMessages(
    tales_align(sub, residue_col = "dom_code", repeat_sims = FALSE)))
  expect_equal(as.data.frame(a), as.data.frame(b))
})


test_that("repeat_sims = \"rvd\" is refused for a repeat-code alignment", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:4], ]
  expect_error(
    suppressMessages(tales_align(sub, residue_col = "dom_code", repeat_sims = "rvd")),
    class = "tantale_error_msa_sim_table")
})

test_that("the default RVD alignment is unchanged by this feature", {
  # opt-in, so an existing call must give exactly what it gave before
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:4], ]
  a <- suppressWarnings(suppressMessages(tales_align(sub, residue_col = "rvd")))
  b <- suppressWarnings(suppressMessages(
    tales_align(sub, residue_col = "rvd", repeat_sims = FALSE)))
  expect_equal(as.data.frame(a), as.data.frame(b))
})
