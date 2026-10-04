# Errors raised by the package carry classed conditions, so they can be
# asserted by class rather than by matching message text. See ledger 9.5:
# before the cli conversion many of these were bare stop() calls whose
# condition message was the empty string.

fixture_tales <- function() {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  tales_quietly(out$tale_parts)
}

test_that("every error the package raises carries the tantale_error parent", {
  # a representative sample across the different subsystems
  expect_error(tales(data.frame(nope = 1)), class = "tantale_error")
  expect_error(pairwise_distances(data.frame(id1 = "a")), class = "tantale_error")
  expect_error(tales_compare_distal(data.frame(a = 1)), class = "tantale_error")
})


#### tales_compare_distal() guards ####

test_that("an unknown aln_method is rejected by name", {
  x <- fixture_tales()
  # This guard was unreachable before: logger::log_errors() && stop(...) short
  # circuited, so the user saw a complaint about calling handlers instead.
  expect_error(tales_compare_distal(x, aln_method = "not_a_method"),
               class = "tantale_error_aln_method")
  expect_error(tales_compare_distal(x, aln_method = "not_a_method"), "not_a_method")
})

test_that("tales_compare_distal() refuses a tales with neither aa_seq nor dna_seq", {
  # dropping aa_seq alone is no longer an error: dna_seq is translated instead
  x <- fixture_tales()
  x$aa_seq <- NULL
  expect_warning(suppressMessages(tales_compare_distal(x)),
                 class = "tantale_warning_translated_aa")
  x$dna_seq <- NULL
  expect_error(tales_compare_distal(x), class = "tantale_error_compare_no_aa")
})


#### the RVD conversion internals repaired in ledger 2 ####

rvd_align_fixture <- function() {
  matrix(c("NI", "HD", NA, "NI", NA, "NG"), nrow = 2, byrow = TRUE,
         dimnames = list(c("a", "b"), NULL))
}

test_that(".rvd_to_match_align() runs against the internal rvdSimDf", {
  # Its default argument used to be tantale::rvdSimDf, which errors because
  # rvdSimDf lives in sysdata.rda -- the function could never run (ledger 2).
  m <- rvd_align_fixture()
  m[is.na(m)] <- "NG"
  out <- tantale:::.rvd_to_match_align(m)
  expect_true(is.matrix(out))
  expect_equal(dim(out), dim(m))
  expect_type(out, "double")
})


#### Empty-input guards ####

test_that("projecting a tales with nothing left to render errors", {
  x <- fixture_tales()
  # keep only the terminus sentinels, then ask for repeats only
  anchors <- x[x$rvd %in% tales_anchor_codes(), ]
  expect_error(tales_rvd_strings(anchors, repeats_only = TRUE),
               class = "tantale_error_projection_empty")
})

test_that("aligning an empty tales errors", {
  x <- fixture_tales()
  expect_error(tales_align(x[0, ], residue_col = "rvd"),
               class = "tantale_error_msa_empty")
})

test_that("a part with no amino acid sequence is reported with its array", {
  x <- fixture_tales()
  bad <- unique(x$array_id)[1]
  x$aa_seq[x$array_id == bad][1] <- NA_character_
  expect_error(suppressWarnings(tales_compare_distal(x)),
               class = "tantale_error_parts_no_aa")
  # the offending array must be named, not just counted
  expect_error(suppressWarnings(tales_compare_distal(x)), bad, fixed = TRUE)
})

test_that("an empty-string aa_seq is reported too, not just NA", {
  # The guard covers is.na() | == "", but the list of affected arrays used to
  # collect only the NA ones, so an empty string produced a message naming
  # nothing at all.
  x <- fixture_tales()
  bad <- unique(x$array_id)[1]
  x$aa_seq[x$array_id == bad][1] <- ""
  expect_error(suppressWarnings(tales_compare_distal(x)), bad, fixed = TRUE)
})


#### Specific classes on the last generic-only errors (ledger §9.5) ####

test_that(".split_list() rejects input it cannot split", {
  expect_error(.split_list(list(c("NI", "HD"))),
               class = "tantale_error_bad_argument")
  expect_error(.split_list(1:3), class = "tantale_error_bad_argument")
})

test_that("correct_tales() asks for corrected_path before running anything", {
  subj <- system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                      package = "tantale", mustWork = TRUE)
  expect_error(correct_tales(subj), class = "tantale_error_missing_output")
})

test_that("correct_tales() names a missing input file", {
  missing <- file.path(tempdir(), "no_such_file.fa")
  expect_error(correct_tales(missing), class = "tantale_error_missing_file")
  expect_error(correct_tales(missing), "no_such_file.fa", fixed = TRUE)
})

test_that(".rvds_from_annotale_file() rejects a file AnnoTALE did not name", {
  expect_error(.rvds_from_annotale_file("rvds.fasta"),
               class = "tantale_error_annotale_file")
})
