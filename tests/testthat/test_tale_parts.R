
# The err_* fixtures are copies of example_output with one inconsistency
# planted in ROI_00004 (data-raw/make_telltale_test_fixtures.R).

test_that("an array without an N-terminus part is kept, renumbered from its first repeat", {
  expect_warning(talParts <- .tale_parts(test_path("data_for_tests", "err_missing_nterm")),
                 class = "tantale_warning_terminus_absent")
  roi4 <- talParts[talParts$array_id == "ROI_00004", ]
  expect_false("N-terminus" %in% roi4$domain_type)
  expect_identical(roi4$domain_type[1], "repeat")
  expect_identical(sort(roi4$position_in_array), seq_len(nrow(roi4)))
  expect_identical(roi4$domain_type[which.max(roi4$position_in_array)], "C-terminus")
})

test_that("RVDs that do not correspond to the repeat parts are an error", {
  expect_error(.tale_parts(test_path("data_for_tests", "err_array_count")),
               class = "tantale_error_parts_inconsistent")
})

test_that("an array whose protein and DNA parts disagree is left out", {
  expect_warning(talParts <- .tale_parts(test_path("data_for_tests", "err_missing_dna")),
                 class = "tantale_warning_parts_inconsistent")
  expect_false("ROI_00004" %in% talParts$array_id)
  expect_setequal(unique(talParts$array_id), c("ROI_00001", "ROI_00002", "ROI_00003"))
})

# a writable copy of example_output with its array report edited by `edit`
edited_output <- function(edit) {
  dir <- file.path(tempfile(), "example_output")
  dir.create(dirname(dir))
  file.copy(test_path("data_for_tests", "example_output"), dirname(dir), recursive = TRUE)
  report <- readr::read_tsv(file.path(dir, "array_report.tsv"), show_col_types = FALSE)
  readr::write_tsv(edit(report), file.path(dir, "array_report.tsv"))
  dir
}

test_that("terminus codes follow the protein-profile columns of the array report", {
  talParts <- .tale_parts(test_path("data_for_tests", "example_output"))
  expect_setequal(talParts$rvd[talParts$domain_type == "N-terminus"], "NTERM")
  expect_setequal(talParts$rvd[talParts$domain_type == "C-terminus"], "CTERM")

  dir <- edited_output(function(r) {
    r$cterm_aa_hit[r$array_id == "ROI_00002"] <- FALSE
    r
  })
  talParts <- .tale_parts(dir)
  cterm <- talParts[talParts$domain_type == "C-terminus", ]
  expect_identical(cterm$rvd[cterm$array_id == "ROI_00002"], "XXXXX")
  expect_setequal(cterm$rvd[cterm$array_id != "ROI_00002"], "CTERM")
})

test_that("a directory from an older tell_tales() is an error", {
  dir <- edited_output(function(r) dplyr::select(r, -nterm_aa_hit, -cterm_aa_hit))
  expect_error(.tale_parts(dir), class = "tantale_error_telltale_outdated")
})

test_that("repeat RVDs are AnnoTALE's", {
  talParts <- .tale_parts(test_path("data_for_tests", "example_output"))
  roi1 <- talParts[talParts$array_id == "ROI_00001" & talParts$domain_type == "repeat", ]
  annotale <- .split_list(test_path("data_for_tests", "example_output", "annotale",
                                    "ROI_00001", "TALE_RVDs.fasta"))[[1]]
  expect_identical(roi1$rvd[order(roi1$position_in_array)], annotale)
})

test_that("tale_parts is of expected dims",
          {talParts <- .tale_parts(test_path("data_for_tests", "example_output"))
          expect_identical(dim(talParts), c(96L, 9L))}
)

test_that("tale_parts emits the canonical snake_case column vocabulary",
          {talParts <- .tale_parts(test_path("data_for_tests", "example_output"))
          expect_setequal(names(talParts),
                          c("array_id", "domain_type", "position_in_crd", "dna_seq",
                            "source_directory", "position_in_array", "aa_seq",
                            "rvd", "seqnames"))
          expect_false(any(grepl("[A-Z]", names(talParts))))}
)
