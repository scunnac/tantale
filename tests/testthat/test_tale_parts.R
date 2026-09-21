
test_that("tale_parts warns if AnnoTALE did not report on a N-Term seq which is present in the rvd seq file from telltale",
          {expect_warning(.tale_parts(test_path("data_for_tests", "err_missing_nterm")))}
)

test_that("tale_parts errors if AnnoTALE did not report on a TALE which is present in the rvd seq file from telltale",
          {expect_warning(.tale_parts(test_path("data_for_tests", "err_array_count")))}
)

test_that("tale_parts warns if AnnoTALE prot and dna files have different number of sequence parts",
          {expect_warning(.tale_parts(test_path("data_for_tests", "err_missing_dna")))}
)

test_that("tale_parts is of expected dims",
          {talParts <- .tale_parts(test_path("data_for_tests", "example_output"))
          expect_true(identical(dim(talParts), c(96L, 9L)))}
)

test_that("tale_parts emits the canonical snake_case column vocabulary",
          {talParts <- .tale_parts(test_path("data_for_tests", "example_output"))
          expect_setequal(names(talParts),
                          c("array_id", "domain_type", "position_in_crd", "dna_seq",
                            "source_directory", "position_in_array", "aa_seq",
                            "rvd", "seqnames"))
          expect_false(any(grepl("[A-Z]", names(talParts))))}
)
