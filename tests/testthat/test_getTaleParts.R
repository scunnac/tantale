
test_that("getTaleParts warns if AnnoTALE did not report on a N-Term seq which is present in the rvd seq file from telltale",
          {expect_warning(getTaleParts(test_path("data_for_tests", "tellTaleErrorOutputMissingN-term")))}
)

test_that("getTaleParts errors if AnnoTALE did not report on a TALE which is present in the rvd seq file from telltale",
          {expect_warning(getTaleParts(test_path("data_for_tests", "tellTaleErrorOutputInconsistentArrayNumber")))}
)

test_that("getTaleParts warns if AnnoTALE prot and dna files have different number of sequence parts",
          {expect_warning(getTaleParts(test_path("data_for_tests", "tellTaleErrorMissingAnnotaleDnaDomain")))}
)

test_that("getTaleParts is of expected dims",
          {talParts <- getTaleParts(test_path("data_for_tests", "tellTaleExampleOutput"))
          expect_true(identical(dim(talParts), c(96L, 9L)))}
)
