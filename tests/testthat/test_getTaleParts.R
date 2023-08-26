
test_that("getTaleParts errors if AnnoTALE did not report on a N-Term seq which is present in the rvd seq file from telltale",
          {expect_error(getTaleParts(test_path("data_for_tests", "tellTaleErrorOutputMissingN-term")))}
          )

test_that("getTaleParts errors if AnnoTALE did not report on a TALE which is present in the rvd seq file from telltale",
          {expect_error(getTaleParts(test_path("data_for_tests", "tellTaleErrorOutputInconsistentArrayNumber")))}
          )

test_that("getTaleParts errors if AnnoTALE prot and dna files have different number of components",
          {expect_error(getTaleParts(test_path("data_for_tests", "tellTaleErrorMissingAnnotaleDnaDomain")))}
)
