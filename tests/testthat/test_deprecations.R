# The superseded functions still work, but warn and point at their replacement.
# Their behaviour is tested against the internal implementations (.tale_parts(),
# .split_list(), .build_repeat_msa()) in the respective test files, so that those
# assertions cannot be satisfied by the deprecation warning alone.

test_that("tale_parts() is deprecated in favour of tales_from_telltale()", {
  w <- NULL
  suppressWarnings(withCallingHandlers(
    tale_parts(test_path("data_for_tests", "tellTaleExampleOutput")),
    deprecatedWarning = function(c) w <<- conditionMessage(c)
  ))
  expect_match(w, "tales_from_telltale")
})

test_that("split_list() is deprecated in favour of as_tales()", {
  w <- NULL
  suppressWarnings(withCallingHandlers(
    split_list(test_path("data_for_tests", "tellTaleExampleOutput", "rvdSequences.fas"),
               sep = "-"),
    deprecatedWarning = function(c) w <<- conditionMessage(c)
  ))
  expect_match(w, "as_tales")
})

test_that("build_repeat_msa() is deprecated in favour of tales_align()", {
  w <- NULL
  suppressWarnings(withCallingHandlers(
    build_repeat_msa(
      input_seqs = system.file("extdata", "small_Out_CodedRepeats.fa",
                               package = "tantale", mustWork = TRUE),
      sep = " "
    ),
    deprecatedWarning = function(c) w <<- conditionMessage(c)
  ))
  expect_match(w, "tales_align")
})

test_that("the deprecated wrappers still return what they always did", {
  parts <- suppressWarnings(
    tale_parts(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
  expect_identical(dim(parts), c(96L, 9L))

  seqs <- suppressWarnings(
    split_list(test_path("data_for_tests", "tellTaleExampleOutput", "rvdSequences.fas"),
               sep = "-")
  )
  expect_type(seqs, "list")
  expect_setequal(lengths(seqs), c(28L, 16L, 28L, 24L))
})
