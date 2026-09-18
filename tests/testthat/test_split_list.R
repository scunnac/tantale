

test_that("split_list output warnings when last element in 'vectorized' sequence is empty", {
  expect_warning(.split_list(test_path("data_for_tests", "Out_CodedRepeats.fa"), sep = " "))
})

test_that("split_list output list of expected shape with with a fasta file as input and sep ' '", {
  suppressWarnings(l <- .split_list(test_path("data_for_tests", "Out_CodedRepeats.fa"), sep = " "))
  expect_setequal(lapply(l, length), c(24, 28, 16, 28, 20, 24, 21, 18, 19, 24, 28, 16, 28,
                     20, 23, 21, 14, 19, 20, 24, 21, 18, 24, 28, 16, 28, 19))
})

test_that("split_list output list of expected shape with a 'BStringSet' object as input and sep '-'", {
  l <- .split_list(
    Biostrings::readBStringSet(test_path("data_for_tests", "tellTaleExampleOutput", "rvd_sequences.fas")),
    sep = "-")
  expect_setequal(lapply(l, length), c(28, 16, 28, 24))
})

test_that("split_list output list of expected shape with a fasta file as input and sep '-'", {
  l <- .split_list(
    test_path("data_for_tests", "tellTaleExampleOutput", "rvd_sequences.fas"),
    sep = "-")
  expect_setequal(lapply(l, length), c(28, 16, 28, 24))
})

test_that("split_list output list of expected shape with list of strings as input and sep '-'", {
  l <- .split_list(
    Biostrings::readBStringSet(test_path("data_for_tests", "tellTaleExampleOutput", "rvd_sequences.fas")) %>%
      as.character() %>% as.list(),
    sep = "-")
  expect_setequal(lapply(l, length), c(28, 16, 28, 24))
})

test_that("split_list output an error if strings is not of expected type", {
  l <- .split_list(
    Biostrings::readBStringSet(test_path("data_for_tests", "tellTaleExampleOutput", "rvd_sequences.fas")) %>%
      as.character() %>% as.list(),  sep = "-")
  expect_error(.split_list(l))
})



