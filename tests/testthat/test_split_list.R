

test_that("toListOfSplitedStr output warnings when last element in 'vectorized' sequence is empty", {
  expect_warning(toListOfSplitedStr(test_path("data_for_tests", "Out_CodedRepeats.fa"), sep = " "))
})

test_that("toListOfSplitedStr output list of expected shape with with a fasta file as input and sep ' '", {
  suppressWarnings(l <- toListOfSplitedStr(test_path("data_for_tests", "Out_CodedRepeats.fa"), sep = " "))
  expect_setequal(lapply(l, length), c(24, 28, 16, 28, 20, 24, 21, 18, 19, 24, 28, 16, 28,
                     20, 23, 21, 14, 19, 20, 24, 21, 18, 24, 28, 16, 28, 19))
})

test_that("toListOfSplitedStr output list of expected shape with a 'BStringSet' object as input and sep '-'", {
  l <- toListOfSplitedStr(
    Biostrings::readBStringSet(test_path("data_for_tests", "tellTaleExampleOutput", "rvdSequences.fas")),
    sep = "-")
  expect_setequal(lapply(l, length), c(28, 16, 28, 24))
})

test_that("toListOfSplitedStr output list of expected shape with a fasta file as input and sep '-'", {
  l <- toListOfSplitedStr(
    test_path("data_for_tests", "tellTaleExampleOutput", "rvdSequences.fas"),
    sep = "-")
  expect_setequal(lapply(l, length), c(28, 16, 28, 24))
})

test_that("toListOfSplitedStr output list of expected shape with list of strings as input and sep '-'", {
  l <- toListOfSplitedStr(
    Biostrings::readBStringSet(test_path("data_for_tests", "tellTaleExampleOutput", "rvdSequences.fas")) %>%
      as.character() %>% as.list(),
    sep = "-")
  expect_setequal(lapply(l, length), c(28, 16, 28, 24))
})

test_that("toListOfSplitedStr output an error if atomicStrings is not of expected type", {
  l <- toListOfSplitedStr(
    Biostrings::readBStringSet(test_path("data_for_tests", "tellTaleExampleOutput", "rvdSequences.fas")) %>%
      as.character() %>% as.list(),  sep = "-")
  expect_error(toListOfSplitedStr(l))
})



