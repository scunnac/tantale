
test_that(".distalPairwiseAlign2 output a tibble with the expected dims", {
  taleParts <- getTaleParts(system.file("extdata", "tellTaleExampleOutput", package = "tantale", mustWork = T)) %>%
    dplyr::mutate(partId = paste(arrayID, position, sep = "_"))
  partAaStringSet <- Biostrings::AAStringSet(taleParts$aaSeq)
  names(partAaStringSet) <- taleParts$partId
  pairAlignScores <- .distalPairwiseAlign2(partAaStringSet, ncores = 1, condaBinPath = "/home/cunnac/bin/miniconda3/condabin/conda")
  expect_true(identical(dim(pairAlignScores), c(9216L,18L)))
})

