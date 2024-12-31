
taleParts <- getTaleParts(test_path("data_for_tests", "tellTaleExampleOutput")) %>%
  dplyr::mutate(partId = paste(arrayID, positionInArray, sep = "_"))
partAaStringSet <-  Biostrings::AAStringSet(taleParts$aaSeq)
names(partAaStringSet) <- taleParts$partId

test_that(".distalPairwiseAlign output a tibble with the expected dims", {
  pairAlignScores <- .distalPairwiseAlign(partAaStringSet, ncores = 4)
  expect_true(identical(dim(pairAlignScores), c(9216L,5L)))
})

test_that(".distalPairwiseAlign2 output a tibble with the expected dims", {
  pairAlignScores <- .distalPairwiseAlign2(partAaStringSet, condaBinPath = "auto")
  expect_true(identical(dim(pairAlignScores), c(9216L,19L)))
})

test_that(".distalPairwiseAlign3 output a tibble with the expected dims", {
  pairAlignScores <- .distalPairwiseAlign3(partAaStringSet)
  expect_true(identical(dim(pairAlignScores), c(9216L,3L)))
})





# pairAlignScores$raw %>%  hist(breaks = 100)
# pairAlignScores$Dissim %>%  hist(breaks = 100)
# pairAlignScores %>% dplyr::filter(Dissim < 1000, Dissim > 30)
# pairAlignScores %>% dplyr::filter(Dissim < 10, Dissim >= 0)
# pairAlignScores %>% dplyr::filter(Dissim < 20, Dissim > 10)
# 
# pairAlignScores %<>% dplyr::mutate(Sim = 100/(1+exp(-1*-0.9*(Dissim-3))))
# pairAlignScores$Sim %>%  hist(breaks = 100)
# skimr::skim(pairAlignScores$Dissim)
# ggplot2::ggplot(pairAlignScores, mapping = ggplot2::aes(x= Dissim, y= Sim)) +
#   ggplot2::geom_point(alpha = 0.1)



