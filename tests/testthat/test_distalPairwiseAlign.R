
tale_parts <- tale_parts(test_path("data_for_tests", "tellTaleExampleOutput")) %>%
  dplyr::mutate(partId = paste(arrayID, positionInArray, sep = "_"))
part_aa_set <-  Biostrings::AAStringSet(tale_parts$aaSeq)
names(part_aa_set) <- tale_parts$partId

test_that(".pairwise_align_biostrings output a tibble with the expected dims", {
  pair_align_scores <- .pairwise_align_biostrings(part_aa_set, ncores = 4)
  expect_true(identical(dim(pair_align_scores), c(9216L,5L)))
})

test_that(".pairwise_align_mmseq2 output a tibble with the expected dims", {
  pair_align_scores <- .pairwise_align_mmseq2(part_aa_set, conda_bin = "auto")
  expect_true(identical(dim(pair_align_scores), c(9216L,19L)))
})

test_that(".pairwise_align_decipher output a tibble with the expected dims", {
  pair_align_scores <- .pairwise_align_decipher(part_aa_set)
  expect_true(identical(dim(pair_align_scores), c(9216L,3L)))
})





# pair_align_scores$raw %>%  hist(breaks = 100)
# pair_align_scores$Dissim %>%  hist(breaks = 100)
# pair_align_scores %>% dplyr::filter(Dissim < 1000, Dissim > 30)
# pair_align_scores %>% dplyr::filter(Dissim < 10, Dissim >= 0)
# pair_align_scores %>% dplyr::filter(Dissim < 20, Dissim > 10)
# 
# pair_align_scores %<>% dplyr::mutate(Sim = 100/(1+exp(-1*-0.9*(Dissim-3))))
# pair_align_scores$Sim %>%  hist(breaks = 100)
# skimr::skim(pair_align_scores$Dissim)
# ggplot2::ggplot(pair_align_scores, mapping = ggplot2::aes(x= Dissim, y= Sim)) +
#   ggplot2::geom_point(alpha = 0.1)



