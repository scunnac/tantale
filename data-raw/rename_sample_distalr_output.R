# One-off conversion of tests/testthat/data_for_tests/sampleDistalrOutput.rds
# to the current column names, when tales() and pairwise_distances() stopped
# renaming legacy spellings on the way in (ledger §37). The fixture itself
# predates data-raw/ and has no generator; this records what was changed.
#
# Run with:  Rscript data-raw/rename_sample_distalr_output.R

path <- "tests/testthat/data_for_tests/sampleDistalrOutput.rds"
d <- readRDS(path)
renames <- c(TAL1 = "id1", TAL2 = "id2", RepU1 = "id1", RepU2 = "id2",
             Sim = "sim", Dissim = "dissim", arlemScore = "arlem_score",
             maxLength = "max_length", normArlemScore = "norm_arlem_score")
for (el in c("repeat.similarity", "tal.similarity")) {
  hit <- names(d[[el]]) %in% names(renames)
  names(d[[el]])[hit] <- renames[names(d[[el]])[hit]]
}
saveRDS(d, path)
