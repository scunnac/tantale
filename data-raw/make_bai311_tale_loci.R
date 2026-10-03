# Builds tests/testthat/data_for_tests/bai311_tale_loci.fa: the TALE loci of
# the BAI3-1-1 example genome, each with 3 kb on either side, so that
# test_correct_tales.R can run TALEcorrection without the whole genome, which
# no longer ships with the package (ledger §50). Loci from tell_tales() on the
# uncorrected genome (cterm_min_score = 300); overlapping windows merged.
# Run from the package root after tantale_setup(install = TRUE).

genome <- Biostrings::readDNAStringSet(tantale::tantale_genome("BAI3-1-1"))
loci <- IRanges::IRanges(
  start = c(206203, 266994, 1736263, 2202289, 2206727, 2209940, 4182275, 4256992),
  end   = c(209612, 269893, 1739877, 2206618, 2209831, 2214255, 4186201, 4260502)
)
windows <- IRanges::reduce(IRanges::resize(loci, IRanges::width(loci) + 6000L,
                                           fix = "center"))
excerpt <- Biostrings::subseq(rep(genome["contig_1"], length(windows)),
                              start = IRanges::start(windows),
                              end = IRanges::end(windows))
names(excerpt) <- sprintf("contig_1_%d_%d", IRanges::start(windows), IRanges::end(windows))
Biostrings::writeXStringSet(excerpt, "tests/testthat/data_for_tests/bai311_tale_loci.fa")
