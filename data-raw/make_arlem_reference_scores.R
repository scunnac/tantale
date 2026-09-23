# Reference scores from the ARLEM binary, for tests/testthat/test_arlem_r.R.
#
# tantale used to run a bundled ARLEM 1.0 executable (Linux x86-64 only,
# licence unclear for redistribution; ledger 28). It was replaced by an R
# implementation of the same model (R/arlem.R, ledger 33) and removed from
# the package. This script records what the binary answered, while it still
# ran, so the R implementation stays tested against the binary's own
# output. Re-running it needs that executable, e.g. from the v0.1.9553
# history bundle or inst/tools/arlem/arlem in a checkout before ledger 33.
#
# Run from the package root with ARLEM_BIN pointing at the executable:
#   ARLEM_BIN=inst/tools/arlem/arlem Rscript data-raw/make_arlem_reference_scores.R

devtools::load_all(".", quiet = TRUE)

arlem_bin <- Sys.getenv("ARLEM_BIN")
stopifnot(nzchar(arlem_bin), file.exists(arlem_bin))

write_cost_file <- function(cost, dup, ins) {
  k <- nrow(cost)
  m <- matrix(format(cost), k, k)
  m[lower.tri(m)] <- ""
  diag(m) <- ""
  rows <- apply(m, 1, function(r) trimws(paste(r[r != ""], collapse = " ")))
  path <- tempfile()
  writeLines(c(paste("# Type no.", k),
               paste("# Types", paste(seq_len(k), collapse = " ")),
               paste("# Indel align", ins), paste("# Indel hist", ins),
               paste("# Dup", dup), "# matrix", rows[-k]), path)
  path
}

# The binary's raw pairwise scores, as ARLEM reports them (0-based, id1 < id2)
binary_scores <- function(maps, cost, dup, ins, insert) {
  seqs <- tempfile()
  writeLines(paste0(">a", seq_along(maps), "\n", maps), seqs)
  out <- system2(arlem_bin, c("-f", seqs, "-cfile", write_cost_file(cost, dup, ins),
                              "-align", if (insert) "-insert"),
                 stdout = TRUE, stderr = TRUE)
  stopifnot(is.null(attr(out, "status")))
  hits <- regmatches(out, regexec("Score of aligning Seq:([0-9]+), Seq:([0-9]+) =([0-9]+)", out))
  hits <- do.call(rbind, lapply(hits[lengths(hits) > 0], function(x) as.numeric(x[2:4])))
  stopifnot(nrow(hits) == choose(length(maps), 2))
  tibble::tibble(id1 = hits[, 1], id2 = hits[, 2], arlem_score = hits[, 3])
}

random_case <- function(insert) {
  k <- sample(2:30, 1)
  cost <- ceiling(as.matrix(stats::dist(matrix(stats::runif(k * 3, 0, 60), k))))
  dimnames(cost) <- list(seq_len(k), seq_len(k))
  maps <- replicate(sample(2:5, 1), {
    units <- character()
    while (length(units) < 30) {
      units <- c(units, rep(as.character(sample(k, 1)), sample(c(1, 1, 1, 2, 3, 5), 1)))
    }
    paste(units[seq_len(sample(1:30, 1))], collapse = " ")
  })
  dup <- sample(c(0, 1, 5, 10, 30), 1)
  ins <- sample(c(1, 10, 25, 100, 1000), 1)
  list(maps = maps, cost = cost, dup = dup, ins = ins, insert = insert,
       scores = binary_scores(maps, cost, dup, ins, insert))
}

set.seed(20260923)
cases <- c(lapply(1:25, function(i) random_case(insert = TRUE)),
           lapply(1:10, function(i) random_case(insert = FALSE)))

# Real TALE arrays with tantale's own costs, as tales_tale_distances() ran them
d <- readRDS("tests/testthat/data_for_tests/sampleDistalrOutput.rds")
coded <- .coded_seq_set(tibble::as_tibble(d$tale_parts))
cost <- suppressMessages(.arlem_cost_matrix(domain_distances(d$repeat.similarity)))
real <- list(maps = stats::setNames(as.character(coded), names(coded)),
             cost = cost, dup = 10, ins = 10, insert = TRUE,
             scores = binary_scores(as.character(coded), cost, 10, 10, TRUE))

saveRDS(list(random = cases, real = real),
        "tests/testthat/data_for_tests/arlem_reference_scores.rds",
        compress = "xz")
