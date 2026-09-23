# Retired: .run_arlem() and .arlem_cost_file()
#
# Superseded by .arlem_scores_r() (R/arlem.R), an R implementation of
# ARLEM's alignment model (Abouelhoda, Giegerich, Behzadi & Steyaert, APBC
# 2008; https://doi.org/10.1142/S0219720009004060). These two functions
# wrote ARLEM's cost file and ran the ARLEM 1.0 executable that used to ship
# as inst/tools/arlem/arlem, parsing its stdout. The R implementation gives
# identical scores (checked on ~6000 random pairs and on real TALE arrays;
# tests/testthat/data_for_tests/arlem_reference_scores.rds holds the
# binary's own answers, made by data-raw/make_arlem_reference_scores.R).
#
# Why the executable went: it ran only on Linux x86-64, it printed
# "unauthorized commercial usage and distribution of this program is
# prohibited", and no source or conda build could be found (ledger 28).
# It was deleted from the package rather than moved here, on the
# maintainer's instruction, so this code no longer runs as is: the default
# `arlem` path points at a file that no longer exists. It is kept as the
# record of how the executable was driven (options, cost-file format,
# output parsing).
#
# Moved out of R/distalr.R on 2026-09-23. See dev/restructuring-notes.md
# sections 28 and 33.

#' Write ARLEM's cost file
#'
#' @param mat The matrix from `.arlem_cost_matrix()`.
#' @return The path of the cfile.
#' @noRd
.arlem_cost_file <- function(mat) {
  header <- c(glue::glue("# Type no. ", nrow(mat)),
              glue::glue("# Types ", paste(seq_len(nrow(mat)), collapse = " ")),
              glue::glue("# Indel align ", .arlem_indel_cost),
              glue::glue("# Indel hist ", .arlem_indel_cost),
              glue::glue("# Dup ", .arlem_dup_cost), "# matrix")

  mat <- matrix(format(mat), ncol = ncol(mat),
                dimnames = list(rownames(mat), colnames(mat)))
  mat[lower.tri(mat)] <- ""
  diag(mat) <- ""
  lines <- apply(mat, 1, function(row) gsub("^[ ]+", "", paste(row, collapse = " ")))
  lines <- lines[seq_len(length(lines) - 1L)]   # drop the empty last row
  
  cfile <- tempfile()
  writeLines(cfile, text = c(header, lines))
  cfile
}

#' Run ARLEM over the coded strings and parse its stdout
#' @noRd
.run_arlem <- function(coded, cfile,
                       arlem = system.file("tools", "arlem", "arlem",
                                           package = "tantale", mustWork = TRUE)) {
  seqfile <- tempfile(fileext = ".fasta")
  outfile <- tempfile(fileext = ".txt")
  errfile <- tempfile(fileext = ".txt")
  Biostrings::writeXStringSet(coded, filepath = seqfile, format = "fasta",
                              width = 20000L)
  cmd <- glue::glue("{shQuote(arlem)} -f {shQuote(seqfile)} -cfile {shQuote(cfile)} -align -insert -showalign > {shQuote(outfile)}")
  cli::cli_inform("Running ARLEM version 1.0 : ")
  cli::cli_inform("Copyright by Mohamed I. Abouelhoda")
  cli::cli_inform("Plz. cite Abouelhoda, Giegerich, Behzadi, and Steyaert")
  .tantale_exec(cmd, stderr_file = errfile, what = "ARLEM")
  raw <- readLines(outfile, warn = FALSE)

  scores <- grep("Score of aligning Seq:", raw, value = TRUE)
  # A binary that cannot run here (another architecture, a lost execute bit)
  # can still exit 0 from the shell's point of view on some systems; the
  # score count is the check that does not depend on that.
  expected <- choose(length(coded), 2)
  if (length(scores) != expected) {
    cli::cli_abort(
      c("ARLEM did not report a score for every pair of arrays.",
        "x" = "Expected {expected} score{?s}, got {length(scores)}.",
        "i" = "Binary: {.file {arlem}}, platform {.val {R.version$platform}}.",
        "i" = "The bundled ARLEM is a Linux x86-64 executable."),
      class = c("tantale_error_arlem_failed", "tantale_error"))
  }
  scores <- gsub("Score of aligning Seq:([0-9]+), Seq:([0-9]+) =([0-9]+)",
                 "\\1|\\2|\\3", scores)
  scores <- strsplit(scores, split = "\\|")
  scores <- tibble::as_tibble(
    do.call(rbind, lapply(scores, function(s) t(as.matrix(as.numeric(s))))),
    .name_repair = "minimal")
  colnames(scores) <- c("id1", "id2", "arlem_score")
  .arlem_scores_long(scores, length(coded))
}
