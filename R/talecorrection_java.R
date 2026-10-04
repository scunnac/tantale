

#' Correct TALE ORFs in error-prone sequences
#'
#' @description
#' 
#' A much faster alternative to running \code{\link[tantale:tell_tales]{tell_tales}}
#' in correction mode on error-prone sequences such as ONT-assembled genomes:
#' a wrapper around the Java TALE correction tool from this GitHub
#' \href{https://github.com/Jstacs/Jstacs/tree/master/projects/talecorrect}{page}.
#' It takes an input fasta file and writes a file with corrected indels in TALE coding sequences
#' using an approach described in the Erkes et al. \href{https://doi.org/10.1186/s12864-023-09228-1}{paper}.
#' \code{\link[tantale:tell_tales]{tell_tales}} can subsequently be run on the corrected
#' sequences in no correction mode.
#'
#' @param uncorrected_path Path to the input sequence file
#' @param corrected_path Path of the output file. Required: tantale does not
#'   write to the working directory unasked.
#' @param hmm_path Path to the folder containing the profile HMM files. The
#' default points to the ones built from Xanthomonas oryzae pv. oryzae (Xoo)
#' templates, in the tools [tantale_setup()] downloads; Xoc (X. oryzae pv.
#' oryzicola) ones are in the sibling folder \code{Xoc}.
#' Please see the GitHub
#' \href{https://github.com/Jstacs/Jstacs/tree/master/projects/talecorrect}{page}
#' for instructions on building custom profiles.
#' @param return_corrections Specify \code{TRUE} if you want the list of executed
#' operations on the sequence as a tibble.
#' @param conda_bin Path to your Conda binary file if you need to specify a
#'   path different from the one that is automatically searched by the
#'   reticulate package functions.
#' @param java_args A single string of options for the Java virtual
#'   machine, placed before \code{-jar}, such as \code{"-Xmx8G"} to raise
#'   its memory limit. The default, \code{""}, leaves Java's own defaults.
#'   TALEcorrection itself takes no options beyond its inputs.
#' 
#' @return A tibble if \code{return_corrections} is \code{TRUE} or the path to the
#' corrected sequences file.
#' 
#' @export
#' @family TALE discovery
#' @examples
#' \donttest{
#' # Needs nhmmer and a Java runtime.
#' subj <- system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
#'                     package = "tantale")
#' out_fa <- tempfile(fileext = ".fa")
#' correct_tales(uncorrected_path = subj, corrected_path = out_fa)
#' }
correct_tales <- function(uncorrected_path ,
                     corrected_path,
                     hmm_path = .tantale_tool("talecorrection_hmm"),
                     return_corrections = FALSE,
                     conda_bin = "auto",
                     java_args = "") {
  .check_jar_args(NULL, java_args, fn = "correct_tales")
  pathToTALECorrection <- .tantale_tool("talecorrection")
  outputFolder <- tempfile(pattern = "correct_tales")
  dir.exists(outputFolder) || dir.create(outputFolder, recursive = TRUE)
  # Keys match TALEcorrection.jar's own flags (n=/c=/r=), not domain
  # initials: its help text names c= "C-terminus nHMMER File" and r=
  # "Repeats nHMMER File" -- confirmed by running it, not assumed from the
  # flag letters, which read the other way round at a glance.
  domains <- c(N = "N-terminus.10bpRepeat1", C = "C-terminus", R = "repeat")
  if (!fs::file_exists(uncorrected_path)) {
    cli::cli_abort("No sequence file at {.file {uncorrected_path}}.",
                   class = c("tantale_error_missing_file", "tantale_error"))
  }
  # Checked here: the path is otherwise first read after the searches and
  # the correction have run.
  if (missing(corrected_path)) {
    cli::cli_abort(
      c("{.arg corrected_path} is required.",
        "i" = "Give the file to write the corrected sequences to."),
      class = c("tantale_error_missing_output", "tantale_error"))
  }
  
  #### run nHMMER ####
  # Resolved to an absolute path, which matters here more than elsewhere:
  # these are three commands in one string, and under `conda run` only the
  # first of a compound string ran inside the environment (the rest picked
  # up /usr/bin/nhmmer, HMMER 3.4 against the pinned 3.3.2). Joined with
  # "&&" so that the exit status checked is that of any failing search.
  nhmmer <- shQuote(.tantale_bin("nhmmer", conda_bin = conda_bin))
  nhmmerCmd <- paste(glue::glue(
    "{nhmmer} {shQuote(file.path(hmm_path, paste0(domains, '.hmm')))} {shQuote(uncorrected_path)}",
    " > {shQuote(file.path(outputFolder, paste0('out_nhmmer.', domains, '.txt')))}",
    .sep = ""), collapse = " && ")
  envReady <- !as.logical(.create_tantale_env(conda_bin = conda_bin))
  if (envReady) {
    cli::cli_inform("Running nHMMER")
    .tantale_exec(nhmmerCmd, what = "nHMMER")
  } else {
    .abort_no_env("nHMMER")
  }
  
  #### run TALEcorrection ####
  cli::cli_inform("Performing TALEs cds correction on provided sequences.")
  hmmerOut <- function(d) shQuote(file.path(outputFolder, paste0("out_nhmmer.", d, ".txt")))
  talecorCmd <- glue::glue("java {java_args} -jar {shQuote(pathToTALECorrection)} correct s={shQuote(uncorrected_path)}",
                  "n={hmmerOut(domains[\"N\"])} r={hmmerOut(domains[\"R\"])}",
                  "c={hmmerOut(domains[\"C\"])} outdir={shQuote(outputFolder)}", .sep = " ")
  # Its standard output was never shown; its standard error still is.
  .tantale_exec(paste(talecorCmd, "> /dev/null"), what = "TALEcorrection")
  correctionsTble <- readr::read_tsv(file = file.path(outputFolder, "substitionList.tsv"),
                                     show_col_types = FALSE) %>%
    dplyr::rename(posInOriginSeq = `position in uncorrected sequences`)
  
  file.copy(file.path(outputFolder, "correctedTALEs.fa"), file.path(corrected_path),
            overwrite = FALSE)
  unlink(outputFolder, recursive = TRUE) 
  
  if (return_corrections) {
    return(correctionsTble)
  } else {
    return(file.path(corrected_path))
  }
}







