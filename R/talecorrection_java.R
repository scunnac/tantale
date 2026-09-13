

#' Correct TALE ORFs in error-prone sequences
#'
#' @description
#' 
#' As a way faster alternative to run \code{\link[tantale:tell_tales]{tell_tales}}
#' in correction mode on error prone sequences such as ONT assembled genomes,
#' we provide a wrapper around the java binaries from this gitHub
#' \href{https://github.com/Jstacs/Jstacs/tree/master/projects/talecorrect}{page}.
#' It takes an input fasta file and output a file with corrected indels in TALE coding sequences
#' using an approach described in the Erkes et al. \href{https://doi.org/10.1186/s12864-023-09228-1}{paper}.
#' \code{\link[tantale:tell_tales]{tell_tales}} can subsequently be run on the corrected
#' sequences in no correction mode.
#'
#' @param uncorrected_path Path to the input sequence file
#' @param corrected_path Path of the ouput file
#' @param hmm_path Path the folder containning the profile HMM files. The default
#' value points to the ones build from Xanthomonas oryzae pv. oryzae templates.
#' Xox ones are also available in the parent directory. Please see the gitHub
#' \href{https://github.com/Jstacs/Jstacs/tree/master/projects/talecorrect}{page}
#' for instructions on building custom profiles.
#' @param return_corrections Specify \code{TRUE} if you want the list of executed
#' operations on the sequence as a tibble.
#' @param conda_bin Path to your Conda binary file if you need to specify a
#'   path different from the one that is automatically searched by the
#'   reticulate package functions.
#' 
#' 
#' @return A tibble if \code{return_corrections} is \code{TRUE} or the path to the
#' corrected sequences file.
#' 
#' @export
#' @family TALE discovery
correct_tales <- function(uncorrected_path ,
                     corrected_path = file.path(getwd(), "correctedTALEs.fa"),
                     hmm_path = system.file("tools", "talecorrect", "HMMs", "Xoo", package = "tantale", mustWork = T),
                     return_corrections = FALSE,
                     conda_bin = "auto") {
  
  pathToTALECorrection <- system.file("tools", "talecorrect", "TALEcorrection.jar", package = "tantale", mustWork = T)
  outputFolder <- tempfile(pattern = "correct_tales")
  dir.exists(outputFolder) || dir.create(outputFolder, recursive = TRUE)
  domains <- c(N = "N-terminus.10bpRepeat1", C = "repeat", R = "C-terminus")
  if (!fs::file_exists(uncorrected_path)) {
    cli::cli_abort("The provided input file does not exists", class = c("tantale_error"))
  }
  
  #### run nHMMER ####
  nhmmerCmd <- paste(glue::glue("nhmmer {hmm_path}/{domains}.hmm {uncorrected_path} > {outputFolder}/out_nhmmer.{domains}.txt",
                       .sep = "; "), collapse = "; ")
  envReady <- !as.logical(.create_tantale_env(conda_bin = conda_bin))
  if (envReady) {
    cli::cli_inform("Running nHMMER")
    res <- .run_in_conda(env_name = "tantale",
                            conda_bin = conda_bin,
                            command = nhmmerCmd
    )
    if (res) {
      cli::cli_warn("The following nHMMER commands failed:")
      cli::cli_abort("{nhmmerCmd}", class = c("tantale_error"))
    }
  } else {
    stop("Could not create the tantale conda environment on your machine to run nHMMER...")
  }
  
  #### run TALEcorrection ####
  cli::cli_inform("Performing TALEs cds correction on provided sequences.")
  talecorCmd <- glue::glue("java -jar {pathToTALECorrection} correct s={uncorrected_path}",
                  "n={outputFolder}/out_nhmmer.{domains[\"N\"]}.txt r={outputFolder}/out_nhmmer.{domains[\"R\"]}.txt",
                  "c={outputFolder}/out_nhmmer.{domains[\"C\"]}.txt outdir={outputFolder}", .sep = " ")
  res <- try(system(command = talecorCmd, intern = TRUE))
  if (class(res) == "try-error") {
    cli::cli_warn("The following TALEcorrection commands failed:")
    cli::cli_abort("{talecorCmd}", class = c("tantale_error"))
  }
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







