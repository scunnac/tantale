

#' Correct TALE ORFs in error-prone sequences
#'
#' @description
#' 
#' As a way faster alternative to run \code{\link[tantale:tellTale]{tellTale}}
#' in correction mode on error prone sequences such as ONT assembled genomes,
#' we provide a wrapper around the java binaries from this gitHub
#' \href{https://github.com/Jstacs/Jstacs/tree/master/projects/talecorrect}{page}.
#' It takes an input fasta file and output a file with corrected indels in TALE coding sequences
#' using an approach described in the Erkes et al. \href{https://doi.org/10.1186/s12864-023-09228-1}{paper}.
#' \code{\link[tantale:tellTale]{tellTale}} can subsequently be run on the corrected
#' sequences in no correction mode.
#'
#' @param tellTaleOutDir Path to a \code{\link[tantale:tellTale]{tellTale}} run
#'   output directory
#' @param uncorrectedAssemblyPath Path to the input sequence file
#' @param correctedAssemblyPath Path of the ouput file
#' @param pathToHMMs Path the folder containning the profile HMM files. The default
#' value points to the ones build from Xanthomonas oryzae pv. oryzae templates.
#' Xox ones are also available in the parent directory. Please see the gitHub
#' \href{https://github.com/Jstacs/Jstacs/tree/master/projects/talecorrect}{page}
#' for instructions on building custom profiles.
#' @param returnCorrectionsTble Specify \code{TRUE} if you want the list of executed
#' operations on the sequence as a tibble.
#' @param condaBinPath Path to your Conda binary file if you need to specify a
#'   path different from the one that is automatically searched by the
#'   reticulate package functions.
#' 
#' 
#' @return A tibble if \code{returnCorrectionsTble} is \code{TRUE} or the path to the
#' corrected sequences file.
#' 
#' @export
correcTales <- function(uncorrectedAssemblyPath ,
                     correctedAssemblyPath = file.path(getwd(), "correctedTALEs.fa"),
                     pathToHMMs = system.file("tools", "talecorrect", "HMMs", "Xoo", package = "tantale", mustWork = T),
                     returnCorrectionsTble = FALSE,
                     condaBinPath = "auto") {
  
  pathToTALECorrection <- system.file("tools", "talecorrect", "TALEcorrection.jar", package = "tantale", mustWork = T)
  outputFolder <- tempfile(pattern = "correcTales")
  dir.exists(outputFolder) || dir.create(outputFolder, recursive = TRUE)
  domains <- c(N = "N-terminus.10bpRepeat1", C = "repeat", R = "C-terminus")
  invisible(fs::file_exists(uncorrectedAssemblyPath)) || logger::log_error("The provided input file does not exists") & stop()
  
  #### run nHMMER ####
  nhmmerCmd <- paste(g("nhmmer {pathToHMMs}/{domains}.hmm {uncorrectedAssemblyPath} > {outputFolder}/out_nhmmer.{domains}.txt",
                       .sep = "; "), collapse = "; ")
  envReady <- !as.logical(createTantaleEnv(condaBinPath = condaBinPath))
  if (envReady) {
    logger::log_info("Running nHMMER")
    logger::log_debug("Invoking nHMMER using the following command:\n {nhmmerCmd}")
    res <- systemInCondaEnv(envName = "tantale",
                            condaBinPath = condaBinPath,
                            command = nhmmerCmd
    )
    if (res) {
      logger::log_error("The following nHMMER commands failed:")
      logger::log_error("{nhmmerCmd}")
      stop()
    }
  } else {
    stop("Could not create the tantale conda environment on your machine to run nHMMER...")
  }
  
  #### run TALEcorrection ####
  logger::log_info("Performing TALEs cds correction on provided sequences.")
  talecorCmd <- g("java -jar {pathToTALECorrection} correct s={uncorrectedAssemblyPath}",
                  "n={outputFolder}/out_nhmmer.{domains[\"N\"]}.txt r={outputFolder}/out_nhmmer.{domains[\"R\"]}.txt",
                  "c={outputFolder}/out_nhmmer.{domains[\"C\"]}.txt outdir={outputFolder}", .sep = " ")
  res <- try(system(command = talecorCmd, intern = TRUE))
  if (class(res) == "try-error") {
    logger::log_error("The following TALEcorrection commands failed:")
    logger::log_error("{talecorCmd}")
    stop()
  }
  correctionsTble <- readr::read_tsv(file = file.path(outputFolder, "substitionList.tsv"),
                                     show_col_types = FALSE) %>%
    dplyr::rename(posInOriginSeq = `position in uncorrected sequences`)
  
  file.copy(file.path(outputFolder, "correctedTALEs.fa"), file.path(correctedAssemblyPath),
            overwrite = FALSE)
  unlink(outputFolder, recursive = TRUE) 
  
  if (returnCorrectionsTble) {
    return(correctionsTble)
  } else {
    logger::log_debug("Here is the list of fixes made to the provided sequence:")
    logger::skip_formatter(as.character(knitr::kable(correctionsTble))) %>% logger::log_debug()
    return(file.path(correctedAssemblyPath))
  }
}







