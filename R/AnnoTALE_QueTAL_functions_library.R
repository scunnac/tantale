
##### Utility functions ####
## -convert TALE RVD sequences into different formats
## - run external TALE clustering tools (FuncTALE and AnnoTALE)




.reformat_array_report <- function(f, name_prefix = NULL) {
  # Trying to guess prefix if it is not provided
  strain <- name_prefix
  if(is.null(name_prefix)) {
    #strain <- regmatches(dirname(f), regexpr("MAI[0-9]{1,2}", dirname(f)))
    strain <- sub("^.*Sebra/(.+)_[0-9].+$", "\\1", dirname(f))
  }
  if(identical(length(strain), 0L)) {strain <- "NA"}

  # Fetch the content of an arrayReportFile
  TALEs <- subset(read.delim(f), selectedForAssembly)

  # Reformating TALE arrays seq of RVDs to comply with functal requirements
  pattern <- "^(BBB-)*([^(ZZZ)]*)(-ZZZ)*$"
  isFullLength <- grepl(pattern, TALEs$SeqOfRVD)
  fullLengthTALEs <- TALEs[isFullLength, c("arrayID", "SeqOfRVD")] # Keeping only full length RVD arrays
  fullLengthTALEs$SeqOfRVD <- gsub(pattern, "\\2", fullLengthTALEs$SeqOfRVD) # Remonving start and end flags on sequences
  TALENames <- paste(strain, fullLengthTALEs$arrayID, sep = "x")

  paste0(">", TALENames, "\t", fullLengthTALEs$SeqOfRVD) # QueTAL formatted TALEs

}


.annotale_to_quetal_rvd <- function(input_file, output_file = "RVDSeqs.QueTal.fasta") {
  # Need a AnnoTALE "TALE_RVDs.fasta" - like RVD file
  TALERVDSeqs <- Biostrings::readBStringSet(filepath = input_file)
  TALERVDSeqs <- as.character(TALERVDSeqs)
  writeLines(text = paste0(">", names(TALERVDSeqs), "\t", TALERVDSeqs), con = output_file) # functal formatted TALEs writen in text file
}


.quetal_to_annotale_rvd <- function(input_file, output_file = "RVDSeqs.AnnoTALE.fasta") {
  # Need a file with RVD sequences in the QueTal specific format (>TaleA\tNN-NH-N*)
  TALERVDSeqs <- read.table(input_file, header = FALSE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)
  # Construct an XStringSet
  TALERVDSeqsBS <- Biostrings::BStringSet(x=TALERVDSeqs[,2])
  names(TALERVDSeqsBS) <- TALERVDSeqs[,1]
  # Remove heading '>' sign
  names(TALERVDSeqsBS) <- gsub(pattern = "^>", replacement = "", x = names(TALERVDSeqsBS), perl = TRUE)
  # Make the names prettier
  names(TALERVDSeqsBS) <- gsub(pattern = "(MAI\\d{1,3}).*TALE(\\d{1,3}).*$",
                               replacement = "Tal\\2-\\1", x = names(TALERVDSeqsBS), perl = TRUE)
  # Write the sequences to disc in fasta format
  Biostrings::writeXStringSet(TALERVDSeqsBS, filepath = output_file, format="fasta")
}



#### run external TALE RVD sequences inferrence and clustering tools ####

#' Runs the "predict" and "analyze" steps of AnnoTALE on a fasta file.
#'
#' A R wrapper around the
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/26876161}{AnnoTALE} 'AnnoTALE.jar
#' predict' and 'AnnoTALE.jar analyze' shell calls. The whole AnnoTALE workflow
#' can be completed by a subsequent call to the \code{\link{run_annotale_build}}
#' function.
#'
#' @param fasta_file Path to a fasta file containing DNA (?) sequences to be
#'   analyzed for TALE content.
#' @param output_dir Directory where output will be written (created if does not
#'   exist).
#' @param prefix A scalar character vector containing a prefix that will be
#'   appended to TALE names by AnnoTALE. If not supplied, the function will try
#'   to guess the prefix from the input file name.
#' @param annotale_jar Path to the AnnoTALE jar file if you want to use another
#'   version than the one provided with tantale.
#' @return Returns invisibly the exit code of the shell call to the last
#'   AnnoTALE step (ie '0' if successful).
#' @export
#' @family external TALE tools
run_annotale_predict <- function(fasta_file,
                            output_dir = getwd(),
                            prefix = NULL,
                            annotale_jar = system.file("tools", "AnnoTALEcli-1.5.jar", package = "tantale", mustWork = T)
                            ) {
  # Define output dirs for the various stages of AnnoTALE
  stopifnot(dir.exists(output_dir) || dir.create(path = output_dir, showWarnings = TRUE, recursive = TRUE, mode = "775"))
  predict_dir <- file.path(output_dir, "Predict")
  analyze_dir <- file.path(output_dir, "Analyze")
  # Define a prefix for TALEs (the strain or assembly ID) derived from the genome file name.
  if( is.null(prefix) ) {
  prefix <- gsub(pattern = "^(.*)\\.(fasta|fa|fas)$" , replacement  = "\\1", basename(fasta_file), perl = TRUE)
  }
  # Run the "predict" stage of AnnoTALE
  comPredict <- paste0(
    "java -jar ", annotale_jar,
    " predict",
    " g=", fasta_file,
    " s=", prefix,
    " outdir=", predict_dir
  )
  cli::cli_inform(c("Running AnnoTALE predict for {.val {prefix}}", " " = "{comPredict}"))
  exitPredict <- system(comPredict)
  !exitPredict || stop("##  AnnoTALE predict failed with an error. Aborting...")

  # Run the "analyze" stage of AnnoTALE
  comAnalyze <- paste0(
    "java -jar ", annotale_jar,
    " analyze ",
    " t=", shQuote(list.files(predict_dir, pattern = "^TALE_DNA_sequences_", full.names = TRUE)),
    " outdir=", shQuote(analyze_dir)
  )
  cli::cli_inform(c("Running AnnoTALE analyze for {.val {prefix}}", " " = "{comAnalyze}"))
  exitAnalyze <- system(comAnalyze)
  return(invisible(exitAnalyze))
}






#' Run the "build" stage of AnnoTALE.
#'
#' A R wrapper around the \href{https://www.ncbi.nlm.nih.gov/pubmed/26876161}{AnnoTALE} 'AnnoTALE.jar build' program.
#' It usually takes is input from the file generated by the \code{\link{run_annotale_predict}}
#' function.
#'
#' @param fasta_file Path to a fasta file containing TALE sequences as
#'   returned by AnnoTALE (?) to be classified into groups.
#' @param output_dir Directory where output will be written (created if does not
#'   exist).
#' @param annotale_jar Path to the AnnoTALE jar file if you want to use another
#'   version than the one provided with tantale.
#' @return Returns invisibly the exit code of the shell call to Annotale (ie '0' if successful).
#' @export
#' @family external TALE tools
run_annotale_build <- function(fasta_file,
                          output_dir = getwd(),
                          annotale_jar = system.file("tools", "AnnoTALEcli-1.5.jar", package = "tantale", mustWork = T)
                          ) {
  if(! dir.exists(output_dir)) dir.create(path = output_dir, showWarnings = TRUE, recursive = TRUE, mode = "775")
  comBuild <- paste0(
    "java -Xms512M -Xmx6G -jar ", annotale_jar,
    " build ",
    " t=", shQuote(fasta_file),
    " outdir=", shQuote(output_dir)
  )
  cli::cli_inform(c("Running AnnoTALE build", " " = "{comBuild}"))
  exitBuild <- system(comBuild)
  return(invisible(exitBuild))
}



#' Run functal from QueTAL to build a phylogenetic tree of TALE RVD sequences.
#'
#' A R wrapper around the \href{https://doi.org/10.3389/fpls.2015.00545}{QueTAL} 'functal' perl script.
#'
#' @param tal_file Path to a QueTAL-formatted file of TALE RVD sequences.
#' @param tree_format Tree layout passed to functal's `-n` option (default `"fan"`).
#' @param output_prefix Prefix used for functal's output file names.
#' @param output_dir Directory where output will be copied (default: current working directory).
#' @param functal_path Path to the functal perl script if you want to use another
#'   version than the one provided with tantale.
#' @param conda_bin Path to your Conda binary file if you need to specify a
#'   non-standard location, otherwise leave to "auto".
#' @return Returns invisibly the exit code of the shell call to functal (ie '0' if successful).
#' @export
#' @family external TALE tools
functal <- function(tal_file,
                    tree_format = "fan",
                    output_prefix = "FuncTALE",
                    output_dir = getwd(),
                    functal_path = system.file("tools", "QueTAL_v1.1", "FuncTAL", "FuncTAL_v.1.1.pl", package = "tantale", mustWork = T),
                    conda_bin = "auto") {
  # Running functal from QueTAL_v1.1
  # A few observations:
  # Refuse to use another output directory than the "Ouputs" one in the program folder
  # Cannot invoke the program from another working directory than the one where the pl script is located
  # Crashes when provided the CDS of the TALES from Hinda's Malian strains
  # So here is a caller function to get around these issues:

  functal_dir <- dirname(functal_path)

  # Assembling the command to be run
  # -I functal_dir is required for perl to find Statistics.pm, which ships
  # alongside the script rather than as an installed module. List::MoreUtils
  # and Bio::Perl come from the 'tantale' conda environment.
  functal_cmd <- paste("perl", "-I", functal_dir, functal_path,
                      "-n", tree_format,
                      tal_file,
                      output_prefix)

  # Run the command inside the tantale conda environment
  cli::cli_inform(c("Running functal", " " = "{functal_cmd}"))
  envReady <- !as.logical(.create_tantale_env(conda_bin = conda_bin))
  if (envReady) {
    exitCom <- .run_in_conda(env_name = "tantale",
                                conda_bin = conda_bin,
                                command = functal_cmd,
                                cwd = functal_dir)
  } else {
    stop("Could not create the tantale conda environment on your machine to run functal...")
  }
  if (exitCom != 0) {
    stop("functal failed (perl exit code ", exitCom, "). See console output above for details.")
  }

  # Transferring the ouput to the output dir and deleting it in the functal "Outputs" directory
  functal_output_files <-
    list.files(file.path(functal_dir, "Outputs"), full.names = TRUE)
  file.copy(
    from = functal_output_files,
    to = output_dir,
    overwrite = FALSE,
    recursive = FALSE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  unlink(functal_output_files, recursive = TRUE, force = FALSE)
  return(invisible(exitCom))
}
