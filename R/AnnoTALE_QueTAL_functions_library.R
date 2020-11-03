
##### Utility functions ####
## -convert TALE RVD sequences into different formats
## - run external TALE clustering tools (FuncTALE and AnnoTALE)
## - plot pretty alignments between RVD seq and EBE



reformatTALEsFromArrayReportFile <- function(f, namePrefix = NULL) {
  # Trying to guess prefix if it is not provided
  strain <- namePrefix
  if(is.null(namePrefix)) {
    #strain <- regmatches(dirname(f), regexpr("MAI[0-9]{1,2}", dirname(f)))
    strain <- sub("^.*Sebra/(.+)_[0-9].+$", "\\1", dirname(f))
  }
  if(identical(length(strain), 0L)) {strain <- "NA"}

  # Fetch the content of an arrayReportFile
  TALEs <- subset(read.delim(f), selectedForAssembly)

  # Reformating TALE arrays seq of RVDs to comply with FuncTAL requirements
  pattern <- "^(BBB-)*([^(ZZZ)]*)(-ZZZ)*$"
  isFullLength <- grepl(pattern, TALEs$SeqOfRVD)
  fullLengthTALEs <- TALEs[isFullLength, c("arrayID", "SeqOfRVD")] # Keeping only full length RVD arrays
  fullLengthTALEs$SeqOfRVD <- gsub(pattern, "\\2", fullLengthTALEs$SeqOfRVD) # Remonving start and end flags on sequences
  TALENames <- paste(strain, fullLengthTALEs$arrayID, sep = "x")

  paste0(">", TALENames, "\t", fullLengthTALEs$SeqOfRVD) # QueTAL formatted TALEs

}


AnnoTALE2QueTALRVD <- function(inputFile, outputFile = "RVDSeqs.QueTal.fasta") {
  # Need a AnnoTALE "TALE_RVDs.fasta" - like RVD file
  TALERVDSeqs <- Biostrings::readBStringSet(filepath = inputFile)
  TALERVDSeqs <- as.character(TALERVDSeqs)
  writeLines(text = paste0(">", names(TALERVDSeqs), "\t", TALERVDSeqs), con = outputFile) # FuncTAL formatted TALEs writen in text file
}


QueTALRVD2AnnoTALE <- function(inputFile, outputFile = "RVDSeqs.AnnoTALE.fasta") {
  # Need a file with RVD sequences in the QueTal specific format (>TaleA\tNN-NH-N*)
  TALERVDSeqs <- read.table(inputFile, header = FALSE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)
  # Construct an XStringSet
  TALERVDSeqsBS <- Biostrings::BStringSet(x=TALERVDSeqs[,2])
  names(TALERVDSeqsBS) <- TALERVDSeqs[,1]
  # Remove heading '>' sign
  names(TALERVDSeqsBS) <- gsub(pattern = "^>", replacement = "", x = names(TALERVDSeqsBS), perl = TRUE)
  # Make the names prettier
  names(TALERVDSeqsBS) <- gsub(pattern = "(MAI\\d{1,3}).*TALE(\\d{1,3}).*$",
                               replacement = "Tal\\2-\\1", x = names(TALERVDSeqsBS), perl = TRUE)
  # Write the sequences to disc in fasta format
  Biostrings::writeXStringSet(TALERVDSeqsBS, filepath = outputFile, format="fasta")
}



#### run external TALE RVD sequences inferrence and clustering tools ####

#' Runs the "predict" and "analyze" steps of AnnoTALE on a fasta file.
#'
#' A R wrapper around the
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/26876161}{AnnoTALE} 'AnnoTALE.jar
#' predict' and 'AnnoTALE.jar analyze' shell calls. The whole AnnoTALE workflow
#' can be completed by a subsequent call to the \code{\link{buildAnnoTALE}}
#' function.
#'
#' @param inputFastaFile Path to a fasta file containing DNA (?) sequences to be
#'   analyzed for TALE content.
#' @param outputDir Directory where output will be written (created if does not
#'   exist).
#' @param prefix A scalar character vector containing a prefix that will be
#'   appended to TALE names by AnnoTALE. If not supplied, the function will try
#'   to guess the prefix from the input file name.
#' @param annoTALE Path to the AnnoTALE jar file if you want to use another
#'   version than the one provided with tantale.
#' @return Returns invisibly the edit code of the shell call to the last
#'   AnnoTALE step (ie '0' if successful).
#' @export
analyzeAnnoTALE <- function(inputFastaFile,
                            outputDir = getwd(),
                            prefix = NULL,
                            annoTALE = system.file("tools", "AnnoTALEcli-1.4.1.jar", package = "tantale", mustWork = T)
                            ) {
  # Define output dirs for the various stages of annoTALE
  stopifnot(dir.exists(outputDir) || dir.create(path = outputDir, showWarnings = TRUE, recursive = TRUE, mode = "775"))
  annoTALEPredictDir <- file.path(outputDir, "Predict")
  annoTALEAnalyzeDir <- file.path(outputDir, "Analyze")
  # Define a prefix for TALEs (the strain or assembly ID) derived from the genome file name.
  if( is.null(prefix) ) {
  prefix <- gsub(pattern = "^(.*)\\.(fasta|fa|fas)$" , replacement  = "\\1", basename(inputFastaFile), perl = TRUE)
  }
  # Run the "predict" stage of annoTALE
  comPredict <- paste0(
    "java -jar ", annoTALE,
    " predict",
    " g=", inputFastaFile,
    " s=", prefix,
    " outdir=", annoTALEPredictDir
  )
  cat("##  Now running annoTALE predict for", prefix, "using the following command:\n##  ",  comPredict, "\n")
  exitPredict <- system(comPredict)
  ! exitPredict || stop("##  annoTALE predict failed with an error. Aborting...")

  # Run the "analyze" stage of annoTALE
  comAnalyze <- paste0(
    "java -jar ", annoTALE,
    " analyze ",
    " t=", shQuote(list.files(annoTALEPredictDir, pattern = "^TALE_DNA_sequences_", full.names = TRUE)),
    " outdir=", shQuote(annoTALEAnalyzeDir)
  )
  cat("##  Now running annoTALE analyze for", prefix, "using the following command:\n##  ",  comAnalyze, "\n")
  exitAnalyze <- system(comAnalyze)
  return(invisible(exitAnalyze))
}






#' Run the "build" stage of AnnoTALE.
#'
#' A R wrapper around the \href{https://www.ncbi.nlm.nih.gov/pubmed/26876161}{AnnoTALE} 'AnnoTALE.jar build' program.
#' It usually takes is input from the file generated by the \code{\link{analyzeAnnoTALE}}
#' function.
#'
#' @param TALESeqsFastaFile Path to a fasta file containing TALE sequences as
#'   returned by AnnoTALE (?) to be classified into groups.
#' @param outputDir Directory where output will be written (created if does not
#'   exist).
#' @param annoTALE Path to the AnnoTALE jar file if you want to use another
#'   version than the one provided with tantale.
#' @return Returns invisibly the exit code of the shell call to Annotale (ie '0' if successful).
#' @export
buildAnnoTALE <- function(TALESeqsFastaFile,
                          outputDir = getwd(),
                          annoTALE = system.file("tools", "AnnoTALEcli-1.4.1.jar", package = "tantale", mustWork = T)
                          ) {
  if(! dir.exists(outputDir)) dir.create(path = outputDir, showWarnings = TRUE, recursive = TRUE, mode = "775")
  comBuild <- paste0(
    "java -Xms512M -Xmx6G -jar ", annoTALE,
    " build ",
    " t=", shQuote(TALESeqsFastaFile),
    " outdir=", shQuote(outputDir)
  )
  cat("Now running annoTALE build using the following command:\n",  comBuild, "\n")
  exitBuild <- system(comBuild)
  return(invisible(exitBuild))
}



#' @export
FuncTAL <- function(TALfile,
                    treeFormat = "fan",
                    outputPrefix = "FuncTALE",
                    outputDir = "/home/cunnac/Documents",
                    FunctTAL = system.file("tools", "QueTAL_v1.1", "FuncTAL", "FuncTAL_v.1.1.pl", package = "tantale", mustWork = T)) {
  # Running FuncTAL from QueTAL_v1.1
  # A few observations:
  # Refuse to use another output directory than the "Ouputs" one in the program folder
  # Cannot invoke the program from another working directory than the one where the pl script is located
  # Crashes when provided the CDS of the TALES from Hinda's Malian strains
  # So here is a caller function to get around these issues:

  FunctTALDir <- dirname(FunctTAL)

  # Assembling the command to be run
  goToFuncTALDir <- paste("cd", FunctTALDir)
  runFuncTAL <- paste(FunctTAL,
                      "-n", treeFormat,
                      TALfile,
                      outputPrefix)
  com <- paste(goToFuncTALDir, runFuncTAL, sep = ";")

  # Run the command with system
  cat("Now running FunctTAL using the following command:\n",
      com,
      "\n\n")
  exitCom <- system(com)

  # Transferring the ouput to the output dir and deleting it in the FuncTAL "Outputs" directory
  FuncTALEOuputFiles <-
    list.files(file.path(FunctTALDir, "Outputs"), full.names = TRUE)
  file.copy(
    from = FuncTALEOuputFiles,
    to = outputDir,
    overwrite = FALSE,
    recursive = FALSE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  unlink(FuncTALEOuputFiles, recursive = TRUE, force = FALSE)
}
