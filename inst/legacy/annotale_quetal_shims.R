# Retired: AnnoTALE <-> QueTAL RVD-format shims
#
# Thin file-format adapters between two external tools incompatible on-disk
# conventions for the same content -- AnnoTALE writes RVDs as standard FASTA,
# QueTAL/FuncTAL wants one line per TALE with a tab between id and RVD string.
# No biological or computational transformation happens in any of them.
#
# Retired because they are dead as a set: they fed the FuncTAL branch, which is
# itself dormant. Kept rather than deleted, since the format conventions they
# encode are the only record of how the two tools disagree.
#
# Moved out of R/ on 2026-09-13. See dev/restructuring-notes.md section 3.

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
