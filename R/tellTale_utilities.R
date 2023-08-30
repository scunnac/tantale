

getHMMER <- function() {
  pathOfHmmerBinsDir <- system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = TRUE)
  return(pathOfHmmerBinsDir)
}


checkHMMER <- function(hmmerpath) {
  cmd <- file.path(hmmerpath, "hmmsearch -h | grep \"^#\"")
  if (system(command = cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)) {
    stop("HMMER is not in PATH. Follow instructions at http://hmmer.org/documentation.html to install it.")
  } else {
    out <- system(command = cmd,intern = TRUE)
    cat(out[2:3], sep = "\n")
  }
}


writeHMMFile <- function(hmmerpath = NULL, alignmentFile, HMMOutFile) {
  if (is.null(hmmerpath)) hmmerpath <- getHMMER()
  checkHMMER(hmmerpath)
  buildCmd <- paste(file.path(hmmerpath,"hmmbuild"),
                    HMMOutFile,
                    alignmentFile,
                    sep = " ")
  commandOut <- system(command = buildCmd, ignore.stderr = FALSE, intern = TRUE)
  return(commandOut)
}


runHmmerSearch <- function(hmmerpath = NULL, subjectFile, hmmFile, searchTblOutFile, humReadableOutFile) {
  if (is.null(hmmerpath)) hmmerpath <- getHMMER()
  checkHMMER(hmmerpath)
  searchCmd <- paste(file.path(hmmerpath, "hmmsearch"),
                     "--tblout",
                     searchTblOutFile,
                     hmmFile,
                     subjectFile,
                     ">",
                     humReadableOutFile,
                     sep = " "
  )
  system(command = searchCmd, ignore.stderr = FALSE, intern = TRUE)
}


runNhmmerSearch <-  function(hmmerpath = NULL, subjectFile, hmmFile, searchTblOutFile, humReadableOutFile) {
  if (is.null(hmmerpath)) hmmerpath <- getHMMER()
  checkHMMER(hmmerpath)
  searchCmd <- paste(file.path(hmmerpath, "nhmmer"),
                     "--tblout",
                     searchTblOutFile,
                     hmmFile,
                     subjectFile,
                     ">",
                     humReadableOutFile,
                     sep = " "
  )
  system(command = searchCmd, ignore.stderr = FALSE, intern = TRUE)
}


runHmmalign <- function(hmmerpath = NULL, hmmFile, seqsFile, alignOutFile) {
  if (is.null(hmmerpath)) hmmerpath <- getHMMER()
  checkHMMER(hmmerpath)
  alignCmd <- paste(file.path(hmmerpath, "hmmalign"),
                    "--outformat Phylip", #Stockholm, SELEX, Clustal, Phylip, Pfam, A2M, PSIBLAST.
                    "--trim",
                    hmmFile,
                    seqsFile,
                    ">", alignOutFile,
                    sep = " "
  )
  system(command = alignCmd, ignore.stderr = FALSE, intern = TRUE)
}


correction_tibble <- function(indels) {
  indelsTble <- lapply(indels, function(lst) {
    info <- tibble::tibble()
    colnames(info) <- c("variable","value")
    for (talOrfID in 1:length(lst)) {
      element <- lst[talOrfID]
      if (length(unlist(element)) == 0L) {
        next()
      } else {
        info %<>% dplyr::bind_rows(tibble::tibble(variable = names(element), value = as.numeric(unlist(element))))
      }
    }
    return(info)
  }
  ) %>% dplyr::bind_rows(.id = "Seq")
  return(indelsTble)
} 


#' annotale output 
#' @exportClass annout
#' @import Biostrings
annout <- setClass(
  # Set the name for the class
  Class = "annout",
  
  # Define the slots
  slots = c(
    domainsReport = "data.frame"
  ),
  
  contains = "AAStringSet",
  
  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  validity = function(object) {
    val <- is.data.frame(object@domainsReport)
    return(val)
  }
)


hitsReportToGFF <- function(f = "hitsReport.csv") {
  # Convert the info contained in a HitReport file into a GFF file for display by
  # a genome viewer.
  # The f parameter corresponds to the path to a hitsReport file.
  # Read the file as a data.frame
  hitsReport <- read.delim(f)
  # Create a GenomicRange that will be converted.
  hitsGR <- GenomicRanges::makeGRangesFromDataFrame(hitsReport, keep.extra.columns=TRUE)
  # Write a gff3 file to disk with this info.
  rtracklayer::export.gff3(hitsGR,
                           con = file.path(dirname(f), paste0(sub("\\..*$", "", basename(f)), ".gff"))
  )
  
}

## !! THIS SHOULD BE MADE OBSOLETE AND CODE USING IT SHOULD BE MODIFIED
extractSeqsfromHits <- function(nhmmerTabularOutputSelect, DNAsequences){
  repeatSeqsSetList <- mapply(
    function(hitID, start, end, strand, subjectID, sequences) {
      seq <- XVector::subseq(sequences[subjectID], start, end)
      if (strand == "-") {seq <- Biostrings::reverseComplement(seq)}
      names(seq) <- hitID
      return(seq)
    },
    hitID = nhmmerTabularOutputSelect$hitID,
    start = nhmmerTabularOutputSelect$start,
    end = nhmmerTabularOutputSelect$end,
    strand = nhmmerTabularOutputSelect$strand,
    subjectID = nhmmerTabularOutputSelect$target_name,
    MoreArgs = list(sequences = DNAsequences),
    USE.NAMES = FALSE)
  do.call(c, repeatSeqsSetList)
}

