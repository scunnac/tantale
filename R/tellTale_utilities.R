

.get_hmmer <- function() {
  pathOfHmmerBinsDir <- system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = TRUE)
  return(pathOfHmmerBinsDir)
}


.check_hmmer <- function(hmmer_path) {
  cmd <- file.path(hmmer_path, "hmmsearch -h | grep \"^#\"")
  if (system(command = cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)) {
    stop("HMMER is not in PATH. Follow instructions at http://hmmer.org/documentation.html to install it.")
  } else {
    out <- system(command = cmd,intern = TRUE)
    cli::cli_inform(gsub("^#[ ]?", "", out[2:3]))
  }
}






.run_nhmmer_search <- function(hmmer_path = NULL, subject_file, hmm_file, search_out_file, readable_out_file) {
  if (is.null(hmmer_path)) hmmer_path <- .get_hmmer()
  .check_hmmer(hmmer_path)
  searchCmd <- paste(file.path(hmmer_path, "nhmmer"),
                     "--tblout",
                     search_out_file,
                     hmm_file,
                     subject_file,
                     ">",
                     readable_out_file,
                     sep = " "
  )
  system(command = searchCmd, ignore.stderr = FALSE, intern = TRUE)
}




.correction_tibble <- function(indels) {
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
#' @importFrom methods setClass
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


.hits_report_to_gff <- function(f = "hitsReport.csv") {
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


