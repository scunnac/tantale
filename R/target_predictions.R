

##### Compute TALE targets predictions #####


#' Run TALE target predictions on DNA sequence(s) using PrediTale.
#'
#'
#' A R wrapper around the
#' \href{https://www.jstacs.de/index.php/PrediTALE}{PrediTale} 'PrediTALE.jar
#' preditale' module. Takes a list of TALE RVD sequences and a fasta file of DNA
#' sequences and runs PrediTale.
#'
#' @param rvdSeqs Tale RVD sequences are supplied as either a fasta file (atomic
#'   character vector) with Tale info (name) in title and sequences of RVD as a
#'   space or '-' separeted string or as a Biostrings XStringSet with sequences
#'   of RVD similarly formated. See the
#'   \href{https://www.jstacs.de/index.php/PrediTALE}{PrediTale} man page for
#'   how to encode RVDs present on aberrant repeats.
#' @param subjDnaSeqFile Expects a character vector specifying the path to the
#'   fasta file holding subject DNA sequence(s).
#' @param optParam An atomic character vector specifying optionnal parameters
#'   for PrediTALE.jar preditale (eg "Strand=\"forward strand\"").
#' @param outDir Expects a character vector specifying the path to an output
#'   directory. If not supplied, output files will be temporary.
#' @param predictorPath If you want to use another version of "PrediTALE.jar"
#'   than the one supplied with tantale, specify its path here.
#' @return A tibble with the EBE predictions. \strong{Note that column names
#'   have been modified} relative to the column names found in the originale
#'   programs's output in order to homogenize column names across TALE target
#'   prediction programs in tantale
#' @export
preditale <- function(rvdSeqs, subjDnaSeqFile, optParam = "", outDir = NULL,
                      predictorPath = system.file("tools", "PrediTALE.jar", package = "tantale", mustWork = T)) {
  # Checking input args
  if (class(rvdSeqs) == "character") {
    rvdSeqsTest <- Biostrings::readBStringSet(rvdSeqs)
    rvdSeqsFile <- rvdSeqs
  } else if (class(rvdSeqs) == "BStringSet") {
    rvdSeqsFile <- tempfile()
    Biostrings::writeXStringSet(rvdSeqs, rvdSeqsFile)
  } else {
    stop("##  Something is wrong with the value provided for rvdSeqs. It must be either\n",
         "##  the path to a fasta file containing strings of RVD sequences (space/dash-separated\n",
         "##  rvd) or a Biostrings XStringSet object.")
  }
  if(!file.exists(subjDnaSeqFile)) stop("Unable to find the set of target DNA sequences (fasta file) at the specified location. Please verify the file exists")
  if (is.null(outDir)) {
    outDir <- tempfile(pattern = "preditale_")
    dir.create(outDir, recursive = TRUE)
  }
  if (length(f <- list.files(path = outDir, pattern = "^Predicted_binding.*tsv$", full.names = TRUE)) != 0) {
    stop("The output directory '", outDir, "' already contains files that are possibly previous results of the preditale function. ",
    "Cannot proceed. Please remove the following files:\n", paste("-", f, sep = " ", collapse = "\n"))
  }
  # Assembling preditale command
  cmd <- glue::glue("java -Xms512M -Xmx2G -jar {predictorPath} preditale {optParam} TALEs={rvdSeqsFile} s={subjDnaSeqFile} outdir={outDir}")
  glue::glue("## Invoking Preditale using the following command:\n", stringr::str_wrap(cmd, 80), "\n")
  # Running Preditale
  system(command = cmd)
  # Parsing output
  predFiles <- list.files(path = outDir, pattern = "^Predicted_binding.*tsv$", full.names = TRUE)
  predictions <- lapply(predFiles, function(f) {
    suppressMessages(readr::read_tsv(f, show_col_types = FALSE))
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::rename(subjSeqId = `# Seq-ID`,
                  start = Position,
                  strand = Strand,
                  score = Score,
                  ebeSeq = Sequence,
                  pval = `Approx. p-value`,
                  rvds = RVDs,
                  taleId = TALE) %>%
    dplyr::mutate(start = start + 1) %>%
    dplyr::mutate(end = nchar(ebeSeq) + (start - 1)) %>%
    dplyr::relocate(end, .after = start) %>%
    dplyr::mutate_if(is.factor, as.character)
  # Returning
  return(predictions)
}

# java -Xms512M -Xmx2G -jar /home/cunnac/Lab-Related/MyScripts/xanthopore-scripts/tantale/inst/tools/PrediTALE.jar preditale \
# TALEs=/home/cunnac/Lab-Related/MyScripts/xanthopore-scripts/tantale/inst/extdata/cladeIII_sweet_targeting_control_TALEs.fa \
# s=/home/cunnac/Lab-Related/MyScripts/xanthopore-scripts/tantale/inst/extdata/cladeIII_sweet_promoters.fasta \
# outdir=$(pwd)

#' Run TALE target predictions on DNA sequence(s) using Talvez
#'
#'
#' A R wrapper around the
#' \href{https://doi.org/10.1371/journal.pone.0068464}{Talvez} predictor perl
#' script. Takes a list of TALE RVD sequences and a fasta file of DNA sequences
#' and runs Talvez
#'
#' Note that this talvez wrapper, uses the RVD - Nucleotide specificity matrices
#' used with talvez, it is not possible to use custom ones.
#' Note also that the talvez script is run in a conda environment providing the
#' necessary dependencies. This environment will be created automatically if
#' necessary.
#'
#' @param rvdSeqs Tale RVD sequences are supplied as either a fasta file (atomic
#'   character vector) with Tale info (name) in title and sequences of RVD as a
#'   space or '-' separeted string or as a Biostrings XStringSet with sequences
#'   of RVD similarly formated.
#' @param subjDnaSeqFile Expects a character vector specifying the path to the
#'   fasta file holding subject DNA sequence(s). Talvez forbids to have
#'   sequences in the file wrapped at a fixed width. The function uses
#'   Biostrings to unwrap them but if you use subject sequences longer than
#'   20kb, this will fail and you are advised to unwrap your sequences before
#'   hand.
#' @param optParam An atomic character vector specifying optionnal parameters
#'   for the Talvez script (eg "-t 0 -l 19"). \strong{These may not include} the
#'   '-e' and '-z' options specifying the matrix files.
#' @param outDir Expects a character vector specifying the path to an output
#'   directory. If not supplied, output files will be temporary.
#' @param talvezDir If you want to use another version of Talvez than the one
#'   supplied with tantale, specify the path of the directory containing the
#'   necessary files here.
#' @param condaBinPath Path to your Conda binary file if you need to specify
#'   a path different from the one that is automatically searched by the
#'   reticulate package functions.
#' @return A tibble with the EBE predictions. \strong{Note that column names
#'   have been modified} relative to the column names found in the originale
#'   programs's output in order to homogenize column names across TALE target
#'   prediction programs in tantale.
#' @export
talvez <- function(rvdSeqs, subjDnaSeqFile, optParam = "-t 0 -l 19", outDir = NULL,
                   talvezDir = system.file("tools", "TALVEZ_3.2", package = "tantale", mustWork = T),
                   condaBinPath = "auto") {

  # Checking input args
  if (class(rvdSeqs) == "character") {
    rvdSeqs <- Biostrings::readBStringSet(rvdSeqs)
  } else if (class(rvdSeqs) == "BStringSet") {
    rvdSeqs <- rvdSeqs
  } else {
    stop("##  Something is wrong with the value provided for rvdSeqs. It must be either\n",
         "##  the path to a fasta file containing strings of RVD sequences (space/dash-separated\n",
         "##  rvd) or a Biostrings XStringSet object.")
  }
  if (!file.exists(subjDnaSeqFile)) stop("Unable to find the set of target DNA sequences (fasta file)",
                                        " at the specified location. Please verify the file exists")

  # Creating a temporary output dir to run everything inside it
  tempOutDir <- tempfile(pattern = "talvez_")
  dir.create(tempOutDir, recursive = TRUE)

  # Copying a reformatted copy of subject DNA seq fasta file to tempOutDir
  # This is necessary if only to make sure that fixed width formatting of
  # sequences is loosened otherwise talvez complains.
  Biostrings::writeXStringSet(x = Biostrings::readDNAStringSet(filepath = subjDnaSeqFile),
                              filepath = file.path(tempOutDir, basename(subjDnaSeqFile)),
                              format = "fasta", width = 20000L)

  # Copying the content of the talvez scripts dir in the temp tempOutDir
  # because it considerably facilitate interaction with talvez.
  if(!all(
    file.copy(from = list.files(talvezDir, full.names = TRUE, include.dirs = TRUE),
              to = tempOutDir, overwrite = TRUE, recursive = TRUE)
  )) stop("Unable to copy talvez scripts to temporary location...")

  # Formatting rvd sequences to fit the talvez format and write to tempfile
  rvdSeqsFileForTv <- tempfile(pattern = "rvdSeqsTalvez_", tmpdir = tempOutDir, fileext = ".tsv")
  writeLines(text = paste(">", names(rvdSeqs), "\t", as.character(rvdSeqs), sep = ""),
             con = rvdSeqsFileForTv)
  
  # Assembling and running talvez command
  envReady <- !as.logical(createTantaleEnv(condaBinPath = condaBinPath))
  if (envReady) {
    cmd <- glue::glue(#"cd {tempOutDir};",
                      "perl TALVEZ_3.2.pl {optParam} -e mat1 -z mat2 {basename(rvdSeqsFileForTv)} {basename(subjDnaSeqFile)}")
    logger::log_info("Invoking Talvez using the following command:\n {stringr::str_wrap(cmd, 80)}")
    res <- systemInCondaEnv(envName = "tantale",
                            condaBinPath = condaBinPath,
                            cwd = tempOutDir,
                            command = cmd)
  } else {
    stop("Could not create the tantale conda environment on your machine to run Talvez...")
  }

  # Parsing and reformating output
  predictions <- suppressMessages(
    readr::read_tsv(file.path(tempOutDir, "output_complete"), show_col_types = FALSE)
    ) %>%
  dplyr::rename(subjSeqId = `SEQ_ID`,
                start = TALBS_start,
                end = TALBS_end,
                strand = EBEstrand,
                score = SCORE,
                ebeSeq = TALBS_sequence,
                rank = RANK,
                rvds = TAL_SEQ,
                taleId = TAL_ID) %>%
    dplyr::mutate(
      subjSeqId = gsub(pattern = "^>", "", subjSeqId),
      taleId = gsub(pattern = "^>", "", taleId),
      strand = gsub(pattern = "strand$", "", strand),
      start = start + 1 # THIS SEEEMS TO BE NECESSARY TO EXTRACT EXACT EBE FROM SUBJ SEQ...
                  ) %>%
    dplyr::select(c(1,2,3,4,5,7,8,9, 10))
  # copy output files to outDir if it is not null
  if (!is.null(outDir)) {
    list.files(tempOutDir, all.files = TRUE)
    file.copy(from = file.path(tempOutDir, c("output_complete", "tmp")), to = outDir, overwrite = TRUE, recursive = TRUE)
    # list.files(tempOutDir, all.files = TRUE)
  }
  # returning
  return(predictions)
}




# cd /home/cunnac/Lab-Related/MyScripts/xanthopore-scripts/tantale/inst/tools/TALVEZ_3.2
# perl /home/cunnac/Lab-Related/MyScripts/xanthopore-scripts/tantale/inst/tools/TALVEZ_3.2/TALVEZ_3.2.pl \
# -t 100 -l 19  \
# -e /home/cunnac/Lab-Related/MyScripts/xanthopore-scripts/tantale/inst/tools/TALVEZ_3.2/mat1 \
# -z /home/cunnac/Lab-Related/MyScripts/xanthopore-scripts/tantale/inst/tools/TALVEZ_3.2/mat2 \
# TALCONTROL_RVDs.txt \
# TALCONTROL_Promoters.fasta


##### Displaying TALE RVD sequences - predicted target DNA sequences alignemnts #####


computeRVDSeqEBESeqMatchQualityString <- function(RVDSeq, EBESeq, RVDNucAssocMat = rvdToNtAssocMat) {
  # Function that return a vector of numeric scores reflecting how good is the match between RVD and nucleotide at each successive position
  # Expects sequences as as a scalar string. They are split on "-" for RVD seqs and "" for EBE DNA sequences.
  # For each RVD
  # For the corresponding nucleotide at this position on the EBE
  # Look in the Nucl-RVD association matrix for the score in the corresponding cell
  # How does this score compare to scores of other nucleotides for this RVD.
  # Come up with some kind of a scoring function for these comparaison (best = 3, worst = 1, intermediate = 2)
  # Attribute a score refecting how well the RVD match the nucleotide.
  # Return a vector of match quality scores of lenght equal to the number of RVDs - Nucleotide pairs in input sequences
  RVDSeqVector <- unlist(stringr::str_split(RVDSeq, pattern = "-"))
  EBESeqVector <- unlist(stringr::str_split(EBESeq, pattern = ""))
  if(length(RVDSeqVector) != length(EBESeqVector)) stop("Number of elements in RVD and DNA sequences are not equal.")

  mapply(FUN = function(RSV, ESV, m) {
    RVDScores <- RVDNucAssocMat[rownames(RVDNucAssocMat) == RSV, ]
    if (nrow(RVDScores) == 0) RVDScores <- RVDNucAssocMat[rownames(RVDNucAssocMat) == "XX", ] # in case the RVD under consideration has no dedicated row in the matrix
    RVDScoresRanks <- rank(RVDScores, ties.method = "min") #
    RVDScoresRanksInverted <- max(RVDScoresRanks) + 1 - RVDScoresRanks
    if (RVDScoresRanksInverted[ESV] == min(RVDScoresRanks)) {3}
    else if (RVDScoresRanksInverted[ESV] == max(RVDScoresRanks)) {1}
    else 2
  }, RSV = RVDSeqVector, ESV = EBESeqVector, MoreArgs = list(m = RVDNucAssocMat))

}



#' Plot TALE RVD sequences along a potential DNA target region
#'
#'This function enable visual inspection of TALE target predictions results that are located \strong{within} a specified subject DNA sequence region
#'
#'RVDs sequences predicted to target an EBE on the sense strand of the DNA sequence are plotted on top of the double stranded DNA sequence in parallel to its cognate EBE which is highlighted on the corresponding strand of the DNA sequence. Those preicted to target an EBE on the opposite strand are displayed below.
#'
#'Individual RVDs are printed inside colored boxes. The color of the boxes indicate to which degree the RVD is predicted to have affinity with the corresponding nucleotide on the DNA sequence at that position relatively to other nucleotides. The RVDs labelled "OO" correspond to the first non-canonical repeat also called repeat zero in TALE protein squences.
#'
#'Numeric values inside the boxes located immediately to the right of the TALE labels reflect prediction scores.
#'
#'
#' @param predResults A tibble of prediction results obtained with \code{\link[tantale:preditale]{preditale}} or \code{\link[tantale:talvez]{talvez}} or a custom table in this format.
#' @param subjDnaSeqFile The fasta file of subject DNA sequences that was used to predict DNA binding elements.
#' @param filterRange A length one genomic ranges in the form of a properly formatted character string (eg. "chr2:56-125") or an atomic GenomicRanges object. This argument specify the DNA region that will be plotted together with predicted binding TALEs RVD sequences whose predicted EBE lies \strong{entirely whithin}.
#'
#' @return Returns a ggplot object that can be further altered using ggplot2 package functions.
#' @export
plotTaleTargetPred <- function(predResults, subjDnaSeqFile, filterRange) {
  ######### Check and parse arguments
  subjDnaSeqs <- Biostrings::readDNAStringSet(subjDnaSeqFile)
  predsGr <- predResults %>% #dplyr::filter(strand == "-") %>%
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "subjSeqId", keep.extra.columns = TRUE)
  if (!all(as.character(BSgenome::getSeq(subjDnaSeqs, predsGr)) == predsGr$ebeSeq)) stop(
    "EBE sequences in the target predictions table did not match those extracted from the subjDnaSeqs!\n",
    "Verify that the content of the objects supplied as parameters are consistent.")

  if (length(filterRange) != 1L) stop(
    "The range used for filtering the displayed region must be of lenght one."
  )
  filterRange <- as(filterRange, "GRanges")


  ############ Prepare subjDnaSeqs for plotting
  yposSenseStrd <- 0.2
  yposAntisenseStrd <- -0.2

  tidySubjSeqs <- lapply(1:length(subjDnaSeqs),
                         function(i) {
                           oneSeq <- subjDnaSeqs[i]
                           dplyr::bind_rows(`+` = tidy_biostrings_msa(oneSeq),
                                            `-` = tidy_biostrings_msa(Biostrings::complement(oneSeq)),
                                            .id = "strand") %>%
                             dplyr::as_tibble() %>%
                             dplyr::mutate(yPos = dplyr::if_else(strand == "+",
                                                                 yposSenseStrd,
                                                                 yposAntisenseStrd)
                                           ) %>%
                             dplyr::rename(seqnames = name, xPos = position)
                         }
  ) %>% dplyr::bind_rows()



  grSubjSeqs <- tidySubjSeqs %>%
    GenomicRanges::makeGRangesFromDataFrame(start.field = "xPos",
                                            end.field = "xPos",
                                            keep.extra.columns = TRUE)

  ######### FILTER WHAT WILL BE DISPLAYED
  if (!is.null(filterRange)) {
    filteredPreds <- IRanges::subsetByOverlaps(predsGr, filterRange, type = "within") %>%
      as.data.frame(optional = TRUE, stringsAsFactors = FALSE) %>%
      dplyr::as_tibble(.name_repair = "minimal") %>%
      dplyr::rename(subjSeqId = seqnames) %>%
      dplyr::mutate(width = NULL) #%>% print(n = Inf)

    relevantTidySubjSeqs <- grSubjSeqs %>%
      IRanges::subsetByOverlaps(filterRange, type = "within", ignore.strand = TRUE) %>%
      as.data.frame(optional = TRUE, stringsAsFactors = FALSE) %>%
      dplyr::as_tibble(.name_repair = "minimal") %>%
      dplyr::rename(subjSeqId = seqnames, xPos = start) %>%
      dplyr::mutate(width = NULL, start = NULL, end = NULL, strand = NULL) #%>% print(n = Inf)
  } else {
    filteredPreds <- predsGr
    relevantTidySubjSeqs <- tidySubjSeqs
  }

  ############ Prepare preds for plotting
  # Append RVD OO
  predsForPlot <- filteredPreds %>% dplyr::mutate(rvds = paste("OO", rvds, sep = "-"))
  # Create yPos
  predsForPlot %<>% dplyr::group_by(subjSeqId, strand) %>%
    dplyr::mutate(
      yPos = dplyr::if_else(strand == "+",
                            rank(start, ties.method = "random"),
                            -1L*(rank(start, ties.method = "random"))
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(subjSeqId, yPos) #%>% print(n = Inf)

  # Create rvd and xPos
  predsForPlot %<>% dplyr::group_by_all() %>%
    dplyr::group_modify( ~{
      tibble::tibble(
        rvd = if(.y$strand == "+") {unlist(stringr::str_split(.y$rvds, pattern = "-"))}
        else {sapply(unlist(stringr::str_split(.y$rvds, pattern = "-")), rev)},
        xPos = if(.y$strand == "+") .y$start:.y$end else .y$end:.y$start,
        rvd2ntMatchScore = computeRVDSeqEBESeqMatchQualityString(RVDSeq = .y$rvds, EBESeq = .y$ebeSeq)
      )
    }) %>%
    dplyr::mutate(rvd = sapply(stringr::str_split(string = rvd, pattern = ""), paste, collapse = "\n"))

  # Create a derived tibble to add prediction scores and EBE box to the plot
  predsForScoreAndEbe <- predsForPlot %>% dplyr::group_by(subjSeqId, taleId, ebeSeq, score, strand, yPos) %>%
    dplyr::summarise(n = dplyr::n(), start = min(start), end = max(end)) %>%
    dplyr::mutate(
      yPosOnSeq = dplyr::if_else(strand == "+",
                                 yposSenseStrd,
                                 yposAntisenseStrd)
      )



  ############ Assembling a ggplot object
  # COLOR SCALE FOR NT
  ntColScale <- biovizBase::getBioColor("DNA_BASES_N")
  #ntColScale %>% pals::pal.safe()
  p <- ggplot2::ggplot(data = predsForPlot, ggplot2::aes(x = xPos, y = yPos)) +

    ggplot2::geom_tile (data = predsForScoreAndEbe, #<--------- EBE highlight box
                        mapping = ggplot2::aes(x = (start + (end - start) / 2),
                                               y = yPosOnSeq,
                                               width = (end - start +1)
                        ),
                        height = 0.4,
                        colour = "grey30",
                        size = 0.2,
                        fill = "khaki3", #  "yellowgreen",
                        alpha = 0.2,
    ) +

    ggplot2::geom_text(data = relevantTidySubjSeqs, #<--------- double stranded DNA sequence
              mapping = ggplot2::aes(label = character, colour = character),
              size = 3,
              fontface = "bold",
              show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(name = NULL, values = ntColScale) +
    
    ggplot2::geom_label(ggplot2::aes(label = rvd, fill = factor(rvd2ntMatchScore)), #<--------- RVDs
              size = 2.8,
              fontface = "bold",
              color = "white",
              label.padding = ggplot2::unit(0.1, "lines"),
              label.size = 0) +
    ggplot2::scale_fill_viridis_d(name = "RVD-DNA base match score",
                         option = "magma",
                         end = 0.5,
                         alpha = 0.7,
                         direction = -1) +

    ggplot2::geom_label(data = predsForScoreAndEbe, #<--------- Score
                        mapping = ggplot2::aes(label = sprintf("%05.2f",score)),
                        x = min(relevantTidySubjSeqs$xPos)+ 1,
                        fill = "white",
                        size = 2.8,
                        fontface = "bold",
                        color = "black",
                        label.padding = ggplot2::unit(0.1, "lines"),
                        label.size = 0.2) +

    ggplot2::scale_x_continuous( #<--------- x axis
      name = paste0("Position on subject DNA sequence: '", unique(as.character(filteredPreds$subjSeqId)), "'"),
      na.value = 0,
      breaks = function(r) seq(round(floor(r[1]), -1), r[2], by = 10),
      minor_breaks = function(r) seq(round(floor(r[1]), -1), r[2], by = 1),
      expand = ggplot2::expansion(mult = 0, add = c(1,1)),
      position = "top") +

    ggplot2::scale_y_continuous( #<--------- y axis
      name = NULL,
      breaks = unique(c(predsForPlot$yPos, relevantTidySubjSeqs$yPos)),
      minor_breaks = NULL,
      labels = function(brks) {
        sapply(brks, function(brk){
          if (brk == yposSenseStrd) return("5'")
          if (brk == yposAntisenseStrd) return("3'")
          return(predsForPlot$taleId[match(brk, predsForPlot$yPos)])
        }, simplify = TRUE)
      },
      expand = ggplot2::expansion(add = 0.5)
    ) +

    ggplot2::theme_linedraw() + #<--------- theme parameters
    ggplot2::theme(legend.position = "bottom",
          text = ggplot2::element_text(face = "bold",size = 11),
          panel.border = ggplot2::element_blank())

  return(p)
}




##' Convert msa file/object to tidy data frame.
##'
##'
##' @title tidy_msa
##' @param msa multiple sequence alignment file or sequence object in 
##' DNAStringSet, RNAStringSet, AAStringSet, BStringSet, DNAMultipleAlignment, 
##' RNAMultipleAlignment, AAMultipleAlignment, DNAbin or AAbin
##' @param start start position to extract subset of alignment
##' @param end end position to extract subset of alignemnt
##' @author Modified from Guangchuang Yu
##' @noRd
tidy_biostrings_msa <- function(msa, start = NULL, end = NULL) {
  aln <- msa
  alnmat <- lapply(seq_along(aln), function(i) {
    ##Preventing function collisions
    base::strsplit(as.character(aln[[i]]), '')[[1]]
  }) %>% do.call('rbind', .)
  ## for DNAbin and AAbin
  alndf <- as.data.frame(alnmat, stringsAsFactors = FALSE)
  
  if(unique(names(aln)) %>% length == length(aln)) {
    alndf$name = names(aln)
  }else{
    stop("Sequences must have unique names")
  }
  cn = colnames(alndf)
  cn <- cn[!cn %in% "name"]
  df <- gather(alndf, "position", "character", cn)
  
  y <- df
  y$position = as.numeric(sub("V", "", y$position))
  y$character = toupper(y$character)
  
  y$name = factor(y$name, levels=rev(names(aln)))
  
  
  if (is.null(start)) start <- min(y$position)
  if (is.null(end)) end <- max(y$position)
  
  y <- y[y$position >=start & y$position <= end, ]
  
  return(y)
}









