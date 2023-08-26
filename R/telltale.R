
# subjectFile = system.file("extdata", "bai3_sample_tal_regions.fasta", package = "tantale", mustWork = T)
# outputDir = tempdir(check = TRUE)
# hmmFilesDir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T)
# hmmerpath = system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = T)
# talArrayCorrection = TRUE
# refForTalArrayCorrection = system.file("extdata", "decipher_ref_tales_aa.fa.gz", package = "tantale", mustWork = T)
# frameShiftCorrection = -11
# TALE_NtermDNAHitMinScore = 300
# repeatDNAHitMinScore = 20
# TALE_CtermDNAHitMinScore = 200
# minDomainHitsPerSubjSeq = 4
# mergeHits = TRUE
# minGapWidth = 35
# taleArrayStartAnchorCode = "NTERM"
# taleArrayEndAnchorCode = "CTERM"
# appendExtremityCodes = TRUE
# rvdSep = "-"
# extendedLength = 300
# ... = NULL

# subjectFile = subjectFile
# outputDir = file.path("/home/cunnac/TEMP", gsub("(\\.fasta)|(\\.fa)|(\\.fna)|(\\.fsa)", "", basename(subjectFile)))
# hmmFilesDir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T)
# hmmerpath = system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = T)
# talArrayCorrection = FALSE
# refForTalArrayCorrection = system.file("extdata", "decipher_ref_tales_aa.fa.gz", package = "tantale", mustWork = T)
# frameShiftCorrection = -11
# TALE_NtermDNAHitMinScore = 300
# repeatDNAHitMinScore = 20
# TALE_CtermDNAHitMinScore = 200
# minDomainHitsPerSubjSeq = 4
# mergeHits = TRUE
# minGapWidth = 35
# taleArrayStartAnchorCode = "NTERM"
# taleArrayEndAnchorCode = "CTERM"
# appendExtremityCodes = TRUE
# rvdSep = "-"
# extendedLength = 300
# ... = NULL







#' Search and report on the features of TALE protein domains potentially encoded
#' in subject DNA sequences
#'
#' \code{tellTale} has been primarily written to report on 'corrected' TALE RVD
#' sequences in indels prone, noisy DNA sequences (suboptimally polished genomes
#' assembly, raw reads of long read sequencing technologies [eg PacBio, ONT])
#' that would otherwise be missed by conventional tools (eg AnnoTALE).
#'
#' The approach is first to use \href{http://hmmer.org/}{HMMER} to find and
#' categorize regions in the input DNA sequence that are related to the coding
#' sequence of canonical TALE protein domains (N-Term, repeats, C-term). Hits
#' that are (nearly [see the minGapWidth parameter]) adjacent are grouped in
#' "taleArrays" which are considered as potential tal genes.
#'
#'
#' If the \code{talArrayCorrection} parameter is turned off, the longest
#' predicted open reading frame (+extendedLength) for each talArray is fed to
#' \href{http://www.jstacs.de/index.php/AnnoTALE}{AnnoTALE} to detect TALE
#' domains in the predicted translation product. The Results should hence be
#' very similar to what would be obtained with AnnoTALE, plus many additional
#' informative output files such as tabular reports.
#'
#' If \code{talArrayCorrection} is turned on, these talearrays are passed to the
#' \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}} function that
#' attemps to 'correct' potential frameshifts in the taleArray sequences. This
#' conveniently removes many artefactual indels but bear in mind that this may
#' also \strong{erroneously} 'correct' genuine frame shifts which can be highly
#' relevant especially for truncTALEs or iTALES. The resulting 'corrected'
#' taleArray open reading frames are then passed to AnnoTALE.
#'
#' Note that occasionally, when a putative open reading frame does not encode a
#' canonical TALE protein (early frame shift, incomplete ORF, etc...), the
#' "analyze" module of AnnoTALE outputs DNA parts but no protein parts and/or
#' RVD sequence. This should be detected and reported in the tellTale log.
#'
#'
#' @param subjectFile Fasta file with DNA sequence(s) to be searched for the
#'   presence of TALE coding sequences (CDS).
#' @param outputDir Path of the output directory. If not specified, results will
#'   be written to current working folder.
#' @param hmmFilesDir Specify the path to a folder holding the hmmfiles if you
#'   do not want to use the ones provided with tantale.
#' @param TALE_NtermDNAHitMinScore Minimal nhmmer score cut_off value to
#'   consider the hit as genuine
#' @param repeatDNAHitMinScore Minimal nhmmer score cut_off value to consider
#'   the hit as genuine
#' @param TALE_CtermDNAHitMinScore Minimal nhmmer score cut_off value to
#'   consider the hit as genuine
#' @param minDomainHitsPerSubjSeq Minimum number of nhmmer hits for a subject
#'   sequence to be reported as having TALE diagnostic regions. This is a way to
#'   simplify output a little by getting ride of uninformative sequences
#' @param mergeHits Perform overlapping hits merging per domain type. Should not
#'   be modified.
#' @param minGapWidth Minimum gap in base pairs between two tale domain hits for
#'   them to be considered distinct. If the length of the gap is below this
#'   value, domains are considered "contiguous" and grouped in the same array.
#' @param appendExtremityCodes Set this to \code{FALSE} if you do not want the
#'   N- and C-TREM anchor codes in the output sequences of RVD
#' @param rvdSep Symbol acting as a separator in RVD sequences
#' @param hmmerpath Specify the path to a directory holding the HMMER executable
#'   if you do not want to use the ones provided with tantale.
#' @param extendedLength number of nucleotides to extend in 3'-end at the tal
#'   ORF prediction stage.
#' @param talArrayCorrection True or False
#' @param refForTalArrayCorrection Reference AA sequences for tal array
#'   predicted ORF correction if you do not want to use the ones provided with
#'   tantale.
#' @param frameShiftCorrection This is an internal parameter of the
#'   \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}} function. The
#'   default is 11 and fiddle with this at your own risk...
#' @param ... Additional parameters for the
#'   \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}} function.
#' @return This functions has only side effects (writing files, mostly).
#'   However, if everything ran smoothly, it will invisibly return the path of
#'   the directory where output files were written.
#'
#'
#'   List of output files:
#'   \itemize{
#'   \item allRanges.gff: gff file of all Tal arrays detected by HMMer
#'   \item arrayReport.tsv: report of all Tal arrays. In the arrayReport.tsv,
#'   column \emph{predicted_dels_count}/\emph{predicted_ins_count} shows the
#'   number of putative deletions/insertions in the raw sequences that have been
#'   corrected in the corrected sequences with the
#'   function \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}}.
#'   \item hitsReport.tsv: report of all hits detected by HMMer
#'   \item hitsReport.gff: gff file of all hits detected by HMMer
#'   \item domainsReport.tsv: report of all Tal amino acid domains detected by AnnoTALE analyze
#'   \item putativeTalOrf.fasta: Tal putative ORFs
#'   \item pseudoTalCds.fasta: pseudo Tal CDS, putative Tal array ORFs detected by HMMer for whch
#'    AnnoTALE analyze failed to find RVD(s).
#'   \item rvdSequences.fas: Sequence of RVDs (separated by rvdSep) predicted to be encoded in the Tal array
#'    ORFs by AnnoTALE. Note that if appendExtremityCodes is \code{TRUE} (by default),
#'    the N- and C-TREM anchor codes will be appended at the beginning and end of the sequences
#'    if the corresponding domain coding sequence was wound by HMMer at the DNA level.
#'    If no such HMMer hits were found, the "XXXXX" string will be appended
#'    to denote that AA sequences outside of the RVD array are likely to be atypical.
#'   \item C-terminusAAAlignment.html: protein alignment of all C-termini
#'   \item C-terminusDNAAlignment.html: DNA alignment of all C-termini
#'   \item N-terminusAAAlignment.html: protein alignment of all N-termini
#'   \item N-terminusDNAAlignment.html: DNA alignment of all N-termini
#'   \item TALE_CDS_all_diagnostic_regions_hmmfile.out: HMMER profile used for tale cds search.
#'   \item hmmerSearchOut.txt: ignore
#'   \item nhmmerHumanReadableOutputOfLastRun.txt: primary HMMER output file.
#'   \item tellTale.log: a log file
#'   \item annotale folder: folder containing result of AnnoTALE analyze for all Tal arrays
#'   \item CorrectionAlignmentAA folder: folder containing protein alignment of Tal array detected by HMMer and corrected Tal array if \code{talArrayCorrection} = TRUE
#'   \item CorrectionAlignmentDNA folder: folder containing DNA alignment of Tal array detected by HMMer and corrected Tal array if \code{talArrayCorrection = TRUE}
#'   }
#' @export
tellTale <- function(
  subjectFile,
  outputDir = getwd(),
  hmmFilesDir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T),
  TALE_NtermDNAHitMinScore = 300,
  repeatDNAHitMinScore = 20,
  TALE_CtermDNAHitMinScore = 200,
  minDomainHitsPerSubjSeq = 4,
  mergeHits = TRUE,
  minGapWidth = 35,
  appendExtremityCodes = TRUE,
  rvdSep = "-",
  hmmerpath = system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = T),
  extendedLength = 300,
  talArrayCorrection = FALSE,
  refForTalArrayCorrection = system.file("extdata", "decipher_ref_tales_aa.fa.gz", package = "tantale", mustWork = T),
  frameShiftCorrection = -11,
  ...
) {

  # @param taleArrayStartAnchorCode This scalar character vector will symbolize a
  #   TALE N-TERM CDS hit in the RVD sequence
  # @param taleArrayEndAnchorCode This scalar character vector will symbolize a
  #   TALE C-TERM CDS hit in the RVD sequence
  taleArrayStartAnchorCode <- "NTERM"
  taleArrayEndAnchorCode <- "CTERM"
  taleArrayAtypicalExtremityCode <- "XXXXX"

  #### TODO ####
  # Add an ooptional argument that olds the circularity status of molecules in genome
  # update the seqinfo objects accrodingly. This may solve some issues if a tale is located
  # at the junction of extremities in a circular molecule.
  # Could also implement an autotmated mecanisms "findCircular" that would
  # find a temr (eg 'circular') in the sequence title and act accrodingly.
  
  
  ####   Paths of output files   ####
  dir.create(outputDir, recursive = T, mode = "755", showWarnings = FALSE)
  ## Path of the directories where DECIPHER alignments will be written
  alignmentDNADir <- file.path(outputDir, "CorrectionAlignmentDNA")
  dir.create(alignmentDNADir, showWarnings = F)
  alignmentAADir <- file.path(outputDir, "CorrectionAlignmentAA")
  dir.create(alignmentAADir, showWarnings = F)
  # fasta of tals orfs that have rvds
  putatieOrfOfTaleWithRvdFile <- file.path(outputDir, "putativeTalOrf.fasta")
  # fasta of tals orfs that were not predicted to contain rvds
  pseudoTalFile <- file.path(outputDir, "pseudoTalCds.fasta")
  # annotale output directory
  annotaleMainDir <- file.path(outputDir, "annotale")# tempfile(pattern = "annotale_", tmpdir = outputDir)
  dir.create(annotaleMainDir)
  ## Tabular file reporting on individual TALE domain hits
  hitsReportFile <- file.path(outputDir, "hitsReport.tsv")
  domainsReportFile <- file.path(outputDir, "domainsReport.tsv")
  ## Tabular file reporting on putative TALEs (contiguous arrays of domain hits)
  arrayReportFile <- file.path(outputDir, "arrayReport.tsv")
  ## Gff file with all the identified domains and arrays and their associated data
  affRangesGffFile <- file.path(outputDir, "allRanges.gff")
  ## A fasta file of the selected seq of RVDs without the - separator
  seqsOfRVDFile <- file.path(outputDir, "rvdSequences.fas")
  ## A fasta file with array ORFs DNA sequences
  #arrayOrfsSeqFile <- file.path(outputDir, "arrayOrfs.fas")
  ## A text file where logging info and some general analysis measures are written
  analysisLogFile <- file.path(outputDir, "tellTale.log")
  
  
  ####   Checks for parameters and other things   ####

  ## Deal with spaces in sequence names because this messes up parsing of HMMER output
  logger::log_info("HMMER is very picky about forbiden characters in sequence name. Renaming sequences in {subjectFile}.")
  originalSeqs <- Biostrings::readDNAStringSet(filepath = subjectFile)
  Rsamtools::indexFa(subjectFile)
  originalSeqInfo <- Rsamtools::seqinfo(Rsamtools::FaFile(subjectFile))
  originalSeqlevels <- names(originalSeqs)
  foolproofSeqlevels <- paste0("seq", 1:length(originalSeqlevels))
  names(originalSeqlevels) <- foolproofSeqlevels
  names(originalSeqs) <- foolproofSeqlevels
  logger::log_info("Original seq names : {glue::glue_collapse(originalSeqlevels, sep = ' ; ')}.")
  logger::log_info("Dummy seq names : {glue::glue_collapse(names(originalSeqlevels), sep = ' ; ')}.")
  subjectFile <- tempfile()
  Biostrings::writeXStringSet(originalSeqs, filepath = subjectFile)
  
  ####   Full paths of input HMM files for TALE domains (DNA and AA)   ####
  
  ## TODO come up with a mechanism for the user to be able to provide the FULL PATH
  ## to custom hmm !!! hmmFilesDir parameter is useless unless custom hmm are named
  ## as specified below
  TALE_NtermDNAHMMFile <- file.path(hmmFilesDir, "Xo_TALE_Nterm_CDS_profile.hmm")
  repeatDNAHMMFile <- file.path(hmmFilesDir, "Xo_TALE_repeat_CDS_profile.hmm")
  TALE_CtermDNAHMMFile <- file.path(hmmFilesDir, "Xo_TALE_Cterm_CDS_profile.hmm")
  repeatAAHMMFile <- file.path(hmmFilesDir, "Xo_TALE_repeat_AA_profile.hmm")
  
  DNAHMMFiles <- c(TALE_NtermDNAHMMFile, repeatDNAHMMFile, TALE_CtermDNAHMMFile)
  mergedDNAHMMFile <- file.path(outputDir, "TALE_CDS_all_diagnostic_regions_hmmfile.out")
  
  
  
  ####   Concatenate HMM files for TALE DNA motifs and parse HMM names   ####
  hmmslines <- plyr::llply(DNAHMMFiles, function(x) txt <- readLines(con = x))
  names(DNAHMMFiles) <- plyr::llply(hmmslines, function(x) {
    hmmName <- grep("NAME", x, perl = TRUE, value = TRUE)
    hmmName <- unlist(strsplit(hmmName, split = "\\s+"))
    if (length(hmmName) != 2) stop("One or several profile ",
                                   "HMM have a name with spaces. ",
                                   "Please remove them in the file at the Tag 'NAME'")
    hmmName <- hmmName[2]
  }
  )
  DNAHMMNames <- names(DNAHMMFiles)
  TALE_NtermDNAHMMName <- names(DNAHMMFiles)[1]
  repeatDNAHMMName <- names(DNAHMMFiles)[2]
  TALE_CtermDNAHMMName <- names(DNAHMMFiles)[3]
  
  writeLines(text = unlist(hmmslines), con = mergedDNAHMMFile)
  
  ####   Perform TALE domain CDS search with HMMER  #####
  searchOutFile <- file.path(outputDir, "hmmerSearchOut.txt")
  
  runNhmmerSearch(hmmerpath = hmmerpath,
                  subjectFile = subjectFile,
                  hmmFile = mergedDNAHMMFile,
                  searchTblOutFile = searchOutFile,
                  humReadableOutFile = file.path(outputDir, "nhmmerHumanReadableOutputOfLastRun.txt"))
  
  ####   Load, process, filter TALE domain CDS HMMER hit results    ####
  ## Loading search tabular output file
  nhmmerTabularOutput <- try(read.table(searchOutFile), silent = TRUE)
  if (class(nhmmerTabularOutput) == "try-error") {
    warning("NhmmerSearch found no TALE cds hit in ", subjectFile , " Exitting...")
    return(invisible(outputDir))
  }

  colnames(nhmmerTabularOutput) <- c("target_name", "accession", "query_name", "accession", "hmmfrom", "hmm_to", "alifrom",
                                     "ali_to", "envfrom", "env_to", "sq_len", "strand", "Evalue", "score", "bias", "description_of_target")
  
  ## filtering results differentially depending on the query HMM
  nhmmerTabularOutput <- subset(nhmmerTabularOutput,
                                query_name == TALE_NtermDNAHMMName & score >= TALE_NtermDNAHitMinScore |
                                  query_name == repeatDNAHMMName & score >= repeatDNAHitMinScore |
                                  query_name == TALE_CtermDNAHMMName & score >= TALE_CtermDNAHitMinScore
  )
  nhmmerTabularOutput <- droplevels(nhmmerTabularOutput)
  if(nrow(nhmmerTabularOutput) == 0L) {
    logger::log_warn("No record remains after filtering NhmmerSearch hits based on score. Exitting...")
    return(invisible(outputDir))
  }
  ## Add a hitID column
  nhmmerTabularOutput$hitID <- paste("DOM", sprintf("%05.0f", 1:nrow(nhmmerTabularOutput)), sep="_")
  
  ## Trick to re-order positions in an increasing order to satisfy IRanges() in preparation of creating a GRanges
  nhmmerTabularOutput[,c("start", "end")] <- plyr::adply(.data = nhmmerTabularOutput[,c("envfrom", "env_to")],
                                                         .margins = 1, .fun = c(min, max))[,-(1:2)]
  rownames(nhmmerTabularOutput) <- nhmmerTabularOutput$hitID
  
  #####   Storing all info about individual TALE domains in a GenomicRanges object   ####
  
  ## Filter out target DNA sequences that have too few repeat CDSs
  ## NB: for the sake of consistency  it would be better just to filter out
  ## from any further consideration the ARRAYS shorter than a certain value (say 5).
  ## WHAT DO WE DO ABOUT THAT?
  temp_df <- plyr::ddply(nhmmerTabularOutput[,-20], ~ target_name + sq_len, nrow) # I do not know why but it fails to work if I leave the RVD column (#20)
  nhmmerTabularOutput <- subset(nhmmerTabularOutput, target_name %in% temp_df[temp_df$V1 > minDomainHitsPerSubjSeq, "target_name"])
  nhmmerTabularOutput <- droplevels(nhmmerTabularOutput)
  

  ## Creating a GRanges object from nhmmerOutput
  nhmmerOutputGR <- GenomicRanges::makeGRangesFromDataFrame(
    df = nhmmerTabularOutput, keep.extra.columns = TRUE,
    seqnames.field = "target_name"
  )
  names(nhmmerOutputGR) <- nhmmerOutputGR$hitID
  
  ## Updating seqinfo with original seqinfo from the sequences before renaming
  nhmmerOutputGR <- GenomeInfoDb::renameSeqlevels(nhmmerOutputGR, value = originalSeqlevels)
  GenomeInfoDb::seqinfo(nhmmerOutputGR, pruning.mode = "coarse") <- originalSeqInfo[GenomeInfoDb::seqlevels(nhmmerOutputGR)]
  
  nhmmerOutputGRBeforeMerge <- nhmmerOutputGR
  
  #####   Domain-wise merge of overlapping hits  #####
  if (mergeHits) {
    ## Split the ranges by domain type
    nhmrOutGrByDomain <- GenomicRanges::split(nhmmerOutputGR, f = nhmmerOutputGR$query_name)
    ## Perform merge
    reducedOlapGr <- lapply(nhmrOutGrByDomain, function(gr) {
      reducedOlapGr <- as.data.frame(GenomicRanges::findOverlaps(gr, gr,
                                                                 minoverlap = 2,
                                                                 type = "any",
                                                                 ignore.strand=FALSE,
                                                                 select = "all")) %>%
        dplyr::group_by(queryHits) %>%
        dplyr::group_map({
          ~ GenomicRanges::reduce(gr[as.numeric(.x$subjectHits)], with.revmap = FALSE)
        }) %>%
        plyranges::bind_ranges() %>%
        unique()
      formerIDs <- as.data.frame(GenomicRanges::findOverlaps(gr, reducedOlapGr,
                                                             minoverlap = 2,
                                                             type = "within",
                                                             ignore.strand=FALSE,
                                                             select = "all")) %>%
        dplyr::group_by(subjectHits) %>%
        dplyr::group_map({
          ~ paste0(gr[as.numeric(.x$queryHits)]$hitID, collapse = "|")
        }) %>%
        unlist()
      reducedOlapGr$nhmmerHitID <- formerIDs
      return(reducedOlapGr)
    }) %>%
      plyranges::bind_ranges(.id = "query_name")
    reducedOlapGr$hitID <- paste("MDOM", sprintf("%05.0f", 1:length(reducedOlapGr)), sep="_")
    names(reducedOlapGr) <- reducedOlapGr$hitID
    
    ## From there on, use the merged hits GRanges
    nhmmerOutputGR <- reducedOlapGr
  } else {
    nhmmerOutputGR <- nhmmerOutputGRBeforeMerge
  }
  
  
  #####   Extract the DNA sequences of the hits and record in GRanges mcols   #####
  ## Load in R the DNA sequences that are queried for TALE CDS
  subjectDNASequences <- Biostrings::readDNAStringSet(filepath = subjectFile)
  names(subjectDNASequences) <- originalSeqlevels[match(names(subjectDNASequences), names(originalSeqlevels))]
  
  ## Extract the DNA sequences of the hits
  hitSeqs <- BSgenome::getSeq(subjectDNASequences, nhmmerOutputGR)
  S4Vectors::mcols(nhmmerOutputGR) %<>% cbind(
    data.frame(
      "seq" = as.character(hitSeqs),
      "codon_count" = Biostrings::nchar(hitSeqs) %/% 3,
      "frameshift_count" = Biostrings::nchar(hitSeqs) %% 3,
      row.names = NULL,
      check.rows = T)
  )
  
  
  
  
  #####   Group (nearly) adjacent hits in "TALE array" regions   #####
  
  ## Group "contiguous" hits (repeats or other regions) in a GRangesList
  ## Use reduce to obtain the list of regions (arrays of domains for the time being) that span "contiguous" hits
  ## Here contiguous is defined as hits that are less than minGapWidth bp appart
  arraysGR <- GenomicRanges::reduce(nhmmerOutputGR,
                                    drop.empty.ranges=FALSE,
                                    min.gapwidth= minGapWidth,
                                    with.revmap=TRUE,
                                    ignore.strand=FALSE)
  revmap <- S4Vectors::mcols(arraysGR)$revmap  # an IntegerList
  
  
  ## Use the mapping from reduced to original ranges to group the original ranges by reduced range:
  hitsByArraysLst <- BiocGenerics::relist(nhmmerOutputGR[unlist(revmap)], revmap)
  names(hitsByArraysLst) <- paste("ROI", sprintf("%05.0f", 1:length(hitsByArraysLst)), sep = "_")
  names(arraysGR) <- names(hitsByArraysLst)
  # hitsByArraysLst <- BiocGenerics::sort(hitsByArraysLst) # just to make sure...
  
  ## Make sure that hits do not overlap for some weird reason
  doHitsOverlap <- !GenomicRanges::isDisjoint(hitsByArraysLst)
  if (any(doHitsOverlap)) {
    warning("It appears that some hmmer hits actually overlap.\n It is thus possible that the inferred sequences of RVDs have artefactual insertions.\n")
    warning(paste0("Please check the hits in the following RegionsOfInterest:", "\n",
                   paste(names(doHitsOverlap)[doHitsOverlap], collapse = "\n"), "\n")
    )
  }
  
  
  #####   Storing all info about domain arrays in a central object   ####
  ## Populate metadata about the elements of the list of arrays
  
  
  S4Vectors::mcols(hitsByArraysLst) <- S4Vectors::DataFrame(
    arrayID = names(hitsByArraysLst),
    OriginalSubjectName = sapply(hitsByArraysLst,
                                 function(x) unique(as.character(GenomicRanges::seqnames(x)))),
    Start = BiocGenerics::start(arraysGR),
    End = BiocGenerics::end(arraysGR),
    Strand = BiocGenerics::strand(arraysGR),
    NumberOfHits = S4Vectors::elementNROWS(hitsByArraysLst),
    ArraySeq = BSgenome::getSeq(subjectDNASequences, arraysGR),
    AllDomains = sapply(hitsByArraysLst,
                        function(x) {
                          all(
                            c(TALE_NtermDNAHMMName, repeatDNAHMMName,
                              TALE_CtermDNAHMMName) %in% as.character(x$query_name)
                          )
                        }
    )
  )
  

  #####   Extend DNA Tal arrays   #####
  ## Extract the genomic sequence of arrays +-bp on the borders
  completeArraysGR <- arraysGR
  extdCompleteArraysGR <- GenomicRanges::resize(completeArraysGR,
                                                width = GenomicRanges::width(completeArraysGR) + extendedLength,
                                                fix = "start", ignore.strand = FALSE) %>%
    GenomicRanges::trim(use.names = TRUE)
  extdCompleteArraysSeqs <- BSgenome::getSeq(subjectDNASequences, extdCompleteArraysGR)
  
  

  if (!talArrayCorrection) {
    #### Get ORFs from uncorrected Tal arrays if frame shifts correction is OFF ####
    orfs <- systemPipeR::predORF(x = extdCompleteArraysSeqs,
                                 n = 1, type = "gr", mode = "ORF", strand = "sense")
    fullTalOrf <- BSgenome::getSeq(extdCompleteArraysSeqs, orfs)
    names(fullTalOrf) <- as.character(GenomicRanges::seqnames(orfs))
    TalOrfForAnnoTALE <- fullTalOrf
  } else { 
    #### Correct Tal arrays frame shifts if requested  ####
    ####   Run CorrectFrameshifts   ####
    logger::log_info("Correcting putative TALE coding sequences. Be patient, this may take a LONG time...")
    AAref <- Biostrings::readAAStringSet(refForTalArrayCorrection, seek.first.rec = TRUE, use.names = TRUE)
    rawArraySeq <- extdCompleteArraysSeqs
    ArrayCorrection <- DECIPHER::CorrectFrameshifts(rawArraySeq,
                                                    AAref, type = "both",
                                                    maxComparisons = length(AAref),
                                                    frameShift = frameShiftCorrection, ...)
    logger::log_info("Correction of putative TALE coding sequences is done!")
    corrExtdCompleteArraysSeqs <- ArrayCorrection$sequences
    
    ####   Correction stats   ####
    deletions_count <- correction_tible(ArrayCorrection$indels) %>%
      dplyr::group_by(Seq) %>%
      dplyr::count(variable, name = "predicted_dels_count") %>%
      dplyr::filter(variable == "deletions") %>%
      dplyr::select(-variable)
    insertions_count <- correction_tible(ArrayCorrection$indels) %>%
      dplyr::group_by(Seq) %>%
      dplyr::count(variable, name = "predicted_ins_count") %>%
      dplyr::filter(variable == "insertions") %>%
      dplyr::select(-variable)
    
    
    S4Vectors::mcols(hitsByArraysLst) <- merge(S4Vectors::mcols(hitsByArraysLst), 
                                               dplyr::full_join(insertions_count, deletions_count, by = "Seq"), 
                                               by.x = "arrayID", by.y = "Seq", all.x = T) 
    
    S4Vectors::mcols(hitsByArraysLst)[c("predicted_dels_count", "predicted_ins_count")] %<>% apply(., 2, function(v) ifelse(is.na(v), 0, v))
    
    #### Multiple alignenment of orginal vs corrected vs corrected+N/C subtituted sequences  ####
    
    for (n in names(corrExtdCompleteArraysSeqs)[vcountPattern("N", corrExtdCompleteArraysSeqs) > 0]) {
      logger::log_warn("After correction, {n} sequence contains 'N's which will be substituted by 'C's in order",
                       "to run AnnoTALE analyze for RVDs prediction.")
    }
    substCorrExtdCompleteArraysSeqs <- Biostrings::chartr("N", "C", corrExtdCompleteArraysSeqs)
    

    
    for (n in names(corrExtdCompleteArraysSeqs)) {
      rawSeq <- rawArraySeq[n]
      names(rawSeq) <- paste0("raw_", n)
      correctedSeq <- corrExtdCompleteArraysSeqs[n]
      names(correctedSeq) <- paste0("corrected_", n)
      substitutedSeq <- substCorrExtdCompleteArraysSeqs[n]
      names(substitutedSeq) <- paste0("forAnnoTALE", n)
      seqToAlign <- c(rawSeq, correctedSeq, substitutedSeq)
      alignedSeqs <- DECIPHER::AlignSeqs(seqToAlign)
      DECIPHER::BrowseSeqs(alignedSeqs,
                           htmlFile = file.path(alignmentDNADir, glue::glue("CorrectionAlignmentDNA_{n}.html")),
                           openURL = F, colWidth = 120)
      
      seqToAlignTranslated <- Biostrings::translate(seqToAlign, no.init.codon = T, if.fuzzy.codon = "solve")
      alignedSeqsTranslated <- DECIPHER::AlignSeqs(seqToAlignTranslated)
      DECIPHER::BrowseSeqs(alignedSeqsTranslated,
                           htmlFile = file.path(alignmentAADir, glue::glue("CorrectionAlignmentAA_{n}.html")),
                           openURL = F, colWidth = 120)
      
    }
    
    ####   TalOrfForAnnoTALE   ####
    orfs <- systemPipeR::predORF(x = substCorrExtdCompleteArraysSeqs,
                                 n = 1, type = "gr", mode = "ORF", strand = "sense")
    TalOrfForAnnoTALE <- BSgenome::getSeq(substCorrExtdCompleteArraysSeqs, orfs)
    names(TalOrfForAnnoTALE) <- as.character(GenomicRanges::seqnames(orfs))
    
    #### Get ORFs from corrected arrays sequences that have 'Ns' substituted by 'Cs' ####
    
    fullTalOrf <- TalOrfForAnnoTALE
    
    # The idea here was to convert back the Ns that were substituted for AnnoTALE in order to output
    # an unsubstituted orf.
    # I am not sure the code below properly does that and retrospectively,
    # it may be a problem for downstream bioinformatic analyses to have Ns in the sequence...
    #
    # insPosition <- vmatchPattern("N", corrExtdCompleteArraysSeqs)
    # for (i in names(insPosition)) {
    #   if (length(width(insPosition[[i]])) == 0) next()
    #   insPositionAfCorr <- insPosition[[i]][BiocGenerics::start(insPosition[[i]]) < width(fullTalOrf[i])]
    #   fullTalOrf[i] <- replaceAt(fullTalOrf[i], insPositionAfCorr, value = "N")
    # }
    
    
  }
  
  
  
  #### annoTALE analyze on tal ORFs  ####
  
  # Shall we also run the predict stage of annotale? May be it will do a better job at
  # finding orf and/or filtering out "pseudo tales" because some times analyse output a RVD from
  # a domain that does not look like a repeat....
  
  
  AnnoTALEanalyze <- function(inputFastaFile,
                              outputDir = getwd(),
                              prefix = NULL,
                              annoTALE = system.file("tools", "AnnoTALEcli-1.4.1.jar",
                                                     package = "tantale", mustWork = T)
                              ) {
    # Define output dirs for the various stages of annoTALE
    stopifnot(dir.exists(outputDir) || dir.create(path = outputDir, showWarnings = TRUE,
                                                  recursive = TRUE, mode = "775"))
    # Define a prefix for TALEs (assembly ID) derived from the genome file name.
    if (is.null(prefix)) {
      prefix <- gsub(pattern = "^(.*)\\.(fasta|fa|fas)$" ,
                     replacement  = "\\1", basename(inputFastaFile),
                     perl = TRUE)
    }
    # Run the "analyze" stage of annoTALE
    comAnalyze <- paste0(
      "java -jar ", annoTALE,
      " analyze ",
      " t=", inputFastaFile,
      " outdir=", outputDir
    )
    cat("##  Now running annoTALE analyze for", prefix, "using the following command:\n##  ",  comAnalyze, "\n")
    exitAnalyze <- system(comAnalyze, ignore.stdout = TRUE, ignore.stderr = TRUE)
    return(invisible(exitAnalyze))
  }
  
  ## run annotale for tal putative orfs
  
  annoTaleMessages <- character()
  
  ##!!! NOTE: correction could easily be parallelized in the "sapply" below...
  annoTaleOut <- sapply(names(TalOrfForAnnoTALE), function(talOrfID) {
    AnnotaleDir <- file.path(annotaleMainDir, talOrfID)
    dir.create(AnnotaleDir)
    TalOrf <- TalOrfForAnnoTALE[talOrfID]
    correctedTalOrfFile <- file.path(AnnotaleDir, "putativeTalOrf.fasta")
    Biostrings::writeXStringSet(TalOrf, correctedTalOrfFile)
    # Run Annotale on corrected ORF
    checkAnnoTale <- try(AnnoTALEanalyze(correctedTalOrfFile, AnnotaleDir), silent = TRUE)
    # Get Annotale output with conditions handling
    dna_parts_files <- file.path(AnnotaleDir, "TALE_DNA_parts.fasta")
    prot_parts_files <- file.path(AnnotaleDir, "TALE_Protein_parts.fasta")
    annoTaleRVD_file <- file.path(AnnotaleDir, "TALE_RVDs.fasta")
    seqOfRVDs <- try(Biostrings::readAAStringSet(annoTaleRVD_file,
                                                 seek.first.rec = T,
                                                 use.names = T),
                     silent = TRUE)
    prot_parts <- try(Biostrings::readAAStringSet(prot_parts_files), silent = TRUE)
    
    if (any(
      inherits(checkAnnoTale, "try-error"), # in case annotale does not work 
      if (inherits(prot_parts, "try-error")) { # in case annotale does not return a prot_parts file or if it is empty.
        file.exists(prot_parts_files) && file.remove(prot_parts_files)
        TRUE
      } else {
        if (file.exists(prot_parts_files) && length(prot_parts) == 0L) {
          file.remove(prot_parts_files) # should also return TRUE
        }
      },
      if (inherits(seqOfRVDs, "try-error")) { # in case annotale works but cannot find rvds or rvd seq file is empty.
        file.exists(annoTaleRVD_file) && file.remove(annoTaleRVD_file)
        TRUE 
      } else {
        if (Biostrings::width(seqOfRVDs) == 0) file.remove(annoTaleRVD_file) # should also return TRUE
      }
    )) {
      annoTaleMessages <<- c(annoTaleMessages,
                            (m <- glue::glue("Annotale failed to parse TALE domains for {talOrfID}.")))
      logger::log_warn(m)
      return(annout(Biostrings::AAStringSet(), domainsReport = data.frame()))
    }
    names(seqOfRVDs) <- talOrfID
    
    ## domains report
    stops <- Biostrings::vcountPattern("*", prot_parts)
    domainsReport <- tibble::tibble("arrayID" = talOrfID,
                                "seqnames" = S4Vectors::mcols(hitsByArraysLst)$OriginalSubjectName[S4Vectors::mcols(hitsByArraysLst)$arrayID == talOrfID],
                                "query_name" = gsub("(.+\\: )|( \\d+)", "", names(prot_parts)),
                                "codon_count" = width(prot_parts) - stops
                                )
    
    annotale_output <- annout(seqOfRVDs, domainsReport = domainsReport)
    return(annotale_output)
  }, simplify = T, USE.NAMES = F)
  
  seqsOfRVDs <- unlist(Biostrings::AAStringSetList(annoTaleOut))
  domainsReport <- do.call(rbind, lapply(annoTaleOut, function(x) x@domainsReport))
  
  
  #### TODO  #####
  # The exact content of the files below needs to be reassesed and 
  # we need to determine if this is really what we want.
  # save tals orfs that have rvds
  Biostrings::writeXStringSet(fullTalOrf[names(fullTalOrf) %in% names(seqsOfRVDs)], putatieOrfOfTaleWithRvdFile)
  
  # save tals that DO NOT have rvds
  Biostrings::writeXStringSet(extdCompleteArraysSeqs[!names(extdCompleteArraysSeqs) %in% names(seqsOfRVDs)], pseudoTalFile)
  
  aberrantRepeat <- sapply(seqsOfRVDs, function(s) {
    ifelse(length(s) > 0, grepl("[a-z]", s), NA)
  })
  
  seqsOfRVDs <- gsub("\\-", rvdSep, seqsOfRVDs) %>% Biostrings::AAStringSet()
  
  if (appendExtremityCodes) {
    # This is necessary for other tantale utilities that can operate on 'full' domains sequences, ie downstream of distal, for TALE  domains sequences alignments.
    for (s in names(seqsOfRVDs)) {
      hitsByArray <- hitsByArraysLst[[s]]
      seqsOfRVDs[s] <- paste(ifelse(TALE_NtermDNAHMMName %in% as.character(hitsByArray$query_name),
                                    taleArrayStartAnchorCode, taleArrayAtypicalExtremityCode), 
                             seqsOfRVDs[s], 
                             sep = rvdSep)
      seqsOfRVDs[s] <- paste(seqsOfRVDs[s], 
                             ifelse(TALE_CtermDNAHMMName %in% as.character(hitsByArray$query_name),
                                    taleArrayEndAnchorCode, taleArrayAtypicalExtremityCode), 
                             sep = rvdSep)
    }
  }
  
  S4Vectors::mcols(hitsByArraysLst) <- merge(S4Vectors::mcols(hitsByArraysLst),
                                             data.frame(SeqOfRVD = seqsOfRVDs,
                                                        aberrantRepeat = aberrantRepeat,
                                                        arrayID = names(seqsOfRVDs)
                                                        ),
                                             by = "arrayID", 
                                             all.x = T)
  S4Vectors::mcols(hitsByArraysLst)$SeqOfRVD[is.na(S4Vectors::mcols(hitsByArraysLst)$SeqOfRVD)] <- ""
  
  #### Align N-term and C-term ####
  dnaPartFiles <- list.files(annotaleMainDir, "TALE_DNA_parts.fasta", recursive = T, full.names = T)
  for (part in c("N-terminus", "C-terminus")) {
    allpart <- sapply(dnaPartFiles, function(p) {
      allpart <- Biostrings::readDNAStringSet(p, seek.first.rec = T)
      onepart <- allpart[grepl(part, names(allpart))]
      talRoi <- basename(dirname(p))
      names(onepart) <- talRoi
      return(onepart)
    }, simplify = "array", USE.NAMES = F) %>%
      Biostrings::DNAStringSetList() %>%
      unlist()
    dnaAlignment <- DECIPHER::AlignSeqs(allpart)
    DECIPHER::BrowseSeqs(dnaAlignment,
                         htmlFile = file.path(outputDir, glue::glue("{part}DNAAlignment.html")),
                         openURL = F, colWidth = 120)
  }
  
  aaPartFiles <- list.files(annotaleMainDir, "TALE_Protein_parts.fasta", recursive = T, full.names = T)
  endsAA <- sapply(c("N-terminus", "C-terminus"), function(part) {
    allpart <- sapply(aaPartFiles, function(p) {
      allpart <- Biostrings::readAAStringSet(p, seek.first.rec = T)
      onepart <- allpart[grepl(part, names(allpart))]
      talRoi <- basename(dirname(p))
      names(onepart) <- talRoi
      return(onepart)
    }, simplify = "array", USE.NAMES = F) %>%
      Biostrings::AAStringSetList() %>%
      unlist()
    dnaAlignment <- DECIPHER::AlignSeqs(allpart)
    DECIPHER::BrowseSeqs(dnaAlignment,
                         htmlFile = file.path(outputDir, glue::glue("{part}AAAlignment.html")),
                         openURL = F, colWidth = 120)
    return(allpart)
  }, USE.NAMES = T)
  
  
  #### Protein length of N-term and C-term ####
  
  endsAAlength <- lapply(names(endsAA), function(e) {
    stringset <- endsAA[e] %>% Biostrings::AAStringSetList(., use.names = F) %>% unlist()
    df <- data.frame(names(stringset), BiocGenerics::width(stringset))
    colnames(df) <- c("arrayID", paste0(e, "AAlength"))
    return(df)
  })
  
  
  S4Vectors::mcols(hitsByArraysLst) <- merge(S4Vectors::mcols(hitsByArraysLst),
                                             do.call(merge, endsAAlength),
                                             by = "arrayID", 
                                             all.x = T)
  
  #### 
  
  
  ## Merge with arrays metadata in mcols(hitsByArraysLst)
  moreInfo <- merge(
    S4Vectors::mcols(hitsByArraysLst),
    data.frame(arrayID = names(fullTalOrf),
               LongestOrfLength = Biostrings::nchar(fullTalOrf),
               OrfCovOverArrayLength = round(100 * Biostrings::nchar(fullTalOrf)/GenomicRanges::width(extdCompleteArraysSeqs[names(fullTalOrf)])),
               LongestORFSeq = fullTalOrf),
    by = "arrayID", all.x = TRUE, sort = FALSE)
  rownames(moreInfo) <- moreInfo$arrayID
  S4Vectors::mcols(hitsByArraysLst) <- moreInfo[rownames(S4Vectors::mcols(hitsByArraysLst)),]
  #str(S4Vectors::mcols(hitsByArraysLst))
  

  ####   IS THIS STILL USEFULL?   ####
  ####   Look at gaps between HitDomains on the same subject sequence to detect potential missed repeats   #####
  arraysBySeqlevelLst <- split(arraysGR, GenomicRanges::seqnames(arraysGR)) # group  arraysGR by seqlevel
  
  gaplengthBetweenHitDomains <- sapply(arraysBySeqlevelLst, function(x) {
    dists <- t(as.data.frame(GenomicRanges::distanceToNearest(x)))[3,]
  })
  gaplengthBetweenHitDomains <- unlist(gaplengthBetweenHitDomains)
  gaplengthBetweenHitDomains <- gaplengthBetweenHitDomains[!is.na(gaplengthBetweenHitDomains)]
  gaplengthBetweenHitDomainsbelow500 <- gaplengthBetweenHitDomains[gaplengthBetweenHitDomains <= 500]
  quartilesGapLength <- quantile(gaplengthBetweenHitDomainsbelow500,  probs = c(0.25, 0.50, 0.75))
  ##ggplot(data.frame(gapSize = gaplengthBetweenHitDomainsbelow500), aes(x=gapSize)) + geom_histogram(binwidth=10)
  
  
  ####   Write tabulated reports   #####
  ## Report on the hmmer hits. both  as a  tab-delimited  and a gff
  hitsReport <- lapply(as.list(hitsByArraysLst, use.names = TRUE),
                       function(gr) {
                         tibble::as_tibble(as.data.frame(gr))
                       }
  ) %>% dplyr::bind_rows(.id = "arrayID")
  readr::write_tsv(x = hitsReport, file = hitsReportFile)
  hitsReportToGFF(hitsReportFile) # saving to gff format
  
  readr::write_tsv(x = domainsReport, file = domainsReportFile)
  
  ##   Report with info on arrays, including the seq of RVD
  arrayReport <- as.data.frame(
    S4Vectors::mcols(hitsByArraysLst)[
      order(
        S4Vectors::mcols(hitsByArraysLst)$OriginalSubjectName,
        S4Vectors::mcols(hitsByArraysLst)$NumberOfHits
      ),]
  )
  readr::write_tsv(x = arrayReport, file = arrayReportFile)
  
  
  ## Write a gff with all collated
  allGR <- c(unlist(hitsByArraysLst),
             GenomicRanges::makeGRangesFromDataFrame(
               S4Vectors::mcols(hitsByArraysLst),
               seqnames.field= "OriginalSubjectName",
               keep.extra.columns= TRUE),
             if (exists("reducedOlapGr")) nhmmerOutputGRBeforeMerge else NULL
  )
  rtracklayer::export.gff3(allGR, affRangesGffFile)
  
  ####   Extract various DNA/AA sequences of interest and save to files   #####
  
  
  ## Write a fasta file of the seq of RVDs
  seqsOfRVDs <- Biostrings::BStringSet(S4Vectors::mcols(hitsByArraysLst)$SeqOfRVD)
  names(seqsOfRVDs) <- S4Vectors::mcols(hitsByArraysLst)$arrayID
  seqsOfRVDs <- seqsOfRVDs[!width(seqsOfRVDs) == 0]
  Biostrings::writeXStringSet(x = seqsOfRVDs, seqsOfRVDFile)
  
  ####   Generate info messages and log file about the analysis   #####
  
  ## counts of appearance of each RVD type (excluding N- and C- terms symbols) for the log file
  #RVDtbl <- table(subset(unlist(hitsByArraysLst), query_name ==repeatDNAHMMName, drop = TRUE)$RVD)
  ## Total count of repeat CDS after filtering for uniformative subject seqs for the log file
  numberOfRepeatHitsAfterFiltering <- length(subset(unlist(hitsByArraysLst), query_name == repeatDNAHMMName))
  ## Distribution of the number of hits per array
  countsHitsByArrayDistri <- summary(S4Vectors::mcols(hitsByArraysLst)$NumberOfHits)
  ## Number of domains in arrays that display all domain types
  # completeArrayLengths <- subset(S4Vectors::mcols(hitsByArraysLst), AllDomains)$NumberOfHits
  
  
  ## might have been cleaner with a glue approach
  txt <- c(
    "#****************************************",
    "#**   tellTale analysis done     **",
    
    paste("Current date:", date(), sep = "\t"),
    "#_________Provided I/O parameters __________",
    paste("File of subject DNA sequences:", subjectFile, sep = "\t"),
    paste("TALE N-term CDS region detection HMM file:", TALE_NtermDNAHMMFile, sep = "\t"),
    paste("TALE repeat unit CDS detection HMM file:", repeatDNAHMMFile, sep = "\t"),
    paste("TALE C-term CDS region detection HMM file:", TALE_CtermDNAHMMFile, sep = "\t"),
    paste("Output directory:", outputDir, sep = "\t"),
    
    "#____________Other parameters________________",
    paste(Quote(TALE_NtermDNAHitMinScore),":", TALE_NtermDNAHitMinScore, sep = "\t"),
    paste(Quote(repeatDNAHitMinScore),":", repeatDNAHitMinScore, sep = "\t"),
    paste(Quote(TALE_CtermDNAHitMinScore),":", TALE_CtermDNAHitMinScore, sep = "\t"),
    paste(Quote(minDomainHitsPerSubjSeq),":", minDomainHitsPerSubjSeq, sep = "\t"),
    paste(Quote(mergeHits),":", mergeHits, sep = "\t"),
    # paste(Quote(repMsaMethod),":", repMsaMethod, sep = "\t"),
    paste(Quote(minGapWidth),":", minGapWidth, sep = "\t"),
    paste(Quote(taleArrayStartAnchorCode),":", taleArrayStartAnchorCode, sep = "\t"),
    paste(Quote(taleArrayEndAnchorCode),":", taleArrayEndAnchorCode, sep = "\t"),
    paste(Quote(extendedLength),":", extendedLength, sep = "\t"),
    paste(Quote(talArrayCorrection),":", talArrayCorrection, sep = "\t"),
    paste(Quote(refForTalArrayCorrection),":", refForTalArrayCorrection, sep = "\t"),
    paste(Quote(frameShiftCorrection),":", frameShiftCorrection, sep = "\t"),
    
    "#__________Summary measures of TALE search outcome__________",
    paste("Number of analysed subject sequences :", length(subjectDNASequences), sep = "\t"),
    paste("Total number of TALE repeat DNA coding sequence motif hits found with the nhmmer approach:",
          numberOfRepeatHitsAfterFiltering, sep = "\t"),
    #paste("Total number of repeat HMM hits on the corresponding set of translated DNA hits:", sum(RVDtbl), sep = "\t"),
    
    paste("Total number of subject seqs with TALE motif hits after low hit number filtering:",
          length(GenomeInfoDb::seqlevels(arraysGR)), sep = "\t"),
    paste("Total number of distinct regions (repeat arrays) with adjacent TALE motifs :", nrow(arrayReport), sep = "\t"),
    paste("Total number of 'complete' arrays (with both N- and C-term flanking motifs):",
          sum(S4Vectors::mcols(hitsByArraysLst)$AllDomains),	sep = "\t"),
    
    #paste("Total number of distinct types of RVD:", nrow(RVDtbl), sep = "\t"),
    
    paste("Minimum array length (number of TALE domain hits):", min(arrayReport$NumberOfHits), sep = "\t"),
    paste("Maximum array length:", max(arrayReport$NumberOfHits), sep = "\t"),
    paste("Median array length:", median(arrayReport$NumberOfHits), sep = "\t"),
    # paste("Length of the longest 'complete' array:", max(completeArrayLengths),	sep = "\t"),
    # paste("Length of the shortest 'complete' array:", min(completeArrayLengths),	sep = "\t"),
    
    #message("Distribution of the number of TALE domain hits per array:\n")
    #message(paste(names(countsHitsByArrayDistri), countsHitsByArrayDistri, sep = "\t", collapse = "\n"))
    
    paste("Number of gaps of size below 500nt between TALE motifs arrays:", length(gaplengthBetweenHitDomainsbelow500)/2, sep = "\t"),
    
    paste("First quartile of size of gaps (below 500nt) between TALE motifs arrays:", quartilesGapLength[1], sep = "\t"),
    paste("Median size of gaps (below 500nt) between TALE motifs arrays:", quartilesGapLength[2], sep = "\t"),
    paste("Upper quartile of size of gaps (below 500nt) between TALE motifs arrays:", quartilesGapLength[3], sep = "\t"),
    
    "#__________Noteworthy AnnoTale issues__________",
    paste("#", annoTaleMessages),
    
    "#*************************\n"
  )
  
  message(paste(txt, collapse = "\n"))
  logf <- file(analysisLogFile, open = "w")
  writeLines(text = txt, con = logf)
  close(logf)
  return(invisible(outputDir))
}

#' This function name is deprecated and will ultimately be removed.
#' It corresponds to the \link{tellTale} which should be used instead.
#'
#' @export
tellTale2 <- tellTale


