correction_tible <- function(indels) {
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
  "annout",
  
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

# subjectFile = "/home/baotram/tal/xanthopore-scripts/tantale/test/telltale_vs_annotale/genomes/BAI3-1-1.fa"
# outputDir = "/home/baotram/tal/xanthopore-scripts/tantale/test/telltale2/bai311_correct"
# hmmFilesDir = "/home/baotram/tal/xanthopore-scripts/tantale/inst/extdata/hmmProfile"
# hmmerpath = "/home/baotram/tal/xanthopore-scripts/tantale/inst/tools/hmmer-3.3/bin/"
# talArrayCorrection = TRUE
# refForTalArrayCorrection = "/home/baotram/tal/xanthopore-scripts/talomes_analysis/correctframshift_ref/AA_ref.fa.gz"
# frameShiftCorrection = -11
# minRatioOfGapForColMasking = 0.8
# TALE_NtermDNAHitMinScore = 300
# repeatDNAHitMinScore = 20
# TALE_CtermDNAHitMinScore = 200
# minDomainHitsPerSubjSeq = 4
# mergeHits = TRUE
# # repMsaMethod = "decipher"
# minGapWidth = 35
# minDomainHitsPerArrayForAssembl  = 5
# taleArrayStartAnchorCode = "N-TERM"
# taleArrayEndAnchorCode = "C-TERM"
# appendExtremityCodes = TRUE
# rvdSep = "-"
# extendedLength = 300

#' Search and report on the features of TALE protein domains potentially encoded
#' in subject DNA sequences
#'
#' \code{tellTale} has been primarily written to report on 'corrected' TALE RVD
#' sequences in noisy DNA sequences (suboptimally polished genomes assembly, raw
#' reads of long read sequencing technologies [eg PacBio, ONT]) that would
#' otherwise be missed by conventional tools (eg AnnoTALE).
#'
#' It works but is far from optimal for 'corrected' repeat CDS and  \strong{it
#' currently tends to remove 'unconventional' portions of repeat CDS}. This may
#' be problematic for downstream analysis with DisTAL, especially for 'aberrant'
#' repeats.
#'
#' N- and C-terminal domains CDS are not currently 'corrected'.
#'
#' It will work also if the subject DNA sequences are of high quality
#' \strong{but may 'correct' genuine frame shifts in TALE repeat array CDS}.
#'
#' Will output a number of files including tabular reports on the TALE arrays
#' (ie full length or partial tal gene CDS) and on the identified coding
#' sequences for TALE domains, including repeat CDS.
#'
#' @param subjectFile Fasta file with DNA sequence(s) to be searched for the
#'   presence of TALE coding sequences (CDS).
#' @param outputDir Path of the output directory. If not specified, resluts will
#'   be written to working folder.
#' @param hmmFilesDir Specify the path to a folder holding the hmmfiles if you
#'   do not want to use the ones provided with tantale.
#' @param minRatioOfGapForColMasking Columns of the tale repeat CDS alignment
#'   that contain a gap in a fraction of sequence higher than this value (betwen
#'   0 and 1) will be masked from the alignment when translating the DNA
#'   sequences to protein.
#' @param TALE_NtermDNAHitMinScore Minimal nhmmer score cut_off value to
#'   consider the hit as genuine
#' @param repeatDNAHitMinScore Minimal nhmmer score cut_off value to consider
#'   the hit as genuine
#' @param TALE_CtermDNAHitMinScore Minimal nhmmer score cut_off value to
#'   consider the hit as genuine
#' @param minDomainHitsPerSubjSeq Minimum number of nhmmer hits for a subject
#'   sequence to be reported as having TALE diagnostic regions. This is a way to
#'   simplify output a little by getting ride of uninformative sequences
#' @param mergeHits Perform overlapping hits merging per domain type.
#' @param minGapWidth Minimum gap in base pairs between two tale domain hits for
#'   them to be considered distinct. If the length of the gap is below this
#'   value, domains are considered "contiguous" and grouped in the same array.
#' @param minDomainHitsPerArrayForAssembl DEPRECATED argument. Used to speficy
#'   the Minimum number of repeat in an array for its seq of RVD to be
#'   considered for assembly. This is a way to get ride of sequences that are
#'   too short reasonably be of any help for assembly
#' @param taleArrayStartAnchorCode This scalar character vector will symbolize a
#'   TALE N-TERM CDS hit in the RVD sequence
#' @param taleArrayEndAnchorCode This scalar character vector will symbolize a
#'   TALE C-TERM CDS hit in the RVD sequence
#' @param appendExtremityCodes Set this to \code{FALSE} if you do not want the N- and C-TREM anchor codes in the output sequences of RVD
#' @param rvdSep Symbol acting as a separator in RVD sequences
#' @param hmmerpath Specify the path to a directory holding the HMMER
#'   executable if you do not want to use the ones provided with tantale.
#' @param extendedLength number of nucleotides to extend in 3'-end
#' @param talArrayCorrection True or False
#' @param refForTalArrayCorrection reference AA sequences for tal array correction
#' @param frameShiftCorrection default = 11
#' @param ... \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}}
#' @return This functions has only side effects (writing files, mostly).
#'   However, if everything ran smoothly, it will invisibly return the path of
#'   the directory where output files were written.
#'   In the arrayReport.tsv, column \emph{predicted_dels_count}/\emph{predicted_ins_count} shows the number of putative deletions/insertions in the raw sequences that have been corrected in the corrected sequences by deletion/insertion with the function \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}}.
#'   LIst of output files:
#'   \itemize{
#'   \item allRanges.gff: gff file of all Tal arrays detected by HMMer
#'   \item arrayReport.tsv: report of all Tal arrays
#'   \item hitsReport.tsv: report of all hits detected by HMMer
#'   \item hitsReport.gff: gff file of all hits detected by HMMer
#'   \item domainsReport.tsv: report of all Tal domains detected by AnnoTALE analyze
#'   \item putativeTalOrf.fasta: Tal putative ORFs
#'   \item pseudoTalCds.fasta: pseudo Tal CDS, Tal arrays detected by HMMer but failed in AnnoTALE analyze
#'   \item C-terminusAAAlignment.html: protein alignment of all C-termini
#'   \item C-terminusDNAAlignment.html: DNA alignment of all C-termini
#'   \item N-terminusAAAlignment.html: protein alignment of all N-termini
#'   \item N-terminusDNAAlignment.html: DNA alignment of all N-termini
#'   \item TALE_CDS_all_diagnostic_regions_hmmfile.out: 
#'   \item hmmerSearchOut.txt: 
#'   \item nhmmerHumanReadableOutputOfLastRun.txt: 
#'   \item tellTale.log: 
#'   \item temp_annotale folder: folder containing result of AnnoTALE analyze for all Tal arrays
#'   \item CorrectionAlignmentAA folder: folder containing protein alignment of Tal array detected by HMMer and corrected Tal array if \code{talArrayCorrection} = TRUE
#'   \item CorrectionAlignmentDNA folder: folder containing DNA alignment of Tal array detected by HMMer and corrected Tal array if \code{talArrayCorrection} = TRUE
#' }
#' @export
tellTale2 <- function(
  subjectFile,
  outputDir = getwd(),
  hmmFilesDir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T),
  minRatioOfGapForColMasking = 0.8,
  TALE_NtermDNAHitMinScore = 300,
  repeatDNAHitMinScore = 20,
  TALE_CtermDNAHitMinScore = 200,
  minDomainHitsPerSubjSeq = 4,
  mergeHits = TRUE,
  # repMsaMethod = "decipher",
  minGapWidth = 35,
  minDomainHitsPerArrayForAssembl  = 5,
  taleArrayStartAnchorCode = "NTERM",
  taleArrayEndAnchorCode = "CTERM",
  appendExtremityCodes = TRUE,
  rvdSep = "-",
  hmmerpath = system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = T),
  extendedLength = 300,
  talArrayCorrection = TRUE,
  refForTalArrayCorrection = system.file("extdata", "AA_ref.fa.gz", package = "tantale", mustWork = T),
  frameShiftCorrection = -11,
  ...
) {
  
  dir.create(outputDir, recursive = T, mode = "755")
  
  #####   Paths of output files   #####
  ## Sequences of hmmer hits extracted from subject sequence using hmmer hit positions
  repeatDNASeqsFile <- file.path(outputDir, "hmmerRepeatHitsDNASeqs.fas")
  ## hmmerAlign alignment of repeat CDS
  repeatsAlignOutFile <- file.path(outputDir, "repeatsAlignOut.txt")
  ## sequences of the 'corrected/cleaned' repeat CDS
  correctedRepeatsAlignOutFile <- file.path(outputDir, "hmmerCleanedAlignOut.fas")
  ## sequences of the 'corrected/cleaned' repeat Translations after removing stop codon translations (*)
  translatedCorrectedRepeatsAlignOutFile <- file.path(outputDir, "translatedCleanRepeatsAlignment.pep")
  ## output of hmmer search on the repeat AA sequences
  repeatAaSearchOutFile <- file.path(outputDir, "hmmerAASearchOut.txt")
  ## hmmerAlign alignment of the repeat AA sequences
  repeatsAaAlignOutFile <- file.path(outputDir, "hmmerAAAlignOut.txt")
  ## fasta file of central repeat domain "corrected" CDS for each array
  correctedArrayDnaSeqsFile <- file.path(outputDir, "correctedRepeatDNASeqs.fas")
  ## fasta file of translation of products of central repeat domain "corrected" CDS for each array
  correctedArrayAaSeqsFile <- file.path(outputDir, "correctedRepeatAASeqs.fas")
  
  
  ## Fasta (.pep) file where sequences of translated TALE N-Term CDS domains are written
  fullNtermAAseqFile <- file.path(outputDir, "fullNtermAAseq.pep")
  ## Fasta (.pep) file where sequences of translated TALE C-Term CDS domains are written
  fullCtermAAseqFile <- file.path(outputDir, "fullCtermAAseq.pep")
  
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
  arrayOrfsSeqFile <- file.path(outputDir, "arrayOrfs.fas")
  ## A text file where logging info and some general analysis measures are written
  analysisLogFile <- file.path(outputDir, "tellTale.log")
  
  
  #####   Checks for parameters and other things   #####
  if (minDomainHitsPerSubjSeq > minDomainHitsPerArrayForAssembl) {
    warning("The value of the minDomainHits parameter for subject sequence selection\n",
            "is more stringent than the value of the minNumberOfDomainHitsForAssembly!?\n",
            "Is this really what you want to do?")
  }
  
  
  
  
  ## Deal with spaces in sequence names because this messes up parsing of HMMER output
  originalSeqs <- Biostrings::readDNAStringSet(filepath = subjectFile)
  originalSeqlevels <- names(originalSeqs)
  foolproofSeqlevels <- paste0("seq", 1:length(originalSeqlevels))
  names(originalSeqlevels) <- foolproofSeqlevels
  names(originalSeqs) <- foolproofSeqlevels
  subjectFile <- tempfile()
  Biostrings::writeXStringSet(originalSeqs, filepath = subjectFile)
  
  #####   Full paths of input HMM files for TALE domains (DNA and AA)   ####
  
  ## TODO come up with a mechanism for the user to be able to provide the FULL PATH
  ## to custom hmm !!! hmmFilesDir parameter is useless unless custom hmm are named
  ## as specified below
  TALE_NtermDNAHMMFile <- file.path(hmmFilesDir, "Xo_TALE_Nterm_CDS_profile.hmm")
  repeatDNAHMMFile <- file.path(hmmFilesDir, "Xo_TALE_repeat_CDS_profile.hmm")
  TALE_CtermDNAHMMFile <- file.path(hmmFilesDir, "Xo_TALE_Cterm_CDS_profile.hmm")
  repeatAAHMMFile <- file.path(hmmFilesDir, "Xo_TALE_repeat_AA_profile.hmm")
  
  DNAHMMFiles <- c(TALE_NtermDNAHMMFile, repeatDNAHMMFile, TALE_CtermDNAHMMFile)
  mergedDNAHMMFile <- file.path(outputDir, "TALE_CDS_all_diagnostic_regions_hmmfile.out")
  
  
  
  #####   Concatenate HMM files for TALE DNA motifs and parse HMM names   #####
  hmmslines <- plyr::llply(DNAHMMFiles, function(x) txt <- readLines(con = x))
  names(DNAHMMFiles) <- plyr::llply(hmmslines, function(x) {
    hmmName <- grep("NAME", x, perl = TRUE, value = TRUE)
    hmmName <- unlist(strsplit(hmmName, split = "\\s+"))
    if (length(hmmName) != 2) stop("One or several profile ",
                                   "HMM have a name with with spaces. ",
                                   "Please remove them in the file at the Tag 'NAME'")
    hmmName <- hmmName[2]
  }
  )
  DNAHMMNames <- names(DNAHMMFiles)
  TALE_NtermDNAHMMName <- names(DNAHMMFiles)[1]
  repeatDNAHMMName <- names(DNAHMMFiles)[2]
  TALE_CtermDNAHMMName <- names(DNAHMMFiles)[3]
  
  writeLines(text = unlist(hmmslines), con = mergedDNAHMMFile)
  
  #####   Perform TALE domain CDS search with HMMER  #####
  searchOutFile <- file.path(outputDir, "hmmerSearchOut.txt")
  
  runNhmmerSearch(hmmerpath = hmmerpath,
                  subjectFile = subjectFile,
                  hmmFile = mergedDNAHMMFile,
                  searchTblOutFile = searchOutFile,
                  humReadableOutFile = file.path(outputDir, "nhmmerHumanReadableOutputOfLastRun.txt"))
  
  #####   Load, process, filter TALE domain CDS HMMER hit results    #####
  ## Loading search tabular output file
  nhmmerTabularOutput <- read.table(searchOutFile)
  
  colnames(nhmmerTabularOutput) <- c("target_name", "accession", "query_name", "accession", "hmmfrom", "hmm_to", "alifrom",
                                     "ali_to", "envfrom", "env_to", "sq_len", "strand", "Evalue", "score", "bias", "description_of_target")
  
  ## filtering results differentially depending on the query HMM
  nhmmerTabularOutput <- subset(nhmmerTabularOutput,
                                query_name == TALE_NtermDNAHMMName & score >= TALE_NtermDNAHitMinScore |
                                  query_name == repeatDNAHMMName & score >= repeatDNAHitMinScore |
                                  query_name == TALE_CtermDNAHMMName & score >= TALE_CtermDNAHitMinScore
  )
  
  nhmmerTabularOutput <- droplevels(nhmmerTabularOutput)
  
  ## Add a hitID column
  nhmmerTabularOutput$hitID <- paste("DOM", sprintf("%05.0f", 1:nrow(nhmmerTabularOutput)), sep="_")
  
  ## Trick to re-order positions in an increasing order to satisfy IRanges() in preparation of creating a GRanges
  nhmmerTabularOutput[,c("start", "end")] <- plyr::adply(.data = nhmmerTabularOutput[,c("envfrom", "env_to")], .margins = 1, .fun = c(min, max))[,-(1:2)]
  # nhmmerTabularOutput[,c("start", "end")] <- plyr::adply(.data = nhmmerTabularOutput[,c("alifrom", "ali_to")], .margins = 1, .fun = c(min, max))[,-(1:2)]
  
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
    seqnames.field= "target_name"
  )
  names(nhmmerOutputGR) <- nhmmerOutputGR$hitID
  GenomeInfoDb::seqinfo(nhmmerOutputGR) <- GenomeInfoDb::Seqinfo(
    seqnames = as.character(temp_df$target_name),
    seqlengths=temp_df$sq_len,
    isCircular=NA,
    genome=NA
  )
  GenomeInfoDb::seqlevels(nhmmerOutputGR) <- originalSeqlevels[match(GenomeInfoDb::seqlevels(nhmmerOutputGR), names(originalSeqlevels))]
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
  stitchSeqs <- function(grl, seqType, sep = "-", onlyRepeats = FALSE) {
    BiocGenerics::sapply(grl, function(x) {
      if (onlyRepeats) {
        gr <- BiocGenerics::subset(x = x, query_name == repeatDNAHMMName)
      } else {
        gr <- x
      }
      v <- as.character(S4Vectors::mcols(gr)[, seqType])
      if (all(as.character(BiocGenerics::strand(gr)) == "+")) {
        seq <- paste0(v, collapse = sep)
      } else if (all(as.character(BiocGenerics::strand(gr)) == "-")) {
        seq <- paste0(rev(v), collapse = sep)
      } else {
        stop("The current array has domains on different strands. Something is very wrong!")
      }
    }, simplify = TRUE, USE.NAMES = TRUE)
  }
  
  
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
  
  
  
  #str(S4Vectors::mcols(hitsByArraysLst))
  #####   Extend DNA Tal arrays   #####
  ## Extract the genomic sequence of arrays +-bp on the borders
  ## run systemPipeR::predORF()
  ## ask if the longest orf GRanges on the plus strand (perfectly) matches with the complete array GRanges
  ## Record this info in the arrayReport object
  ## Export the dna seq of the longest orf as the tal sequences
  completeArraysGR <- arraysGR #subset(arraysGR, S4Vectors::mcols(hitsByArraysLst)$AllDomains) #
  extdCompleteArraysGR <- GenomicRanges::resize(completeArraysGR,
                                                width = GenomicRanges::width(completeArraysGR) + extendedLength,
                                                fix="start", ignore.strand=FALSE) #%>%
  # GenomicRanges::resize(., width = GenomicRanges::width(.) + 3,
  #                       fix="end", ignore.strand=FALSE)
  extdCompleteArraysSeqs <- BSgenome::getSeq(subjectDNASequences, extdCompleteArraysGR)
  
  
  # correct tal arrays if requested
  if (!talArrayCorrection) {
    #### Get ORFs from uncorrected Tal arrays ####
    # fullTalOrf <- extdCompleteArraysSeqs
    # TalOrfForAnnoTALE <- extdCompleteArraysSeqs
    orfs <- systemPipeR::predORF(x = extdCompleteArraysSeqs,
                                 n=1, type = "gr", mode = "ORF", strand = "sense")
    fullTalOrf <- BSgenome::getSeq(extdCompleteArraysSeqs, orfs)
    names(fullTalOrf) <- as.character(GenomicRanges::seqnames(orfs))
    TalOrfForAnnoTALE <- fullTalOrf
  } else { 
    #### Correct Tal arrays ####
    AAref <- Biostrings::readAAStringSet(refForTalArrayCorrection, seek.first.rec = TRUE, use.names = TRUE)
    rawArraySeq <- extdCompleteArraysSeqs
    ArrayCorrection <- DECIPHER::CorrectFrameshifts(rawArraySeq, AAref, type = "both", maxComparisons = length(AAref), frameShift = frameShiftCorrection, ...)
    corrExtdCompleteArraysSeqs <- ArrayCorrection$sequences
    
    #### Correction stats ####
    deletions_count <- correction_tible(ArrayCorrection$indels) %>% dplyr::group_by(Seq) %>% dplyr::count(variable, name = "predicted_dels_count") %>% dplyr::filter(variable == "deletions") %>% dplyr::select(-variable)
    insertions_count <- correction_tible(ArrayCorrection$indels) %>% dplyr::group_by(Seq) %>% dplyr::count(variable, name = "predicted_ins_count") %>% dplyr::filter(variable == "insertions") %>% dplyr::select(-variable)
    
    
    S4Vectors::mcols(hitsByArraysLst) <- merge(S4Vectors::mcols(hitsByArraysLst), 
                                               dplyr::full_join(insertions_count, deletions_count, by = "Seq"), 
                                               by.x = "arrayID", by.y = "Seq", all.x = T) 
    
    S4Vectors::mcols(hitsByArraysLst)[c("predicted_dels_count", "predicted_ins_count")] %<>% apply(., 2, function(v) ifelse(is.na(v), 0, v))
    
    #### Correction alignments ####
    for (n in names(corrExtdCompleteArraysSeqs)[vcountPattern("N", corrExtdCompleteArraysSeqs) > 0]) {
      warning(glue::glue("Sequence {n} contains 'N's and will be substituted by 'C's in order to run AnnoTALE analyze for RVDs prediction."))
    }
    TalOrfForAnnoTALE <- Biostrings::chartr("N", "C", corrExtdCompleteArraysSeqs)
    
    alignmentDNADir <- file.path(outputDir, "CorrectionAlignmentDNA")
    dir.create(alignmentDNADir, showWarnings = F)
    
    
    alignmentAADir <- file.path(outputDir, "CorrectionAlignmentAA")
    dir.create(alignmentAADir, showWarnings = F)
    
    # correctedTalOrf <- corrExtdCompleteArraysSeqs[sapply(names(corrExtdCompleteArraysSeqs), function(n) n %in% c(insertions_count$Seq, deletions_count$Seq))]
    for (n in names(corrExtdCompleteArraysSeqs)) {
      rawSeq <- rawArraySeq[n]
      names(rawSeq) <- paste0("raw_", n)
      correctedSeq <- corrExtdCompleteArraysSeqs[n]
      names(correctedSeq) <- paste0("corrected_", n)
      substitutedSeq <- TalOrfForAnnoTALE[n]
      names(substitutedSeq) <- paste0("forAnnoTALE", n)
      seqToAlign <- c(rawSeq, correctedSeq, substitutedSeq)
      alignedSeqs <- DECIPHER::AlignSeqs(seqToAlign)
      DECIPHER::BrowseSeqs(alignedSeqs, htmlFile = file.path(alignmentDNADir, glue::glue("CorrectionAlignmentDNA_{n}.html")), openURL = F, colWidth = 120)
      
      seqToAlignTranslated <- Biostrings::translate(seqToAlign, no.init.codon = T, if.fuzzy.codon = "solve")
      alignedSeqsTranslated <- DECIPHER::AlignSeqs(seqToAlignTranslated)
      DECIPHER::BrowseSeqs(alignedSeqsTranslated, htmlFile = file.path(alignmentAADir, glue::glue("CorrectionAlignmentAA_{n}.html")), openURL = F, colWidth = 120)
      
    }
    
    #### detect ORFs from corrected arrays that have been substituted 'N' by 'C' ####
    orfs <- systemPipeR::predORF(x = TalOrfForAnnoTALE,
                                 n=1, type = "gr", mode = "ORF", strand = "sense")
    TalOrfForAnnoTALE <- BSgenome::getSeq(TalOrfForAnnoTALE, orfs)
    names(TalOrfForAnnoTALE) <- as.character(GenomicRanges::seqnames(orfs))
    
    fullTalOrf <- TalOrfForAnnoTALE
    insPosition <- vmatchPattern("N", corrExtdCompleteArraysSeqs)
    for (i in names(insPosition)) {
      if (length(width(insPosition[[i]])) == 0) next()
      insPositionAfCorr <- insPosition[[i]][BiocGenerics::start(insPosition[[i]]) < width(fullTalOrf[i])]
      fullTalOrf[i] <- replaceAt(fullTalOrf[i], insPositionAfCorr, value = "N")
    }
  }
  
  
  
  #### Get RVD seqs from Tal ORFs by annoTALE analyze ####
  
  AnnoTALEanalyze <- function(inputFastaFile,
                              outputDir = getwd(),
                              prefix = NULL,
                              annoTALE = system.file("tools", "AnnoTALEcli-1.4.1.jar", package = "tantale", mustWork = T)
  ) {
    # Define output dirs for the various stages of annoTALE
    stopifnot(dir.exists(outputDir) || dir.create(path = outputDir, showWarnings = TRUE, recursive = TRUE, mode = "775"))
    # Define a prefix for TALEs (the strain or assembly ID) derived from the genome file name.
    if( is.null(prefix) ) {
      prefix <- gsub(pattern = "^(.*)\\.(fasta|fa|fas)$" , replacement  = "\\1", basename(inputFastaFile), perl = TRUE)
    }
    # Run the "analyze" stage of annoTALE
    comAnalyze <- paste0(
      "java -jar ", annoTALE,
      " analyze ",
      " t=", inputFastaFile,
      " outdir=", outputDir
    )
    cat("##  Now running annoTALE analyze for", prefix, "using the following command:\n##  ",  comAnalyze, "\n")
    exitAnalyze <- system(comAnalyze)
    return(invisible(exitAnalyze))
  }
  
  # tempdir for annotale
  tmpAnnotaleDir <- tempfile(pattern = "temp_annotale_", tmpdir = outputDir)
  dir.create(tmpAnnotaleDir)
  
  # run annotale for tal arrays on tempdir
  
  # setClass(
  #   # Set the name for the class
  #   "annout",
  # 
  #   # Define the slots
  #   slots = c(
  #     domainsReport = "data.frame"
  #   ),
  # 
  #   contains = "AAStringSet",
  # 
  #   # Make a function that can test to see if the data is consistent.
  #   # This is not called if you have an initialize function defined!
  #   validity = function(object) {
  #     val <- is.data.frame(object@domainsReport)
  #     return(val)
  #   }
  # )
  
  annoTaleOut <- sapply(names(TalOrfForAnnoTALE), function(talOrfID) {
    AnnotaleDir <- file.path(tmpAnnotaleDir, talOrfID)
    dir.create(AnnotaleDir)
    
    TalOrf <- TalOrfForAnnoTALE[talOrfID]
    
    correctedTalOrfFile <- file.path(AnnotaleDir, "putativeTalOrf.fasta")
    Biostrings::writeXStringSet(TalOrf, correctedTalOrfFile)
    
    checkAnnoTale <- try(AnnoTALEanalyze(correctedTalOrfFile, AnnotaleDir))
    if (!is.null(attr(checkAnnoTale, "condition"))) return(annout(Biostrings::AAStringSet(), domainsReport = data.frame())) # in case annotale does not work
    annoTaleRVD <- file.path(AnnotaleDir, "TALE_RVDs.fasta")
    seqOfRVDs <- try(Biostrings::readAAStringSet(annoTaleRVD, seek.first.rec = T, use.names = T))
    if (!is.null(attr(seqOfRVDs, "condition"))) return(annout(Biostrings::AAStringSet(), domainsReport = data.frame())) # in case annotale works but cannot find rvds
    if (Biostrings::width(seqOfRVDs) == 0) return(annout(Biostrings::AAStringSet(), domainsReport = data.frame())) # in case annotale does not complain about rvds but does not output any rvds
    names(seqOfRVDs) <- talOrfID
    
    ## domains report
    parts_files <- list.files(AnnotaleDir, "TALE_Protein_parts.fasta", recursive = T, full.names = T)
    parts <- Biostrings::readAAStringSet(parts_files)
    if (Biostrings::width(parts) == 0) {
      file.remove(parts_files)
      return(annout(Biostrings::AAStringSet(), domainsReport = data.frame()))
      } # in case annotale works but protein parts file is empty
    stops <- Biostrings::vcountPattern("*", parts)
    domainsReport <- data.frame("arrayID" = talOrfID,
                                "seqnames" = S4Vectors::mcols(hitsByArraysLst)$OriginalSubjectName[S4Vectors::mcols(hitsByArraysLst)$arrayID == talOrfID],
                                "query_name" = gsub("(.+\\: )|( \\d+)", "", names(parts)),
                                "codon_count" = width(parts) - stops
                                )
    
    annotale_output <- annout(seqOfRVDs, domainsReport = domainsReport)
    return(annotale_output)
  }, simplify = T, USE.NAMES = F)
  
  seqsOfRVDs <- unlist(Biostrings::AAStringSetList(annoTaleOut))
  domainsReport <- do.call(rbind, lapply(annoTaleOut, function(x) x@domainsReport))
  
  # save tals that have rvds
  correctedTalArrayFile <- file.path(outputDir, "putativeTalOrf.fasta")
  Biostrings::writeXStringSet(fullTalOrf[names(fullTalOrf) %in% names(seqsOfRVDs)], correctedTalArrayFile)
  
  # save tals that DO NOT have rvds
  pseudoTalFile <- file.path(outputDir, "pseudoTalCds.fasta")
  Biostrings::writeXStringSet(extdCompleteArraysSeqs[!names(extdCompleteArraysSeqs) %in% names(seqsOfRVDs)], pseudoTalFile)
  
  # unlink(tmpAnnotaleDir, recursive = T)
  aberrantRepeat <- sapply(seqsOfRVDs, function(s) {
    ifelse(length(s) > 0, grepl("[a-z]", s), NA)
  })
  
  seqsOfRVDs <- gsub("\\-", rvdSep, seqsOfRVDs) %>% Biostrings::AAStringSet()
  
  if (appendExtremityCodes) {
    for (s in names(seqsOfRVDs)) {
      hitsByArray <- hitsByArraysLst[[s]]
      seqsOfRVDs[s] <- paste(ifelse(TALE_NtermDNAHMMName %in% as.character(hitsByArray$query_name), taleArrayStartAnchorCode, "XXXXX"), 
                             seqsOfRVDs[s], 
                             sep = rvdSep)
      seqsOfRVDs[s] <- paste(seqsOfRVDs[s], 
                             ifelse(TALE_CtermDNAHMMName %in% as.character(hitsByArray$query_name), taleArrayEndAnchorCode, "XXXXX"), 
                             sep = rvdSep)
    }
  }
  
  S4Vectors::mcols(hitsByArraysLst) <- merge(S4Vectors::mcols(hitsByArraysLst),
                                             data.frame(SeqOfRVD = seqsOfRVDs, aberrantRepeat = aberrantRepeat, arrayID = names(seqsOfRVDs)),
                                             by = "arrayID", 
                                             all.x = T)
  S4Vectors::mcols(hitsByArraysLst)$SeqOfRVD[is.na(S4Vectors::mcols(hitsByArraysLst)$SeqOfRVD)] <- ""
  
  #### Align N-term and C-term ####
  dnaPartFiles <- list.files(tmpAnnotaleDir, "TALE_DNA_parts.fasta", recursive = T, full.names = T)
  for (part in c("N-terminus", "C-terminus")) {
    allpart <- sapply(dnaPartFiles, function(p) {
      allpart <- Biostrings::readDNAStringSet(p, seek.first.rec = T)
      onepart <- allpart[grepl(part, names(allpart))]
      talRoi <- basename(dirname(p))
      names(onepart) <- talRoi
      return(onepart)
    }, simplify = "array", USE.NAMES = F) %>% Biostrings::DNAStringSetList() %>% unlist()
    dnaAlignment <- DECIPHER::AlignSeqs(allpart)
    DECIPHER::BrowseSeqs(dnaAlignment, htmlFile = file.path(outputDir, glue::glue("{part}DNAAlignment.html")), openURL = F, colWidth = 120)
  }
  
  aaPartFiles <- list.files(tmpAnnotaleDir, "TALE_Protein_parts.fasta", recursive = T, full.names = T)
  endsAA <- sapply(c("N-terminus", "C-terminus"), function(part) {
    allpart <- sapply(aaPartFiles, function(p) {
      allpart <- Biostrings::readAAStringSet(p, seek.first.rec = T)
      onepart <- allpart[grepl(part, names(allpart))]
      talRoi <- basename(dirname(p))
      names(onepart) <- talRoi
      return(onepart)
    }, simplify = "array", USE.NAMES = F) %>% Biostrings::AAStringSetList() %>% unlist()
    dnaAlignment <- DECIPHER::AlignSeqs(allpart)
    DECIPHER::BrowseSeqs(dnaAlignment, htmlFile = file.path(outputDir, glue::glue("{part}AAAlignment.html")), openURL = F, colWidth = 120)
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
  #####   Look at gaps between HitDomains on the same subject sequence to detect potential missed repeats   #####
  arraysBySeqlevelLst <- split(arraysGR, GenomicRanges::seqnames(arraysGR)) # group  arraysGR by seqlevel
  
  gaplengthBetweenHitDomains <- sapply(arraysBySeqlevelLst, function(x) {
    dists <- t(as.data.frame(GenomicRanges::distanceToNearest(x)))[3,]
  })
  gaplengthBetweenHitDomains <- unlist(gaplengthBetweenHitDomains)
  gaplengthBetweenHitDomains <- gaplengthBetweenHitDomains[!is.na(gaplengthBetweenHitDomains)]
  gaplengthBetweenHitDomainsbelow500 <- gaplengthBetweenHitDomains[gaplengthBetweenHitDomains <= 500]
  quartilesGapLength <- quantile(gaplengthBetweenHitDomainsbelow500,  probs = c(0.25, 0.50, 0.75))
  ##ggplot(data.frame(gapSize = gaplengthBetweenHitDomainsbelow500), aes(x=gapSize)) + geom_histogram(binwidth=10)
  
  
  #####   Write tabulated reports   #####
  ## Report on the hmmer hits. both  as a  tab-delimited  and a gff
  hitsReport <- lapply(as.list(hitsByArraysLst, use.names = TRUE),
                       function(gr) {
                         tibble::as_tibble(as.data.frame(gr))
                       }
  ) %>% dplyr::bind_rows(.id = "arrayID")
  readr::write_tsv(x = hitsReport, path = hitsReportFile)
  hitsReportToGFF(hitsReportFile) # saving to gff format
  
  readr::write_tsv(x = domainsReport, path = domainsReportFile)
  
  ##   Report with info on arrays, including the seq of RVD
  arrayReport <- as.data.frame(
    S4Vectors::mcols(hitsByArraysLst)[
      order(
        S4Vectors::mcols(hitsByArraysLst)$OriginalSubjectName,
        S4Vectors::mcols(hitsByArraysLst)$NumberOfHits
      ),]
  )
  readr::write_tsv(x = arrayReport, path = arrayReportFile)
  
  
  ## Write a gff with all collated
  allGR <- c(unlist(hitsByArraysLst),
             GenomicRanges::makeGRangesFromDataFrame(
               S4Vectors::mcols(hitsByArraysLst),
               seqnames.field= "OriginalSubjectName",
               keep.extra.columns= TRUE),
             if (exists("reducedOlapGr")) nhmmerOutputGRBeforeMerge else NULL
  )
  rtracklayer::export.gff3(allGR, affRangesGffFile)
  
  #####   Extract various DNA/AA sequences of interest and save to files   #####
  
  
  ## Write a fasta file of the seq of RVDs
  seqsOfRVDs <- Biostrings::BStringSet(S4Vectors::mcols(hitsByArraysLst)$SeqOfRVD)
  names(seqsOfRVDs) <- S4Vectors::mcols(hitsByArraysLst)$arrayID
  seqsOfRVDs <- seqsOfRVDs[!width(seqsOfRVDs) == 0]
  # seqsOfRVDs <- chartr("X", "U", seqsOfRVDs)
  # availableAA <- AA_ALPHABET[!AA_ALPHABET %in% uniqueLetters(seqsOfRVDs)]
  Biostrings::writeXStringSet(x = seqsOfRVDs, seqsOfRVDFile)
  
  #####   Generate info messages and log file about the analysis   #####
  
  ## counts of appearance of each RVD type (excluding N- and C- terms symbols) for the log file
  RVDtbl <- table(subset(unlist(hitsByArraysLst), query_name ==repeatDNAHMMName, drop = TRUE)$RVD)
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
    paste(Quote(minRatioOfGapForColMasking),":", minRatioOfGapForColMasking, sep = "\t"),
    paste(Quote(minDomainHitsPerSubjSeq),":", minDomainHitsPerSubjSeq, sep = "\t"),
    paste(Quote(mergeHits),":", mergeHits, sep = "\t"),
    # paste(Quote(repMsaMethod),":", repMsaMethod, sep = "\t"),
    paste(Quote(minGapWidth),":", minGapWidth, sep = "\t"),
    paste(Quote(minDomainHitsPerArrayForAssembl),":", minDomainHitsPerArrayForAssembl, sep = "\t"),
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
    paste("Total number of repeat HMM hits on the corresponding set of translated DNA hits:",
          sum(RVDtbl), sep = "\t"),
    
    paste("Total number of subject seqs with TALE motif hits after low hit number filtering:",
          length(GenomeInfoDb::seqlevels(arraysGR)), sep = "\t"),
    paste("Total number of distinct regions (repeat arrays) with adjacent TALE motifs :", nrow(arrayReport), sep = "\t"),
    paste("Total number of 'complete' arrays (with both N- and C-term flanking motifs):",
          sum(S4Vectors::mcols(hitsByArraysLst)$AllDomains),	sep = "\t"),
    
    paste("Total number of distinct types of RVD:", nrow(RVDtbl), sep = "\t"),
    
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
    
    "#*************************\n"
  )
  
  message(paste(txt, collapse ="\n"))
  logf <- file(analysisLogFile, open = "w")
  writeLines(text = txt, con = logf)
  close(logf)
  return(invisible(outputDir))
}
