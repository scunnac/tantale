
#' Search and report on the features of TALE protein domains potentially encoded
#' in subject DNA sequences
#'
#' This function is maintained in the package temporarily and should not be used
#' for other purpose than curiosity. You are better off using \link{tellTale}.
#'
#' \code{tellTaleLegacy}, the predecessor of the current \link{tellTale}} function
#' and has been primarily written to report on 'corrected' TALE RVD
#' sequences in noisy DNA sequences (suboptimally polished genomes assembly, raw
#' reads of long read sequencing technologies [eg PacBio, ONT]) that would
#' otherwise be missed by conventional tools (eg AnnoTALE).
#'
#' It works but is far from optimal for 'corrected' repeat CDS and \strong{it
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
#' @param repMsaMethod Value is a character string being either "decipher" or "hmmalign".
#'   This parameter defines which method is used to compute tale repeat CDS
#'   multiple alignments as a prerequisite for repeat CDS indel correction.
#'   hmmaling tends to remove a couple of nucleotide at the extremities of the
#'   repeats when they do not match the hmm profile. In contrast, the DECIPHER
#'   package function tends to leave spurious sequences at the end of repeat
#'   regions (notably for the last half repeat). In both cases, this can bias
#'   subsequent tale analysis based on DisTALE that relies on repeat sequences.
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
#' @return This functions has only side effects (writing files, mostly).
#'   However, if everything ran smoothly, it will invisibly return the path of
#'   the directory where output files were written.
#'
#' @export
tellTaleLegacy <- function(
  subjectFile,
  outputDir = getwd(),
  hmmFilesDir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T),
  minRatioOfGapForColMasking = 0.8,
  TALE_NtermDNAHitMinScore = 300,
  repeatDNAHitMinScore = 20,
  TALE_CtermDNAHitMinScore = 200,
  minDomainHitsPerSubjSeq = 4,
  mergeHits = TRUE,
  repMsaMethod = "decipher",
  minGapWidth = 35,
  minDomainHitsPerArrayForAssembl  = 5,
  taleArrayStartAnchorCode = "N-TERM",
  taleArrayEndAnchorCode = "C-TERM",
  appendExtremityCodes = TRUE,
  rvdSep = " ",
  hmmerpath = system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = T)
) {

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

  #####   Creating multiple alignments of the TALE repeat CDS hits with hmmalign or DECIPHER::AlignTranslation   #####
  ## Extract the DNA sequences of the repeat hits
  repeatSeqsSet <- Biostrings::DNAStringSet(subset(nhmmerOutputGR, query_name == repeatDNAHMMName)$seq)
  names(repeatSeqsSet) <- subset(nhmmerOutputGR, query_name == repeatDNAHMMName)$hitID

    multipleAlignment <- switch(repMsaMethod,
    hmmalign = {
      #for unclear reasons, Hmmalign removes a few residues at the
      # extremities of the seq HOW CAN WE FIX THAT? -> tried to make the hmm profile
      # for repeats more diverse... We will see how it goes...
      # Alternatively one could try multialigners in the msa or muscle package to
      # see how they compare to hmmalign
      # Alternatives to hmmaling:
      # Tried muscle -> not good
      # Tried the functions in the amazing DECIPHER package -> potentially promizing

      ## Save in a fasta file
      Biostrings::writeXStringSet(repeatSeqsSet, repeatDNASeqsFile)
      ## Run hmmalign
      runHmmalign(hmmerpath = hmmerpath,
                  hmmFile = repeatDNAHMMFile,
                  seqsFile = repeatDNASeqsFile,
                  alignOutFile = repeatsAlignOutFile
      )
      ## Read Multiple alignment of TALE repeat unit DNA seqs in R
      multipleAlignment <- Biostrings::readDNAMultipleAlignment(filepath = repeatsAlignOutFile,
                                                                format = "phylip")
    },
    decipher = {
      # !!!!!!!!!!!!!!!
      # DECIPHER::CorrectFrameshifts(repeatSeqsSet)
      # This function if it works well may allow to completely revise the paradigm
      # used in tellTale to correct repeat sequences: we would no longer need to rely on masking
      # As a bonus we may also be able to correct N and C term CDS...
      # To use it we need to have the AA sequences of the reference
      # domains used to build the hmm profiles.
      # In addition we can use DECIPHER::BrowseSeqs to write to output html files of
      # alignments for direct inspection!!! (use sf::path_ext_remove())
      DNA <- DECIPHER::AlignTranslation(repeatSeqsSet) # align the translation then reverse translate
      multipleAlignment <- Biostrings::DNAMultipleAlignment(DNA)
      #!!!! DECIPHER::CorrectFrameshifts(repeatSeqsSet)
      #!!!! DECIPHER::BrowseSeqs(DNA, highlight=1) # view the alignment
    }
  )

  ## Write output html files of alignments for direct inspection
  DECIPHER::BrowseSeqs(as(multipleAlignment, "DNAStringSet"),
                       htmlFile = paste0(gsub("^(.*)\\..*$", "\\1", repeatsAlignOutFile), ".html"),
                       openURL = FALSE)

  #######   Remove artefactual insertions in multiple alignment   #####
  ## Mask columns with artefactual nucleotide insertions relative to HMM
  ## (hopefully caused by sequencing errors and not genuine polymorphism)
  cleanMultipleAlignment <- Biostrings::maskGaps(multipleAlignment,
                                                 min.fraction = minRatioOfGapForColMasking,
                                                 min.block.width = 1
                                                 )

  #####   Compute counts of indels relative to a "canonical" repeat CDS   ####
  ## WOULD BE NICE TO KEEP TRACK OF THE NUMBER OF GENUINE IN/DEL FOR EACH REPEAT
  ## DEVIATING FROM THE CANONICAL REPEAT to include in final report. For each
  ## repeat sequ in the alignment compute the number of positions with a
  ## nucleotide (rather than a gap) in the masked columns This may be done with
  ## maskedncol(x): Returns the number of masked aligned characters in x? This
  ## will give the number of insertions for each seqs => store in a variable and
  ## add to the hit report df

  ## count number of insertion
  insCol <- cleanMultipleAlignment
  Biostrings::colmask(insCol, append = "replace", invert = TRUE) <- Biostrings::colmask(insCol)
  insColStringSet <- as(insCol, "DNAStringSet")
  insColNo <- BiocGenerics::width(insColStringSet) - Biostrings::letterFrequency(insColStringSet, letters = "-")

  ## count number of deletion
  # REMEMBER that in tests, sometimes the number of deletion is overestimated
  # because for unclear reasons, Hmmalign removes a few residues at the
  # extremities of the seq HOW CAN WE FIX THAT? -> tried to make the hmm profile
  # for repeats more diverse... We will see how it goes...
  #
  # Alternatively one could try multialigners in the msa or muscle package to
  # see how they compare to hmmalign

  delColStringSet <- as(cleanMultipleAlignment, "DNAStringSet")
  delColNo <- Biostrings::letterFrequency(delColStringSet, letters = "-")
  indelReport <- data.frame(hitID = names(as.character(multipleAlignment)),
                            corrected_ins_count = as.vector(insColNo),
                            corrected_del_count = as.vector(delColNo), #
                            stringsAsFactors = FALSE,
                            check.names = FALSE)

  #####   Translate the repeat CDS alignment and write to file   #####
  cleanDNASeqs <- as(cleanMultipleAlignment, "DNAStringSet")
  cleanDNASeqs <- Biostrings::chartr("-", "N", cleanDNASeqs) # Substitute gaps for N
  Biostrings::writeXStringSet(cleanDNASeqs, correctedRepeatsAlignOutFile) # Save aligned TALE repeat unit DNA sequences to file

  translatedCleanAlign <- Biostrings::translate(cleanDNASeqs, if.fuzzy.codon = "solve", no.init.codon = T) # Translate with Biostrings function

  ## !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!!
  ## Implement a mechanism to keep track of repeats that contain an inframe stop codon
  ## !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!!


  translatedCleanAlign <- Biostrings::chartr("*", "X", translatedCleanAlign) # It looks like if stop codon are there, hmmalign do not align them with other residues and create gaps

  Biostrings::writeXStringSet(translatedCleanAlign, translatedCorrectedRepeatsAlignOutFile) # Save TALE repeat unit infered AA sequences to file
  DECIPHER::BrowseSeqs(as(translatedCleanAlign, "AAStringSet"),
                       htmlFile = paste0(gsub("^(.*)\\..*$", "\\1", translatedCorrectedRepeatsAlignOutFile), ".html"),
                       openURL = FALSE)


  #####    Perform RVD AA pattern search on cleaned TALE repeat AA seqs   #####
  ## As opposed to extracting positions 12-13 from AArepeatsAlignOutFile which seems like a reasonable and simple
  ## thing to do, this extra step of running hmmer on the translated candidate repeat sequences
  ## was included in the pipeline in order to guard against
  ## situations where, for some reasons, some DNA sequences in the degaped msa had their frame shifted
  ## or if their translation product did not look like the AA sequence of a TALE repeat.
  ## Not very clearly worded...
  runHmmerSearch(hmmerpath = hmmerpath, subjectFile = translatedCorrectedRepeatsAlignOutFile,
                 hmmFile = repeatAAHMMFile,
                 searchTblOutFile = repeatAaSearchOutFile,
                 humReadableOutFile = file.path(outputDir, "hmmsearchHumanReadableOutputOfLastRun.txt"))

  #####   Creating multiple repeat AA alignments with hmmalign   ###
  runHmmalign(hmmerpath = hmmerpath, hmmFile = repeatAAHMMFile,
              seqsFile = translatedCorrectedRepeatsAlignOutFile,
              alignOutFile = repeatsAaAlignOutFile)


  #####   Extracting RVDs   #####
  ## Read Multiple alignment of TALE repeat unit protein seqs in R
  multipleProtAlignment <- Biostrings::readAAMultipleAlignment(filepath = repeatsAaAlignOutFile, format = "phylip") #"fasta" (the default), stockholm, or "clustal" "phylip"
  DECIPHER::BrowseSeqs(as(multipleProtAlignment, "AAStringSet"),
                       htmlFile = paste0(gsub("^(.*)\\..*$", "\\1", repeatsAaAlignOutFile), ".html"),
                       openURL = FALSE)


  RVDAlignment <- multipleProtAlignment
  Biostrings::colmask(RVDAlignment) <- IRanges::IRanges(start=1,end=11) # Mask all AA before the first RVD residue
  Biostrings::colmask(RVDAlignment, append = "union") <- IRanges::IRanges(start=14, end = ncol(multipleProtAlignment)) # Mask all AA after the second RVD residue
  RVDs <- as(RVDAlignment, "AAStringSet") # Get the set of RVDs
  RVDs <- Biostrings::chartr("-", "X", RVDs) # Recode the "gaps" or incomplete repeats sequences with * rather than - which messes up the sequences of RVDs where "-" is used as a separator
  RVDs <- as.character(RVDs)

  ## Could be done:
  ## susbsitute spurious/artefactual RVDs for "XX"
  ## Procedure: generate a list of all the descirbed RVDs and convert all others to "XX"
  ## Append a column for each the substituted and the original vector of RVDs
  ## CONS: Would masks genuine novel RVD...

  #####   Storing info about RVDs and potential artefactual indels in repeat CDS   ####
  repeatsReport <-
    merge(
      data.frame(
        "hitID" = names(RVDs),
        "RVD" = RVDs,
        row.names = NULL, stringsAsFactors = FALSE,
        check.rows = T),
    indelReport,
    by = "hitID", all = TRUE, sort = FALSE) %>%
    merge(
      data.frame(
      hitID = names(cleanDNASeqs),
      corrected_seq = as.character(cleanDNASeqs)
      ),
    by = "hitID", all = TRUE, sort = FALSE) %>%
    merge(
      data.frame(
      hitID = names(translatedCleanAlign),
      corrected_AA_seq = as.character(translatedCleanAlign)
      ),
    by = "hitID", all = TRUE, sort = FALSE
    )

  repeatsReport %<>%
    merge(x =  S4Vectors::mcols(nhmmerOutputGR),
          y = .,
          by = "hitID",
          all.x = TRUE, sort = FALSE
    )
  row.names(repeatsReport) <- repeatsReport$hitID
  S4Vectors::mcols(nhmmerOutputGR) <- repeatsReport[names(nhmmerOutputGR),]

  ## Assign an anchor code to the RVD column for N-Term and C-Term motif hits
  S4Vectors::mcols(nhmmerOutputGR)[nhmmerOutputGR$query_name == TALE_NtermDNAHMMName, "RVD"] <- taleArrayStartAnchorCode
  S4Vectors::mcols(nhmmerOutputGR)[nhmmerOutputGR$query_name == TALE_CtermDNAHMMName, "RVD"] <- taleArrayEndAnchorCode
  #nhmmerOutputGR$RVD <- as.factor(nhmmerOutputGR$RVD)





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
    SeqOfRVD =   stitchSeqs(hitsByArraysLst, "RVD", sep = rvdSep, onlyRepeats = !appendExtremityCodes),
    ArraySeq = BSgenome::getSeq(subjectDNASequences, arraysGR),
    CorrectedRepArrayDnaSeq = Biostrings::DNAStringSet(
      stitchSeqs(hitsByArraysLst, "corrected_seq", sep = "", onlyRepeats = TRUE)
      ),
    CorrectedRepArrayAaSeq = Biostrings::AAStringSet(
      stitchSeqs(hitsByArraysLst, "corrected_AA_seq", sep = "", onlyRepeats = TRUE)
      ),
    SelectedForAssembly = S4Vectors::elementNROWS(hitsByArraysLst) >= minDomainHitsPerArrayForAssembl,
    AllDomains = sapply(hitsByArraysLst,
                           function(x) {
                             all(
                               c(TALE_NtermDNAHMMName, repeatDNAHMMName,
                                 TALE_CtermDNAHMMName) %in% x$query_name
                             )
                           }
    )
  )
  #str(S4Vectors::mcols(hitsByArraysLst))
  #####   Detects longest ORF in (complete) arrays   #####
  ## Extract the genomic sequence of ("full length") arrays +- 5bp on the borders
  ## run systemPipeR::predORF()
  ## ask if the longest orf GRanges on the plus strand (perfectly) matches with the complete array GRanges
  ## Record this info in the arrayReport object
  ## Export the dna seq of the longest orf as the tal sequences
  completeArraysGR <- arraysGR #subset(arraysGR, S4Vectors::mcols(hitsByArraysLst)$AllDomains) #
  extdCompleteArraysGR <- GenomicRanges::resize(completeArraysGR,
                        width = GenomicRanges::width(completeArraysGR) + 200,
                        fix="start", ignore.strand=FALSE) %>%
    GenomicRanges::resize(., width = GenomicRanges::width(.) + 3,
                          fix="end", ignore.strand=FALSE)
  extdCompleteArraysSeqs <- BSgenome::getSeq(subjectDNASequences, extdCompleteArraysGR)
  orfs <- systemPipeR::predORF(x = extdCompleteArraysSeqs,
                               n=1, type = "gr", mode = "ORF", strand = "sense")


  fullTalCds <- BSgenome::getSeq(extdCompleteArraysSeqs, orfs)
  names(fullTalCds) <- as.character(GenomicRanges::seqnames(orfs))

  ## Merge with arrays metadata in mcols(hitsByArraysLst)
  moreInfo <- merge(
    S4Vectors::mcols(hitsByArraysLst),
    data.frame(arrayID = names(fullTalCds),
               LongestOrfLength = Biostrings::nchar(fullTalCds),
               OrfCovOverArrayLength = round(100 * Biostrings::nchar(fullTalCds)/GenomicRanges::width(completeArraysGR[names(fullTalCds)])),
               LongestORFSeq = fullTalCds),
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

  ## extract N-term CDS from hmmer output hits and translate them
  NtermHitsByArray <- lapply(hitsByArraysLst, function(x) {
    x[x$query_name == TALE_NtermDNAHMMName]}) # extract N-term hits
  NtermHitsByArray <-as(NtermHitsByArray, "CompressedGRangesList")
  fullNtermCDS <- range(NtermHitsByArray) # remove gaps, overlaps to get full length N-term
  fullNtermCDS  <- as.data.frame(fullNtermCDS )
  colnames(fullNtermCDS) <- c("group", "hitID", "target_name", "start", "end", "width", "strand")
  NtermSeqsSet <- extractSeqsfromHits(fullNtermCDS, subjectDNASequences) # extract N-term hit DNA sequences
  fullNtermAAseq <- Biostrings::translate(NtermSeqsSet, if.fuzzy.codon = "solve")
  fullNtermAAseq <- Biostrings::chartr("*", "X", fullNtermAAseq)
  Biostrings::writeXStringSet(fullNtermAAseq, fullNtermAAseqFile)

  ## extract C-term CDS from hmmer output hits and translate them
  CtermHitsByArray <- lapply(hitsByArraysLst, function(x) {
    x[x$query_name == TALE_CtermDNAHMMName]})
  CtermHitsByArray <-as(CtermHitsByArray, "CompressedGRangesList")

  fullCtermCDS <- range(CtermHitsByArray) # remove gaps, overlaps to get full length C-term
  fullCtermCDS <- as.data.frame(fullCtermCDS)
  colnames(fullCtermCDS) <- c("group", "hitID", "target_name", "start", "end", "width", "strand")
  CtermSeqsSet <- extractSeqsfromHits(fullCtermCDS, subjectDNASequences)
  fullCtermAAseq <- Biostrings::translate(CtermSeqsSet, if.fuzzy.codon = "solve", no.init.codon = TRUE)
  fullCtermAAseq <- Biostrings::chartr("*", "X", fullCtermAAseq)
  Biostrings::writeXStringSet(fullCtermAAseq, fullCtermAAseqFile)

  ## Write a fasta file of central domain "corrected" CDS for each array
  cleanDNAarraySeqs <- Biostrings::DNAStringSet(S4Vectors::mcols(hitsByArraysLst)$CorrectedRepArrayDnaSeq)
  names(cleanDNAarraySeqs) <- S4Vectors::mcols(hitsByArraysLst)$arrayID
  Biostrings::writeXStringSet(cleanDNAarraySeqs, correctedArrayDnaSeqsFile)

  ## Write a fasta file for protein repeats sequences
  cleanAAarraySeqs <- Biostrings::AAStringSet(S4Vectors::mcols(hitsByArraysLst)$CorrectedRepArrayAaSeq)
  names(cleanAAarraySeqs) <- S4Vectors::mcols(hitsByArraysLst)$arrayID
  Biostrings::writeXStringSet(cleanAAarraySeqs, correctedArrayAaSeqsFile)

  ## Write a fasta file of the seq of RVDs
  seqsOfRVDs <- Biostrings::BStringSet(S4Vectors::mcols(hitsByArraysLst)$SeqOfRVD)
  names(seqsOfRVDs) <- S4Vectors::mcols(hitsByArraysLst)$arrayID
  # seqsOfRVDs <- chartr("X", "U", seqsOfRVDs)
  # availableAA <- AA_ALPHABET[!AA_ALPHABET %in% uniqueLetters(seqsOfRVDs)]
  Biostrings::writeXStringSet(x = seqsOfRVDs, seqsOfRVDFile)

  ## Write a fasta file of the ORF of the TALE arrays
  arraysWithOrf <- subset(S4Vectors::mcols(hitsByArraysLst), !is.na(LongestORFSeq))
  arrayOrfsDnaSeqs <- Biostrings::DNAStringSet(arraysWithOrf$LongestORFSeq)
  names(arrayOrfsDnaSeqs) <- arraysWithOrf$arrayID
  Biostrings::writeXStringSet(arrayOrfsDnaSeqs, arrayOrfsSeqFile)

  #####   Generate info messages and log file about the analysis   #####

  ## counts of appearance of each RVD type (excluding N- and C- terms symbols) for the log file
  RVDtbl <- table(subset(unlist(hitsByArraysLst), query_name ==repeatDNAHMMName, drop = TRUE)$RVD)
  ## Total count of repeat CDS after filtering for uniformative subject seqs for the log file
  numberOfRepeatHitsAfterFiltering <- length(subset(unlist(hitsByArraysLst), query_name == repeatDNAHMMName))
  ## Distribution of the number of hits per array
  countsHitsByArrayDistri <- summary(S4Vectors::mcols(hitsByArraysLst)$NumberOfHits)
  ## Number of domains in arrays that display all domain types
  completeArrayLengths <- subset(S4Vectors::mcols(hitsByArraysLst), AllDomains)$NumberOfHits


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
    paste(Quote(repMsaMethod),":", repMsaMethod, sep = "\t"),
    paste(Quote(minGapWidth),":", minGapWidth, sep = "\t"),
    paste(Quote(minDomainHitsPerArrayForAssembl),":", minDomainHitsPerArrayForAssembl, sep = "\t"),
    paste(Quote(taleArrayStartAnchorCode),":", taleArrayStartAnchorCode, sep = "\t"),
    paste(Quote(taleArrayEndAnchorCode),":", taleArrayEndAnchorCode, sep = "\t"),

    "#__________Summary measures of TALE search outcome__________",
    paste("Number of analysed subject sequences :", length(subjectDNASequences), sep = "\t"),
    paste("Total number of TALE repeat DNA coding sequence motif hits found with the nhmmer approach:",
          numberOfRepeatHitsAfterFiltering, sep = "\t"),
    paste("Number of column masked in the repeat unit hits alignment:",
          Biostrings::maskedncol(cleanMultipleAlignment), sep = "\t"),
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
    paste("Length of the longest 'complete' array:", max(completeArrayLengths),	sep = "\t"),
    paste("Length of the shortest 'complete' array:", min(completeArrayLengths),	sep = "\t"),

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
