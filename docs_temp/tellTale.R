#!/usr/bin/Rscript --vanilla

#### Loading required packages #####
suppressPackageStartupMessages(library("optparse"))

#### COMMAND LINE ARGUMENTS PARSING ######
option_list <- list(
    make_option(c("-S", "--subjectFile"),
                type = "character",
                default = "/home/baotram/tal/hmmprofile/control_Xo_genomes_for_annotale/PXO99A.fasta",
                help = "Path to the subject file."),
    make_option(c("-O", "--outputDir"),
                type = "character",
                default = "/home/baotram/tal/rvdArrayReader_test/tellTale_test/PXO99A",
                help = "Path to the output directory."),
    make_option(c("-H", "--hmmFileDir"),
                type = "character",
                default = "/home/baotram/tal/taleRepeatArraySolverPipeline_v1/inputHMMFiles",
                help = "Path to the folder containing input HMM files."),
    make_option(c("-R", "--minRatioOfGapForColMasking"),
                type = "numeric",
                default = 0.8,
                help = "Minimum ratio of a gap appeared in columns of the table tale repeat CDS aligment that allows masking (between 0 and 1). The value is 0.8 by default."),
    make_option(c("-N", "--TALE_NtermDNAHitMinScore"),
                type = "numeric",
                default = 40,
                help = "nhmer score cut-off value for N-terminus hits. The value is 40 by default."),
    make_option(c("-r", "--repeatDNAHitMinScore"),
                type = "numeric",
                default = 20,
                help = "nhmer score cut-off value for repeat DNA hits. The value is 20 by default."),
    make_option(c("-C", "--TALE_CtermDNAHitMinScore"),
                type = "numeric",
                default = 30,
                help = "nhmer score cut-off value for C-terminus hits. The value is 30 by default."),
    make_option(c("-D", "--minDomainHitsPerSubjSeq"),
                type = "numeric",
                default = 4,
                help = "Minimum number of nhmmer hits for a subject sequence to be reported as having TALE diagnostic regions. The value is 4 by default."),
    make_option(c("-W", "--minGapWidth"),
                type = "numeric",
                default = 35,
                help = "Minimum gap between two tale domain hits for them to be considered 'contiguous' and grouped in the same array. The value is 35 by default."),
    make_option(c("-A", "--minDomainHitsPerArrayForAssembl"),
                type = "numeric",
                default = 5,
                help = "Minimum number of repeat in an array for its seq of RVD to be considered for assembly. The value is 5 by default."),
    make_option(c("-B", "--taleArrayStartAnchorCode"),
                type = "character",
                default = "BBB",
                help = "Code for N-terminus. The code is 'BBB' by default"),
    make_option(c("-Z", "--taleArrayEndAnchorCode"),
                type = "character",
                default = "ZZZ",
                help = "Code for C-terminus. The code is 'ZZZ' by default."),
    make_option(c("-P", "--hmmerpath"),
                type = "character",
                default = "/home/baotram/miniconda3/envs/telltales/bin/",
                help = "Pathway to execute hmmer.")
)

##### RETRIEVEING PARAMS from optparse #####
myArgs <- parse_args(
  OptionParser(usage = "%prog [-S SubjectFile] [-O OutputDir] [-H HMMFileDir] [options]", option_list = option_list,
               description = "Run the RVD Reader ... something")
)


##### FUNCTION DEFINITION #####

# A function that checks that hmmer is accessible"
checkHMMER <- function(hmmerpath) {
  cmd <- paste0(hmmerpath, "hmmsearch -h | grep \"^#\"")
  if (system(command = cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)) {
    stop("HMMER is not in PATH. Follow instructions at http://hmmer.org/documentation.html to install it.")
  } else {
     out <- system(command = cmd,intern = TRUE)
      cat(out[2:3], sep = "\n")
    }
}


writeHMMFile <- function(hmmerpath, alignmentFile, HMMOutFile) {
  buildCmd <- paste(paste0(hmmerpath,"hmmbuild"),
                    HMMOutFile,
                    alignmentFile,
                    sep = " ")
  commandOut <- system(command = buildCmd, ignore.stderr = FALSE, intern = TRUE)
  return(commandOut)
}


runHmmerSearch <- function(hmmerpath, subjectFile, hmmFile, searchTblOutFile, humReadableOutFile) {
  searchCmd <- paste(paste0(hmmerpath, "hmmsearch"),
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


runNhmmerSearch <-  function(hmmerpath, subjectFile, hmmFile, searchTblOutFile, humReadableOutFile) {
  searchCmd <- paste(paste0(hmmerpath, "nhmmer"),
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



runHmmalign <- function(hmmerpath, hmmFile, seqsFile, alignOutFile) {
  alignCmd <- paste(paste0(hmmerpath, "hmmalign"),
                    "--outformat Phylip", #Stockholm, SELEX, Clustal, Phylip, Pfam, A2M, PSIBLAST.
                    "--trim",
                    hmmFile,
                    seqsFile,
                    ">", alignOutFile,
                    sep = " "
  )
  system(command = alignCmd, ignore.stderr = FALSE, intern = TRUE)
}

hitsReportToGFF <- function(f = "hitsReport.csv") {
  # Convert the info contained in a HitReport file into a GFF file for display by
  # a genome viewer.
  # The f parameter corresponds to the path to a hitsReport file.

  # Read the file as a data.frame
  hitsReport <- read.delim(f)

  # make sure genome coordinates are in increasing order.
  isIncreasing <- hitsReport$alifrom <= hitsReport$ali_to
  st <- ifelse(isIncreasing, hitsReport$alifrom, hitsReport$ali_to)
  en <- ifelse(!isIncreasing, hitsReport$alifrom, hitsReport$ali_to)

  # Create a GenomicRange that will be converted.
  hitsGR <-
    GenomicRanges::GRanges(
      seqnames = hitsReport$target_name,
      ranges = IRanges(start = st, end = en, names = hitsReport$hitID),
      strand = hitsReport$strand,
      arrayID = hitsReport$arrayID,
      TALEDomainType = hitsReport$query_name,
      hmmfrom = hitsReport$hmmfrom,
      hmm_to = hitsReport$hmm_to,
      Evalue = hitsReport$Evalue,
      RVD = hitsReport$RVD
    )

  # Write a gff3 file to disk with this info.
  rtracklayer::export.gff3(hitsGR,
                           con = file.path(dirname(f), paste0(sub("\\..*$", "", basename(f)), ".gff"))
  )

}

extractSeqsfromHits <- function(nhmmerTabularOutputSelect, DNAsequences){
  repeatSeqsSetList <- mapply(function(hitID, start, end, strand, subjectID, sequences) {
    seq <- subseq(sequences[subjectID], start, end)
    if (strand == "-") {seq <- reverseComplement(seq)}
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


rvdArrayReader <- function(
  subjectFile,
  outputDir,
  hmmFilesDir,
  minRatioOfGapForColMasking, # columns of the tale repeat CDS alignment that contain a gap in a fraction of sequence higher than this value (betwen 0 and 1) will be masked from the alignment when translating the DNA sequences to protein.
  TALE_NtermDNAHitMinScore, # nhmmer score cut_off value
  repeatDNAHitMinScore, # nhmmer score cut_off value
  TALE_CtermDNAHitMinScore, # nhmmer score cut_off value
  minDomainHitsPerSubjSeq, # Minimum number of nhmmer hits for a subject sequence to be reported as having TALE diagnostic regions. This is a way to simplify output a little by getting ride of uninformative sequences
  minGapWidth, # minimum gap between two tale domain hits for them to be considered "contiguous" and grouped in the same array.
  minDomainHitsPerArrayForAssembl, # Minimum number of repeat in an array for its seq of RVD to be considered for assembly. This is a way to get ride of sequences that are too short reasonably be of any help for assembly
  taleArrayStartAnchorCode,
  taleArrayEndAnchorCode,
  hmmerpath
) {
###############################################################################
# rvdArrayReader
# Detecting and translating the sequences of RVD in
# TALE proteins repeat domains encoded in noisy DNA sequences
###############################################################################
suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(library(XVector))

# source(file = file.path("RVD_seq_elucidation_library.R"))

# Parameter consistency check for the user.
if (minDomainHitsPerSubjSeq > minDomainHitsPerArrayForAssembl) {
  warning("The value of the minDomainHits parameter for subject sequence selection
						is more stringent than the value of the minNumberOfDomainHitsForAssembl!?
						Is this really what you want to do?")
}


##########################################################
## Full paths of input HMM files for TALE domains (DNA and AA)
##########################################################

TALE_NtermDNAHMMFile <- file.path(hmmFilesDir, "TALE_N-term_CDS_diagnostic_region_hmmfile.out")

repeatDNAHMMFile <- file.path(hmmFilesDir, "TalC_RVDs_CDS_hmmfile.out")

TALE_CtermDNAHMMFile <- file.path(hmmFilesDir, "TALE_C-term_CDS_diagnostic_region_hmmfile.out")


DNAHMMFiles <- c(TALE_NtermDNAHMMFile, repeatDNAHMMFile, TALE_CtermDNAHMMFile)

mergedDNAHMMFile <- file.path(outputDir, "TALE_CDS_all_diagnostic_regions_hmmfile.out")

repeatAAHMMFile <- file.path(hmmFilesDir, "ALL_CONTROL_TAL_CDS_RVD_AA_hmmfile.out")


##########################################################
## Concatenate HMM files for TALE DNA motifs and parse HMM names
##########################################################

hmmslines <- llply(DNAHMMFiles, function(x) txt <- readLines(con = x))

names(DNAHMMFiles) <- llply(hmmslines, function(x) {
  hmmName <- grep("NAME", x, perl = TRUE, value = TRUE)
  hmmName <- unlist(strsplit(hmmName, split = "\\s+"))
  if (length(hmmName) != 2) stop("One or several profile HMM have a name with withe spaces. Please remove them in the file at the Tag 'NAME'")
  hmmName <- hmmName[2]

}
)

DNAHMMNames <- names(DNAHMMFiles)
TALE_NtermDNAHMMName <- names(DNAHMMFiles)[1]
repeatDNAHMMName <- names(DNAHMMFiles)[2]
TALE_CtermDNAHMMName <- names(DNAHMMFiles)[3]

writeLines(text = unlist(hmmslines), con = mergedDNAHMMFile)


##########################################################
## Perform TALE repeat unit DNA pattern search
##########################################################
## Use nhmmer from HMMER for searching a file with TALE CDS for RVDs

searchOutFile <- file.path(outputDir, "hmmerSearchOut.txt")

runNhmmerSearch(hmmerpath = hmmerpath,
                subjectFile = subjectFile,
                hmmFile = mergedDNAHMMFile,
                searchTblOutFile = searchOutFile,
                humReadableOutFile = file.path(outputDir, "nhmmerHumanReadableOutputOfLastRun.txt"))


##########################################################
## Load, process, filter TALE repeat unit CDS HMMER hit results
##########################################################

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

####!!!!!!!!!NOT TESTED¡¡¡¡¡¡¡¡¡
## sorting data frame
##	newOrdering <- with(nhmmerOutput,
##			order(target_name, strand, start)
##	)
##	nhmmerOutput <- nhmmerOutput[newOrdering,]
####!!!!!!!!!NOT TESTED¡¡¡¡¡¡¡¡¡


## Add a TALE hitID column
nhmmerTabularOutput$hitID <- paste("DOM", sprintf("%05.0f", 1:nrow(nhmmerTabularOutput)), sep="_")

## Trick to re-order positions in an increasing order to satisfy IRanges() in preparation of creating a GRanges
nhmmerTabularOutput[,c("start", "end")] <- adply(.data = nhmmerTabularOutput[,c("alifrom", "ali_to")], .margins = 1, .fun = c(min, max))[,-(1:2)]
rownames(nhmmerTabularOutput) <- nhmmerTabularOutput$hitID


##########################################################
## Extracting the sequences of the TALE repeat CDS hits and save in a fasta file
##########################################################

## Load in R the DNA sequences that are queried for RVD CDS
subjectDNASequences <- readDNAStringSet(filepath = subjectFile)

## Select only hits corresponding to a repeat CDS in the nhmmer tabular output
nhmmerTabularOutputRepeatsOnly <- subset(nhmmerTabularOutput, query_name == repeatDNAHMMName)

## Extract the DNA sequences of the repeat hits
repeatSeqsSet <- extractSeqsfromHits(nhmmerTabularOutputRepeatsOnly, subjectDNASequences)

## Save in a fasta file
repeatDNASeqsFile <- file.path(outputDir, "hmmerRepeatHitsDNASeqs.fas")
writeXStringSet(repeatSeqsSet, repeatDNASeqsFile)


##########################################################
## Creating multiple alignments with hmmalign
##########################################################

alignOutFile <- file.path(outputDir, "hmmerAlignOut.txt")

runHmmalign(hmmerpath = hmmerpath, hmmFile = repeatDNAHMMFile, seqsFile = repeatDNASeqsFile, alignOutFile = alignOutFile)


##########################################################
## Remove artefactual insertions in multiple alignment
##########################################################

## Read Multiple alignment of TALE repeat unit DNA seqs in R
multipleAlignment <- readDNAMultipleAlignment(filepath = alignOutFile, format = "phylip") #"fasta" (the default), stockholm, or "clustal" "phylip"

## Mask columns with artefactual nucleotide insertions relative to HMM
## (hopefully caused by sequencing errors and not genuine polymorphism)
cleanMultipleAlignment <- maskGaps(multipleAlignment, min.fraction = minRatioOfGapForColMasking, min.block.width = 1)

## WOULD BE NICE TO KEEP TRACK OF THE NUMBER OF GENUINE IN/DEL FOR EACH REPEAT DEVIATING FROM THE CANONICAL REPEAT to include in final report
## For each repeat sequ in the alignment compute the number
## of positions with a nucleotide (rather than a gap) in the masked columns
## This may be done with maskedncol(x): Returns the number of masked aligned characters in x?
## This will give the number of insertions for each seqs => store in a variable and add to the hit report df

## count number of insertion
insCol <- cleanMultipleAlignment
colmask(insCol, append = "replace", invert = TRUE) <- colmask(insCol)
insColStringSet <- as(insCol, "DNAStringSet")
insColNo <- width(insColStringSet) - letterFrequency(insColStringSet, letters = "-")

## count number of deletion
delColStringSet <- as(cleanMultipleAlignment, "DNAStringSet")
delColNo <- letterFrequency(delColStringSet, letters = "-")

# indelReport <- as.data.frame(cbind(names(as.character(multipleAlignment)), insColNo, delColNo))
indelReport <- as.data.frame(cbind(names(as.character(multipleAlignment)), insColNo, delColNo))
colnames(indelReport) <- c("hitID", "number_of_insertion", "number_of_deletion")


##########################################################
## Translate the resulting alignment and write to file
##########################################################

cleanDNASeqs <- as(cleanMultipleAlignment, "DNAStringSet") # Convert DNAMultipleAlignment object into a DNAStringSet

## Now lets get the number of dels:
## count the number of "-" in each sequence  => store in a variable and add to the hit report df

cleanDNASeqs <- chartr("-", "N", cleanDNASeqs) # Substitute gaps for N



cleanedAlignOutFile <- file.path(outputDir, "hmmerCleanedAlignOut.fas")
writeXStringSet(cleanDNASeqs, cleanedAlignOutFile) # Save aligned TALE repeat unit DNA sequences to file

translatedCleanAlign <- translate(cleanDNASeqs, if.fuzzy.codon = "solve") # Translate with Biostrings function

## Implement a mechanism to keep track of repeats that contain a stop codon

translatedCleanAlign <- chartr("*", "X", translatedCleanAlign) # It looks like if stop codon are there, hmmalign do not align them with other residues and create gaps

translatedCleanAlignOutFile <- file.path(outputDir, "translatedCleanRepeatsAlignment.pep")
writeXStringSet(translatedCleanAlign, translatedCleanAlignOutFile) # Save TALE repeat unit infered AA sequences to file




##########################################################
## Perform RVD AA pattern search on cleaned TALE repeat unit DNA sequences after translation
##########################################################
## As opposed to extracting positions 12-13 from AAAlignOutFile which seems like a reasonable and simple
## thing to do, this extra step of running hmmer on the translated candidate repeat sequences
## was included in the pipeline in order to guard against (filter out)
## situations where, for some reasons, some DNA sequences in the degaped msa had their frame shifted
## or if their translation product did not look like the AA sequence of a TALE repeat.
## Not very clearly worded...

repeatAAsearchOutFile <- file.path(outputDir, "hmmerAASearchOut.txt")

runHmmerSearch(hmmerpath = hmmerpath, subjectFile = translatedCleanAlignOutFile,
               hmmFile = repeatAAHMMFile,
               searchTblOutFile = repeatAAsearchOutFile,
               humReadableOutFile = file.path(outputDir, "hmmsearchHumanReadableOutputOfLastRun.txt"))

##########################################################
## Creating multiple alignments with hmmalign
##########################################################

AAHitsAlignOutFile <- file.path(outputDir, "hmmerAAAlignOut.txt")

runHmmalign(hmmerpath = hmmerpath, hmmFile = repeatAAHMMFile,
            seqsFile = translatedCleanAlignOutFile,
            alignOutFile = AAHitsAlignOutFile)

##########################################################
## Extract exact RVDs
##########################################################

## Read Multiple alignment of TALE repeat unit protein seqs in R
multipleProtAlignment <- readAAMultipleAlignment(filepath = AAHitsAlignOutFile, format = "phylip") #"fasta" (the default), stockholm, or "clustal" "phylip"


RVDAlignment <- multipleProtAlignment
colmask(RVDAlignment) <- IRanges(start=1,end=11) # Mask all AA before the first RVD residue
colmask(RVDAlignment, append = "union") <- IRanges(start=14, end = ncol(multipleProtAlignment)) # Mask all AA after the second RVD residue
RVDs <- as(RVDAlignment, "AAStringSet") # Get the set of RVDs
RVDs <- chartr("-", "X", RVDs) # Recode the "gaps" or incomplete repeats sequences with X rather than - which messes up the sequences of RVDs where "-" is used as a separator
RVDs <- as.character(RVDs)

## Could be done:
## susbsitute spurious/artefactual RVDs for "XX"
## Procedure: generate a list of all the descirbed RVDs and convert all others to "XX"
## Append a column for each the substituted and the original vector of RVDs


## Add the infered diresidues to the dataframe
nhmmerTabularOutput <- merge(nhmmerTabularOutput, cbind(RVD = RVDs, hitID = names(RVDs)), by = "hitID", all.x = TRUE, sort = FALSE)
nhmmerTabularOutput$RVD <- as.character(nhmmerTabularOutput$RVD)

## Assign an anchor code to the RVD column for N-Term and C-Term motif hits
nhmmerTabularOutput[nhmmerTabularOutput$query_name == TALE_NtermDNAHMMName, "RVD"] <- taleArrayStartAnchorCode
nhmmerTabularOutput[nhmmerTabularOutput$query_name == TALE_CtermDNAHMMName, "RVD"] <- taleArrayEndAnchorCode
nhmmerTabularOutput$RVD <- as.factor(nhmmerTabularOutput$RVD)

## counts of appearance of each RVD type (excluding anchors)
RVDtbl <- table(as.vector(subset(nhmmerTabularOutput, query_name ==repeatDNAHMMName , select = RVD, drop = TRUE)))

## Store this info in a variable to be used when building the log file
numberOfRepeatHitsBeforeFiltering <- nrow(subset(nhmmerTabularOutput, query_name == repeatDNAHMMName))


## outputing some info in stder
##message(paste("The hmmer approach found a total of",
##				nrow(subset(nhmmerTabularOutput, query_name == repeatDNAHMMName)),
##						"repeats DNA coding sequences in",
##			length(subjectDNASequences), "subject sequences.\n",
##			"Performing a search with an amino acid repeat HMM on this set of translated DNA hits retrieved a total of",
##			sum(RVDtbl), "RVDs that distribute as follow:\n", sep =" ")
##
##a_ply(as.data.frame(RVDtbl), 1, .fun = function(x) message(paste(x$Var1, x$Freq, sep = " => ")))



##########################################################
## Group (nearly) adjacent hits in "TALE array" regions
##########################################################

## Filter out target DNA sequences that have too few repeat CDSs
## NB: for the sake of consistency  it would be better just to filter out from any further consideration the ARRAYS shorter than a certain value (say 5). Could be done on line 351
temp_df <- ddply(nhmmerTabularOutput[,-20], ~ target_name + sq_len, nrow) # I do not know why but it fails to work if I leave the RVD column (#20)
nhmmerTabularOutput <- subset(nhmmerTabularOutput, target_name %in% temp_df[temp_df$V1 > minDomainHitsPerSubjSeq, "target_name"])
nhmmerTabularOutput <- droplevels(nhmmerTabularOutput)
temp_df <- ddply(nhmmerTabularOutput[,-20], ~ target_name + sq_len, nrow)

## Creating a GRanges object from nhmmerOutput
nhmmerOutputGR <- makeGRangesFromDataFrame(df = nhmmerTabularOutput, keep.extra.columns = TRUE, seqnames.field= "target_name")

## Assembling a Seqinfo object representing the set of subject sequences for the GR of hmmer hits
seqinfo(nhmmerOutputGR) <- Seqinfo(seqnames = as.character(temp_df$target_name), seqlengths=temp_df$sq_len, isCircular=NA, genome=NA)

## Group "contiguous" hits (repeats or other regions) in a GRangesList
## Use reduce to obtain the list of regions (arrays of repeats for the time being) that span "contiguous" hits
## Here contiguous is defined as hits that are less than 34bp appart
putativeArraysGR <- reduce(nhmmerOutputGR, drop.empty.ranges=FALSE, min.gapwidth= minGapWidth,
                           with.revmap=TRUE, ignore.strand=FALSE)
revmap <- mcols(putativeArraysGR)$revmap  # an IntegerList


## Use the mapping from reduced to original ranges to group the original ranges by reduced range:
hitsByArraysLst <- relist(nhmmerOutputGR[unlist(revmap)], revmap)
names(hitsByArraysLst) <- paste("ROI", sprintf("%05.0f", 1:length(hitsByArraysLst)), sep = "_")
hitsByArraysLst <- sort(hitsByArraysLst) # just to make sure...

### Or use it to split the DataFrame of original metadata columns by
### reduced range:
##relist(mcols(gr)[unlist(revmap), ], revmap)  # a SplitDataFrameList


##########################################################
## Build reports on TALE repeat sequences including their sequence of RVDs
##########################################################

## Populate metadata about the elements of the list of arrays
mcols(hitsByArraysLst) <- DataFrame(
  arrayID = names(hitsByArraysLst),
  OriginalSubjectName = sapply(hitsByArraysLst, function(x) as.character(seqnames(x))[1]),
  Start = start(putativeArraysGR),
  End = end(putativeArraysGR),
  Strand = strand(putativeArraysGR),
  NumberOfHits = elementNROWS(hitsByArraysLst),
  SeqOfRVD = sapply(hitsByArraysLst, function(x) {
    seq <- as.character(mcols(x)[, "RVD"])
    if (as.character(strand(x))[1] == "-") seq <- rev(seq)
    seq <- paste(seq, collapse = "-")
    return(seq)
  }),
  selectedForAssembly = elementNROWS(hitsByArraysLst) >= minDomainHitsPerArrayForAssembl
)

## Make sure that hits do not overlap for some weird reason
doHitsOverlap <- !isDisjoint(hitsByArraysLst)
if (any(doHitsOverlap)) {
  warning("It appears that some hmmer hits actually overlap.\n It is thus possible that the inferred sequences of RVDs have artefactual insertions.\n")
  warning(paste0("Please check the hits in the following RegionsOfInterest:", "\n",
                 paste(names(doHitsOverlap)[doHitsOverlap], collapse = "\n"), "\n")
  )

}



##########################################################
## Extract N-term and C-term CDS then translate them
##########################################################

## extract N-term CDS from hmmer output hits and translate them
NtermHitsByArray <- lapply(hitsByArraysLst, function(x) {
  x[x$query_name == TALE_NtermDNAHMMName]}) # extract N-term hits

NtermHitsByArray <-as(NtermHitsByArray, "CompressedGRangesList")

fullNtermCDS <- range(NtermHitsByArray) # remove gaps, overlaps to get full length N-term
fullNtermCDS  <- as.data.frame(fullNtermCDS )
colnames(fullNtermCDS) <- c("group", "hitID", "target_name", "start", "end", "width", "strand")
NtermSeqsSet <- extractSeqsfromHits(fullNtermCDS, subjectDNASequences) # extract N-term hit DNA sequences
fullNtermAAseq <- translate(NtermSeqsSet, if.fuzzy.codon = "solve")

fullNtermAAseq <- chartr("*", "X", fullNtermAAseq)
fullNtermAAseqFile <- file.path(outputDir, "fullNtermAAseq.pep")
writeXStringSet(fullNtermAAseq, fullNtermAAseqFile)


## extract C-term CDS from hmmer output hits and translate them
CtermHitsByArray <- lapply(hitsByArraysLst, function(x) {
  x[x$query_name == TALE_CtermDNAHMMName]})
CtermHitsByArray <-as(CtermHitsByArray, "CompressedGRangesList")

fullCtermCDS <- range(CtermHitsByArray) # remove gaps, overlaps to get full length C-term
fullCtermCDS <- as.data.frame(fullCtermCDS)
colnames(fullCtermCDS) <- c("group", "hitID", "target_name", "start", "end", "width", "strand")
CtermSeqsSet <- extractSeqsfromHits(fullCtermCDS, subjectDNASequences)
fullCtermAAseq <- translate(CtermSeqsSet, if.fuzzy.codon = "solve", no.init.codon = TRUE)

fullCtermAAseq <- chartr("*", "X", fullCtermAAseq)
fullCtermAAseqFile <- file.path(outputDir, "fullCtermAAseq.pep")
writeXStringSet(fullCtermAAseq, fullCtermAAseqFile)


## Look at gaps between HitDomains on the same subject sequence to
## detect potential missed repeats
arraysBySeqlevelLst <- split(putativeArraysGR, seqnames(putativeArraysGR)) # group  putativeArraysGR by seqlevel

gaplengthBetweenHitDomains <- sapply(arraysBySeqlevelLst, function(x) {
  dists <- t(as.data.frame(distanceToNearest(x)))[3,]
})
gaplengthBetweenHitDomains <- unlist(gaplengthBetweenHitDomains)
gaplengthBetweenHitDomains <- gaplengthBetweenHitDomains[!is.na(gaplengthBetweenHitDomains)]
gaplengthBetweenHitDomainsbelow500 <- gaplengthBetweenHitDomains[gaplengthBetweenHitDomains <= 500]
quartilesGapLength <- quantile(gaplengthBetweenHitDomainsbelow500,  probs = c(0.25, 0.50, 0.75))
##ggplot(data.frame(gapSize = gaplengthBetweenHitDomainsbelow500), aes(x=gapSize)) + geom_histogram(binwidth=10)

#################
## Write a fasta file of central domain "corrected" CDS for each array.
## This is buildt by concatenating in the correct order the "corrected" repeat sequences fore each array.

cleanDNAarraySeqsList <- lapply(hitsByArraysLst, function(x) {
  if (all(as.character(x@strand) == "+")) {
    hitSeqByArray <- lapply(x$hitID, function(y) toString(cleanDNASeqs[names(cleanDNASeqs) == y]))
    do.call(paste0, hitSeqByArray)
  } else {
    hitSeqByArray <- lapply(x$hitID, function(y) toString(reverse(cleanDNASeqs[names(cleanDNASeqs) == y])))
    reverse(do.call(paste0, hitSeqByArray))
  }
})

cleanDNAarraySeqs <- DNAStringSet(unlist(cleanDNAarraySeqsList))
names(cleanDNAarraySeqs) <- paste(names(subjectDNASequences), names(cleanDNAarraySeqs), sep = "|")

cleanDNAarraySeqsFile <- file.path(outputDir, "correctedRepeatDNASeqs.fas")
writeXStringSet(cleanDNAarraySeqs, cleanDNAarraySeqsFile)

#################
## Write a fasta file for protein repeats sequences

cleanAAarraySeqs <- translate(cleanDNAarraySeqs, if.fuzzy.codon = "solve", no.init.codon = TRUE)
names(cleanAAarraySeqs) <- paste(names(subjectDNASequences), names(cleanAAarraySeqs), sep = "|")

cleanAAarraySeqsFile <- file.path(outputDir, "correctedRepeatAASeqs.fas")
writeXStringSet(cleanAAarraySeqs, cleanAAarraySeqsFile)


#################
## Write a report on the hmmer hits. both  as a  tab-delimited  and a gff
lstOfDf <- mapply(function(x, arrayID, hitsSequences) {
  cbind(
    target_name = as.character(seqnames(x)),
    strand = as.character(strand(x)),
    arrayID, as.data.frame(mcols(x))
  )
},
x = hitsByArraysLst,
arrayID = names(hitsByArraysLst),
SIMPLIFY = FALSE
)
hitsReport <- do.call(rbind, lstOfDf)

repeatSeqs <- repeatSeqsSet[unlist(subset(hitsReport, query_name == repeatDNAHMMName, select = hitID))]

hitsReport <- merge(hitsReport, cbind(hitID = names(repeatSeqs), hitSequence = as.character(repeatSeqs)), all.x = TRUE, sort = FALSE)

hitsReport <- join(hitsReport, indelReport, by = "hitID")

hitsReportFile <- file.path(outputDir, "hitsReport.csv")
write.table(hitsReport,
            file = hitsReportFile,
            sep = "\t",
            row.names = FALSE)


hitsReportToGFF(hitsReportFile) # saving to gff format
rtracklayer::export.gff3(hitsByArraysLst, file.path(outputDir, "hitsReportExperimental.gff"))  # saving to gff format  using a non tested way.



#################
## Write a table report with info on arrays, including the seq of RVD

arrayReportFile <- file.path(outputDir, "arrayReport.txt")

arrayReport <- as.data.frame(mcols(hitsByArraysLst)[order(mcols(hitsByArraysLst)$OriginalSubjectName, mcols(hitsByArraysLst)$NumberOfHits),])

write.table(arrayReport,
            file = arrayReportFile,
            sep = "\t",
            row.names = FALSE)

##ggplot(arrayReport, aes(x=NumberOfHits)) + geom_histogram(binwidth=1)


#################
## Write a fasta file of the selected seq of RVDs without the - separator
seqsOfRVDFile <- file.path(outputDir, "arraysSeqOfRVDsForAssembly.fas")

isSelectedArray <- arrayReport$selectedForAssembly

seqsOfRVDs <- BStringSet(gsub("-", replacement = "", arrayReport$SeqOfRVD[isSelectedArray], ignore.case = FALSE, fixed = TRUE))
names(seqsOfRVDs) <- arrayReport$arrayID[isSelectedArray]
##seqsOfRVDs <- chartr("X", "U", seqsOfRVDs)
## availableAA <- AA_ALPHABET[!AA_ALPHABET %in% uniqueLetters(seqsOfRVDs)]

writeXStringSet(x = seqsOfRVDs, seqsOfRVDFile)

#################
## Generate info messages and log file about the analysis:

stats <- summary(arrayReport$NumberOfHits)

isCompleteArray <- grepl(paste("^", taleArrayStartAnchorCode, ".*", taleArrayEndAnchorCode, "$", sep = ""),
                         arrayReport$SeqOfRVD, perl = TRUE)

completeArrayLengths <- arrayReport$NumberOfHits[isCompleteArray]

## might have been cleaner with a sprintf approach
txt <- c(
  "#****************************************",
  "#**   rvdArrayReader analysis done     **",

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
  paste(Quote(minGapWidth),":", minGapWidth, sep = "\t"),
  paste(Quote(minDomainHitsPerArrayForAssembl),":", minDomainHitsPerArrayForAssembl, sep = "\t"),
  paste(Quote(taleArrayStartAnchorCode),":", taleArrayStartAnchorCode, sep = "\t"),
  paste(Quote(taleArrayEndAnchorCode),":", taleArrayEndAnchorCode, sep = "\t"),

  "#__________Summary measures of TALE search outcome__________",
  paste("Number of analysed subject sequences :", length(subjectDNASequences), sep = "\t"),
  paste("Total number of TALE repeat DNA coding sequence motif hits found with the nhmmer approach:",
        numberOfRepeatHitsBeforeFiltering, sep = "\t"),
  paste("Number of column masked in the repeat unit hits alignment:", maskedncol(cleanMultipleAlignment), sep = "\t"),
  paste("Total number of repeat HMM hits on the corresponding set of translated DNA hits:", sum(RVDtbl), sep = "\t"),

  paste("Total number of subject seqs with TALE motif hits after low hit number filtering:", length(seqlevels(putativeArraysGR)), sep = "\t"),
  paste("Total number of distinct regions (repeat arrays) with adjacent TALE motifs :", nrow(arrayReport), sep = "\t"),
  paste("Total number of arrays selected for SPA assembly:", sum(isSelectedArray),sep = "\t"),
  paste("Total number of 'complete' arrays (with both N- and C-term flanking motifs):", sum(isCompleteArray),	sep = "\t"),

  paste("Total number of distinct types of RVD:", nrow(RVDtbl), sep = "\t"),

  paste("Minimum array length (number of TALE motif hits):", min(arrayReport$NumberOfHits), sep = "\t"),
  paste("Maximum array length:", max(arrayReport$NumberOfHits), sep = "\t"),
  paste("Median array length:", median(arrayReport$NumberOfHits), sep = "\t"),
  paste("Length of the longest 'complete' array:", max(completeArrayLengths),	sep = "\t"),
  paste("Length of the shortest 'complete' array:", min(completeArrayLengths),	sep = "\t"),

  #message("Distribution of the number of TALE motif hits per array:\n")
  #message(paste(names(stats), stats, sep = "\t", collapse = "\n"))

  paste("Number of gaps of size below 500nt between TALE motifs arrays:", length(gaplengthBetweenHitDomainsbelow500)/2, sep = "\t"),

  paste("First quartile of size of gaps (below 500nt) between TALE motifs arrays:", quartilesGapLength[1], sep = "\t"),
  paste("Median size of gaps (below 500nt) between TALE motifs arrays:", quartilesGapLength[2], sep = "\t"),
  paste("Upper quartile of size of gaps (below 500nt) between TALE motifs arrays:", quartilesGapLength[3], sep = "\t"),

  "#*************************\n"
)




message(paste(txt, collapse ="\n"))

analysisLogFile <- file.path(outputDir, "rvdArrayReader.log")
logf <- file(analysisLogFile, open = "w")
writeLines(text = txt, con = logf)
close(logf)
}

####FUNCTION CALLING####
rvdArrayReader(subjectFile= myArgs$subjectFile,
               outputDir = myArgs$outputDir,
               hmmFilesDir = myArgs$hmmFileDir,
               minRatioOfGapForColMasking = myArgs$minRatioOfGapForColMasking, # columns of the tale repeat CDS alignment that contain a gap in a fraction of sequence higher than this value (betwen 0 and 1) will be masked from the alignment when translating the DNA sequences to protein.
               TALE_NtermDNAHitMinScore = myArgs$TALE_NtermDNAHitMinScore, # nhmmer score cut_off value
               repeatDNAHitMinScore = myArgs$repeatDNAHitMinScore, # nhmmer score cut_off value
               TALE_CtermDNAHitMinScore = myArgs$TALE_CtermDNAHitMinScore, # nhmmer score cut_off value
               minDomainHitsPerSubjSeq = myArgs$minDomainHitsPerSubjSeq, # Minimum number of nhmmer hits for a subject sequence to be reported as having TALE diagnostic regions. This is a way to simplify output a little by getting ride of uninformative sequences
               minGapWidth = myArgs$minGapWidth, # minimum gap between two tale domain hits for them to be considered "contiguous" and grouped in the same array.
               minDomainHitsPerArrayForAssembl = myArgs$minDomainHitsPerArrayForAssembl, # Minimum number of repeat in an array for its seq of RVD to be considered for assembly. This is a way to get ride of sequences that are too short reasonably be of any help for assembly
               taleArrayStartAnchorCode = myArgs$taleArrayStartAnchorCode,
               taleArrayEndAnchorCode = myArgs$taleArrayEndAnchorCode,
               hmmerpath = myArgs$hmmerpath)


# specify signal sendt to calling shell environment in case of successful execution
quit(save = "no", status = 0, runLast = FALSE)
