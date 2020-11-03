

#######################################################
###             Parallelized functions             ####
###   for running TALE target prediction scripts   ####
#######################################################

# Reformat geonomes and promoteromes fasta files so that
# Talvez and Talgetter can use them for EBE search
reformatFastaForEBEPredictions <- function(file) {
# Reformat geonomes and promoteromes fasta files so that
# Talvez and Talgetter can use them for EBE search
# Return value is the full name of the reformated file
# Author: sebastien
	###############################################################################
# Goals:
# - remove trailling EOL in sequence lines that Talvez do not like
# - remove nucleotides that are not GAT or C because Talgeter complains
	###############################################################################

	library(Biostrings)

# Initiate file names and locations
	initialFileName <- file
# initialFileName <- file.path(getwd(), "candidateRLKgenomicSequences.fa")

	tempFileName <- tempfile(pattern = "ReformatedFasta", tmpdir = getwd(), fileext = ".fa")
	message(paste("Writting to temp file:", tempFileName))
	reformatedFileName <- file.path(getwd(), paste("reformated", basename(initialFileName), sep = "_"))

# Loading the temp sequence file in R
# This will automatically REMOVE invalid one-letter sequence codes
	temp <- readDNAStringSet(initialFileName)
	table(width(temp))
	fasta.seqlengths(initialFileName, seqtype="DNA")
# Writting the modified sequences to disk in a fasta file
	writeXStringSet(temp, tempFileName, append=FALSE, format="fasta")

# Removing line feeds at the end of sequence lines belonging to the same sequence (>)
# and writting a properly formated fasta file
	cmnd <- paste("perl /media/cunnac/DONNEES/CUNNAC/Lab-Related/MyScripts/Utils/format_fasta2.pl", tempFileName, ">", reformatedFileName)
	system(command = cmnd, intern = FALSE)

# Delete the temp file
	unlink(tempFileName)
	message("\nReformating complete!\n")
	return(reformatedFileName)
}



pTalgetter <- function(TALs, cl, predictorPath, promotersFile) {

	if(!file.exists(predictorPath)) stop("Unable to find the Talgetter program at the specified location. Please verify the file exists")
	if(!file.exists(promotersFile)) stop("Unable to find the set of target DNA sequences (fasta file) at the specified location. Please verify the file exists")

	talgetterCommand <- paste("java -Xms512M -Xmx2G -jar ",
			predictorPath, " input=", promotersFile,
			" rvd=", TALs$RVDseq,
			" top=200",
			" model=TALgetter13",
			" 2>> pTalgetter_StdErr.txt",
			sep = "") # assemble a vector of talgetter commands to run predictions for individual TALs
	message(paste(talgetterCommand, collapse = "\n"))

	allTALPreds <- clusterMap(cl = cl, fun = function(comm, talKEY) {
				TALPred <- system(command = comm, ignore.stderr = FALSE, intern = TRUE) # run the predictions in parallel and collect stdoutput

				TALPred <- do.call(rbind, strsplit(TALPred, "\t", fixed = TRUE)) # split lines of output on tab and merge them row-wise in a matrix
				# tests wether the considered TAL has given no predicitons and output a row of 7 "" cells...
				if (!(length(TALPred) > 1 && TALPred[1,1] %in% "# ID")) {
					TALPred <- c("# ID", "Position", "Distance to end", "Sequence", "Matches", "Score", "p-value", "E-value", "talKEY")
				} else {
					TALPred <- cbind(TALPred, talKEY) # append the talKey to its cognate predicitons
				}
			},
			comm = talgetterCommand, talKEY = TALs$talKEY,
			RECYCLE = FALSE
	)
	allTALPreds <- do.call(rbind, allTALPreds) # agregate the predictions for individual TAL in a single data frame
	colnames(allTALPreds) <-  c("# ID", "Position", "Distance to end", "Sequence", "Matches", "Score", "p-value", "E-value", "talKEY") # add column names

	allTALPreds <- allTALPreds[!allTALPreds[,1] %in% "# ID", ]# get ride of the column names
	return(allTALPreds)
}


#java -Xms512M -Xmx2G -jar /mnt/DONNEES/CUNNAC/Lab-Related/MyScripts/TALETargets_Interact/Version_3-0/PredictorScripts/TALgetter/TALgetter.jar \
#input=/mnt/DONNEES/CUNNAC/Lab-Related/Exp_Projects/Xo-Rice_genomics/OsIndicaASM465v1_GenePromoterome500bp_reformatedTalvez.fa \
#rvd="NS-NG-NS-HD-NI-NG-NN-NG-HD-NI-NN-N*-NI-NN-HD-NG-NI-NN-N*-HD-NN-NG" \
#top=200 > TALgetter.out



pTalvez <- function(TALs, cl, predictorPath, promotersFile) {

	if(!file.exists(predictorPath)) stop("Unable to find the Talgetter program at the specified location. Please verify the file exists")
	if(!file.exists(promotersFile)) stop("Unable to find the set of target DNA sequences (fasta file) at the specified location. Please verify the file exists")

# Split the data frame of all TAL queries into a list of TAL queries data.frames based on the number of
# available cores for parallel processing.
	TALsubsetsList <- lapply(clusterSplit(cl, seq_along(TALs$talKEY)), function(x) {TALs[x,]})

# Run Talvez predictions for individual TAL queries sub data.frames in parallel
	allTALPreds <- parLapply(cl = cl, TALsubsetsList, function(TALsDF, predictorPath, promotersFile) {

				# format TALs query appropriately and save in a temp file
				TALRVDFile <- tempfile(pattern = "TALSubset", tmpdir = getwd(), fileext = ".rvd")
				writeLines(text = paste(">", TALsDF$talKEY, "\t", TALsDF$RVDseq, sep = ""),
						con = TALRVDFile)

				# talvez temp output file and talvez script file location
				talvezOutFile <- paste(TALRVDFile, "_complete.txt", sep = "")
				talvezDir <- dirname(predictorPath)

				# assemble the command for Talvez query
				talvezCommand <- paste("perl ",
						predictorPath,
						" -t 0 -l 19",
						" -e ", file.path(talvezDir, "mat1"),
						" -z ", file.path(talvezDir, "mat2"),
						" ",
						TALRVDFile, " ",
						promotersFile, " ",
						TALRVDFile,
						" 2>> pTalvez_StdErr.txt",
						sep = "")
				message(talvezCommand)

				# run predictions
				system(command = talvezCommand, ignore.stderr = FALSE, intern = TRUE)

				# read temp prediciton file and import as a df
				TALPred <- read.delim(file = talvezOutFile)

				## Remove temporary files and dir:
				# tmp of Talvez
				unlink(file.path(getwd(),"tmp"), recursive = TRUE)
				# TALRVDFile and TALRVDFile output: file.path(getwd(), paste(TALRVDFile, "_complete.txt", sep = ""))
				unlink(c(talvezOutFile, TALRVDFile), recursive = FALSE)

				# return the df of predictions for this TAL subset
				return(TALPred)
			}, predictorPath = predictorPath, promotersFile = promotersFile
	)
	allTALPreds <- do.call(rbind, allTALPreds) # agregate the predictions for individual TAL subsets in a single data frame
	return(allTALPreds)
}


#perl /mnt/DONNEES/CUNNAC/Lab-Related/Exp_Projects/Xo-Rice_genomics/EBEPredictorPrograms/TALVEZ/TALVEZ_3.3.pl \
#-t 100 -l 19  \
#-e /mnt/DONNEES/CUNNAC/Lab-Related/Exp_Projects/Xo-Rice_genomics/EBEPredictorPrograms/TALVEZ/mat1 \
#-z /mnt/DONNEES/CUNNAC/Lab-Related/Exp_Projects/Xo-Rice_genomics/EBEPredictorPrograms/TALVEZ/mat2 \
#/mnt/DONNEES/CUNNAC/Lab-Related/Exp_Projects/Xo-Rice_genomics/EBEPredictorPrograms/TALVEZ/TALCONTROL_RVDs.txt \
#/mnt/DONNEES/CUNNAC/Lab-Related/Exp_Projects/Xo-Rice_genomics/temp.fa \





predictEBEs <- function(TALs, queryDNASeqFile, talgetterPath, talvezPath, cl, ...) {
	##################################################################
	##                Run TALgetter predictions                     ##
	##################################################################
	talgetterOutFile <- paste(basename(queryDNASeqFile), "_", basename(talgetterPath), ".txt", sep = "")
	tempTg <- pTalgetter(TALs = TALs, cl = cl, predictorPath = talgetterPath, promotersFile = queryDNASeqFile)
	write.csv(tempTg, file = talgetterOutFile)

	##################################################################
	##                   Run Talvez predictions                     ##
	##################################################################
	# NB: would need to modify talvez code to include a reformating of the
	# fasta file if it does not comply with talvez the requirements.
	# it would be best to include the code of the
	# format_fasta2.pl script.
	talvezOutFile <- paste(basename(queryDNASeqFile), "_", basename(talvezPath), ".txt", sep = "")
	tempTv <- pTalvez(TALs = TALs, cl = cl, predictorPath = talvezPath, promotersFile = queryDNASeqFile)
	write.csv(tempTv, file = talvezOutFile)

	##################################################################
	##  Parsing predictions  in appropriate format for DB storage   ##
	##################################################################
	# could have used a lapply(), I guess...
	parsedTv <- ParseTALvezHitLikeFile(my.file = talvezOutFile,
			TALs.info = TALs,
			cl = cl, ...)

	parsedTg <- ParseTALvezHitLikeFile(my.file = talgetterOutFile,
			TALs.info = TALs,
			cl = cl, ...)
	predictions <- rbind(parsedTv, parsedTg)

	##################################################################
	##                         Clean up                             ##
	##################################################################
	unlink(c(talvezOutFile, talgetterOutFile), recursive = FALSE)
	return(predictions)
}







###################################
###         ParseEBEHits        ###
###################################



ParseTALvezHitLikeFile <- function(my.file, TALs.info = NULL, predictType,
		DNAType, infoOnSeqs = chrominfo, cl = cl,
		firstEBE_ID = 0, rerank = Inf, rankCutoff = NULL, ...) {

	# Guess the nature of the prediction algorithm from the file name
	if (missing(predictType) || is.null(predictType)) {
		if (grepl(pattern = "TAlvez", basename(my.file), ignore.case = TRUE)) {predictType <- "Talvez"} else {
			if (grepl(pattern = "Talgetter", basename(my.file), ignore.case = TRUE)) {predictType <- "Talgetter"} else {
				stop("The predictor program is not specified or is unrecognized. Cannot proceed with parsing the predictions file.")
			}
		}
	}

	# Make sure the supplied value of predictType is appropiate
	if (!(length(predictType) == 1L || predictType %in% c("Talvez", "Talgetter"))) {
		stop("The supplied value of the predictor program is not recognized.  Cannot proceed with parsing the predictions file.")
	}
	print(predictType)

	# Make sure that the infoOnSeqs object is of class "SeqInfo"
	if (missing(infoOnSeqs)) {
		if (DNAType == "genome") {stop("Unable to process genomic predictions without info on molecule names and length provided in the infoOnSeqs argument.")}
	} else {
		if (DNAType == "genome" & !is(infoOnSeqs, "Seqinfo")) stop("The provided value of the infoOnSeqs parameter is not of class 'SeqInfo'")
	}


	#####################################################################################
	# Load the content of the TALE target predictor output file and returns a data.frame.
	#####################################################################################
	EBEs.raw <- read.csv(my.file, header = TRUE,
			check.names = FALSE, stringsAsFactors = TRUE, strip.white = TRUE)

	# Standard names for a EBE prediction table
	standardEBETableColNames <- c(
			"TALkey",
			"SeqID",
			"SCORE",
			"SEQUENCE",
			"PREDICTYPE",
			"StartPos",
			"EndPos",
			"Rank",
			"DistFromEnd",
			"p-value")

	# Deal with empty prediction files
	if(nrow(EBEs.raw) == 0) {
		EBEs.raw <- matrix(ncol = length(standardEBETableColNames), nrow = 0)
		colnames(EBEs.raw) <- standardEBETableColNames
		message("Returning an empty prediction table !!!")
		return(EBEs.raw)
	}

	# Append a "PREDICTYPE" column that record the nature of the prediction algorithm
	EBEs.raw <- cbind(EBEs.raw, PREDICTYPE = predictType)

	# Parse Talvez prediction files
	if (identical(predictType, "Talvez")) {
		EBEs.raw <- EBEs.raw[,c("TAL_ID", "SEQ_ID", "SCORE", "SEQ", "PREDICTYPE",
						"INI", "END", "RANK", "DISTANCE_FROM_END")] # keep  usefull columns only

		# Add a p-vaue column
		EBEs.raw <- cbind(EBEs.raw, "p-value" = NA)
		EBEs.raw[,"p-value"] <- as.numeric(EBEs.raw[,"p-value"])

		# Rename the columns to match the format of the table in the DB
		colnames(EBEs.raw) <- standardEBETableColNames
	}

	# Parse Talgetter prediction files
	if (identical(predictType, "Talgetter")) {
		# keep  usefull columns only
		EBEs.raw <- EBEs.raw[,c("# ID", "Position",	"Distance to end", "Sequence", "Score", "p-value",	"talKEY", "PREDICTYPE")]
		# add "fake" columns for those that are missing relative to the DB table.
		EBEs.raw <- cbind(EBEs.raw, EndPos = integer(length = 1))
		EBEs.raw <- cbind(EBEs.raw, Rank = integer(length = 1))

		# and re-order to match the format of the table in the DB
		EBEs.raw <- EBEs.raw[, c("talKEY",
						"# ID",
						"Score",
						"Sequence",
						"PREDICTYPE",
						"Position",
						"EndPos",
						"Rank",
						"Distance to end",
						"p-value")]

		# Rename the columns to match the format of the table in the DB
		colnames(EBEs.raw) <- standardEBETableColNames

		# Add 1 to the start position because Talgetter uses zero for origin of numbering
		EBEs.raw$StartPos <- EBEs.raw$StartPos+1

		# Calculate the EndPos values and update the data frame accordingly
		EBEs.raw$EndPos <- EBEs.raw$StartPos + nchar(as.vector(EBEs.raw$SEQUENCE)) - 1

	}

	EBEs.raw$TALkey <- sub(">", "", EBEs.raw$TALkey) # remove leading > if present

	# Whenever possible, make sure that the "TAL ids" levels in the predictor output file are found
	# in the talKey field of the TALs.info table
	if (is.null(TALs.info)) {
		message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n
						BEWARE, no TALs.info argument value provided.\n
						'TAL ids' levels in the TALvez output file will not be checked to have a match\n
						in the talKEYS levels of the TALs.info table!\n
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
	} else {
		if (!all(levels(EBEs.raw$TALkey) %in% TALs.info$talKEY)) {
			warning("BEWARE, some 'TAL ids' levels in the TALvez output file were not found in the Desc field of the TALs.info table!\n")}
	}

	# Add a numeric EBE_ID column
	# Not very smart because this is not handeled consistently for genome and promoter predictions
	# For the moment the firstEBE_ID argument allows to create values that are not in the database
	# by providing an origin number for starting the numbering.
	EBEs.raw <- cbind(EBE_ID = 1:nrow(EBEs.raw), EBEs.raw)

	# Convert distance to end to a positive value
	EBEs.raw$DistFromEnd <- abs(EBEs.raw$DistFromEnd)

	# if rerank is not NULL, change the value of the Rank field to the rank of their score
	# in the set of EBEs for a given TAL that are located at most at the rerank value distance
	# from the promoter 3' end. For those EBEs that are located further upstream, the rank is
	# set to 60000
	# After this reranking operation is done, the rankCutoff parameter allows to keep only the EBEs
	# that rank the same or better than the provided rankCutoff cutoff.
	if (is.numeric(rerank)) {
		indexes <- tapply(rownames(EBEs.raw), EBEs.raw$TALkey, function(x) x) # for each TALkey, fetch the corresponding predictions (raws) and create a list of rownames vectors
		EBEs.raw <- do.call(rbind,
				parLapply(cl, indexes, function(x, EBEs.raw, rankCutoff, rerank) {
							EBEs <- EBEs.raw[x,]
							EBErerankIN <- EBEs[abs(EBEs$DistFromEnd) <= rerank, ] # select EBEs that are at most rerank base pairs from the end (5' end) of the sequence
							EBErerankOUT <- EBEs[!abs(EBEs$DistFromEnd) <= rerank, ] # select the complement of the EBEs set
							if(nrow(EBErerankOUT) != 0) {EBErerankOUT$Rank <- Inf} # Assign an infinite rank to the EBEs that are too far away from the end of the sequence

							# Perform the actual reranking of the EBE hits based on decreasing score.
							ranking <- numeric()
							ranking[order(EBErerankIN$SCORE, decreasing = TRUE)] <- 1:nrow(EBErerankIN)
							EBErerankIN$Rank <- ranking

							if (nrow(EBErerankOUT) == 0) {EBEs <- EBErerankIN} else {EBEs <- rbind(EBErerankIN, EBErerankOUT)} # merge the two EBEs sets (IN and OUT)
							if (is.numeric(rankCutoff)) {
								EBEs <- EBErerankIN[EBErerankIN$Rank <= rankCutoff, ]  # keep only EBEs with rank equal or below the cut off
							}
							return(EBEs)
						}, EBEs.raw = EBEs.raw, rankCutoff =rankCutoff, rerank = rerank)
		)

	} else {
		if (!is.null(rerank)) {
			message("The provided value for the rerank parameter is not numeric nor NULL and was ignored in the call.\n")}
	}
	EBEs.raw$EBE_ID <- 1:nrow(EBEs.raw)+firstEBE_ID



	#############################PROMOTER ONLY PREDICTIONS###########################################
	if (DNAType == "promoters") {
		# Standard names for a EBE prediction table from promoter sequences
		promoterEBETableColNames <- c(
				"EBE_ID",
				"TALkey",
				"ensembl_gene_id",
				"SCORE",
				"SEQUENCE",
				"PREDICTYPE",
				"StartPos",
				"EndPos",
				"Rank",
				"DistFromEnd",
				"p-value")

		# Rename columns of the EBE.raw df to comply with the table structure in the DB
		colnames(EBEs.raw) <- promoterEBETableColNames

		# Exit function and return the parsed data frame.
		return(EBEs.raw)
	}
	##############################GENOMIC PREDICTIONS####################################################
	if (DNAType == "genome") {

		# Some formalism regarding the naming of the molecules that are fed to the EBE predictors:
		#-----------------------------------------------------------------------------------------
		# For a "whole genome" search, the name of the molecules (generally the plus strand) is appendend with a "_ps" suffix.
		# These molecules are reverse complemented to get the sequence of the complementary strand and named with the name of
		# the initial molecule appended with a "_ms" suffix.

		# First verify that the core names of the seqID in the EBEs.raw data frame are all containned in the infoOnSeqs object
		coreNameOfSequences <- as.factor(gsub("_[pm]s$", "", EBEs.raw$SeqID, perl = TRUE))
		if (!all(levels(coreNameOfSequences) %in% seqnames(infoOnSeqs))) {
			stop("Some of the names of the EBE-associated sequences are not contained in the names attribute of the object passed as infoOnseqs value.")
		}

		vectorOfStrand <- rep("+", nrow(EBEs.raw))

		# Assign a "-" strand values to those rows that have a "_ms" suffix in the SeqID value.
		vectorOfStrand[grepl(pattern = "_ms$", as.vector(EBEs.raw$SeqID), perl = TRUE)] <- "-"

		# Append a "Strand" column that record the polarity of the sequence on which the EBEs were found
		EBEs.raw <- cbind(EBEs.raw, Strand = vectorOfStrand)

		# Get ride of the suffix in the SeqID and keep only the core name of the molecules
		EBEs.raw$SeqID <- coreNameOfSequences

		# Transform coordinates depending on strand
		# Convert relative positions of the minus strand to ascending absolute genome positions
		# using the chromosome sizes defined in the chr.sizes vector
		EBEsOnMinus <- EBEs.raw$Strand == "-" # logical vector of EBEs on the minus strand

		# start = 1 + length of the corresponding chromozome - End postion in the coordinate system
		# where the origin is the first nucleotide of the reverse complemented chromosome sequence
		left <- (1 + seqlengths(infoOnSeqs)[match(EBEs.raw$SeqID[EBEsOnMinus], seqnames(infoOnSeqs))]
					- EBEs.raw$EndPos[EBEsOnMinus])

		# end = 1 + length of the corresponding chromozome - start postion in the coordinate system
		# where the origin is the first nucleotide of the reverse complemented chromosome sequence
		right <- (1 + seqlengths(infoOnSeqs)[match(EBEs.raw$SeqID[EBEsOnMinus], seqnames(infoOnSeqs))]
					- EBEs.raw$StartPos[EBEsOnMinus])

		# Update the data frame with the transformed coordinates
		EBEs.raw$StartPos[EBEsOnMinus] <- left
		EBEs.raw$EndPos[EBEsOnMinus]	<- right

		# Change strand encoding to 1 for + and -1 for -
		levels(EBEs.raw$Strand) <- c(-1, 1)

		# Exit function and return the parsed data frame.
		return(EBEs.raw)

	} else {stop("The value of the DNAType argument is not approriate.")}

}




#################################
###      Parse Query TALS     ###
#################################

TALvezQueryTALsParser <- function(my.file) {
	## Extract the info relative to TALs that were used to interrogate TALvez
	## and return a data.frame

# Import TALvez file with info relative to TALS that were used as a data frame
	TALs.raw <- read.delim(my.file, header = FALSE, sep = "\t",
			check.names = TRUE, stringsAsFactors = TRUE, strip.white = TRUE)

# Assemble all info to have it ready for analysis or export to DB
	colnames(TALs.raw) <- as.character(unlist(TALs.raw[1,]))
	TALs.raw <- TALs.raw[-1,]

#	TALs.info <- data.frame(
#			"talKEY" = 1:nrow(TALs.raw),
#			"gi" = TALs.gi,
#			"gb" = TALS.gb,
#			"Strain" = strains,
#			"Descr" = TALs.raw[,2],
#			"RVDseq" = TALs.raw[,3],
#			stringsAsFactors = FALSE)

	return(TALs.raw)
}




#####################################################################
##     Investigating the composition and RVD-content of TALs       ##
#####################################################################

##################################
## Definition and some methods for a TALRepeatDomain class
setClass("TALRepeatDomain",
		representation(
				names="character",
				strain="character",
				RVDseq="character"),
		validity = function(object) {
			length(object@names) == 1 &&
					length(object@strain) == 1 &&
					all(nchar(object@RVDseq) == 2) &&
					all(object@RVDseq %in% c("H*", "HA", "HD", "HG", "HH", "HI", "HN", "N*",
									"NA", "NC", "ND", "NG", "NH", "NI", "NK", "NN", "NQ", "NS", "NV", "S*", "SN", "SS", "YG",
									"[[", "]]", "NX"
									)
					)
		}
)


TALRepeatDomain <- function(r) {
  ## A constructor for a TALRepeatDomain object
  ## from a df row comming from the data base
  new("TALRepeatDomain", names = r$gb, strain = r$Strain,
      RVDseq = as.vector(unlist(strsplit(r$RVDseq, "-"))))
}


setMethod(names, signature=c("TALRepeatDomain"), definition = function (x) x@names)


getRVD <- function(object, endTags = TRUE) {
# An accessor for the vector of RVD in a TALRepeatDomain object
  x <- object@RVDseq
  if(!endTags) x <- x[!x %in% c("[[", c("]]"))]
  return(x)
}

getStrain <-  function(object) {object@strain}


RVDSeqToXStringSet <- function(object, RVDs = c("H*", "HA", "HD", "HG", "HH", "HI", "HN", "*N", "N*",
				"NA", "NC", "ND", "NG", "NH", "NI", "NK", "NN", "NQ", "NS", "NV", "S*", "SN", "SS", "YG"),
		singleLetterCode = LETTERS[1:24]) {
	# Convert the vector of RVDs in the RVDseq slot of a TALRepeatDomain object into
	# a length one BStringSet object.
	# The diresidues are recoded to single letters using the vectors of RVDs and correpsonding substitute letters,
	# respectively in arguments RVDs and singleLetterCode.
	library(Biostrings)
	if (class(object) == "TALRepeatDomain") {stop("The provided object must be of class 'TALRepeatDomain'.")}
	RVDseq <- getRVD(object)
	RVDseq <- as.factor(RVDseq)
	levels(RVDseq) <- singleLetterCode[RVDs %in% levels(RVDseq)] # substitute the RVD diresidue with a single letter as specified in the arguments
	RVDseqString <- paste(RVDseq, collapse = "")
	RVDseqBStringSet <- BStringSet(x = RVDseqString)
	names(RVDseqBStringSet) <- "bla"
	return(RVDseqBStringSet)
}


# RVDDistMAtrix <-

#################################
## Formerly contained in the
## QueryingBiomart file
#################################


#################################################
# Interaction with Biomart and database feeding #
#################################################

createPlantMart <- function (dataset = "osativa_eg_gene"){
# use biomaRt package to create a Mart connection with the specified plant genome dataset.
# The default plant genome in plant_mart database to be queried is O. sativa (MSU6)
  library("biomaRt")
  return(useMart("plants_mart_12", dataset = dataset))
}



getGRangesContext <- function (mart, gr,
    attributes = c("ensembl_gene_id",
        "chromosome_name", "start_position", "end_position", "strand",
        "description")) {
  ## This function is very efficient for fetching info in biomart for large
  ## numbers of genomic positions (tested with up to 2000) stored in a GenomicRanges object.
  ##  takes a mart object connecting to a plant genome dataset in Bsent to the Biomart dataset
  library("biomaRt")
  library("GenomicRanges")

  # Check if mart is of class mart
  if (!is(mart, "Mart")) {stop("The provided mart connection object is not of class 'Mart'.")}

  # Check if gr is of GenomicRanges class
  if (!is(gr, "GenomicRanges")) {stop("The provided gr object is not of class 'GenomicRanges'.")}

  # vector of selected filters
  filters <- c("chromosomal_region")

  # prepare the loci values as a list
  loci <- paste(as.integer(as.vector(seqnames(gr))),
      start(gr), end(gr), ifelse("+" == (as.vector(strand(gr))), 1, -1), sep = ":")

  # Perform actual query, retrieve and return dataframe
  results <- (getBM(attributes, filters, loci, mart))
  cat("** Query for EBE_ID:\n", paste(names(gr), collapse = "\n"),
      "\n** was completed on ", date(),"\n")

  if(nrow(results) == 0) {return(NULL)} else {return(results)} # Gives the possibility to filter empty results
}



GRanges2DB <- function(mart = genome, gr, query.length = 500,
		log.file = TRUE, DB.con = TALR, table = "TargetGenes",  ...) {
	## Uses in general a large GenomicRanges object
	## Iteratively fetches the info corresponding to positions in Biomart
	## Add the records retrieved from the mart object to the table "table".
	## Relies on the GetGRangesContext function to interogate the mart. Therefore, additional arguments
	## such as biomart attributes vector can be passed in the call to this function.
	## Note that the EBE_IDs are not mapped to the retrieved genes.

	# directs output to "biomartlog.txt" if desired
	if (log.file) sink(file = file.path(getwd(), "biomartlog.txt"), append = TRUE, split = TRUE, type = "output")

	# Check that DB.conn is open
	tryCatch(dbGetInfo(DB.con))

	# Message reporting the time stamp of the function call
	cat("** Sending ",length(gr), "Genomic Positions queries to Biomart on ", date(), "\n")

	# Setting indexes values to iterate over the gr object elements query.length at a time
	# Note that I am re-inventing the wheel : could have been done
	# with a for statement and a j range given by seq(offset, length(gr), by = offset)...
	i <- 1
	j <- i + query.length
	if (j > length(gr)) { j <- length(gr)}

	# Run the query in biomart, retrieve results and append to table in DB
	while (j <= length(gr)) {
		gene.list <- getGRangesContext(mart = mart, gr = gr[i:j])
		if (!is.null(gene.list)) {
			# Convert the Chr names returned by biomart to a managable format:
			if(all(is.numeric(gene.list[,2]) &  gene.list[,2] <= 99 &  gene.list[,2] >=0)) {
				gene.list[,2] <- sprintf("%02d", gene.list[,2])}
			# Populate DB with resulting rows
			dbWriteTable(conn = DB.con, name = table, value = gene.list,
					row.names = FALSE, overwrite = FALSE, append = TRUE)
			##   sqlSave(DB.con, gene.list, tablename = table, append = TRUE,
			## rownames = FALSE, fast = FALSE)
		}

		# Increment appropriately the values of the indexes
		i <- j+1
		if (i > length(gr)) {break}
		if ((i + query.length) > length(gr)) {j <- length(gr)} else {j <- i + query.length}
	}

	# Message reporting the time stamp of the end of the query
	cat("** All Queries completed on ", date(), "\n")

	# close the sink to redirect output to stdout()
	if (log.file) while(sink.number()) {sink()}
}



getAffy <- function(mart, gene.ids,
    attributes = c("ensembl_gene_id", "affy_rice"), ...) {
# Fetch by default rice affymetrix probe set ids using
# a list of ensembl gene ids as queries
# returns only results with a non empty value
  results <- getBM(attributes = attributes,
      filters = "ensembl_gene_id",
      values = unlist(gene.ids, use.names = FALSE),
      mart, uniqueRows = TRUE)
  return(results[apply(results["affy_rice"], 1, nchar) > 0,])

}

getAnnot <- function(mart, gene.ids,
    attributes = c("ensembl_gene_id", "superfamily"), ...) {
# Fetch by default superfamily ids using
# a list of ensembl gene ids as queries
# returns only results with a non empty value
# other portential attributes : "interpro", "interpro_short_description"

  results <- getBM(attributes = attributes,
      filters = "ensembl_gene_id",
      values = unlist(gene.ids, use.names = FALSE),
      mart, uniqueRows = TRUE)
  return(results[!is.na(results[2]),])
}



###################################################################
# Mapping of EBE with target genomics features (all as GRanges)    #
#                   and save in DB                                #
###################################################################

findGRMatchesTypes <- function(query, subject, output="gr", ...) {
  require(GenomicRanges); require(IRanges)
  ## Utility: identify overlaps in range data sets, such as annotation or alignment positions defined
  ## by two GRanges objects. Extends the original types of overlap to include "upstream" and "downstream"
  ## when the findOverlaps match finding engine is used with a maxgap argument different from 0.
  ## Provides various distances between both extremities of query and hit.
  ## Heavily inspired from Thomas Girke's olRanges
  ## Details on usage and use cases of the original function are available here:
  ## http://manuals.bioinformatics.ucr.edu/home/ht-seq#TOC-Analysis-Routines-with-IRanges-Geno


  ## Input check
  if(!(class(query)=="GRanges" & class(subject)=="GRanges")) {
    stop("Query and subject need to be of class GRanges.")
  }
  if(!all(seqlengths(query) == seqlengths(subject))) {
    stop("Query and subject must have identical seqlengths.")
  }
  if(any(is.na(seqlengths(query))) | any(is.na(seqlengths(subject)))) {
    stop("Query and subject seqlengths may not contain NAs values.")
  }

  ## A convenience function: depending on argument extrem returns a Rle vector with the "start" or
  ## "end" of the elements in the GenomicRanges object gr depending on strand.
  extremityPos <- function(gr, extrem){
    ## A convenience function: depending on argument extrem returns a Rle vector with the "start" or
    ## "end" of the elements in the GenomicRanges object gr, depending on their strand.
    ## If "+" or "*", the corresponding value is directly extracted from the range.
    ## If "-", the coordinates of the gr element are transposed in a referential whose origin
    ## is the end of the chromosome or its seqname (origin for numbering is the 5'end of the "-" strand,
    ## equivalent to some kind of reverse operation). See below:
    ##     input gr :     output gr :
    ##     1S    E       1  S    E
    ##     ||    |       |  |    |
    ##     -<*****---    ---*****>-
    ## So in addition to being transposed the positions of the extremities are swapped.

    extremities <- c("start", "end") # allowed values of extrem
    if(!extrem %in% extremities) {stop("the extrem argument must be either 'start' or 'end'.")}
    extremity <- extrem == extremities # the value of extrem converted in a logical vector for indexing.

    ifelse(strand(gr) == "-",
        seqlengths(gr)[as.character(seqnames(gr))] - do.call(extremities[!extremity], list(gr)) + 1,
        do.call(extremities[extremity], list(gr)))
  }

  ## Find overlapping ranges
  olindex <- as.matrix(findOverlaps(query, subject, ...))
  query <- query[olindex[,1]] # a gr object containg the subset of query gr elements
  subject <- subject[olindex[,2]] # a gr object containg the subset of subject gr elements

  ## Data.frame of the strand relative coordinates of query and their hits
  matches.pos <- cbind(
      Qstart = as.integer(extremityPos(query, "start")),
      Qend = as.integer(extremityPos(query, "end")),
      Sstart = as.integer(extremityPos(subject, "start")),
      Send = as.integer(extremityPos(subject, "end"))
  )

  ## Pre-queries for match types
  startup <- matches.pos[,"Qstart"] < matches.pos[,"Sstart"]
  startin <- matches.pos[,"Qstart"] >= matches.pos[,"Sstart"] & matches.pos[,"Qstart"] <= matches.pos[,"Send"]
  startdown <- matches.pos[,"Qstart"] > matches.pos[,"Send"]
  endup <- matches.pos[,"Qend"] < matches.pos[,"Send"]
  endin <- matches.pos[,"Qend"] >= matches.pos[,"Sstart"] & matches.pos[,"Qend"] <=  matches.pos[,"Send"]
  enddown <- matches.pos[,"Qend"] > matches.pos[,"Send"]

  ## Defines the different types of match:

  ## upstream:
  ## Q -----
  ## S       --------------
  upstream <- startup & endup

  ## inside:
  ## Q     -----
  ## S --------------
  inside <- startin & endin

  ## downstream:
  ## Q                -----
  ## S --------------
  downstream <- startdown & enddown

  ## olup:
  ## Q --------------
  ## S       -------------
  olup <- startup & endin

  ## oldown:
  ## Q       -------------
  ## S --------------
  oldown <- startin & enddown

  ## contained:
  ## Q --------------
  ## S     -----
  contained <- startup & enddown

  ## Match types in one vector
  match.type <- rep("", length(matches.pos[,"Qstart"]))

  match.type[upstream] <- "upstream"
  match.type[olup] <- "olup"
  match.type[inside] <- "inside"
  match.type[oldown] <- "oldown"
  match.type[downstream] <- "downstream"
  match.type[contained] <- "contained"


  ## Calculation of various Distances between the extremities of query and their hits
  dist <- cbind(
      SstartMinusQstart = (matches.pos[,"Sstart"] - matches.pos[,"Qstart"]),
      SendMinusQstart = (matches.pos[,"Send"] - matches.pos[,"Qstart"]),
      SendMinusQend = (matches.pos[,"Send"] - matches.pos[,"Qend"]),
      SstartMinusQend = (matches.pos[,"Sstart"] - matches.pos[,"Qend"])
  )


  ## Output type
  oldf <- cbind(data.frame(Qindex=olindex[,1], Sindex=olindex[,2], TypeOfMatch = match.type), dist)
  if(output=="df") {return(oldf)}
  if(output=="gr") {
    elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf[,-1])
    return(query)
  }
  ##  Run matchTypes function
  ##  Sample Data Set
  ##
#gr <- GRanges(seqnames = Rle(rep("Chr1",4)),
#    ranges = IRanges(c(1,2,4,6), end = c(3,7,5,8), names = letters[1:4]),
#    strand = Rle(rep("+",4)),
#    seqlengths = c(Chr1 = 8))
  ##     12345678
  ##     ||||||||
  ## 1 a --
  ## 2 b  ------
  ## 3 c    -
  ## 4 d       --

#query <- gr
#subject <- gr
#matchTypes(query, subject, output = "gr")
#matchTypes(query, subject, output = "df")
#matchTypes(query, subject, output = "df", maxgap = 1)

#strand(gr) <- strand = Rle(rep("-",4))
  ##     87654321
  ##     ||||||||
  ## 1 a --
  ## 2 b  ------
  ## 3 c    -
  ## 4 d       --

#matchTypes(query, subject, output = "df", maxgap = 1)

}


EBEMapping2Gene <- function(query, subject, type = "any" , maxgap = 1000,
    exclude.match.type = NULL,  DB.con = TALR, table = "EBE2Gene") {
  ## The function first map the elements in the GRanges query  that overlap
  ## with elements in the GRanges subject wih parameters type and maxgap
  ## used by the findOverlap function.
  ## The mapping found is written to a DB.con connection in table table.
  ## Note that the definition of "overlap" can be adjusted with maxgap and type parameters

  # Creates a data frame summarizing the matches between query and subject
  matches <- findGRMatchesTypes(query, subject, output = "df", maxgap = maxgap, type = type, select = "all")

  # Replaces indexes by the names of the corresponding elements from GRanges object
  EBEToGene <- cbind(
      data.frame(EBE_ID = as.integer(names(query)[matches$Qindex]),
          ensembl_gene_id = as.character(names(subject)[matches$Sindex]), stringsAsFactors = FALSE),
      matches[-1:-2])


  # Filter out the rows that have TypeOfMatch value = match.type
  EBEToGene <- EBEToGene[!EBEToGene$TypeOfMatch %in% exclude.match.type,]

  # Populate database with the resulting rows of the df
  dbWriteTable(conn = DB.con, name =  table, value = EBEToGene,
		  overwrite = FALSE, append = FALSE, row.names = FALSE)
}



########################################
#    More GenomicRanges methods        #
########################################

GRangeSubsetBy <- function(gr, columnName="type", value) {
  ## For a GRanges object that contain elementMetadata df
  ## Return a subset of the GRanges whose elementMetadata column
  ## specified by columnName matches the value(s) specified by
  ## the value argument
  subset(gr, subset = values(gr)[,columnName] %in% value)
}

# split(x, f, drop=FALSE): Splits x according to f to create a GRangesList object.
# If f is a list-like object then drop is ignored and f is treated as if it was rep(seq_len(length(f)), sapply(f, length))
# so the returned object has the same shape as f (it also receives the names of f).
# Otherwise, if f is not a list-like object, empty list elements are removed from the returned object if drop is TRUE.






lociAsGRange <- function(loci, offset = 0, downstream = TRUE , ID = "EBE_ID", chrom.sizes = os.chr.sizes) {
# This function converts a loci df returned by quering a genomic loci table from DB
# into a GenomicRanges object of the "GenomicRanges" package from Bioconductor.
# Positions in loci are processed as a function of the offset and downstream arguments.
# if offset is >0 the GenomcRanges object contains the loci flanks of size offset
# if donstream is TRUE, the downstream flanks are returned, the upstrea ones otherwise.

  library("GenomicRanges")

  # names(chrom.sizes) <- as.integer(names(chrom.sizes)) # type conversion to comply with GRange
  # chrom.sizes <- unlist(chrom.sizes) # type conversion to comply with GRange

  # Creates the GenomicRanges container
  gr <- GRanges(seqnames = as(as.character(loci$Chr), "Rle"),
      ranges = IRanges(start = loci$StartPos, end = loci$EndPos,
          names = unlist(subset(loci, select = ID), use.names = FALSE)),
      strand = strand(as.integer(loci$Strand)),
      seqlengths = chrom.sizes[names(chrom.sizes) %in% as.character(loci$Chr)] # Include Chromosome length if the Chromosome is in loci
  )
  if (!offset == 0) {gr <- flank(gr, width = offset,  start = !downstream)} # apply offset if necessary
  return(gr)
}

######################################################
##       Interraction with sqlite database         ###
##               DBI and SQL methods               ###
######################################################

dbListFullFields <- function(conn = TALR, table) {
	paste(table, dbListFields(conn = conn, table), sep = ".")
}



fetchDescr <- function(x) {
	# write the gene IDs of interest for querying the DB
	dbWriteTable(conn = TALR, name = "QueryIDs_TEMP",
			value = data.frame(id = x), row.names = FALSE, overwrite = TRUE, append = FALSE)
# get the corresponding descriptions in the DB
	descrip <- dbGetQuery(TALR, "SELECT QueryIDs_TEMP.id, TargetGenes.description
					FROM QueryIDs_TEMP LEFT JOIN TargetGenes ON QueryIDs_TEMP.id=TargetGenes.ensembl_gene_id")
	dbRemoveTable(TALR, "QueryIDs_TEMP") # remove temporary table
	return(descrip)
}


##############################################################
#         Microarray data analysis and storage in DB         #
#                                                            #
##############################################################

# Finding differentially expresssed genes in microarray data and save to DB

diffExpression2DB <- function (fit, contrast.matrix = fit$contrasts,
		min.lfc = log2(2), max.Pval = 0.05,
		DB.con = TALR, table = "ArrayData", return.df = FALSE) {
	## For each comparaison in the provided contrast.matrix,
	## this function save the list of the differentially expressed probesets
	## contained in the fit argument into the specified DB.con table under the name provided by table.
	## It corrects for multiple testing adjustment based on Benjamini and Hochberg's method
	## (false discovery rate).
	## Selected probeSets differ significantly from the hypothesis tested
	## by the constrast as specified by the max.adj.P.Val and min.lfc (minimum
	## log2-fold-change) arguments.
	## The fit argument should be an object of class MArrayLM as produced by lmFit and eBayes

	if(!is(fit, "MArrayLM")) {stop("The fit argument is not of type 'MArrayLM'")}

	global.df <- data.frame()
	for (Contrasts in colnames(contrast.matrix)) {
		# builds a data frame with the results of eBayes
		diff.expr <- topTable(fit, coef = Contrasts, n = Inf, adjust = "BH")
		# filters probesets that meet the criteria specified in the arguments
		# Actually, topTable does just that...
		diff.expr <- subset(diff.expr, subset = abs(logFC) >= min.lfc & P.Value <= max.Pval )
		# Binds the corresponding contrast descriptor to the df if not empty and
		# add a ranking numer based on expression fold values in each contrast
		if(nrow(diff.expr) != 0) {
			ranking <- numeric()
			ranking[order(diff.expr$logFC, decreasing = TRUE)] <- 1:nrow(diff.expr)
			diff.expr <- cbind(diff.expr, logFCRank = ranking, Contrasts)} else {next}
		global.df <- rbind(global.df, diff.expr) # append the diff exp probes for the current contrast to the global data frame.

	}
	# Saves the df in the DB object
	dbWriteTable(conn = DB.con, value = global.df, name = table, row.names = FALSE, append = TRUE, overwrite = FALSE)
	if (return.df) return(global.df)

}


# imports a desing matrix for limma annalysis
importDesign <- function(file = "samples.txt") {
	### !!! NEW VERSION NOT TESTED WITH ALL EXP !!!!
	# imports a desing matrix for limma annalysis
	treatments <- read.delim(file)
	fileNames <- treatments[,1]
	if (ncol(treatments) == 2) {treatments <- as.data.frame(treatments[[-1]])} else {
		treatments <- treatments[,-1]}
	rownames(treatments) <- fileNames
	return(treatments)
}

# Perform the workflow for diff expression with limma using standard parameters
# and write the results to the approriate table in the DB
processRMAToDB <- function(rma = "RMA_processed_data.RData",
		contrasts = c("XooBAI3-H2O", "XooBAI3_Dtalc-H2O", "XooBAI3-XooBAI3_Dtalc"),
		design,
		DB.con = TALR, table = "ArrayData") {
	# Perform the workflow for diff expression with limma using standard parameters
	# and write the results to the approriate table in the DB
# load the processed.data object from an R data file
	load(rma)
	## Make sure the samples in the prelminary design matrix
	## are properly sorted:
#	print(all(sampleNames(processed.data) == rownames(treatments))) # must be TRUE
# fit a linear model to the experimental expression data
	fit <- lmFit(processed.data, design)
# definning a constrast matrix describing
# parameter comparaisons of interest
	contrast.matrix <- do.call(makeContrasts, c(as.list(contrasts),list(levels = design)))

# Computing estimated coefficients and standard errors for
# these contrasts from original linear model fit
	fit <- contrasts.fit(fit, contrast.matrix)
# Uses and empirical Bayes method to moderate standard errors of the estimated log-fold changes
	fit <- eBayes(fit)

# Save differentially expressed probesets into DB for each tested contrast
	diffExpression2DB(fit, min.lfc = log2(2), max.Pval = 1,
			DB.con = DB.con, table = table, return.df = FALSE)
}




