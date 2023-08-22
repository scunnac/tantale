
.getTalePartsFromAFile <- function(fasta) {
  if (grepl("TALE_Protein_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readAAStringSet(fasta)
  if (grepl("TALE_DNA_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readDNAStringSet(fasta)
  if (length(taleStrings) == 0L) {
    tibble::tibble(arrayIDs = character(), domainType = character(),
                   "position" = character(), "string" = character(),
                   "sourceDirectory" = character())
  } else {
    tibble::tibble(arrayID = gsub("(.*): .*", "\\1", names(taleStrings)),
                   domainType = gsub(".*: (.*?)[ ]?[0-9]{0,}$", "\\1", names(taleStrings)),
                   position = 1:length(taleStrings),
                   string = as.character(taleStrings),
                   sourceDirectory = dirname(fasta)
    )
  }
}

#' Fetch Annotale parts from a tellTale output directory.
#'
#' @description
#' 
#' This function get sequences from Annotale "TALE_Protein_parts.fasta" and "TALE_DNA_parts.fasta"
#' files from a specified \code{\link[tantale:tellTale]{tellTale}} run output directory
#' and returns a tibble. Each row describes a domain from a tale array and includes the 'arrayID',
#'  the id of the sequence where this array was found, the 'domainType' (type of domain, repeat,
#'  N-term or C-term), the position of the domain inside the array,
#' the DNA of the corresponding domain and the RVD and amino acid sequences if relevant.
#'
#' @param tellTaleOutDir Path to a \code{\link[tantale:tellTale]{tellTale}} run output directory
#' @return A tibble.
#' @export
getTaleParts <- function(tellTaleOutDir) {
  protPartsFiles <- list.files(tellTaleOutDir, "TALE_Protein_parts.fasta", recursive = T, full.names = T)
  dnaPartsFiles <- list.files(tellTaleOutDir, "TALE_DNA_parts.fasta", recursive = T, full.names = T)
  rvds <- fa2liststr(list.files(tellTaleOutDir, "rvdSequences.fas", recursive = T, full.names = T)) %>%
    lapply(function(x) tibble::tibble(rvd = x, position = 1:length(x) )) %>%
    dplyr::bind_rows(.id = "arrayID")
  
  taleProtString <- lapply(protPartsFiles, .getTalePartsFromAFile) %>% dplyr::bind_rows()
  taleDnaString <- lapply(dnaPartsFiles, .getTalePartsFromAFile) %>% dplyr::bind_rows()
  taleParts <- dplyr::full_join(taleDnaString %>% dplyr::rename(dnaSeq = string),
                                taleProtString %>% dplyr::rename(aaSeq = string),
                                by = c("arrayID", "domainType", "position", "sourceDirectory"))
  
  taleParts %<>% dplyr::group_by(aaSeq) %>%
    dplyr::mutate(domCode = as.character(dplyr::cur_group_id())) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(domCode = dplyr::if_else(is.na(aaSeq), as.character(NA), domCode),
                  aaSeq = gsub("[*]", "", aaSeq)
                  )
  
  taleParts %<>% dplyr::left_join(rvds, by = c("arrayID", "position"))
  
  taleParts %<>% dplyr::left_join(
  readr::read_tsv(list.files(tellTaleOutDir, "hitsReport.tsv", recursive = T, full.names = T),
                  show_col_types = FALSE) %>%
    dplyr::select(arrayID, seqnames) %>%
    dplyr::distinct(), by = "arrayID"
  )
  return(taleParts)
}

identSubMat <- matrix(data = rep(0, times = length(Biostrings::AA_PROTEINOGENIC) ^ 2),
                      nrow = length(Biostrings::AA_PROTEINOGENIC),
                      ncol = length(Biostrings::AA_PROTEINOGENIC),
                      dimnames = list(Biostrings::AA_PROTEINOGENIC, Biostrings::AA_PROTEINOGENIC))
diag(identSubMat) <- 1


.distalPairwiseAlign <- function(partAaStringSet, ncores = 1) {
  bpparam <- BiocParallel::MulticoreParam(ncores, progressbar = TRUE)
  pairAlignScores <- BiocParallel::bplapply(1:length(partAaStringSet), function(i) {
    singleSubAln <- Biostrings::pairwiseAlignment(pattern = partAaStringSet,
                                                  subject = partAaStringSet[i],
                                                  substitutionMatrix = identSubMat, #"BLOSUM62",
                                                  gapOpening = 1, gapExtension = 0.5,
                                                  type = "global", scoreOnly = FALSE)
    tibble::tibble(subj = names(partAaStringSet[i]),
                   pattern = names(alignedPattern(singleSubAln)),
                   score = Biostrings::score(singleSubAln),
                   nedit = Biostrings::nmismatch(singleSubAln)
    )
  }) %>%
    dplyr::bind_rows()
  
  pairAlignScores %<>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      maxLength = max(nchar(partAaStringSet[subj]), nchar(partAaStringSet[pattern])),
      # This is an approximate equivalent of how Alvaro computed dissimilarity in distal
      Dissim = 100 - 100 * (maxLength - score) / maxLength,
      Dissim = ifelse(Dissim < 0, 100, 100 - Dissim)
    ) %>%
    dplyr::ungroup()
  #pairAlignScores$Dissim %>%  hist(breaks = 40)
  return(pairAlignScores)
}

condaBinPath <- "/home/cunnac/bin/miniconda3/condabin/conda"
taleParts <- getTaleParts(system.file("extdata", "tellTaleExampleOutput", package = "tantale", mustWork = T)) %>%
  dplyr::mutate(partId = paste(arrayID, position, sep = "_"))

partAaStringSet <- Biostrings::AAStringSet(taleParts$aaSeq)
names(partAaStringSet) <- taleParts$partId


.distalPairwiseAlign2 <- function(partAaStringSet, ncores = 1, condaBinPath = "auto") {
  outdir <- tempdir(check = TRUE)
  partAaStringSetFile <- file.path(outdir, "taleAsParts.fsa")
  mmseq2DbPath <- file.path(outdir, 'mmseq2DB')
  prefDbPath <- file.path(outdir, 'resultDB_pref')
  alnDbPath <- file.path(outdir, 'resultDB_aln')
  alnTabFile <- file.path(outdir, 'alnRes.tab')
  
  df <- expand.grid(names(partAaStringSet), names(partAaStringSet), stringsAsFactors = FALSE) %>% tibble::as_tibble()
  colnames(df) <- c("query", "target")
  
  Biostrings::writeXStringSet(partAaStringSet, filepath = partAaStringSetFile)
  mmseq2Cmd <- glue::glue("mmseqs createdb {partAaStringSetFile} -v 3 {mmseq2DbPath};",
                          "mmseqs prefilter {mmseq2DbPath} {mmseq2DbPath} {prefDbPath}",
                          "-v 3 --max-seqs 1000000 -s 7.5 --add-self-matches 1",
                          "--cov-mode 0;",
                          "mmseqs align {mmseq2DbPath} {mmseq2DbPath} {prefDbPath} {alnDbPath}",
                          "-v 3 --threads {ncores} --add-self-matches 1 --min-seq-id 0 --cov-mode 0",
                          "-a 1 --alignment-output-mode 0;",
                          "mmseqs convertalis {mmseq2DbPath} {mmseq2DbPath} {alnDbPath} {alnTabFile}",
                          "--format-mode 4",
                          "--format-output query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,mismatch,qcov,tcov",
                          .sep = " ")

  if (!as.logical(createBioPerlEnv(condaBinPath = condaBinPath))) {
    logger::log_info("Invoking mmseq2 using the following command:\n {stringr::str_wrap(mmseq2Cmd, 80)}")
    res <- systemInCondaEnv(envName = "perlforal",
                            condaBinPath = condaBinPath,
                            command = mmseq2Cmd,
                            ignore.stdout = T)
  } else {
    stop("Could not create the Perl infrastructure on your machine to run mmseq2...")
  }
  
  pairAlignScores <- readr::read_tsv(alnTabFile) %>%
    dplyr::mutate(query = as.character(query), target = as.character(target))
  pairAlignScores <- dplyr::left_join(df, pairAlignScores) %>%
    dplyr::rename(pattern = target, subj = query) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Dissim = ifelse(is.na(pident), 100, 100 - pident*(max(qcov, tcov))))
  
  # Check that there is a single pairwise alignment records for all possible parts pairs
  partCombinCounts <- pairAlignScores %>% dplyr::select(pattern, subj) %>%
    dplyr::count() %>%
    dplyr::pull(n)
  if (!all(partCombinCounts == 1L) || length(partCombinCounts) != length(partAaStringSet)^2) {
    stop("There is not a single pairwise alignment records for all possible parts pairs",)
  }
  # pairAlignScores %>% dplyr::group_by(query, target) %>% dplyr::tally() %>%
  #   dplyr::ungroup() %>% dplyr::filter(n != 1)
  # pairAlignScores %>% dplyr::filter(Dissim > 18, Dissim < 80)
  # pairAlignScores$Dissim %>%  hist(breaks = 100)
  
  return(pairAlignScores)
}




#' Emulate DisTal in R
#' @description
#' This is meant to approximate the results of DisTal in R and is very similar to
#' \code{\link[tantale:runDistal]{runDistal}}.
#' It still uses the Arlem binary just like DisTal but performs the rest of the operations
#' with R support. Even with parallel computing it is slower than the original perl code...
#' Computational bottleneck occurs probably at the stage where repeat amino acid sequences are
#' systematically pairwise aligned.
#' 
#' @param taleParts a table of TALE parts as returned by the \code{\link[tantale:getTaleParts]{getTaleParts}} function.
#' @param repeats.cluster.h.cut numeric value to cut the hierarchical clustering tree of the repeat.
#' @param pairwiseAlnMethod Specify the underlying approach for computing pairwise similarities between
#' TALE parts amino acid sequences. Must be either "Biostrings" or "mmseq2" 
#' @param condaBinPath Path to your Conda binary file if you need to specify a
#'   path different from the one that is automatically searched by the
#'   reticulate package functions.
#' @return A list with DisTal output components: 
#' \itemize{
#'   \item repeats.code: a data frame of the unique repeat AA sequences and their numeric codes
#'   \item coded.repeats.str: a list of repeat-coded TALE strings
#'   \item repeat.similarity:  a long, three columns data frame with pairwise similarity scores between repeats
#'   \item tal.similarity: a three columns Tals similarity table with pairwise similarity scores between TALEs
#'   \item tree: Newick format neighbor-joining tree of TALs constructed based on TALEs similarity
#'   \item repeats.cluster: a data frame containing repeat code and repeat clusters.
#' }
#' @export
distalr <- function(taleParts, repeats.cluster.h.cut = 10, ncores = 1,
                    pairwiseAlnMethod = "mmseq2", condaBinPath = "auto") {
  
  #### Reality checks ####
  
  ## Make sure we are dealing only with parts that have defined protein sequences.
  if (any(is.na(taleParts$aaSeq) | taleParts$aaSeq == "")) {
    stop("It seems that some of the provided TALE parts miss an amino acid sequence. Cannot proceed!")
    } 
  
  ## Make sure that arrayID - position combinations are unique
  # in case someone would not have made arrayIDs unique before
  # mixing tale predictions from several genomes...
  arayPosCombinCounts <- taleParts %>%
    dplyr::group_by(arrayID, position) %>%
    dplyr::count() %>%
    dplyr::pull(n)
  if (!all(arayPosCombinCounts == 1L)) {
    stop("Your tale arrays identifers are probably not unique.",
         "\n",
         "Make sure that there is only one part per position, arrayID and DNA sequence name 'seqnames'.")
  }
  
  #### Assemble repeat code strings and write in a file for arlem ####
  repeatStrings <- taleParts %>% dplyr::group_by(arrayID) %>%
    dplyr::arrange(position) %>%
    dplyr::summarise(repeatString = paste(domCode, collapse = " "),
                     posString = paste(position, collapse = " "))
  codesSeqSet <- Biostrings::BStringSet(repeatStrings$repeatString)
  names(codesSeqSet) <- repeatStrings$arrayID
  
  codesSeqsfile <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(codesSeqSet, codesSeqsfile)
  
  #### Prepare Arlem cfile with systematic pairwise distances between 'repeat' units. ####
  taleAaParts <- Biostrings::AAStringSet(taleParts$aaSeq)
  names(taleAaParts) <- taleParts$domCode
  
  # Get pairwise repeat aa sequence dissimilarity scores in a long tibble
  if (pairwiseAlnMethod == "mmseq2") {
    dissimLong <- .distalPairwiseAlign2(partAaStringSet = unique(taleAaParts), ncores = ncores,
                                        condaBinPath = condaBinPath) 
  } else if (pairwiseAlnMethod == "Biostrings") {
    dissimLong <- .distalPairwiseAlign(partAaStringSet = unique(taleAaParts), ncores = ncores) 
  } else {
    stop("'pairwiseAlnMethod' parameter must be either 'Biostrings' or 'mmseq2'")
  }

  dissimLong %<>% dplyr::mutate(Sim = 100 - Dissim,
                                Dissim = round(Dissim) %>% as.character())
  
  # Convert to symetric matrix
  dissimMat <- reshape2::acast(dissimLong, formula = subj ~ pattern, value.var = "Dissim")
  stopifnot(nrow(dissimMat) == ncol(dissimMat))
  # reorder row and colnames because I suspect arlem expect them in increasing order
  dissimMat <- dissimMat[rownames(dissimMat) %>% as.numeric() %>% order(),
                         colnames(dissimMat) %>% as.numeric() %>% order()]
  # Get parameters for arlem
  TypeNo <- glue::glue("# Type no. ", nrow(dissimMat))
  Types <- glue::glue("# Types ", paste(1:nrow(dissimMat), collapse = " "))
  
  # 'erase' lower triangle and diag
  dissimMat[lower.tri(dissimMat)] <- ""
  diag(dissimMat) <- ""
  
  # Convert mat rows to strings of space separated values
  dissimMatLines <- apply(dissimMat, 1, function(row) {
    string <- paste(row, collapse = " ")
    gsub("^[ ]+", "", string)
  })
  # remove empty last line
  dissimMatLines <- dissimMatLines[1:(length(dissimMatLines) - 1)]
  # Add the arlem stuff
  dissimMatLines <- c(TypeNo, Types,
                      "# Indel align 10", "# Indel hist 10", "# Dup 10",
                      "# matrix",
                      dissimMatLines)
  
  # write lines in a temp cfile
  cfile <- tempfile()
  writeLines(cfile, text = dissimMatLines)
  
  #### run arlem with a system call and parse std output ####
  arlemPath <- system.file("tools", "DisTAL1.2_MultipleAlignment","arlem", package = "tantale", mustWork = T)
  arlemCmd <- glue::glue("{arlemPath} -f {codesSeqsfile} -cfile {cfile} -align -insert -showalign")
  
  arlemRes <- grep("Score of aligning Seq:", system(arlemCmd, intern = TRUE), value = TRUE)
  arlemScores <- gsub("Score of aligning Seq:([0-9]*), Seq:([0-9]*) =([0-9]*)", "\\1|\\2|\\3", arlemRes)
  arlemScores <- strsplit(arlemScores, split = "\\|")
  arlemScores <- lapply(arlemScores, function(s) {t(as.matrix(as.numeric(s)))}) %>%
    do.call(rbind, .)
  colnames(arlemScores) <- c("TAL1", "TAL2", "arlemScore")
  
  #### Compute normalized arlem scores and include arrayIDs rather than arlem index ####
  arrayLengths <- taleParts %>% dplyr::group_by(arrayID) %>% dplyr::count()
  
  normArlemScoresTble <- arlemScores %>% tibble::as_tibble() %>%
    dplyr::mutate(
      TAL1 = names(codesSeqSet)[TAL1 + 1],
      TAL2 = names(codesSeqSet)[TAL2 + 1]
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      maxLength = max(arrayLengths$n[arrayLengths$arrayID == TAL1],
                      arrayLengths$n[arrayLengths$arrayID == TAL2]),
      normArlemScore = arlemScore/maxLength,
      Sim = 100 - normArlemScore
    ) %>%
    dplyr::ungroup()
  
  #### return a list of results emulating the return value of tantale::runDistal ####
  outputlist <- list("repeats.code" = taleParts %>%
                       dplyr::group_by(domCode, aaSeq) %>%
                       dplyr::count() %>%
                       dplyr::rename(code = domCode, "AA Seq"  = aaSeq) %>%
                       dplyr::mutate(code = as.integer(code)) %>%
                       dplyr::select(-n) %>%
                       dplyr::ungroup(),
                     "coded.repeats.str" = fa2liststr(codesSeqsfile),
                     "repeat.similarity" = dissimLong %>% dplyr::rename(RepU1 = subj, RepU2 = pattern),
                     "tal.similarity" = normArlemScoresTble,
                     "repeats.cluster" = clusterRep(
                       repeatSimMat = reshape2::acast(dissimLong, formula = subj ~ pattern, value.var = "Sim"),
                       repeats.cluster.h.cut = 10
                     )
  )
  return(outputlist)
}





