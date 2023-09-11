
.getTalePartsFromAFile <- function(fasta) {
  if (grepl("TALE_Protein_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readAAStringSet(fasta)
  if (grepl("TALE_DNA_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readDNAStringSet(fasta)
  if (length(taleStrings) == 0L) {
    logger::log_warn("No part sequence found in: {fasta}. Returning an empty tibble.")
    tbl <- tibble::tibble(arrayIDs = character(),
                   domainType = character(),
                   positionInArray = character(),
                   positionInCrd = character(),
                   string = character(),
                   sourceDirectory = character())
    return(tbl)
  } 
  tbl <- tibble::tibble(arrayID = gsub("(.*): .*", "\\1", names(taleStrings)),
                        domainType = gsub(".*: (.*?)[ ]?[0-9]{0,}$", "\\1", names(taleStrings)),
                        positionInCrd = gsub(".*: repeat[ ]([0-9]{0,})$", "\\1", names(taleStrings)) %>%
                          as.integer() %>%
                          suppressWarnings(),
                        string = as.character(taleStrings) %>% as.vector(),
                        sourceDirectory = dirname(fasta)
  )
  missingTerm <- setdiff(c("N-terminus", "C-terminus"), unique(tbl$domainType))
  if (length(missingTerm) != 0L) {
    logger::log_warn("Array {unique(tbl$arrayID)} is missing a {missingTerm} domain in {fasta}")
    warning()
    missingTerm <- tibble::tibble(arrayID = unique(tbl$arrayID),
                                  domainType = missingTerm,
                                  positionInCrd = NA,
                                  string = NA,
                                  sourceDirectory = dirname(fasta))
    tbl <- dplyr::bind_rows(tbl, missingTerm)
    logger::skip_formatter(as.character(knitr::kable(missingTerm))) %>%
      logger::log_debug()
  }
  tbl %<>% 
    dplyr::rowwise() %>%
    dplyr::mutate(positionInArray = switch(domainType,
                                                 `N-terminus` = 1,
                                                 `repeat` = positionInCrd + 1,
                                                 `C-terminus` = nrow(tbl)
                                                 ))
  return(tbl)
}




.getRvdsFromAnAnnotaleFile <- function(fasta) {
  if (!grepl("TALE_RVDs.fasta", basename(fasta))) {
    logger::log_error("The provided file does not seem to be an AnnoTALE RVDs file: {fasta}")
    stop()
  } else {
    rvdTble <- toListOfSplitedStr(fasta) %>%
      lapply(function(x) tibble::tibble(string = x,
                                        positionInCrd = 1:length(x))
             ) %>%
      dplyr::bind_rows(.id = "arrayID")
    rvdTble <- rvdTble %>% dplyr::mutate(sourceDirectory = dirname(fasta),
                                         domainType = "repeat")
  }
  return(rvdTble)
}


#' Fetch Annotale parts from a tellTale output directory.
#'
#' @description
#'
#' This function get sequences from tellTale "rvdSequences.fas" and AnnoTALE
#' "TALE_Protein_parts.fasta" and "TALE_DNA_parts.fasta" files from a SINGLE
#' \code{\link[tantale:tellTale]{tellTale}} run output directory and returns a
#' tibble. Each row describes a domain from a tale array and includes the
#' 'arrayID', the id of the sequence where this array was found, the
#' 'domainType' (type of domain, repeat, N-term or C-term), the position of the
#' domain inside the array, the DNA of the corresponding domain and the RVD and
#' amino acid sequences if relevant.
#'
#'
#' **IMPORTANT**: telltale MUST have been run with the appendExtremityCodes = TRUE
#'
#' @param tellTaleOutDir Path to a \code{\link[tantale:tellTale]{tellTale}} run
#'   output directory
#' @return A tibble.
#' @export
getTaleParts <- function(tellTaleOutDir) {
  # Get info from telltale output dir
  # !!!! arrayID are assumed to be unique !!!!
  protPartsFiles <- list.files(tellTaleOutDir, "TALE_Protein_parts.fasta", recursive = T, full.names = T)
  dnaPartsFiles <- list.files(tellTaleOutDir, "TALE_DNA_parts.fasta", recursive = T, full.names = T)
  if (tellTaleOutDir %>% dirname() %>% unique() %>% length() != 1L) {
    log_error("The provided path most likely does not correspond to a SINGLE tellTale output directory.")
  }
  # Fetch info from annotale/telltale files with .getTalePartsFromAFile
  taleProtString <- lapply(protPartsFiles, .getTalePartsFromAFile) %>% dplyr::bind_rows()
  taleDnaString <- lapply(dnaPartsFiles, .getTalePartsFromAFile) %>% dplyr::bind_rows()
  #stopifnot(nrow(taleProtString) == nrow(taleDnaString))
  # Join info in a table with one domain per row
  taleParts <- dplyr::full_join(taleDnaString %>% dplyr::rename(dnaSeq = string),
                                taleProtString %>% dplyr::rename(aaSeq = string),
                                by = c("arrayID", "domainType", "positionInArray", "positionInCrd", "sourceDirectory"),
                                relationship = "one-to-one") %>%
    dplyr::mutate(aaSeq = gsub("[*]", "", aaSeq))

  # Get RVDs
  # NOTE: could be easier to get the RVDs directly from AnnoTALE output with
  # .getRvdsFromAnAnnotaleFile() but I currently feel that it is good to
  # be aware of disagreements between AnnoTALE diagnostic on terminal domains presence in AA seqs
  # and nhmmer diagnostic on terminal domains CDS presence on DNA.
  rvds <- toListOfSplitedStr(list.files(path = tellTaleOutDir,
                                        pattern = "rvdSequences.fas",
                                        recursive = F,
                                        full.names = T)
                             ) %>%
    lapply(function(x) tibble::tibble(rvd = x, positionInArray = 1:length(x) )) %>%
    dplyr::bind_rows(.id = "arrayID")  
  anchorCodes <- c("NTERM", "CTERM", "XXXXX")
  
  
  # Some checks on the consistency between parts and rvd sequences
  # if nhmmer did not report on a C-Term CDS, the corresponding domain
  # "CTERM" tag will not be written in the rvd slot of the table.
  arraysConsistency <- dplyr::full_join(taleParts %>% dplyr::count(arrayID, name = "AnnoTALELength"),
                                        rvds %>% dplyr::count(arrayID, name = "rvdFileLength"),
                                        by = dplyr::join_by(arrayID)) %>%
    dplyr::mutate(sameLength = AnnoTALELength == rvdFileLength)
  
  if (any(is.na(arraysConsistency$sameLength))) {
    logger::log_warn("There are mismatches in array IDs between rvd seq file and AnnoTALE parts files:")
    logger::skip_formatter(as.character(knitr::kable(arraysConsistency %>% dplyr::filter(is.na(sameLength))))) %>%
      logger::log_warn()
    warning()
  } else if (!all(arraysConsistency$sameLength, na.rm = TRUE)) {
    logger::log_error("Array lengths are inconsistent between rvd ",
                      "seq file and AnnoTALE parts files:")
    logger::skip_formatter(as.character(knitr::kable(arraysConsistency %>% dplyr::filter(!sameLength)))) %>%
      logger::log_error()
    stop()
  }
  
  # Include RVDs in the talParts tibble
  taleParts <- dplyr::left_join(taleParts, 
                                rvds,
                                by = c("arrayID", "positionInArray"),
                                unmatched = "drop", relationship = "one-to-one")
  # Include seqnames in the talParts tibble
  taleParts %<>% dplyr::left_join(
    readr::read_tsv(list.files(tellTaleOutDir, "hitsReport.tsv", recursive = T, full.names = T),
                    show_col_types = FALSE) %>%
      dplyr::select(arrayID, seqnames) %>%
      dplyr::distinct(),
    by = "arrayID", relationship = "many-to-one"
  )
  # Check talparts
  partsWithMissingAaSeq <- taleParts %>% dplyr::filter(is.na(aaSeq)) %>% dplyr::pull(arrayID) %>% unique()
  partsWithMissingDnaSeq <- taleParts %>% dplyr::filter(is.na(dnaSeq)) %>% dplyr::pull(arrayID) %>% unique()
  partsWithMissingRvdSeq <- taleParts %>% dplyr::filter(is.na(rvd)) %>% dplyr::pull(arrayID) %>% unique()
  if (any(sapply(list(partsWithMissingAaSeq, partsWithMissingDnaSeq, partsWithMissingRvdSeq), length) != 0L)) {
    logger::log_warn("Be aware that the output taleParts tibble has records with missing sequences:")
    logger::log_warn("{unique(c(partsWithMissingAaSeq, partsWithMissingDnaSeq, partsWithMissingRvdSeq))}")
    taleParts %>%
      #dplyr::select(arrayID, domainType, positionInArray, sourceDirectory) %>%
      dplyr::filter(arrayID %in% c(partsWithMissingAaSeq, partsWithMissingDnaSeq, partsWithMissingRvdSeq)) %>%
      knitr::kable() %>% as.character() %>%
      logger::skip_formatter() %>% logger::log_debug()
    warning()
  }
  return(taleParts)
}

identSubMat <- matrix(data = rep(0, times = length(Biostrings::AA_PROTEINOGENIC) ^ 2),
                      nrow = length(Biostrings::AA_PROTEINOGENIC),
                      ncol = length(Biostrings::AA_PROTEINOGENIC),
                      dimnames = list(Biostrings::AA_PROTEINOGENIC, Biostrings::AA_PROTEINOGENIC))
diag(identSubMat) <- 1


.distalPairwiseAlign <- function(partAaStringSet, ncores = 1) {
  bpparam <- BiocParallel::MulticoreParam(ncores, progressbar = TRUE)
  
  if (anyDuplicated(names(partAaStringSet))) {
    stop("Parts in the provided input have duplicated names. Cannot proceeed...")
  }
  
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
      Dissim = ifelse(Dissim < 0, 100, 100 - Dissim),
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-maxLength)
  
  # Check the pairAlignScores tibble
  .checkPairAlignTble(pairAlignScores = pairAlignScores, partAaStringSet = partAaStringSet)
  
  return(pairAlignScores)
}


.distalPairwiseAlign2 <- function(partAaStringSet, ncores = 1, condaBinPath = "auto") {
  outdir <- tempfile(pattern = "distalPairwiseAlign2")
  dir.exists(outdir) || dir.create(outdir, recursive = TRUE)
  partAaStringSetFile <- file.path(outdir, "taleAsParts.fsa")
  mmseq2DbPath <- file.path(outdir, 'mmseq2DB')
  prefDbPath <- file.path(outdir, 'resultDB_pref')
  alnDbPath <- file.path(outdir, 'resultDB_aln')
  alnTabFile <- file.path(outdir, 'alnRes.tab')
  
  df <- expand.grid(names(partAaStringSet),
                    names(partAaStringSet),
                    stringsAsFactors = FALSE
                    ) %>%
    tibble::as_tibble()
  colnames(df) <- c("query", "target")
  if(anyDuplicated(df) != 0) stop("The provided sequences must have unique names.")
  
  Biostrings::writeXStringSet(partAaStringSet, filepath = partAaStringSetFile)
  
  mmseq2createdb <- glue::glue("mmseqs createdb {partAaStringSetFile} {mmseq2DbPath}")

  mmseq2prefilter <- glue::glue("mmseqs prefilter {mmseq2DbPath} {mmseq2DbPath} {prefDbPath}",
                               "-v 3 --threads {max(floor(ncores/2), 1)} --max-seqs 1000 -s 7.5 --add-self-matches 1",
                               "--cov-mode 0", .sep = " ")
  
  mmseq2align <- glue::glue("mmseqs align {mmseq2DbPath} {mmseq2DbPath} {prefDbPath} {alnDbPath}",
                               "-v 3 --threads {ncores} --add-self-matches 1 --min-seq-id 0",
                               "--cov-mode 0 --gap-open aa:11,nucl:5 --gap-extend aa:1,nucl:2",
                               "-a 1 --alignment-mode 3 --alignment-output-mode 0 --seq-id-mode 1",
                               .sep = " ")
  
  mmseq2convertalis <- glue::glue("mmseqs convertalis {mmseq2DbPath} {mmseq2DbPath} {alnDbPath} {alnTabFile}",
                               "--format-mode 4 -v 3",
                               "--format-output query,target,evalue,raw,pident,nident,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qcov,tcov",
                               .sep = " ")
  
  if (!as.logical(createTantaleEnv(condaBinPath = condaBinPath))) {
    logger::log_debug("Invoking mmseq2 using the following command:\n {stringr::str_wrap(mmseq2Cmd, 80)}")
    res <- systemInCondaEnv(envName = "tantale",
                            condaBinPath = condaBinPath,
                            command = mmseq2createdb,
                            ignore.stdout = F)
    res <- systemInCondaEnv(envName = "tantale",
                            condaBinPath = condaBinPath,
                            command = mmseq2prefilter,
                            ignore.stdout = F)
    res <- systemInCondaEnv(envName = "tantale",
                            condaBinPath = condaBinPath,
                            command = mmseq2align,
                            ignore.stdout = F)
    res <- systemInCondaEnv(envName = "tantale",
                            condaBinPath = condaBinPath,
                            command = mmseq2convertalis,
                            ignore.stdout = F)
  } else {
    stop("Could not create the tantale conda environment on your machine to run mmseq2...")
  }
  

  pairAlignScores <- readr::read_tsv(alnTabFile, show_col_types = FALSE) %>%
    dplyr::mutate(target = as.character(target), query = as.character(query)) %>%
    dplyr::group_by(target, query) %>%
    dplyr::slice_max(raw, n = 1, with_ties = FALSE)
  pairAlignScores <- dplyr::left_join(df, pairAlignScores) %>%
    dplyr::rename(subj = target, pattern = query) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Dissim = ifelse(is.na(pident), 100, 100 - pident*min(qcov,tcov))) %>%
    #dplyr::mutate(Dissim = ifelse(is.na(pident), 100, 100 - pident)) %>%
    dplyr::ungroup()

  # Check the pairAlignScores tibble
  .checkPairAlignTble(pairAlignScores = pairAlignScores, partAaStringSet = partAaStringSet)
  unlink(outdir, recursive = TRUE)
  return(pairAlignScores)
}




.distalPairwiseAlign3 <- function(partAaStringSet, ncores = 1) {
  if (anyDuplicated(names(partAaStringSet))) {
    stop("Parts in the provided input have duplicated names. Cannot proceeed...")
  }
  msa <- DECIPHER::AlignSeqs(myXStringSet = partAaStringSet, normPower = 0,
                             processors = ncores, verbose = FALSE) %>%
    DECIPHER::StaggerAlignment(fullLength = TRUE, processors = ncores, verbose = FALSE)
  distMat <- DECIPHER::DistanceMatrix(msa, method = "longest",
                                      includeTerminalGaps = TRUE,
                                      processors = ncores, verbose = FALSE)
  pairAlignScores <- reshape2::melt(as.matrix(distMat)) %>% tibble::as_tibble()
  colnames(pairAlignScores) <- c("pattern", "subj", "Dissim")
  pairAlignScores %<>% dplyr::mutate(pattern = as.character(pattern), subj = as.character(subj))
  pairAlignScores %<>% dplyr::mutate(Dissim = Dissim*100) %>%
    dplyr::ungroup()
  
  # Check the pairAlignScores tibble
  .checkPairAlignTble(pairAlignScores = pairAlignScores, partAaStringSet = partAaStringSet)
  
  return(pairAlignScores)
}





.checkPairAlignTble <- function(pairAlignScores, partAaStringSet) {
  partCombinCounts <- pairAlignScores %>% dplyr::select(pattern, subj) %>%
    dplyr::count(pattern, subj) %>%
    dplyr::pull(n)
  if (!all(partCombinCounts == 1L)) {
    stop("Some alignment pairs have more than one record...",)
  }
  if (length(names(partAaStringSet))^2 != nrow(pairAlignScores)) {
    stop("Some parts pairs are absent from the pairwise parts distance table")
  }
}

#' Report on potential 'pseudo TALEs' in a taleParts object
#' @description
#' NOT TESTED!!!!
#' This displays a compact but information rich view of the TALEs stored in a
#' taleParts object.
#' 
#' @param taleParts a table of TALE parts as returned by the
#' \code{\link[tantale:getTaleParts]{getTaleParts}} function or
#' \code{\link[tantale:distalr]{distalr}}
#' @param sanitize If \code{FALSE}, will return all the arrays with at least one 
#' part with a missing sequence. If \code{TRUE}, will return all the arrays that have
#' no part with a missing sequence.
#' 
#'
#' @return a taleParts object
#' @export
diagnoseTaleParts <- function(taleParts, sanitize = FALSE) {
  # Check talparts
  partsWithMissingAaSeq <- taleParts %>% dplyr::filter(is.na(aaSeq)) %>%
    dplyr::select(arrayID, sourceDirectory) %>%
    distinct()
  partsWithMissingDnaSeq <- taleParts %>% dplyr::filter(is.na(dnaSeq)) %>%
    dplyr::select(arrayID, sourceDirectory) %>%
    distinct()
  partsWithMissingRvdSeq <- taleParts %>% dplyr::filter(is.na(rvd)) %>%
    dplyr::select(arrayID, sourceDirectory) %>%
    distinct()
  problems <- dplyr::bind_rows(partsWithMissingRvdSeq,
                               partsWithMissingDnaSeq,
                               partsWithMissingAaSeq
                               ) %>%
    distinct()
  pseudoTales <- dplyr::left_join(problems, taleParts, relationship = "one-to-many") %>%
    dplyr::arrange(sourceDirectory, arrayID, positionInArray)
  if (nrow(problems) != 0L) {
    logger::log_warn("Be aware that the output taleParts tibble has records with missing sequences")
    warning()
  }
  if (!sanitize) {
    pseudoTales %>% return()
  } else {
    logger::log_info("Returning TALE arrays with no empty sequence parts")
    dplyr::setdiff(taleParts, pseudoTales) %>% return()
  }
}




#' Visualize TALE content in a taleParts object
#' @description
#' NOT TESTED!!!!
#' This displays a compact but information rich view of the TALEs stored in a
#' taleParts object.
#' 
#' @param taleParts a table of TALE parts as returned by the
#' \code{\link[tantale:getTaleParts]{getTaleParts}} function or
#' \code{\link[tantale:distalr]{distalr}}
#'
#' @return The ggplot object
#' @export
plotTaleComposition <- function(taleParts) {
  partsForPlots <- taleParts %>%
    mutate(label = if_else(domainType == "repeat", rvd, ""),
           aaSeqLength = factor(nchar(aaSeq))
    )
  p  <- partsForPlots %>% ggplot(mapping = aes(fill = aaSeqLength,
                                               color = domainType,
                                               label = label,
                                               y = arrayID,
                                               x = positionInArray),
                                 color = isNaAaSeq) +
    scale_color_viridis_d(option = "rocket") +
    scale_fill_discrete() +
    scale_x_continuous(breaks = 1:50, minor_breaks = NULL) +
    geom_point(shape = 21, size = 5, stroke = 0.9) +
    ggnewscale::new_scale_color() +
    ggnewscale::new_scale_fill() +
    geom_text(size = 2.1, color = "white") + 
    facet_grid(seqnames~ ., scales = "free_y", space = "free") +
    labs(title = "Overview of TALE composition by genome") +
    theme_light()
  print(p)
  return(p)
}


# taleParts <- readRDS("/home/cunnac/TEMP/talePartsForDistalr.rds")
# repeats.cluster.h.cut = 10
# ncores = 1
# pairwiseAlnMethod = "DECIPHER"
# condaBinPath = "/home/cunnac/bin/miniconda3/condabin/conda"

#' Emulate DisTal in R
#' @description
#' This is meant to approximate the results of DisTal in R and is very similar to
#' \code{\link[tantale:runDistal]{runDistal}}.
#' It still uses the Arlem binary just like DisTal but performs the rest of the operations
#' with R support and parallelization. Depending on the \code{pairwiseAlnMethod} parameter value it is 
#' much faster than the orginial Perl code and returns similar results. Please take a look at the vignette for
#' tips on how to use it properly.
#' 
#' @param taleParts a table of TALE parts as returned by the \code{\link[tantale:getTaleParts]{getTaleParts}} function.
#' @param repeats.cluster.h.cut numeric value to cut the hierarchical clustering tree of the repeat.
#' @param pairwiseAlnMethod Specify the underlying approach for computing pairwise similarities between
#' TALE parts amino acid sequences. Must be "Biostrings", "mmseq2" or "DECIPHER"
#' @param condaBinPath Path to your Conda binary file if you need to specify a
#'   path different from the one that is automatically searched by the
#'   reticulate package functions.
#' @return A list with DisTal output components: 
#' \itemize{
#'   \item taleParts: the original input tibble with a 'domCode' column corresponding to the unique distal
#'    'code' or label attached to a unique domain sequence. Thus all parts with this sequence will have the same
#'    'domCode'.
#'   \item repeats.code: a data frame of the unique repeat AA sequences and their numeric codes
#'   \item coded.repeats.str: a list of repeat-coded TALE strings
#'   \item repeat.similarity:  a long, three columns data frame with pairwise similarity scores between repeats
#'   \item tal.similarity: a three columns Tals similarity table with pairwise similarity scores between TALEs
#'   \item repeats.cluster: a data frame containing repeat code and repeat clusters.
#' }
#' @export
distalr <- function(taleParts, repeats.cluster.h.cut = 10, ncores = 1,
                    pairwiseAlnMethod = "DECIPHER", condaBinPath = "auto") {
  
  #### Reality checks ####
  
  ## Make sure we are dealing only with parts that have defined protein sequences.
  if (any(is.na(taleParts$aaSeq) | taleParts$aaSeq == "")) {
    logger::log_error("It seems that some of the provided TALE parts miss an amino acid sequence. Cannot proceed!")
    taleParts %>% dplyr::filter(is.na(aaSeq)) %>% 
      knitr::kable() %>% as.character() %>%
      logger::skip_formatter() %>% logger::log_error()
    stop()
  }
  if (any(is.na(taleParts$dnaSeq) | taleParts$dnaSeq == "")) {
    logger::log_warn("It seems that some of the provided TALE parts miss the DNA sequence!")
  } 
  
  ## Make sure that arrayID - position combinations are unique
  # in case someone would not have made arrayIDs unique before
  # mixing tale predictions from several genomes...
  arayPosCombinCounts <- taleParts %>%
    dplyr::group_by(arrayID, positionInArray) %>%
    dplyr::count() %>%
    dplyr::pull(n)
  if (!all(arayPosCombinCounts == 1L)) {
    logger::log_error("Your tale arrays identifers are probably not unique.",
         "\n",
         "Make sure that there is only one part per position per arrayID.")
    stop()
  }
  
  # Assign domain codes
  taleParts %<>% dplyr::group_by(aaSeq) %>%
    dplyr::mutate(domCode = dplyr::cur_group_id() %>% unlist() %>% as.character()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(domCode = dplyr::if_else(is.na(aaSeq), as.character(NA), domCode))
  
  
  #### Assemble repeat code strings and write in a file for arlem ####
  logger::log_info("Assemble repeat code TALE strings and write in a file for ARLEM")
  
  codesSeqsfile <- tempfile(fileext = ".fasta")
  repeatStrings <- taleParts %>%
    dplyr::group_by(arrayID) %>%
    dplyr::arrange(positionInArray) %>%
    dplyr::summarise(repeatString = paste(domCode, collapse = " "),
                     posString = paste(positionInArray, collapse = " "))
  codesSeqSet <- Biostrings::BStringSet(repeatStrings$repeatString)
  names(codesSeqSet) <- repeatStrings$arrayID
  
  # # Must use seqinr because Biostrings wraps sequences in fasta file which messes up Arlem...
  # codesSeqLst <- as.list(repeatStrings$repeatString)
  # names(codesSeqLst) <- repeatStrings$arrayID
  # seqinr::write.fasta(codesSeqLst, names = names(codesSeqLst),
  #                     file.out = codesSeqsfile, as.string = TRUE, nbchar = 10000)
  
  Biostrings::writeXStringSet(x = codesSeqSet,
                              filepath = codesSeqsfile,
                              format = "fasta", width = 20000L)
  
  
  
  #### Compute systematic pairwise dissimilarities (distances) between 'repeat' units. ####
  # Get unique domains sequences
  taleAaParts <- Biostrings::AAStringSet(taleParts$aaSeq)
  names(taleAaParts) <- taleParts$domCode
  uniqueTaleAaParts <- unique(taleAaParts)
  stopifnot(!anyDuplicated(names(uniqueTaleAaParts)))
  stopifnot(!anyDuplicated(names(unique(taleAaParts))))
  
  # Get pairwise repeat aa sequence dissimilarity scores in a long tibble
  logger::log_info("Computing a distance matrix between TALE parts amino acid sequences ",
                   "using: {pairwiseAlnMethod}")
  if (pairwiseAlnMethod == "mmseq2") {
    dissimLong <- .distalPairwiseAlign2(partAaStringSet = uniqueTaleAaParts, ncores = ncores,
                                        condaBinPath = condaBinPath)
    #saveRDS(dissimLong, file = "/home/cunnac/TEMP/dissimLong")
  } else if (pairwiseAlnMethod == "Biostrings") {
    dissimLong <- .distalPairwiseAlign(partAaStringSet = uniqueTaleAaParts, ncores = ncores)
  } else if (pairwiseAlnMethod == "DECIPHER") {
    dissimLong <- .distalPairwiseAlign3(partAaStringSet = uniqueTaleAaParts, ncores = ncores)
  } else {
    logger::log_errors() && stop("'pairwiseAlnMethod' parameter must be either 'Biostrings', 'mmseq2' or 'DECIPHER'")
  }
  dissimLong %<>% dplyr::mutate(Sim = 100 - Dissim)
  # Convert Distance (dissimilarity) measures to Similarity with a four-parameter logistic function
  # pairAlignScores %<>% dplyr::mutate(Sim = 100/(1+exp(-1*-0.9*(Dissim-3))))
  
  
  # Convert to square matrix
  dissimMat <- reshape2::acast(dissimLong, formula = subj ~ pattern, value.var = "Dissim")
  stopifnot(nrow(dissimMat) == ncol(dissimMat))
  # reorder row and colnames because I suspect arlem expect them in increasing order
  dissimMat <- dissimMat[rownames(dissimMat) %>% as.numeric() %>% order(),
                         colnames(dissimMat) %>% as.numeric() %>% order()]
  
  #### Generate an ARLEM cost matrix ####
  if (TRUE) {
    method <- "minkowski"
    logger::log_info("Generate an ARLEM cost matrix which meets triangle inequality criteria by computing ",
                     "the {method} distance between pairwise distance vectors.")
    dissimMat <- as.matrix(stats::dist(dissimMat, method = method, p = 3.5, diag = TRUE, upper = TRUE))
    dissimMat <- dissimMat/max(dissimMat) * 100
  }
  # if (!fossil::tri.ineq(dissimMat)) {
  #   logger::log_error("TALE domains dissimilarity (distance) matrix does not respect the triangle inequality",
  #                     "Arlem will fail. Aborting...")
  #   stop()
  # }
  
  #### Prepare Arlem cfile with systematic pairwise distances between 'repeat' units. ####
  # Get parameters for arlem
  TypeNo <- glue::glue("# Type no. ", nrow(dissimMat))
  Types <- glue::glue("# Types ", paste(1:nrow(dissimMat), collapse = " "))
  # Convert mat to character, beware of the ceiling in conversion...
  dissimMat <- matrix(ceiling(dissimMat) %>% format(),
                      ncol = ncol(dissimMat),
                      dimnames = list(rownames(dissimMat), colnames(dissimMat))
                      )
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
  logger::log_info("Running ARLEM version 1.0 : ")
  logger::log_info("Copyright by Mohamed I. Abouelhoda")
  logger::log_info("Plz. cite Abouelhoda, Giegerich, Behzadi, and Steyaert")
  arlemRawRes <- system(arlemCmd, intern = TRUE)
  logger::log_debug(logger::skip_formatter(arlemRawRes))
  arlemSelfRes <- grep("Processed Seq[.]:", arlemRawRes, value = TRUE) 
  arlemSelfScores <- gsub("Processed Seq[.]: ([0-9]{1,}) Score: ([0-9]{1,}),.*", "\\1|\\2",
                          substring(arlemSelfRes, 1, 35)) %>%
    strsplit(split = "\\|") %>%
    lapply(function(s) t(as.matrix(as.numeric(s)))) %>%
    do.call(rbind, .) %>% tibble::as_tibble(.name_repair = "minimal")
  arlemRes <- grep("Score of aligning Seq:",
                   arlemRawRes, value = TRUE)
  arlemScores <- gsub("Score of aligning Seq:([0-9]+), Seq:([0-9]+) =([0-9]+)", "\\1|\\2|\\3",
                      arlemRes)
  arlemScores <- strsplit(arlemScores, split = "\\|")
  arlemScores <- lapply(arlemScores, function(s) {t(as.matrix(as.numeric(s)))}) %>%
    do.call(rbind, .) %>% tibble::as_tibble(.name_repair = "minimal")
  colnames(arlemScores) <- c("TAL1", "TAL2", "arlemScore")
  # Shaping into matrix to have scores in both directions
  arlemScoresMat <- reshape2::acast(arlemScores, formula = TAL1 ~ TAL2, value.var = "arlemScore")
  arlemScoresMat <- cbind("0" = NA, arlemScoresMat)
  arlemScoresMat <- rbind(arlemScoresMat, "43" = NA)
  # dim(arlemScoresMat)
  # dimnames(arlemScoresMat)
  arlemScores <- stats::as.dist(t(arlemScoresMat), diag = TRUE, upper = TRUE) %>% as.matrix() %>%
    reshape2::melt(value.name = "arlemScore") %>%
    tibble::as_tibble()
  colnames(arlemScores) <- c("TAL1", "TAL2", "arlemScore")
  
  #### Compute normalized arlem scores and include arrayIDs rather than arlem index ####
  arrayLengths <- taleParts %>% dplyr::group_by(arrayID) %>% dplyr::count()
  
  normArlemScoresTble <- arlemScores %>%
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
  
  #### Check features of the Arlem results table
  arraysCount <- codesSeqSet %>% length()
  if (nrow(normArlemScoresTble) != arraysCount^2) {
    allCombs <- expand.grid(names(codesSeqSet), names(codesSeqSet), stringsAsFactors = FALSE) %>% tibble::as_tibble()
    colnames(allCombs) <- c("TAL1", "TAL2")
    absentCombs <- dplyr::left_join(allCombs, normArlemScoresTble) %>%
      dplyr::filter(is.na(Sim))
    cat(knitr::kable(absentCombs), sep = "\n")
    logger::log_error("The TALE similarity table does not have the expected number",
    "of comparisons (number of missing comps: {nrow(absentCombs)})...")
    stop()
  }
  
  
  #### return a list of results emulating the return value of tantale::runDistal ####
  outputlist <- list(taleParts = taleParts,
                     "repeats.code" = taleParts %>%
                       dplyr::group_by(domCode, aaSeq, rvd) %>%
                       dplyr::count() %>%
                       dplyr::rename(code = domCode, "AA Seq"  = aaSeq) %>%
                       dplyr::mutate(code = as.integer(code)) %>%
                       dplyr::select(-n) %>%
                       dplyr::ungroup(),
                     "coded.repeats.str" = toListOfSplitedStr(codesSeqsfile),
                     "repeat.similarity" = dissimLong %>% dplyr::rename(RepU1 = subj, RepU2 = pattern),
                     "tal.similarity" = normArlemScoresTble,
                     "repeats.cluster" = clusterRep(
                       repeatSimMat = reshape2::acast(dissimLong, formula = subj ~ pattern, value.var = "Sim"),
                       repeats.cluster.h.cut = repeats.cluster.h.cut
                     )
  )
  logger::log_info("Finished! Returning a list with the results.")
  return(outputlist)
}




