# roxygen2::roxygenise() to write documentation

#### Supporting functions ####


#'
#' Split strings of TALE sequences (`sep`-separated rvd or distal repeat IDs)
#'
#' @description Load the content of fasta file containing TALE sequences (either RVD or Distal repeat code)
#' and return a list of vectors each one composed of the individual elements of the sequence.
#'
#' @param atomicStrings Either, the path to a fasta file, an AAStringSet or "BStringSet"
#' or or a list. In all cases, the content of these objects is composed of strings of
#' tale sequences (`sep`-separated rvd or distal repeat IDs)
#' @param sep Separator of the elements of the sequence
#'
#' @return A list of named vectors representing the 'splited' sequence.
#'
#' @export
toListOfSplitedStr <- function(atomicStrings, sep = "-") {
  if (length(atomicStrings) > 1 && is.list(atomicStrings)) {
    seqs <- atomicStrings
  } else if (length(atomicStrings) == 1 && is.character(atomicStrings)) {
    stopifnot(fs::file_exists(atomicStrings))
    seqs <- as.character(Biostrings::readBStringSet(atomicStrings), use.names=TRUE)
  } else if (class(atomicStrings) %in% c("AAStringSet", "BStringSet")) {
    seqs <- as.character(atomicStrings, use.names=TRUE)
  } else {
    logger::log_error("Something is wrong with the value provided for atomicStrings.")
    stop()
  }

  seqsAsVectors <- stringr::str_split(seqs, pattern = glue("[{sep}]"))
  seqsAsVectors <- lapply(seqsAsVectors, function(x) { # Remove last residue if it is empty string
    if ( x[length(x)] == "") {
      warning("## Last element in 'vectorized' sequence is empty. It was removed from output.")
      x[-length(x)]
      } else x
  }
  )
  names(seqsAsVectors) <- names(seqs)
  return(seqsAsVectors)
}





formatDistalRepeatDistMat <- function(distalRepeatDistMatFile) {
  # Distal-1.2 repeat distance matrix is 'almost' symetrical but does not contains the diagonal
  # Top triangle of a symetrical :   Symetrical :
  # 1234                              1234
  #  234                              2234
  #   34                              3334
  #    4                              4444
  # It is like this:
  # 234
  # 34
  # 4
  # So, we need to do a few transformations:
  # Loading file content as a matrix
  distalRepeatDist <- as.matrix(
    read.table(distalRepeatDistMatFile,
               sep = " ",
               fill = TRUE,
               blank.lines.skip = TRUE,
               header = FALSE,
               check.names = FALSE,
               comment.char = "#"
    )
  )
  # Getting ride of the last columns that appears because the lines in the file have a final space
  distalRepeatDist <- distalRepeatDist[, -ncol(distalRepeatDist)]
  # Adding a first column
  distalRepeatDist <- cbind(V0 = NA, distalRepeatDist)
  # Adding a last line
  distalRepeatDist <- rbind(distalRepeatDist, NA)

  # Shifting values to the right in rows with NAs
  distalRepeatDist <- t(apply(distalRepeatDist, 1, function(x) {
    row <- x
    c(row[is.na(row)], row[!is.na(row)])
  }))
  # Filling diagonal with 0 values
  diag(distalRepeatDist) <- 0
  # Filling NA values with diagonal symetric values
  newmat <- distalRepeatDist
  for (i in 1:nrow(distalRepeatDist)) {
    for (j in 1:ncol(distalRepeatDist)) {
      if (is.na(distalRepeatDist[i, j])) {
        distalRepeatDist[i, j] <- distalRepeatDist[j, i]
      } else {
        next()
      }
    }
  }
  # Names
  stopifnot(exprs= all.equal(nrow(distalRepeatDist), ncol(distalRepeatDist)))
  rownames(distalRepeatDist) <- 0:(nrow(distalRepeatDist) - 1)
  colnames(distalRepeatDist) <- 0:(ncol(distalRepeatDist) - 1)
  # Similarity rather than dissimilarity (that is what Alvaro does in his scripts)
  distalRepeatSim <- 100 - distalRepeatDist
  # str(distalRepeatSim)
  # image(t(apply(distalRepeatSim, 2, rev)))

  # Output in 'long format'
  distalRepeatSimTable <- reshape2::melt(distalRepeatSim,
                                         as.is = TRUE,
                                         varnames = c("RepU1", "RepU2"),
                                         value.name = "Sim")
  # str(distalRepeatSimTable)
  return(distalRepeatSimTable)
}





#' Generate a mapping between Distal repeat IDs and their cognate RVD.
#'
#' Uses Distal repeat sequences and RVD sequences from a set of TALEs to return
#' the association between repeat ID and RVD.
#'
#' Care must be taken that TALEs in the two sets of sequences have the same name.
#' In addition, the function tries hard to make sure that the two sets of sequences are identical in every ways but the individual 'values' they contain.
#' It is therefore notably important to make sure that the sequences are consistent in whether they include N-term and C-term domains IDs/Tags or not.
#'
#' @param talesRepeatVectors Expects a list of Distal repeat IDs character
#'   vectors. Each \strong{named} element corresponding to a TALE.
#' @param talesRepeatVectors Expects a list of Distal RVDs character
#'   vectors. Each \strong{named} element corresponding to a TALE.
#' @return A two columns repeatID - RVD data frame.
#' @export
getRepeat2RvdMapping <- function(talesRepeatVectors, talesRvdVectors) {
  # Making sure, these objects are indentical in every ways but the actual values of the vectors
  stopifnot(setequal(names(talesRepeatVectors), names(talesRvdVectors)))
  lrep <- sapply(talesRepeatVectors, length)
  lrep <- lrep[order(names(lrep))]
  lrvd <- sapply(talesRvdVectors, length)
  lrvd <- lrvd[order(names(lrvd))]
  stopifnot(names(lrvd) == names(lrep))
  stopifnot(apply(cbind(lrep, lrvd), 1, function(x) x[1] == x[2]))

  # Function to melt the lists
  .l2df <- function(l) {
    dplyr::bind_rows(lapply(l, function(v) data.frame(idx = 1:length(v),
                                                      repeats = v,
                                                      stringsAsFactors = FALSE)
    ),
    .id = "Name")
  }
  talesRepeatDf <- .l2df(talesRepeatVectors)
  talesRvdDf <- .l2df(talesRvdVectors)
  stopifnot(nrow(talesRepeatDf) == nrow(talesRvdDf))
  # Merging to have the repeat ID vs RVDs
  repeat2rvd <- dplyr::full_join(talesRepeatDf, talesRvdDf, by = c("idx" = "idx", "Name" = "Name"))
  colnames(repeat2rvd) <- c("names", "idx", "repeatID", "RVD")
  # Check that for each repeat there is only one corresponding RVD (the converse is NOT true)
  repeat2rvd <- repeat2rvd %>% dplyr::group_by(repeatID, RVD) %>%
    dplyr::summarise(count = dplyr::n())
  check <- repeat2rvd %>%
    dplyr::arrange(repeatID) %>%
    dplyr::group_by(repeatID) %>%
    dplyr::summarise(count = dplyr::n())
  stopifnot(sum(check$count) == length(unique(repeat2rvd$repeatID)))
  # Simplifiy the df
  repeat2rvd <- repeat2rvd %>% dplyr::select(repeatID, RVD) %>% dplyr::ungroup()

  return(repeat2rvd)
}

#' read alignment results of Alvaro's Perl scripts for building TALE groups from DisTal output
#' @description read multiple files (.RVDs/.Dists/.Reps) and output a list whose each element is a dataframe read from each file
read_distal_aligns <- function(file.lists) {
  output <- list()
  for (x in file.lists) {
    tab <- read.delim(x, header = F, row.names = 1, allowEscapes = T, stringsAsFactors = F, na.strings = c("-", "NA"))
    tab <- tab[order(rownames(tab)), ]
    tab <- as.matrix(tab)
    Gname <- gsub("\\..*", "", basename(x))
    rownames(tab) <- gsub(glue::glue("\\|{Gname}"), "", rownames(tab))
    colnames(tab) <- sapply(1:ncol(tab), function(y) paste0("R", y))
    output[[Gname]] <- tab
  }
  return(output)
}


clusterRep <- function(repeatSimMat, repeats.cluster.h.cut) {
  dist_clust <- hclust(as.dist(repeatSimMat))
  dist_cut <- as.data.frame(
    cbind(RepID = dist_clust$labels,
          Rep_clust = cutree(dist_clust, h = repeats.cluster.h.cut)
    )
  )
  dist_cut$Rep_order <- order.dendrogram(as.dendrogram(dist_clust))
  dist_cut <- dist_cut[order(dist_cut$Rep_order),] %>%
    dplyr::as_tibble()
  return(dist_cut)
}









#### MAIN functions ####

#' Run the Distal tool
#' @description Run the Distal tool of the \href{https://doi.org/10.3389/fpls.2015.00545}{QueTal} suite to classify and compare TAL effectors functionally and phylogenetically.
#'
#' @param fasta.file fasta file containing TALE DNA/AA sequences.
#' @param outdir directory to store disTal output.
#' @param repeats.cluster.h.cut numeric value to cut the hierachycal clustering tree of the repeat.
#' @param overwrite logical indicating whether to rerun disTal or only load the existing results.
#' @return A list with DisTal output components: 
#' \itemize{
#'   \item repeats.code: a data frame of the unique repeat AA sequences and there numeric codes
#'   \item coded.repeats.str: a list of repeat-coded TALE strings
#'   \item repeat.similarity:  a long, three columns data frame with pairwise similarity scores between repeats
#'   \item tal.similarity: a three columns Tals similarity table with pairwise similarity scores between TALEs
#'   \item tree: Newick format neighbor-joining tree of TALs constructed based on TALEs similarity
#'   \item repeats.cluster: a data frame containing repeat code and repeat clusters.
#' }
#' @export
runDistal <- function(fasta.file, outdir = NULL, treetype = "p", repeats.cluster.h.cut = 10, overwrite = F) {
  if (is.null(outdir)) {
    outdir <- tempdir(check = TRUE)
  }
  # whether it's necessary to run distal perl script or not
  check_files <- list.files(outdir, "Output(.)*(.mat|.tre|.pdf|.fa|.txt)")
  if (overwrite || length(check_files) < 6) {
    logger::log_info("Running the DisTAL perl program. Be patient this may take a while to complete...",
                     "Pre-existing DisTAL ouput files in {outdir} will be overwritten.")
    fasta.file <- fasta.file
    outdir <- outdir
    sourcecode <- system.file("tools", "DisTAL1.2_MultipleAlignment", package = "tantale", mustWork = T)
    disTal <- file.path(sourcecode, "DisTAL_v1.2_matest_M.pl")
    lib <- file.path(sourcecode, "lib")
    disTalCMD <- paste("perl -I", lib, disTal, "-m T", "-n", treetype, "-o", outdir, shQuote(fasta.file), sep = " ")
    system(disTalCMD, ignore.stdout = T)
  } else {
    logger::log_info("The specified {outdir} already contains all DisTAL output files. ",
                     "The returned results object will be build from their content.")
  }
  
  
  # read repeatscode.txt
  repeatscode_File <- glue::glue("{outdir}/Output_Repeatscode.txt")
  repeatscode <- read.csv(repeatscode_File, sep = "\t", header = F)
  colnames(repeatscode) <- c("code", "AA Seq")
  
  # read codedrepeats.fa
  codedRepeats_File <- glue::glue("{outdir}/Output_CodedRepeats.fa")
  codedRepeats_str <- toListOfSplitedStr(codedRepeats_File)
  
  # read repeatmatrix.mat
  repeatmatrix_File <- glue::glue("{outdir}/Output_Repeatmatrix.mat")
  repeatSim <- formatDistalRepeatDistMat(repeatmatrix_File)
  
  
  
  ## define 'col' based on clustering groups of repeats
  unmelt_repeatSim <-  as.matrix(reshape2::acast(repeatSim, RepU1 ~ RepU2, value.var="Sim"))
  dist_cut <-  clusterRep(repeatSimMat = unmelt_repeatSim,
                          repeats.cluster.h.cut = repeats.cluster.h.cut)
  
  # read tal seqs distance matrix and make sim df
  talMatrix_File <- glue::glue("{outdir}/Output.mat")
  raw_talMatrix <- read.table(talMatrix_File,
                              header = T,
                              sep = "\t",
                              row.names = NULL,
                              check.names = FALSE)
  rownames(raw_talMatrix) <- raw_talMatrix[, 1] # rownames without 'matrix'
  # Removing whitespaces
  rownames(raw_talMatrix) <- stringr::str_trim(rownames(raw_talMatrix))
  colnames(raw_talMatrix) <- stringr::str_trim(colnames(raw_talMatrix))
  # Removing first column with row names and last column with NAs
  raw_talMatrix <- raw_talMatrix[, c(-1, -ncol(raw_talMatrix))]
  raw_talMatrix <- 100 - as.matrix(raw_talMatrix) # Distance to similarity conversion
  # Checking if the matrix is square
  stopifnot(all.equal(ncol(raw_talMatrix), nrow(raw_talMatrix)))
  if (!setequal(colnames(raw_talMatrix), rownames(raw_talMatrix))) {
    logger::log_warn("Sequences of Row and Col names in TALEs distance matrix do not match...")
    logger::log_warn("Set dif Col vs Row names : {setdiff(colnames(raw_talMatrix), rownames(raw_talMatrix))}")
    logger::log_warn("Set dif Row vs Col names : {setdiff(rownames(raw_talMatrix), colnames(raw_talMatrix))}")
    
    }
  talsim <- reshape2::melt(raw_talMatrix)
  colnames(talsim) <-c("TAL1", "TAL2", "Sim")
  
  # read newick tree
  nw_tree_File <-  glue::glue("{outdir}/Output.tre")
  nw_tree <- ape::read.tree(nw_tree_File)
  
  # Assemble return object
  outputlist <- list("repeats.code" = repeatscode %>% tibble::as_tibble(),
                     "coded.repeats.str" = codedRepeats_str,
                     "repeat.similarity" = repeatSim %>% tibble::as_tibble(),
                     "tal.similarity" = talsim %>% tibble::as_tibble(),
                     "tree" = nw_tree,
                     "repeats.cluster" = dist_cut %>% tibble::as_tibble())
  return(outputlist)
}

#' Run Alvaro's Perl scripts for TALEs grouping
#' @description take distal output and classify groups
#' @param path directory containing DisTal output files, or the same object as the output of \code{outdir} for \code{\link[tantale:runDistal]{runDistal}}
#' @param num.groups an integer indicating the number of TALEs groups you want to classify.
#' @param overwrite logical indicating whether to rerun the Perl scripts or only load the existing results.
#' @return A list containing: 
#' \itemize{
#'   \item SeqOfRvdAlignments: a list of matrices of TALES alignment with RVD sequences
#'   \item SeqOfDistancesAlignments: a list of matrices of TALEs alignment with repeat distance
#'   \item repeatUnitsDistanceMatrix:  a matrix of pairwise similarity scores between repeats
#'   \item SeqOfRepsAlignments: a list of matrices of TALEs alignment with repeat codes
#'   \item TALgroups: a data frame of TALEs names and their groups
#' }
#' @export
buildDisTalGroups <- function(path, num.groups, overwrite = F) {
  check_files <- list(
    Aligned_TEV_ALL = file.path(path, "Aligned_TEV_ALL"),
    Coded_Reps_withgroups = file.path(path, glue::glue("Coded_Reps_withgroups_{num.groups}.fa")),
    BigRepDist = file.path(path, "BigRepDist.mat"),
    Multiple_align = file.path(path, "Multiple_align"),
    aligns = file.path(path, "ALIGNS"),
    talgroups = file.path(path, "TALgroups.txt")
  )
  print(check_files)
  if (!all(unlist(lapply(check_files, file.exists))) || overwrite) {
    # get perl scripts and perl lib
    sourcecode <- system.file("tools", "DisTAL1.2_MultipleAlignment", package = "tantale", mustWork = T)
    Analyze_TALs1 <- file.path(sourcecode, "Analyze_TALs1_M.pl")
    Analyze_TALs2 <- file.path(sourcecode, "Analyze_TALs2_M.pl")
    Generate_Bigmat <- file.path(sourcecode, "Generate_Bigmat.pl")
    Alignement_files <- file.path(sourcecode, "Alignement_files_M.pl")
    lib <- file.path(sourcecode, "lib")

    # wrap 'Analyze_TALs1_M.pl'
    mat <- list.files(path, ".mat", full.names = T)
    tal_mat <- mat[!grepl("_Repeatmatrix|BigRepDist", mat)]
    codedRepeats <-  list.files(path, "_CodedRepeats.fa", full.names =T)
    analyze1CMD <- paste("perl -I", lib, Analyze_TALs1, tal_mat, num.groups, codedRepeats, "T", sep = " ")

    # wrap 'Analyze_TALs2_M.pl'
    # Aligned_TEV_ALL <- file.path(path, "Aligned_TEV_ALL")
    # Coded_Reps_withgroups <- list.files(path, "Coded_Reps_withgroups_\\d+", full.names = T)
    analyze2CMD <- paste("perl -I", lib, Analyze_TALs2, check_files$Aligned_TEV_ALL, check_files$Coded_Reps_withgroups, sep = " ")

    # wrap 'Generate_Bigmat.pl'
    repeat_mat <- mat[grepl("_Repeatmatrix", mat)]
    # BigRepDist <- file.path(path, "BigRepDist.mat")
    bigmatCMD <- paste("perl -I", lib, Generate_Bigmat, repeat_mat, check_files$BigRepDist, sep = " ")

    # wrap 'Alignement_files_M.pl'
    Repeatscode <- list.files(path, "_Repeatscode.txt", full.names = T)
    # Multiple_align <- file.path(path, "Multiple_align")
    alignmentCMD <- paste("perl -I", lib, Alignement_files, Repeatscode, check_files$Multiple_align, sep = " ")

    # remove old files
    unlink(check_files, recursive = T)

    # run everything
    print(c(analyze1CMD, analyze2CMD, bigmatCMD, alignmentCMD))
    system(paste(analyze1CMD, analyze2CMD, bigmatCMD, alignmentCMD, sep = " && "), ignore.stdout = T)

    # check if tals classification is done
    if (!file.exists(check_files$talgroups)) {
      rlang::abort(message = glue::glue("Cannot classify {num.groups} groups!"))
    }
  }

  # read outputs
  RVDalignFiles <- list.files(check_files$aligns, ".RVDs", full.names = T)
  RVDgroups <- read_distal_aligns(RVDalignFiles)

  sim_within_group_File <- list.files(check_files$aligns, ".Dists", full.names = T)
  sim_within_group <- read_distal_aligns(sim_within_group_File)
  sim_within_group <- lapply(sim_within_group, function(x) x <- 100 - x)

  # dist_between_rep_File <- file.path(path, "BigRepDist.mat")
  dist_between_rep <- as.matrix(read.table(check_files$BigRepDist, header=TRUE, sep = "\t",row.names = 1, as.is=TRUE))

  RepsalignFile <- list.files(check_files$aligns, ".Reps", full.names = T)
  RepCodegroups <- read_distal_aligns(RepsalignFile)

  # talgroupsFile <- list.files(path, "TALgroups.txt", full.names = T)
  talgroups <- read.table(check_files$talgroups, header = T)
  # talgroups <- readr::read_delim(talgroupsFile," ", escape_double = FALSE, trim_ws = TRUE)

  outputlist <- list("SeqOfRvdAlignments" = RVDgroups, "SeqOfDistancesAlignments" = sim_within_group, "repeatUnitsDistanceMatrix" = dist_between_rep, "SeqOfRepsAlignments" = RepCodegroups, "TALgroups" = talgroups)
  return(outputlist)
  }









#' Heatmap plotting of rvd sequence variants
#' @description The function creates a graphical presentation from a tale annotation table. The output is like a heatmap that presents rvd sequence variants in Tal groups as column and respective strains as rows (or vice versa). It is different from a typical heatmap that it can display more than one value in a cell; for example, if one strain has 2 rvd sequence variants belong to 1 group, it will be displayed by 2 colors in 1 cell.
#' @param tale_annotation a data frame containing at least 3 columns for Tal groups, strain names, and rvd seqs, and 1 row is 1 Tal.
#' @param col "character", column name of \code{tale_annotation} to be displayed as rows in the heatmap (e.g. tal groups).
#' @param row "character", column name of \code{tale_annotation} to be displayed as columns in the heatmap (e.g. strain names).
#' @param value "character", column name for rvdseqs in the \code{tale_annotation}
#' @param truncTaleLab (optional, default = NULL) "character", column name of \code{tale_annotation} labeling the truncTales by TRUE/FALSE value. The truncTales are labeled by "T" in the heatmap cells, but if this argument is called.
#' @param extraCol (optional, default = NULL) "character", column name of \code{tale_annotation} containing other information (e.g. origin). It will be presented in a side bar on the right of the heatmap.
#' @param x.lab,y.lab,title character for x axix, y axis names and title 
#' @param mapcol character vector of colors for the cells.
#' @param mar.side margin of the heatmap for row dendrogram, col dendrogram, rownames, colnames, respectively. (by default, c(5, 5, 3, 3)).
#' @param sepwid numeric value for the width of separator between adjacent cells.
#' @param sepcol character of color for the separator between adjacent cells
#' @param inner_sepcol character of color for the separator between colors within 1 cell if there are more than 1.
#' @param save.path (optional) file path to save the plot, format of the image depends on the file extension. If save.path is NULL, the heatmap will be printed. If save.path is specified, the image file will be created.
#' @export
heatmap_talomes <- function(tale_annotation, col, row, value, truncTaleLab = NULL, extraCol = NULL,
                            x.lab = "TALE Group", y.lab = "Strain", title = "RVD sequences variants",
                            plot.type = "all",
                            mapcol = viridis::viridis(10), mar.side = c(5, 5, 3, 3),
                            sepwid = 5, sepcol = "white", inner_sepcol = "white", save.path = NULL) {

  ## rename colnames of tale annotation
  colnames(tale_annotation)[which(colnames(tale_annotation) == col)] <- "group"
  colnames(tale_annotation)[which(colnames(tale_annotation) == row)] <- "strain"
  colnames(tale_annotation)[which(colnames(tale_annotation) == value)] <- "rvdseq"
  if (!is.null(extraCol)) colnames(tale_annotation)[which(colnames(tale_annotation) == extraCol)] <- "extraCol"
  if (!is.null(truncTaleLab)) colnames(tale_annotation)[which(colnames(tale_annotation) == truncTaleLab)] <- "truncTale"

  tale_annotation %<>% dplyr::mutate(aberrantRepeat = ifelse(grepl("[a-z]", rvdseq), TRUE, FALSE))
  tale_annotation$rvdseq <- toupper(tale_annotation$rvdseq)

  if(is.numeric(tale_annotation$group)) tale_annotation$group <- paste0("G", tale_annotation$group)

  ## convert rvdseq to rvdfac, ranking of abundance
  tale_annotation %<>% dplyr::group_by(group) %>% dplyr::mutate(rvdfac = do.call((function(x) {
    levels(x) <- nlevels(x) + 1 - rank(table(x), ties.method = "first")
    return(as.integer(as.character(x)))
  }), list(as.factor(rvdseq))))

  ## get number of rvdseqs per group per strain
  ## it is used for the layout of heatmap cell
  numAlleles <- reshape2::dcast(tale_annotation, strain ~ group, value.var = "rvdseq", drop = F, fun.aggregate = function(x) length(x))

  rownames(numAlleles) <- numAlleles$strain
  numAlleles <- numAlleles[,-1]

  ## count number of alleles
  ## for column label
  variantsCount <- sapply(sort(unique(tale_annotation$group)), function(g) {
    length(unique(tale_annotation[tale_annotation$group == g,]$rvdseq))
  }, USE.NAMES = T)
  colnames(numAlleles) <- paste0(colnames(numAlleles), " #", variantsCount)

  tale_annotation %<>% dplyr::group_by(group, strain) %>% dplyr::mutate(reprsntRVDfac = ifelse(1 %in% rvdfac, 1, rvdfac[1]))

  ## representative alleles per group per strain
  ## in case a strain has more than 1 rvd seqs in 1 group
  ## find the most abundant rvdseqin that group
  ## if this strain has that rvdseq, take it as representative rvdseq
  ## if not, take randomly 1 rvdseq that strain has
  reprsntAlleles <- reshape2::dcast(tale_annotation, strain ~ group, value.var = "reprsntRVDfac", drop = F, fun.aggregate = function(x) as.integer(x[1]))
  reprsntAlleles <- as.data.frame(lapply(X = reprsntAlleles[,-1], FUN = as.integer))
  rownames(reprsntAlleles) <- rownames(numAlleles)
  colnames(reprsntAlleles) <- colnames(numAlleles)
  reprsntAlleles <- apply(reprsntAlleles, 2, function(x) { # allele '0' if missing
    x[is.na(x)] <- 0
    return(x)})

  ## compute dendrograms with the representative alleles
  Strain_HC <- hclust(ape::dist.gene(x = reprsntAlleles), method = "average")

  strainAlleles <- apply(reprsntAlleles, 1, as.integer)
  rownames(strainAlleles) <- colnames(numAlleles)
  TALE_HC <- hclust(ape::dist.gene(x = strainAlleles), method = "average")


  df <- par(no.readonly = T)
  ## plot sizes
  widleft <- mar.side[1]
  heitop <- mar.side[2]
  widright <- mar.side[3]
  heibot <- mar.side[4]

  if (plot.type == "single") { # plot representative alleles
    codedAlleles1 <- apply(reprsntAlleles, 2, function(x) ifelse(x == 0, NA, x))
    if (hasArg(save.path)) {
      img_format <- gsub(".*\\.", "", basename(save.path))
      img_size <-  list(save.path, width = (widleft + ncol(codedAlleles1) + widright)/2.54, height = (heitop + nrow(codedAlleles1) + heibot)/2.54)
      if (img_format %in% c("bmp", "jpeg", "png", "tiff")) {
        img_size <- c(img_size, units = "in", res = 1440)
      }
      do.call(img_format, img_size)
    }
    gplots::heatmap.2(as.matrix(codedAlleles1),
                      trace = "none",
                      col = mapcol[1:max(codedAlleles1, na.rm = T)],
                      breaks = 0:max(codedAlleles1, na.rm = T),
                      density.info = "none",
                      key = F,
                      xlab = x.lab,
                      ylab = y.lab,
                      margins = c(7, 7),
                      colsep = 0:(ncol(codedAlleles1)-0),
                      rowsep = 0:(nrow(codedAlleles1)-0),
                      sepcolor = sepcol,
                      sepwidth = rep(sepwid/100, 2),
                      Rowv = as.dendrogram(Strain_HC),
                      Colv = as.dendrogram(TALE_HC),
                      main = title,
                      na.color = "grey",
                      lmat = rbind(c(4, 3), c(2, 1)),
                      lhei = c(heitop, nrow(codedAlleles1) + heibot), ##
                      lwid = c(widleft, ncol(codedAlleles1) + widright)
    )

  } else if (plot.type == "all") { # plot all alleles
    uniqueRVD <- numAlleles
    rorder <- order.dendrogram(as.dendrogram(Strain_HC))
    corder <- order.dendrogram(as.dendrogram(TALE_HC))
    uniqueRVD <- uniqueRVD[rev(rorder), corder]


    colmat <- mapcol
    ## plot layout
    nplots <- nrow(uniqueRVD)*ncol(uniqueRVD)
    mainmat <- matrix(1:nplots, nrow = nrow(uniqueRVD))
    left.mat <- matrix(rep(nplots+1, nrow(uniqueRVD)), ncol = 1)
    top.mat <- matrix(c(0, rep(nplots+2, ncol(uniqueRVD))), nrow = 1, byrow = T)
    # right.mat <- matrix(rep(c(rep(0, heitop), (nplots+3):(nrow(uniqueRVD)+nplots+2)),2), ncol = 2, byrow = F)
    # bottom.mat <- matrix(c(rep(0, widleft), (nrow(uniqueRVD)+nplots+3):(nrow(uniqueRVD)+nplots+2+ncol(uniqueRVD)), 0, 0), nrow = 1)
    extcol <- matrix(c(0, rep(nplots+3, nrow(uniqueRVD))), ncol = 1, byrow = F)
    right.mat <- matrix(c(0, rep(nplots+4, nrow(uniqueRVD))), ncol = 1, byrow = F)
    legend.mat <- matrix(c(0, rep(nplots+5, nrow(uniqueRVD))), ncol = 1, byrow = F)
    bottom.mat <- matrix(c(0, rep(nplots+6, ncol(uniqueRVD)), 0, 0, 0), nrow = 1, byrow = T)
    laymat <- rbind(top.mat, cbind(left.mat, mainmat))
    laymat <- rbind(cbind(laymat, extcol, right.mat, legend.mat), bottom.mat)
    titmat <- matrix(c(0, rep(nplots+7, ncol(mainmat)), 0, 0, 0), nrow = 1, byrow = T)
    laymat <- rbind(titmat, laymat)
    on.exit(par(no.readonly = TRUE))

    if (hasArg(save.path)) {
      img_format <- gsub(".*\\.", "", basename(save.path))
      img_size <-  list(save.path, width = sum(widleft, rep(1, ncol(uniqueRVD)), widright)/2.54, height = sum(2, heitop, rep(1, nrow(uniqueRVD)), heibot)/2.54)
      if (img_format %in% c("bmp", "jpeg", "png", "tiff")) {
        img_size <- c(img_size, units = "in", res = 1440)
      }
      do.call(img_format, img_size)
    }

    layout(laymat, widths = c(widleft, rep(1, ncol(uniqueRVD)), ifelse(is.null(extraCol), 0.1, .7), widright, ifelse(is.null(extraCol), 0.1, 3)), heights = c(2, heitop, rep(1, nrow(uniqueRVD)), heibot))
    # layout.show(nplots+7)


    ## plot heatmap
    for (c in 1:ncol(uniqueRVD)) {
      gname <- gsub(" \\#\\d+", "", colnames(uniqueRVD)[c])
      # gname <- gsub("G", "", gname)
      g1 <- tale_annotation[tale_annotation$group == gname,]
      # g1$rvdseq <- as.integer(as.factor(g1$rvdseq))
      for (r in 1:nrow(uniqueRVD)) {
        par(mar = rep(0, 4))
        nelements <- uniqueRVD[r, c]
        if (nelements == 0) {
          image(z = matrix(0), col = "grey", axes = F)
          if (sepwid > 0) box(lwd = sepwid/2, col = sepcol)
        } else {
          g1s1 <- g1[g1$strain == rownames(uniqueRVD[r,]),]
          rvd.factor <- g1s1$rvdfac
          if (is.null(truncTaleLab)) {
            truncTale <- NA
          } else {
            truncTale <- sapply(g1s1$truncTale, function(l) ifelse(l, "T", NA))
          }
          # plot.index <- r + (c-1)*nrow(uniqueRVD)
          color.elements <- colmat[rvd.factor]
          par(mar = rep(0, 4))
          image(z = matrix(1:(nelements), ncol = 1), col = color.elements, axes = F)
          text(seq(0,1, length.out = nelements), 0, labels = truncTale, cex = 1, font = 2, col = "black")
          if (nelements > 1) {
            abline(v = seq(0.5/(nelements-1), 1-.5/(nelements-1), length.out = nelements -1), col = inner_sepcol, lty = 1)
          }
          if (sepwid > 0) {
            par(mar = rep(0, 4))
            box(lwd = sepwid/2, col = sepcol)}
        }

      }
    }

    ## dendrograms
    par(mai = rep(0, 4))
    plot(as.dendrogram(Strain_HC), horiz = TRUE, axes = FALSE, leaflab = "none", yaxs = "i")
    par(mai = rep(0, 4))
    plot(as.dendrogram(TALE_HC), axes = FALSE, leaflab = "none", xaxs = "i")

    if (!is.null(extraCol)) {
      ## extra column
      # extra.bar <- sample(c("Hanoi", "Hatay", "Namdinh", NA), nrow(uniqueRVD), replace = T)
      extra.bar <- sapply(rownames(uniqueRVD), function(s) {
        ext <- unique(tale_annotation$extraCol[tale_annotation$strain == s])
        return(ext)
      })
      rextra.bar <- unique(extra.bar[!is.na(extra.bar)])
      rextra.bar <- data.frame("lab" = sort(rextra.bar), "fac" = 1:length(rextra.bar))
      extra.col <- viridis::viridis(n = nrow(rextra.bar))[sapply(extra.bar, function(e) ifelse(is.na(e), NA, rextra.bar$fac[rextra.bar$lab == e]), simplify = T)]
      extra.col[is.na(extra.col)] <- "gray"
      par(mar = c(0, 0, 0, 0))
      image(z = matrix(1:nrow(uniqueRVD), nrow = 1), col = rev(extra.col), yaxt = "n", xaxt = "n", axes = F)
    } else {
      par(mar = c(0, 0, 0, 0))
      image(z = matrix(1:nrow(uniqueRVD), nrow = 1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    }


    ## rownames
    rownames(uniqueRVD) <- paste0(rownames(uniqueRVD),"  #", rowSums(uniqueRVD, na.rm = T))
    par(mar = c(0, 0.5, 0, 0))
    image(z = matrix(1:nrow(uniqueRVD), nrow = 1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    text(-1, seq(0, 1, length.out = nrow(uniqueRVD)), labels = rev(rownames(uniqueRVD)), font = 1, col = "black", cex = 1.2, adj = 0)
    mtext(side = 4, at = .5, text = "Strain", col = "black", padj = 0, line = -1)

    if (!is.null(extraCol)) {
      ## legend column
      par(mar = c(0,0.5,1,0))
      plot(rep(0, nrow(rextra.bar)), -seq(from = 0, by = .8, length.out = nrow(rextra.bar)), type = "p", pch = 15, col = viridis::viridis(n = nrow(rextra.bar)), axes = F, main = extraCol, xlab = NA, ylab = NA, cex = 4, ylim = c(-nrow(uniqueRVD), 0), xlim = c(0,2))
      text(rep(.3, nrow(rextra.bar)), -seq(from = 0, by = .8, length.out = nrow(rextra.bar)), labels = rextra.bar$lab, font = 1, col = "black", bg = "red", cex = 1.2, adj = 0)
    } else {
      par(mar = c(0,0,0,0))
      image(matrix(0), col = "white", axes = F)
    }


    ## colnames
    par(mar = c(0.5, 0, 0.5, 0))
    image(z = matrix(1:ncol(uniqueRVD), ncol = 1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    text(seq(0, 1, length.out = ncol(uniqueRVD)), 1, labels = colnames(uniqueRVD), font = 1, col = "black", bg = "red", cex = 1.2, srt = 90, adj = 1)
    mtext(side = 1, at = 0.5, text = "TALE Group", col = "black", padj = 0, line = -1)


    ## title
    par(xpd = T, mai = rep(0, 4))
    image(z = matrix(1), col = "white", yaxt = "n", xaxt = "n", axes = F)
    text(0, 1/3, labels = "RVD sequences variants", font = 2, col = "black", cex = 2, pos = 1)


    par(df)
  }
  if (hasArg(save.path)) {
    dev.off()
  }
}

