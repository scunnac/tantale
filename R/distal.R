# roxygen2::roxygenise() to write documentation

#### Some specific supporting functions ####

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
    outdir <- tempfile(pattern = "runDistal")
    dir.exists(outdir) || dir.create(outdir, recursive = TRUE)
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
  colnames(talsim) <- c("TAL1", "TAL2", "Sim")
  
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

  outputlist <- list("SeqOfRvdAlignments" = RVDgroups,
                     "SeqOfDistancesAlignments" = sim_within_group,
                     "repeatUnitsDistanceMatrix" = dist_between_rep,
                     "SeqOfRepsAlignments" = RepCodegroups,
                     "TALgroups" = talgroups
                     )
  return(outputlist)
  }






