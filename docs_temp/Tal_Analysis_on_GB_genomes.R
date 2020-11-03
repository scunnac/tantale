
###############################################
###############################################
## Code used after submission to Genebank of the malian Xoo genomes for Tal-centered downstream analysis
###############################################
###############################################

###############################################
## Loading function definitions to handle RVD seq convertion and run AnnoTALE and QueTal

AnnoQueTALFunctionsFile <- "/media/cunnac/DONNEES/CUNNAC/Lab-Related/Communications/Papers/NewFrontiers2017/GenomesDownstreamAnalysis/TALE_analysis/AnnoTALE_QueTAL_functions_library.R"

source(file = AnnoQueTALFunctionsFile, echo = TRUE)



###############################################
## Running AnnoTALE on the set of fasta genomes submited to GenBank


inputGenomesDir <- "/media/cunnac/DONNEES/CUNNAC/Lab-Related/Communications/Papers/NewFrontiers2017/MalianGenomesAssembly/GeneBank_submission/FASTA_files"
inputGenomeFastaFiles <- list.files(inputGenomesDir,
                                    pattern = "\\.f",
                                    recursive = TRUE,
                                    full.names = TRUE)

TALEAnalysisDir <- "/media/cunnac/DONNEES/CUNNAC/Lab-Related/Communications/Papers/NewFrontiers2017/GenomesDownstreamAnalysis/TALE_analysis"
annoTALEAnalysisDir <- file.path(TALEAnalysisDir, "AnnoTALE")
annoTALEOutputDirs <- file.path(annoTALEAnalysisDir,
                                     gsub("(.*)\\.f.*$", "\\1", basename(inputGenomeFastaFiles)))

invisible(
  mapply(FUN = analyzeAnnoTALE, inputFastaFile = inputGenomeFastaFiles, outputDir = annoTALEOutputDirs)
)


###############################################
## Reading Malian strain gff3 files created by AnnoTALE and create a single annotation table
# The goal is to create an annotation table
# to keep location of tal on the genomes as well as mapping between tal IDs


MalianTALEsAnnoTALEGff3Files <- list.files(annoTALEAnalysisDir, pattern = "\\.gff3$",
                                             full.names = TRUE,
                                             include.dirs = TRUE,
                                             recursive = TRUE)

MalianTALEsAnnotation <- do.call(rbind,
        lapply(
          X = MalianTALEsAnnoTALEGff3Files,
          FUN = function(x)
          {
            taleInfo <- read.delim(x, header = FALSE)[, c(1, 3, 4, 5, 7, 9)]
            taleInfo <- taleInfo[taleInfo$V3 == "CDS",]
            taleInfo$V9 <- gsub("Parent=(.*)$", "\\1", taleInfo$V9)
            colnames(taleInfo) <- c("Genome", "AnnotationElement", "Start", "End", "Strand", "AnnoTALEID")
            taleInfo
          }
        )
)

MalianTALEsAnnotationFile <- file.path(dirname(annoTALEAnalysisDir), "MalianTALEsAnnotation.txt")


###############################################
## Reading Malian strain RVD seq files as well as TALEs from BAI3/MAI1/CFBP1947 and create properly formated files


# Convert reference TALE RVD seqs from the QueTal to the AnnoTALE format
QueTalRefRVDSeqsFile <- "/media/cunnac/DONNEES/CUNNAC/Lab-Related/Communications/Papers/NewFrontiers2017/GenomesDownstreamAnalysis/TALE_analysis/TUStrains_TALEs_RVDSeqs.QueTal.fasta"
AnnoTALERefRVDSeqsFile <- file.path(dirname(QueTalRefRVDSeqsFile), "TUStrains_TALEs_RVDSeqs.AnnoTALE.fasta")
QueTALRVD2AnnoTALE(QueTalRefRVDSeqsFile, AnnoTALERefRVDSeqsFile)

# Concatenate all the "TALE_RVDs.fasta" files for the set of MAlian genomes and save to a fasta file

MalianTALEsRVDSeqAnnoTALEFiles <- list.files(annoTALEAnalysisDir, pattern = "^TALE_RVDs\\.fasta$",
                           full.names = TRUE,
                           include.dirs = TRUE,
                           recursive = TRUE)
allMalianTALEsRVDSeqAnnoTALEFile <- file.path(dirname(QueTalRefRVDSeqsFile), "Malian_TALEs_RVDSeqs.AnnoTALE.fasta")

allMalianTALEsRVDSeqAnnoTALE <- unlist(lapply(MalianTALEsRVDSeqAnnoTALEFiles, FUN = readLines, warn = TRUE, skipNul = TRUE)) # read files content
allMalianTALEsRVDSeqAnnoTALE <- allMalianTALEsRVDSeqAnnoTALE[allMalianTALEsRVDSeqAnnoTALE != ""] # Removing empty lines
allMalianTALEsRVDSeqAnnoTALE <- gsub("(^.*TALE\\d+).*$", "\\1", allMalianTALEsRVDSeqAnnoTALE, perl = TRUE) # Removing unecessary info in sequence titles
writeLines(allMalianTALEsRVDSeqAnnoTALE, allMalianTALEsRVDSeqAnnoTALEFile) # Write resulting seq to disk

# Include Malian genomes TALEs RVD seqs in the annotation table
library(Biostrings)
MalianTALERVDSeqs <- readBStringSet(filepath = allMalianTALEsRVDSeqAnnoTALEFile) # read RVD seqs
MalianTALERVDSeqs <- cbind(names(MalianTALERVDSeqs), as.character(MalianTALERVDSeqs)) # make a df using the bsstringset
colnames(MalianTALERVDSeqs) <-c("AnnoTALEID", "RVDSeq")

MalianTALEsAnnotation <- merge(MalianTALEsAnnotation, MalianTALERVDSeqs) # incorporate RVD seqs in annotation table
write.csv(x = MalianTALEsAnnotation, file = MalianTALEsAnnotationFile, row.names = FALSE) # write table to disk


# Merge reference and malian tales in a single file using the AnnoTALE fasta format
AnnoTALEAllRVDSeqsFile <- file.path(dirname(QueTalRefRVDSeqsFile), "All_TALEs_RVDSeqs.AnnoTALE.fasta")
writeLines(unlist(sapply(X = c(AnnoTALERefRVDSeqsFile, allMalianTALEsRVDSeqAnnoTALEFile), FUN = readLines, USE.NAMES = FALSE)), AnnoTALEAllRVDSeqsFile)

# Convert all tales to the QueTal format and save to file
QueTalAllRVDSeqsFile <- file.path(dirname(QueTalRefRVDSeqsFile), "All_TALEs_RVDSeqs.QueTal.fasta")
AnnoTALE2QueTALRVD(AnnoTALEAllRVDSeqsFile, QueTalAllRVDSeqsFile)


###############################################
## Running AnnoTALE build including ref TALEs from BAI3/MAI1/CFBP1947 to define membership to groups

annoTALEBuildDir <- unique(file.path(dirname(AnnoTALEAllRVDSeqsFile), "AnnoTALEBuild_ALL"))
dir.exists(annoTALEBuildDir) || dir.create(path = annoTALEBuildDir, showWarnings = TRUE, recursive = TRUE, mode = "0775")

TALEDNApartsFile <- file.path(annoTALEBuildDir, "all_TALE_DNA_parts.fasta")

TALEDNAparts <- unlist(sapply(TALEDNAFiles, readLines))
writeLines(TALEDNAparts, TALEDNApartsFile)

buildAnnoTALE(AnnoTALEAllRVDSeqsFile, annoTALEBuildDir)


###############################################
## Plot the Distal tree of African TALEs

library(ape)
library(ggtree)
library(scales)
library(dichromat)


treeFile <- "/media/cunnac/DONNEES/CUNNAC/Lab-Related/Communications/Papers/NewFrontiers2017/GenomesDownstreamAnalysis/TALE_analysis/AfTALs_DisTALresults_ALPQ/Output.tre"

# Read info from file
tree <- read.tree(file = treeFile)
tree <- phangorn::midpoint(tree)
#tree <- reorder(tree)

# Parsing TAL labels
tipLabelsAnnot <- as.data.frame(do.call(rbind, strsplit(tree$tip.label, split = "|", fixed = TRUE)))
names(tipLabelsAnnot) <- c("StrainGroup", "Strain", "TALEGroup")
tipLabelsAnnot <- cbind(tip = tree$tip.label, tipLabelsAnnot)

# Get strain annotation from file
annotFile <- "/media/cunnac/DONNEES/CUNNAC/Lab-Related/Communications/Papers/NewFrontiers2017/GenomesDownstreamAnalysis/phylogenomics/strainAnnotations.txt"
strainsAnnot <- read.delim(annotFile, stringsAsFactors = FALSE)

# Merging all annotations
tipLabelsAnnot <- merge(tipLabelsAnnot, strainsAnnot, by.x = "Strain", by.y = "prettyName", sort = FALSE)
tipLabelsAnnot <- tipLabelsAnnot[, c(2,1, 3:ncol(tipLabelsAnnot))]

tipsByGroup <- by(data = tipLabelsAnnot$tip, INDICES = tipLabelsAnnot$TALEGroup, FUN = function(x) x)


# Building the ggtree tree


myFirstPalette <- c("#000000", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[seq(4, 8, by = 2)]
mySecondPalette <- dichromat_pal("Categorical.12")(12)[seq(4, 12, by = 4)]
myGreyPalette <- rev(RColorBrewer::brewer.pal(9,"Greys"))


myGreyPalette <- rev(RColorBrewer::brewer.pal(3,"Greys"))[c(1,2,3,2,1,2,3,1,3)]


#dichromat_pal("Categorical.12")(12)

p <- ggtree(tree, ladderize = TRUE, layout = 'fan') %<+% tipLabelsAnnot
for (i in seq_along(tipsByGroup)) {
  p <- p + geom_hilight(node= MRCA(p, tipsByGroup[[i]]),
                        fill= myGreyPalette[i],
                        alpha= 0.3,
                        extendto = 1.4
                        ) +
    geom_cladelabel(node= MRCA(p, tipsByGroup[[i]]),
                    offset = 0.011,
                    label= names(tipsByGroup)[i],
                    align=T,
                    angle=0,
                    hjust='center',
                    offset.text = 0.25,
                    fontsize = 4.5,
                    barsize = 0,
                    color = "Black"
                    )
}

p$data$label <- do.call(rbind, strsplit(p$data$label, split = "|", fixed = TRUE))[, 2]
#p$data$label <- gsub("Xooa_", "", chartr("|", "_", p$data$label))
p <- p +  geom_tiplab2(aes(color = country), size = 4.5, offset = 0.01) +
  scale_color_manual(values = myFirstPalette) +
#  scale_color_brewer(type = "qual", palette = "Dark2") +
  theme(legend.position = "top")
p
ggsave(filename = "distalTree.svg", plot = p, width = 10, height = 12)




###############################################
## Building a TALE RVD seq alleles table and display as a heatmap

library(ape)
library(reshape2)
library(RColorBrewer)

MalianTALEsAnnotationFileStable <- file.path(TALEAnalysisDir, "MalianTALEsAnnotationStable.csv")

MalianTALEsAnnotationStable <- read.csv(MalianTALEsAnnotationFileStable)

# Building a cross tabulation table with TALE alleles in cells
TALEAlleleTable <- dcast(MalianTALEsAnnotationStable[ ,c("TALE_group", "Genome", "RVDSeq")],
#                         fill = "Missing",
                         Genome ~ TALE_group, value.var = 'RVDSeq',
                         drop = FALSE)

codedAlleles <- as.data.frame(lapply(X = TALEAlleleTable[,-1], FUN = factor)) # Converting RVD sequences per TAL group to factors
variantsCount <- sapply(codedAlleles, nlevels) # extracting the count of distinct alleles
colnames(codedAlleles) <- paste(colnames(codedAlleles), variantsCount, sep = " #")

codedAlleles <- do.call(cbind, lapply(X = codedAlleles, FUN = function(x) {
  # recoding allele ID to have increasing integer ID as a function of allele abundance
  levels(x) <- nlevels(x) + 1 - rank(table(x), ties.method = "first")
  as.integer(as.character(x))
  }
  )
)
rownames(codedAlleles) <- TALEAlleleTable$Genome


# Performing HC for rows and column of the heat map
Strain_HC <- hclust(dist.gene(x = codedAlleles), method = "average")
TALE_HC <- hclust(dist.gene(x = t(codedAlleles)), method = "average")

# Assembling the heat map figure:
# myPalette <- colorRampPalette(colors=c("yellow","brown","blue", "cyan"))(4)
#myPalette <- brewer.pal(n = max(codedAlleles), name = "Dark2")
myPalette <- c("#000000", "#E69F00", "#56B4E9",
                    "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[seq(1, 8, length.out = max(codedAlleles))]

myPlot <- gplots::heatmap.2(codedAlleles,
                            trace = "none",
                            col = myPalette,
                            density.info = "none",
                            key = FALSE,
                            xlab = "TALE group",
                            ylab = "Strain",
                            margins = c(7, 7),
                            colsep = 0:(ncol(codedAlleles)-0),
                            rowsep = 0:(nrow(codedAlleles)-0),
                            sepcolor="white",
                            sepwidth=c(0.07,0.07),
                            Rowv = as.dendrogram(Strain_HC),
                            Colv = as.dendrogram(TALE_HC),
                            main = "RVD sequences variants"
)


# export as svg:

svg(filename = "RVDseqAllelesTable.svg", width = 7, height = 7 )
# Copy and paste the LATEST version of the code for the heatmap
myPlot <- gplots::heatmap.2(codedAlleles,
                            trace = "none",
                            col = myPalette,
                            density.info = "none",
                            key = FALSE,
                            xlab = "TALE group",
                            ylab = "Strain",
                            margins = c(7, 7),
                            colsep = 0:(ncol(codedAlleles)-0),
                            rowsep = 0:(nrow(codedAlleles)-0),
                            sepcolor="white",
                            sepwidth=c(0.07,0.07),
                            Rowv = as.dendrogram(Strain_HC),
                            Colv = as.dendrogram(TALE_HC),
                            main = "RVD sequences variants"
)
dev.off()


###############################################
## Building a TALE repeat type seq alleles (DisTAL generated) table and display as a heatmap

library(Biostrings)
library(ape)
library(reshape2)
library(RColorBrewer)


disTalString <- readBStringSet("/media/cunnac/DONNEES/CUNNAC/Lab-Related/Communications/Papers/NewFrontiers2017/GenomesDownstreamAnalysis/TALE_analysis/AfTALs_DisTALresults_ALPQ/Output_CodedRepeats.fa")
disTalStringAnnot <- do.call(rbind, strsplit(names(disTalString), split = "|", fixed = TRUE))
colnames(disTalStringAnnot) <- c("StrainGroup", "Strain", "TALEGroup")
disTalStringAnnot <- cbind(as.data.frame(disTalStringAnnot[, -1]), disTalString = as.factor(gsub("^\\s+|\\s+$", "", as.character(disTalString))))

# Building a cross tabulation table with TALE alleles in cells
TALEAlleleTable <- dcast(disTalStringAnnot,
                         Strain ~ TALEGroup,
                         value.var = "disTalString",
                         drop = FALSE)

codedAlleles <- as.data.frame(lapply(X = TALEAlleleTable[,-1], FUN = factor)) # Converting repeat types sequences per TAL group to factors
variantsCount <- sapply(codedAlleles, nlevels) # computing the counts of distinct alleles
colnames(codedAlleles) <- paste(colnames(codedAlleles), variantsCount, sep = " #")

codedAlleles <-
  do.call(cbind, lapply(
    X = codedAlleles,
    FUN = function(x) {
      # recoding allele ID to have increasing integer ID as a function of allele abundance (most abundant is first)
      levels(x) <- nlevels(x) + 1 - rank(table(x), ties.method = "first")
      as.integer(as.character(x))
    }
  ))
rownames(codedAlleles) <- TALEAlleleTable$Strain


# Performing HC for rows and column of the heat map
Strain_HC <- hclust(dist.gene(x = codedAlleles), method = "average")
TALE_HC <- hclust(dist.gene(x = t(codedAlleles)), method = "average")

# Assembling the heat map figure:
# myPalette <- colorRampPalette(colors=c("yellow","brown","blue", "cyan"))(4)
#myPalette <- brewer.pal(n = max(codedAlleles), name = "Dark2")
myPalette <- c("#000000", "#E69F00", "#56B4E9",
               "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[seq(1, 8, length.out = max(codedAlleles))]

myPlot <- gplots::heatmap.2(codedAlleles,
                            trace = "none",
                            col = myPalette,
                            density.info = "none",
                            key = FALSE,
                            xlab = "TALE group",
                            ylab = "Strain",
                            margins = c(7, 7),
                            colsep = 0:(ncol(codedAlleles)-0),
                            rowsep = 0:(nrow(codedAlleles)-0),
                            sepcolor="white",
                            sepwidth=c(0.07,0.07),
                            Rowv = as.dendrogram(Strain_HC),
                            Colv = as.dendrogram(TALE_HC),
                            main = "DisTAL sequences variants"
)

# export as svg:

svg(filename = "DisTALseqAllelesTable.svg", width = 7, height = 7 )
# Copy and paste the LATEST version of the code for the heatmap
myPlot <- gplots::heatmap.2(codedAlleles,
                            trace = "none",
                            col = myPalette,
                            density.info = "none",
                            key = FALSE,
                            xlab = "TALE group",
                            ylab = "Strain",
                            margins = c(7, 7),
                            colsep = 0:(ncol(codedAlleles)-0),
                            rowsep = 0:(nrow(codedAlleles)-0),
                            sepcolor="white",
                            sepwidth=c(0.07,0.07),
                            Rowv = as.dendrogram(Strain_HC),
                            Colv = as.dendrogram(TALE_HC),
                            main = "DisTAL sequences variants"
)

dev.off()


##############################################
## Generating matrix of shared targets ratio based on target predicitons

library(reshape2)

predicitonsFile <- "/media/cunnac/DONNEES/CUNNAC/Lab-Related/Communications/Papers/NewFrontiers2017/GenomesDownstreamAnalysis/TALE_analysis/PredictedTargetsProfiles_ALPQ/AfTALS_vs_Nip500_complete"

predictionsComplete <- read.delim(predicitonsFile, stringsAsFactors = FALSE)

predictions <- subset(predictionsComplete, subset = RANK <= 200, select = c("TAL_ID", "SEQ_ID"))

talIDs <- unique(predictions$TAL_ID)
pairsOfTal <- cbind(rep(talIDs, each = length(talIDs)), rep(talIDs, time = length(talIDs)))

targetOverlapRatios <- lapply(1:nrow(pairsOfTal), function(x, df, pairs) {
  t1targets <- unique(df[df$TAL_ID == pairs[x, 1], "SEQ_ID"])
  t2targets <- unique(df[df$TAL_ID == pairs[x, 2], "SEQ_ID"])
  #alltargets <- unique(c(t1targets, t2targets))
  #length(intersect(t1targets, t2targets))/length(alltargets)*100
  length(intersect(t1targets, t2targets))/length(t1targets)*100 # this value will depend on the orientation of the pair of TALE: the matrix will not be symetrical
}, df = predictions, pairs = pairsOfTal)

targetOverlapRatios <- data.frame(tal1 = pairsOfTal[,1], tal2 = pairsOfTal[,2], overlap = unlist(targetOverlapRatios))
head(targetOverlapRatios)

targetOverlapRatiosMatrix <- acast(targetOverlapRatios, formula = tal1 ~ tal2, value.var = "overlap")



###############################################
## Generating TALE group wise outpout heatmap for the shared target profiles of Af strains

library(ape)
library(reshape2)
library(RColorBrewer)

TargetProfileMatrix <- round(targetOverlapRatiosMatrix, digits = 0)
# TargetProfileMatrix <- read.csv(file = "/media/cunnac/DONNEES/CUNNAC/Lab-Related/Communications/Papers/NewFrontiers2017/GenomesDownstreamAnalysis/TALE_analysis/PredictedTargetsProfiles_ALPQ/OLD/Shared_targets_SEB_MODALPQ.txt",
#                                 row.names = 1,
#                                 check.names = FALSE)

# Parsing row and column names to extract info
colSplit <- do.call(rbind, strsplit(colnames(TargetProfileMatrix), split = "|", fixed = TRUE))
TALEGroups <- unique(colSplit[, 3])
strains <- unique(colSplit[, 2])


# Generating Heat maps of target set overlap pourcentages for each TALE group


for (TALEGroup in TALEGroups[c(-3, -5)]) {

  #TALEGroup <- TALEGroups[-3]

  # Selecting target set overlap pourcentage values for a specified TALE group
  TargetProfileForGroup <- as.matrix(TargetProfileMatrix[grepl(TALEGroup, rownames(TargetProfileMatrix)), grepl(TALEGroup, colnames(TargetProfileMatrix))])

  # Simplifiying row and col labels
  rownames(TargetProfileForGroup) <- gsub(paste0(".+Xooa\\|(.*)\\|", TALEGroup), "\\1", rownames(TargetProfileForGroup))
  colnames(TargetProfileForGroup) <- gsub(paste0(".+Xooa\\|(.*)\\|", TALEGroup), "\\1", colnames(TargetProfileForGroup))

  # Performing HC for rows and column of the heat map
  Strain_HC <- hclust(dist(x = TargetProfileForGroup, method = "euclidean"), method = "average")

  # Assembling the heat map figure:
  # myPalette <- colorRampPalette(colors=c("yellow","brown","blue", "cyan"))(4)
  myPalette <- brewer.pal(n = 10, name = "PiYG") # colorblind friendly colorbrewer palette


  svgFile <- paste0(TALEGroup, "_sharedTargets_heatmap.svg")

  # export heatmap as svg:
  svg(filename = svgFile, width = 7, height = 7 )

  myPlot <- gplots::heatmap.2(TargetProfileForGroup,
                              symm = TRUE,
                              trace = "none",
                              col = myPalette,
                              breaks = seq(from = 0, to = 100, by = 10),
                              density.info = "none",
                              key = TRUE,
                              cellnote = TargetProfileForGroup,
                              notecex = 1.2,
                              notecol = "black",
                              margins = c(7, 7),
                              lwid = c(0.15, 0.85),
                              lhei = c(0.15, 0.7),
                              colsep = 0:(ncol(TargetProfileForGroup) - 0),
                              rowsep = 0:(nrow(TargetProfileForGroup) - 0),
                              sepcolor="white",
                              sepwidth=c(0.07,0.07),
                              Rowv = as.dendrogram(Strain_HC),
                              Colv = as.dendrogram(Strain_HC),
                              dendrogram = "row",
                              xlab = "Strain",
                              ylab = "Strain",
                              main = paste(TALEGroup, "group")
  )
  dev.off()

}

###############################################
## Status of genuine targets of characterised TALE as predicted target of variants in the same group

controlRegulatoryPairs <- data.frame(
  tale = c("TalF", "TalC", "TalB", "TalB"),
  targetGeneName = c("OsSWEET14" , "OsSWEET14", "OsERF123", "OsTFX1"),
  targetLOCID = c("LOC_Os11g31190", "LOC_Os11g31190", "LOC_Os09g39810", "LOC_Os09g29820"),
  stringsAsFactors = FALSE
)

controlsPredictions <- lapply(1:nrow(controlRegulatoryPairs), FUN = function(i) {
    onePair <- controlRegulatoryPairs[i, ]
    df <- subset(predictionsComplete, subset = grepl(onePair$tale, TAL_ID) & grepl(onePair$targetLOCID, SEQ_ID),
                 select = c(TAL_ID, SEQ_ID, SCORE, RANK, TALBS_sequence, Gene_start, Gene_end))
    df
  }
)
controlsPredictions <- do.call(rbind, controlsPredictions)
controlsPredictions <- merge(controlsPredictions, unique(controlRegulatoryPairs[, 2:3]), by.x = "SEQ_ID", by.y = "targetLOCID",
      all.x = FALSE,
      sort = FALSE)
controlsPredictions$TAL_ID <- gsub(">Xooa\\|", "", controlsPredictions$TAL_ID)

controlsPredictions <- cbind(do.call(rbind, strsplit(controlsPredictions$TAL_ID , split = "|", fixed = TRUE))[,2:1],
      controlsPredictions)

colnames(controlsPredictions) <- c("TALE_Group",
                                   "strain",
                                   "target_LOCID",
                                   "TALE_ID",
                                   "score",
                                   "rank",
                                   "EBE_seq",
                                   "gene_start",
                                   "gene_end",
                                   "target_gene_name"
                                   )

write.csv(controlsPredictions, "predictionRankOfControlTargets.csv", row.names = FALSE)


###############################################
## Extracting the regions of talF and talB for 'promoter' and protein multiple alignement.

library(Biostrings)
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(msa)

# Get tal genes positions from the stable annotation table
MalianTALEsAnnotationFileStable <- file.path(TALEAnalysisDir, "MalianTALEsAnnotationStable.csv")
MalianTALEsAnnotationStable <- read.csv(MalianTALEsAnnotationFileStable)

# Convert locations to a GRangesList
talGRsLst <- makeGRangesListFromDataFrame(df = MalianTALEsAnnotationStable,
                                       split.field = "TALE_group",
                                       names.field = "PublicID",
                                       keep.extra.columns = TRUE,
                                       ignore.strand = FALSE,
                                       seqinfo = NULL,
                                       seqnames.field = "Genome",
                                       start.field = "Start",
                                       end.field = c("End"),
                                       strand.field = "Strand",
                                       starts.in.df.are.0based = FALSE)

#TALEGroupSOfInterest <- c("TalB", "TalF")
TALEGroupSOfInterest <- levels(MalianTALEsAnnotationStable$TALE_group)


for (TALEGroupOfInterest in TALEGroupSOfInterest) {
  # Select GRanges corresponding to a TALE group of interest
  talGRsLstOfInterest <- unlist(talGRsLst[names(talGRsLst) %in% TALEGroupOfInterest])

  # Get promoter sequences and reorder them
  talPromotersSeqs <- getPromoterSeq(query = talGRsLstOfInterest, subject = allGenomes, upstream = 200, downstream = 6)
  properStrainsOrder <- c("MAI68", "MAI134", "MAI1", "MAI73", "MAI95", "MAI99", "MAI106", "MAI129", "MAI145", "BAI3", "CFBP1947")
  talPromotersSeqs <- talPromotersSeqs[properStrainsOrder]

  # Perform multiple alignment and print a pretty version of it
  msaLabel <- paste(TALEGroupOfInterest, "promoters", "msa", sep = "_")
  myMsa <- msa(talPromotersSeqs, method = "ClustalW", verbose = FALSE, order = "input")

  msaPrettyPrint(myMsa,
                 file = paste(msaLabel, "pdf", sep = "."),
                 showConsensus = "none",
                 showLogo = "none",
                 paperWidth = 7)
}

###############################################
## Try to generate a Figure illustrating tal genes synteny across genomes

library(ggbio)

unlist(talGRsLst)
names(talGRsLst)
lapply(talGRsLst, mcols)
talGRsLst$TalA


p <- ggplot(talGRsLst) + geom_alignment(aes(color = strand, fill = strand)) + facet_grid(facets = seqnames ~ .)
p
ggsave(filename = "talGenesSyntheny.svg", height = 13, width = 8)



 autoplot(talGRsLst) + facet_grid(facets = seqnames ~ .)


talGR <- makeGRangesFromDataFrame(df = MalianTALEsAnnotationStable,
                                       keep.extra.columns = TRUE,
                                       ignore.strand = FALSE,
                                       seqinfo = NULL,
                                       seqnames.field = "Genome",
                                       start.field = "Start",
                                       end.field = c("End"),
                                       strand.field = "Strand",
                                       starts.in.df.are.0based = FALSE)
names(talGR) <- mcols(talGR)$PublicID


ggplot(talGR) + geom_rect(aes(color = strand, fill = strand)) + facet_grid(facets = seqnames ~ .) + geom_text(aes(label = TALE_group))

ggplot(talGR[!seqnames(talGR) %in% "CFBP1947"]) + geom_arrow(aes(color = strand, fill = strand)) +
  facet_grid(facets = seqnames ~ TALE_group, scales = "free")


###############################################
## Generate colored alignment of a group of TALE RVD seqs based on how they fit to a given EBE


# Read Talvez Nucl-RVD association matrix
RVDNucAssocMatFile <- "/media/cunnac/DONNEES/CUNNAC/Lab-Related/MyScripts/TALETargets_Interact/Version_3-0/PredictorScripts/TALVEZ/mat1"
mat <- read.delim(file = RVDNucAssocMatFile, header = FALSE, na.strings = "", stringsAsFactors = FALSE)
rownames(mat) <- mat[ , 1]
mat[1] <- NULL
colnames(mat) <- c("A", "C", "G", "T")


# Fetch RVD seq of interest in malian Tal annotation table
MalianTALEsAnnotationFileStable <- file.path(TALEAnalysisDir, "MalianTALEsAnnotationStable.csv")
MalianTALEsAnnotationStable <- read.csv(MalianTALEsAnnotationFileStable)


# TalB variants on ERF EBE
strainsOfInterest <- c("MAI1", "BAI3", "MAI68", "MAI134")

talBVariants <- subset(MalianTALEsAnnotationStable,
                       subset = TALE_group == "TalB" & Genome %in% strainsOfInterest,
                       select = c(Genome, TALE_group, RVDSeq))
talBVariantsSeqs <- Biostrings::BStringSet(talBVariants$RVDSeq)
names(talBVariantsSeqs) <- paste(talBVariants$Genome, talBVariants$TALE_group, sep = "_")
talBVariantsSeqs <- as.character(talBVariantsSeqs)

EBE_TalB_ERF <- c(EBE_TalB_ERF = "TGCGATGCGTTTCCCACCTCCCACCTC")

plotRVDSeqsOnEBE(RVDSeqs = talBVariantsSeqs,
                 EBESeq = EBE_TalB_ERF,
                 mat = mat,
                 height = 2.5)



# TalB variants on TFX1 EBE
EBE_TalB_TFX1 <- c(EBE_TalB_TFX1 = "TAAAAGGCCCTCACCAACCCATCGCCT")

plotRVDSeqsOnEBE(RVDSeqs = talBVariantsSeqs,
                 EBESeq = EBE_TalB_TFX1,
                 mat = mat,
                 height = 2.5)


# TalF variants on SWEET14 EBE

strainsOfInterest <- c("MAI1", "BAI3", "MAI68", "CFBP1947")
talFVariants <- subset(MalianTALEsAnnotationStable,
                       subset = TALE_group == "TalF" & Genome %in% strainsOfInterest,
                       select = c(Genome, TALE_group, RVDSeq))
talFVariantsSeqs <- BStringSet(talFVariants$RVDSeq)
names(talFVariantsSeqs) <- paste(talFVariants$Genome, talFVariants$TALE_group, sep = "_")
talFVariantsSeqs <- as.character(talFVariantsSeqs)

EBE_TalF_SWEET14 <- c(EBE_TalF_SWEET14 = "TAAGCTCATCAAGCCTTCA")

plotRVDSeqsOnEBE(RVDSeqs = talFVariantsSeqs,
                 EBESeq = EBE_TalF_SWEET14,
                 mat = mat,
                 height = 2.5)



#c(TalB_MAI1 = "NN-ND-NN-NI-NN-NN-ND-NN-NG-NG-N*-ND-NG-HD-NN-NN-HD-NG-HD-HD-HD-NN-NN-HD-HD-NG")

#c(EBE_TalB_ERF = "TGCGATGCGTTTCCCACCTCCCACCTC")
#RVDSeqsFile <- "/media/cunnac/DONNEES/CUNNAC/Temp/TalBsRVDSeqs.fas"
#RVDSeqs <- readBStringSet(filepath = RVDSeqsFile)
#RVDSeqs <- as.character(RVDSeqs)
