

#### rvdDistMat ####
rvdToNtAssocMatFile <- system.file("tools", "TALVEZ_3.2", "mat1", package = "tantale", mustWork = T)
rvdToNtAssocMat <- read.delim(file = rvdToNtAssocMatFile, header = FALSE, na.strings = "",row.names = 1)

colnames(rvdToNtAssocMat) <- c("A", "C", "G", "T")
rvdSimMat <- cor(t(rvdToNtAssocMat), method = "spearman")
rvdSimDf <- reshape2::melt(rvdSimMat, varnames = c("rvd1", "rvd2"), value.name = "Cor", as.is = TRUE)


#### save sysdata.rda ####
save(list = c("rvdSimDf", "rvdToNtAssocMat"),
     file = file.path(system.file("R", package = "tantale", mustWork = T), "sysdata.rda")
    )
