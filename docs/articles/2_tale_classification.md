# 2. Classifying TALE sequences from genomes with tantale

``` r
library(tantale)
library(ggplot2)
library(tidyverse)
library(Biostrings)
library(parallel)
library(magrittr)
library(DT)
library(gplots)
library(ape)
library(ggtree)
library(corrr)
library(ggcorrplot)
```

Based on the results of the first article “Mining TALE sequences in
genomes with tantale”, this second section will cover steps downstream
of *tal* loci prediction. It will illustrate TALE classification into
groups of related sequences with functions in package `tantale`.  
As a reminder, the genomes used for this article are MAI1, BAI3,
BAI3-1-1 (an error-prone genome with TalC deleted) and PXO86, an
outgroup Asian strain.  

``` r
outdir <- fs::dir_create("~/TEMP/test_tantale") # tempdir(check = TRUE)
load(file.path(outdir, "mining.RData"))
```

## Classify the various TALE sequences into groups of related sequences

To classify TALE groups, we use the function
[`runDistal()`](https://scunnac.github.io/tantale/reference/runDistal.md)
to run the tool DisTal that outputs similarity between TALEs and
similarity between repeats. The TALE similarity table can be taken as
input of
[`groupTales()`](https://scunnac.github.io/tantale/reference/groupTales.md)
for clustering groups, and the repeats similarity can be used for
plotting multiple alignment of Tals within group.

### Run DisTal

We run DisTal for all the putative TALE ORFs (from `telltale` output).

``` r
distal_telltale <- file.path(outdir, "distal")
unlink(distal_telltale, recursive = TRUE)
dir.create(distal_telltale, showWarnings = F)
telltale_cds <- file.path(distal_telltale, "cds.fasta")
unlink(telltale_cds)


distal_telltale_input <- list.files(telltale_outdir,
                                    "putativeTalOrf.fasta",
                                    recursive = T, full.names = T)
distal_telltale_input <- distal_telltale_input[!grepl("ROI", distal_telltale_input, fixed = T)]
for (input in distal_telltale_input) {
  ttRunId <- gsub("\\_.+", "", basename(dirname(input)))
  seqs <- readBStringSet(input, seek.first.rec = T)
  names(seqs) <- paste0(ttRunId, "_", names(seqs))
  seqs <- seqs[width(seqs) != 0]
  #cat(names(seqs))
  writeXStringSet(seqs, telltale_cds, append = T)
}

# run DisTal
system.time({
  distal_output <- runDistal(telltale_cds, outdir = distal_telltale, overwrite = T)
})
#>>    user  system elapsed 
#>> 628.701   0.121 628.859
```

### Run the R implementation of DisTal

First we use the tellTale output directories for each genome and get the
information about TALE parts (ie domains). The final output is a tibble
with a single part per row.

``` r
# distalrOutDir <- file.path(outdir, "distalr")
# unlink(distalrOutDir, recursive = TRUE)
# dir.create(distalrOutDir, showWarnings = F)


taleParts <- lapply(list.dirs(telltale_outdir, recursive = FALSE), function(ttGenomeDir) {
  genomeId <- gsub("_.*correction", "", basename(ttGenomeDir))
  logger::log_debug("Current genome: {genomeId}")
  taleParts <- getTaleParts(ttGenomeDir) %>%
    mutate(arrayID = paste0(genomeId, "_", arrayID))
  return(taleParts)
}) %>%
  bind_rows()

DT::datatable(head(taleParts), caption = "An example of the content of the taleParts tibble.")
```

With this kind of representation it is quite easy to make summary plots
like this one for example:

``` r
partsForPlots <- taleParts %>%
  mutate(label = if_else(domainType == "repeat", rvd, ""),
         aaSeqLength = factor(nchar(aaSeq))
  )
partsForPlots %>% ggplot(mapping = aes(fill = aaSeqLength,
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
```

![](2_tale_classification_files/figure-html/Overview%20of%20TALE%20composition%20by%20genome-1.png)

Now, with this `tibble`, we can also run the R implementation of the
Distal Perl code.

``` r
partsForDistalr <- taleParts

system.time({distalr_bs_output <- distalr(taleParts = partsForDistalr,
                             repeats.cluster.h.cut = 10, ncores = 6,
                             pairwiseAlnMethod = "Biostrings")})
#>>     user   system  elapsed 
#>> 1954.326    2.380  423.562

system.time({distalr_deci_output <- distalr(taleParts = partsForDistalr,
                              repeats.cluster.h.cut = 10, ncores = 6,
                              pairwiseAlnMethod = "DECIPHER")})
#>>    user  system elapsed 
#>>  15.718   0.084   4.875

system.time({distalr_mms_output <- distalr(taleParts = partsForDistalr,
                              repeats.cluster.h.cut = 10, ncores = 6,
                              pairwiseAlnMethod = "mmseq2",
                              condaBinPath = "/home/cunnac/bin/miniconda3/condabin/conda")})
#>> Joining with `by = join_by(query, target)`
#>>    user  system elapsed 
#>>  22.367   5.111  15.790
```

Clearly, the “DECIPHER” method for pairwise alignment of domains aa
sequences is much faster than all the others.

### Contrast the output of the various distal functions

Now, we can take a look at whether the results of the original Distal
script are consistent with our R implementation of the program.

``` r
talSimOutAgregated <- list(distalALV = distal_output$tal.similarity,
     distalDECI = distalr_deci_output$tal.similarity,
     distalBIOS = distalr_bs_output$tal.similarity,
     distalMMS = distalr_mms_output$tal.similarity
) %>% bind_rows(.id = "method")

talSimOutAgregatedLong <- reshape2::dcast(talSimOutAgregated %>% select(method, TAL1, TAL2, Sim),
                                          formula = method ~ TAL1 + TAL2, value.var = "Sim"
                                            )
rownames(talSimOutAgregatedLong) <- talSimOutAgregatedLong$method
talSimOutAgregatedLong$method <- NULL
talSimOutAgregatedLong %<>% as.matrix() %>% t()
rownames(talSimOutAgregatedLong) <- NULL

corr_matrix <- cor(talSimOutAgregatedLong)
knitr::kable(corr_matrix)
```

|            | distalALV | distalBIOS | distalDECI | distalMMS |
|:-----------|----------:|-----------:|-----------:|----------:|
| distalALV  | 1.0000000 |  0.9533716 |  0.9547991 | 0.9501474 |
| distalBIOS | 0.9533716 |  1.0000000 |  0.9961920 | 0.9649442 |
| distalDECI | 0.9547991 |  0.9961920 |  1.0000000 | 0.9732519 |
| distalMMS  | 0.9501474 |  0.9649442 |  0.9732519 | 1.0000000 |

``` r
ggcorrplot(corr_matrix)
```

![](2_tale_classification_files/figure-html/correlation%20of%20TALE%20similarity%20from%20various%20distal%20implementations-1.png)

``` r


talSimOutAgregated <- reshape2::dcast(talSimOutAgregated %>% select(method, TAL1, TAL2, Sim),
                                          formula = TAL1 + TAL2 ~ method, value.var = "Sim"
                                            ) %>% as_tibble()
```

``` r
# talSimOutAgregated %>%
#   select(-distalBIOS) %>%
#   mutate(pointLab = paste(TAL1, TAL2, sep = "|")) %>%
#   ggplot(mapping = aes(x = distalALV, distalDECI)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1)

pairs(talSimOutAgregated[3:6], pch = 20, cex = 0.5, lower.panel = NULL)
```

![](2_tale_classification_files/figure-html/pair%20plot%20of%20TALE%20similarity%20from%20various%20distal%20implementations-1.png)

The correlation of the TALE pairwise similarity scores output by the
various distal flavors are remarkably high.

### Allocate TALEs to groups of related sequences

``` r
k <- 22
groupsALV <- invisible(groupTales(taleSim = distal_output$tal.similarity,
                     plotTree = TRUE, k = k, k_test = 10:40,
                     method = "k-medoids"))
#>> tale tree will not be plotted!
#>> Number of groups is decided based on the provided value of k: 22
```

![](2_tale_classification_files/figure-html/tale%20trees-1.png)

``` r

groupsDECI <- invisible(groupTales(taleSim = distalr_deci_output$tal.similarity,
                     plotTree = TRUE, k = k, k_test = 10:40,
                     method = "k-medoids"))
#>> tale tree will not be plotted!
#>> Number of groups is decided based on the provided value of k: 22
```

![](2_tale_classification_files/figure-html/tale%20trees-2.png)

``` r

groupsBIOS <- invisible(groupTales(taleSim = distalr_bs_output$tal.similarity,
                     plotTree = TRUE, k = k, k_test = 10:40,
                     method = "k-medoids"))
#>> tale tree will not be plotted!
#>> Number of groups is decided based on the provided value of k: 22
```

![](2_tale_classification_files/figure-html/tale%20trees-3.png)

``` r

groupsMMS <- invisible(groupTales(taleSim = distalr_mms_output$tal.similarity,
                     plotTree = TRUE, k = k, k_test = 10:40,
                     method = "k-medoids"))
#>> tale tree will not be plotted!
#>> Number of groups is decided based on the provided value of k: 22
```

![](2_tale_classification_files/figure-html/tale%20trees-4.png)

``` r


grp <- groupTales(taleSim = distal_output$tal.similarity,
                     plotTree = TRUE, k = k,
                     method = "hclust")
```

![](2_tale_classification_files/figure-html/tale%20trees-5.png)

``` r
grp <- groupTales(taleSim = distalr_deci_output$tal.similarity,
                     plotTree = TRUE, k = k,
                     method = "hclust")
```

![](2_tale_classification_files/figure-html/tale%20trees-6.png)

``` r
grp <- groupTales(taleSim = distalr_bs_output$tal.similarity,
                     plotTree = TRUE, k = k,
                     method = "hclust")
```

![](2_tale_classification_files/figure-html/tale%20trees-7.png)

``` r
grp <- groupTales(taleSim = distalr_mms_output$tal.similarity,
                     plotTree = TRUE, k = k,
                     method = "hclust")
```

![](2_tale_classification_files/figure-html/tale%20trees-8.png)

We can compare the hierarchical trees obtained with the original Distal
tree (left) and distalr with DICEPHER (right). Because PXO86 is the only
Asian strain and these TALEs are rather distantly related to African
ones, the part of the trees with PXO86 TALEs are not really meaningful
in this context.

``` r
getTALEsTree <- function(taleSim) {
  distMat <- 100 - reshape2::acast(taleSim, formula = TAL1 ~ TAL2, value.var = "Sim")
  hclust(dist(distMat, method = "euclidean"), method = "ward.D")
}

op <- par()
suppressWarnings(par(cex = 0.5))
ape::cophyloplot(getTALEsTree(distal_output$tal.similarity) %>% as.phylo(),
                 getTALEsTree(distalr_deci_output$tal.similarity) %>% as.phylo(),
                 assoc = matrix(rep(distal_output$tal.similarity$TAL1 %>% unique() %>% as.character(), 2),
                                ncol = 2, byrow = FALSE),
                 length.line = 0, space=200, gap = 10,
                 col = "red")
```

![](2_tale_classification_files/figure-html/tale%20trees%20comparisons%20distal%20vs%20deci-1.png)

``` r
suppressWarnings(par(op))
```

``` r
op <- par()
suppressWarnings(par(cex = 0.5))
ape::cophyloplot(getTALEsTree(distal_output$tal.similarity) %>% as.phylo(),
                 getTALEsTree(distalr_mms_output$tal.similarity) %>% as.phylo(),
                 assoc = matrix(rep(distal_output$tal.similarity$TAL1 %>% unique() %>% as.character(), 2),
                                ncol = 2, byrow = FALSE),
                 length.line = 0, space=200, gap = 10,
                 col = "red")
```

![](2_tale_classification_files/figure-html/tale%20trees%20comparisons%20distal%20vs%20mms-1.png)

``` r
suppressWarnings(par(op))
```

``` r
taleGroups <- groupsDECI %>% as_tibble() %>% arrange(group, name)
knitr::kable(head(taleGroups, n= 12))
```

| name               | group |
|:-------------------|------:|
| BAI3_ROI_00001     |     1 |
| MAI1_ROI_00001     |     1 |
| BAI3-1-1_ROI_00001 |     2 |
| BAI3_ROI_00002     |     2 |
| MAI1_ROI_00002     |     2 |
| BAI3-1-1_ROI_00002 |     3 |
| BAI3_ROI_00003     |     3 |
| MAI1_ROI_00003     |     3 |
| BAI3-1-1_ROI_00003 |     4 |
| BAI3_ROI_00004     |     4 |
| MAI1_ROI_00004     |     4 |
| BAI3-1-1_ROI_00005 |     5 |

Based on the silhouette values,
[`groupTales()`](https://scunnac.github.io/tantale/reference/groupTales.md)
will take the value at the elbow of the curve. In this case we manually
took k = 22.

``` r
save.image(file.path(outdir, "mining.RData"))
```
