# TALEs domains sequences alignment with MAFFT

``` r
library(tidyverse)
library(tantale)
```

This demo will illustrate how to use the `buildRepeatMsa` function to
produce TALE sequence multiple algnments.  

Additionally, using real life distalr output, it will details the
various type of aesthetics that can be used to plot those TALEs msa.

## Producing simple TALE sequence multiple algnments

Here we start with RVD sequences of the African Xoo TALE TalA:

``` r
aln <- buildRepeatMsa(inputSeqs = system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",
                                              package = "tantale", mustWork = TRUE),
                      sep = "-",
                      distalRepeatSims = NULL,
                      mafftOpts = "--localpair --maxiterate 1000 --quiet --reorder --op 0 --ep 5 --thread 1",
                      gapSymbol = "-")
```

The output shown as a simple matrix :

``` r
knitr::kable(aln)
```

|               | 1   | 2   | 3   | 4   | 5   | 6   | 7   | 8   | 9   | 10  | 11  | 12  | 13  | 14  | 15  | 16  | 17  | 18  | 19  | 20  | 21  | 22  | 23  | 24  | 25  | 26  |
|:--------------|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|
| TalA_BAI3     | NN  | NG  | NN  | HD  | HD  | NI  | N\* | NG  | HD  | NI  | NG  | NN  | HD  | NI  | NG  | NI  | NG  | NN  | NG  | HD  | NI  | NI  | NG  | HD  | NN  | NG  |
| TalA_CFBP1947 | NN  | NG  | NN  | HD  | HD  | NI  | N\* | NG  | HD  | NI  | NG  | NN  | HD  | NI  | NG  | NI  | NG  | NN  | NG  | HD  | NI  | NI  | NG  | HD  | NN  | NG  |
| TalA_MAI1     | NN  | N\* | NN  | HD  | HD  | NI  | N\* | NG  | HD  | NI  | NG  | NN  | HD  | NS  | NG  | NI  | NG  | NN  | NG  | HD  | NI  | NI  | NG  | HD  | NN  | NG  |
| TalA_MAI134   | NN  | N\* | NN  | HD  | HD  | NI  | N\* | NG  | HD  | NI  | NG  | NN  | HD  | NS  | NG  | NI  | NG  | NN  | NG  | HD  | NI  | NI  | NG  | HD  | NN  | NG  |
| TalA_MAI106   | NN  | N\* | NN  | HD  | \-  | \-  | \-  | \-  | HD  | NI  | NG  | NN  | HD  | NS  | NG  | NI  | N\* | NN  | NG  | HD  | NI  | NI  | NG  | HD  | NN  | NG  |
| TalA_MAI129   | NN  | N\* | NN  | HD  | \-  | \-  | \-  | \-  | HD  | NI  | NG  | NN  | HD  | NS  | NG  | NI  | N\* | NN  | NG  | HD  | NI  | NI  | NG  | HD  | NN  | NG  |
| TalA_MAI145   | NN  | N\* | NN  | HD  | \-  | \-  | \-  | \-  | HD  | NI  | NG  | NN  | HD  | NS  | NG  | NI  | N\* | NN  | NG  | HD  | NI  | NI  | NG  | HD  | NN  | NG  |
| TalA_MAI68    | NN  | N\* | NN  | HD  | \-  | \-  | \-  | \-  | HD  | NI  | NG  | NN  | HD  | NS  | NG  | NI  | N\* | NN  | NG  | HD  | NI  | NI  | NG  | HD  | NN  | NG  |
| TalA_MAI73    | NN  | N\* | NN  | HD  | \-  | \-  | \-  | \-  | HD  | NI  | NG  | NN  | HD  | NS  | NG  | NI  | N\* | NN  | NG  | HD  | NI  | NI  | NG  | HD  | NN  | NG  |
| TalA_MAI95    | NN  | N\* | NN  | HD  | \-  | \-  | \-  | \-  | HD  | NI  | NG  | NN  | HD  | NS  | NG  | NI  | N\* | NN  | NG  | HD  | NI  | NI  | NG  | HD  | NN  | NG  |
| TalA_MAI99    | NN  | N\* | NN  | HD  | \-  | \-  | \-  | \-  | HD  | NI  | NG  | NN  | HD  | NS  | NG  | NI  | N\* | NN  | NG  | HD  | NI  | NI  | NG  | HD  | NN  | NG  |

Now, it is more evolutionarily relevant to look at Distal
parts/domain/repeat (this is used somehow interchangeably with roughly
the same meaning) code sequences. Needless to say that to obtain TALE
parts code sequences, you need beforehand to run either the original
Distal perl script with `runDistal` or our R implementation `distalr`
with DNA sequences. In this cas we just use a toy example provided with
tantale.

``` r
aln <- buildRepeatMsa(inputSeqs = system.file("extdata", "small_Out_CodedRepeats.fa",
                                              package = "tantale", mustWork = TRUE),
                      sep = " ",
                      distalRepeatSims = NULL,
                      mafftOpts = "--localpair --maxiterate 1000 --quiet --reorder --op 0 --ep 5 --thread 1",
                      gapSymbol = "-")
```

**NOTE**, the difference in the `sep` parameter value.

The output shown as a simple matrix :

``` r
knitr::kable(aln)
```

|                                                            | 1   | 2   | 3   | 4   | 5   | 6   | 7   | 8   | 9   | 10  | 11  | 12  | 13  | 14  | 15  | 16  | 17  | 18  | 19  | 20  | 21  | 22  | 23  | 24  | 25  | 26  | 27  | 28  |
|:-----------------------------------------------------------|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|
| MAI1\|MAI1_chr1\_+\_2218220_2222519_25.5                   | 99  | 48  | 59  | 48  | 54  | 72  | 17  | 40  | 58  | 72  | 17  | 37  | 64  | 31  | 67  | 37  | 62  | 58  | 48  | 37  | 54  | 44  | 24  | 20  | 31  | 63  | 36  | 12  |
| Xoo_BAI3\|Xoo\|BAI3\|Seq1\_+\_2227619_2231921_25.5         | 100 | 48  | 37  | 48  | 54  | 72  | 19  | 41  | 16  | 54  | 62  | 58  | 81  | 31  | 43  | 37  | 44  | 37  | 48  | 75  | 32  | 44  | 21  | 37  | 31  | 63  | 36  | 10  |
| Xoo_CFBP1947\|Xoo\|CFBP1947\|Seq1\_-\_2480624_2484926_25.5 | 100 | 48  | 37  | 48  | 54  | 72  | 19  | 41  | 16  | 54  | 62  | 58  | 81  | 31  | 43  | 37  | 44  | 37  | 48  | 75  | 32  | 44  | 21  | 37  | 31  | 63  | 36  | 10  |

## A detailled illustration of one of the TALE msa plotting functions

First we get the necessary info from object from a previous run:

``` r
distalrOut <- readRDS(file = system.file("extdata", "sampleDistalrOutput.rds",
                                              package = "tantale", mustWork = TRUE))
repeatMsaByGroup <- readRDS(file = system.file("extdata", "sampleRepeatMsaByGroup.rds",
                                              package = "tantale", mustWork = TRUE))

repeatAlign <- repeatMsaByGroup[[4]]
rvdAlign <- convertRepeat2RvdAlign(repeatAlign = repeatAlign,
                                   repeat2RvdMapping = getRepeat2RvdMappingFromDistalr(distalrOut$taleParts))
```

This is what you get if you provide all the arguments:

``` r
p1 <- ggplotTalesMsa(repeatAlign = repeatAlign,
               talsim = distalrOut$tal.similarity,
               repeatSim = distalrOut$repeat.similarity,
               rvdAlign = rvdAlign,
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatSim" # "repeatClust"
)
#>> Registered S3 methods overwritten by 'treeio':
#>>   method              from    
#>>   MRCA.phylo          tidytree
#>>   MRCA.treedata       tidytree
#>>   Nnode.treedata      tidytree
#>>   Ntip.treedata       tidytree
#>>   ancestor.phylo      tidytree
#>>   ancestor.treedata   tidytree
#>>   child.phylo         tidytree
#>>   child.treedata      tidytree
#>>   full_join.phylo     tidytree
#>>   full_join.treedata  tidytree
#>>   groupClade.phylo    tidytree
#>>   groupClade.treedata tidytree
#>>   groupOTU.phylo      tidytree
#>>   groupOTU.treedata   tidytree
#>>   inner_join.phylo    tidytree
#>>   inner_join.treedata tidytree
#>>   is.rooted.treedata  tidytree
#>>   nodeid.phylo        tidytree
#>>   nodeid.treedata     tidytree
#>>   nodelab.phylo       tidytree
#>>   nodelab.treedata    tidytree
#>>   offspring.phylo     tidytree
#>>   offspring.treedata  tidytree
#>>   parent.phylo        tidytree
#>>   parent.treedata     tidytree
#>>   root.treedata       tidytree
#>>   rootnode.phylo      tidytree
#>>   sibling.phylo       tidytree
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa-1.png)

The returned object can be modified :

``` r
class(p1)
#>> [1] "aplot"
str(p1$plotlist, max.level = 1)
#>> List of 2
#>>  $ :List of 9
#>>   ..- attr(*, "class")= chr [1:2] "gg" "ggplot"
#>>  $ :List of 9
#>>   ..- attr(*, "class")= chr [1:3] "ggtree" "gg" "ggplot"
```

Plotting the original alignment

``` r
p1$plotlist[[1]]
```

![](p2_multiple_alignments_files/figure-html/The%20alignment%20part%20of%20the%20aplot-1.png)

Now we can modify the color scales and replot to check this is what we
want.

``` r
cat("modifying the original alignment as a ggplot")
#>> modifying the original alignment as a ggplot
p2 <- p1$plotlist[[1]] + ggplot2::scale_fill_viridis_c() + ggplot2::scale_color_discrete()
#>> Scale for fill is already present.
#>> Adding another scale for fill, which will replace the existing scale.
#>> Scale for colour is already present.
#>> Adding another scale for colour, which will replace the existing scale.
p2
```

![](p2_multiple_alignments_files/figure-html/The%20modified%20alignment%20part%20of%20the%20aplot-1.png)

We can extract the tree part of the object.

``` r
p1$plotlist[[2]]
```

![](p2_multiple_alignments_files/figure-html/The%20tree%20part%20of%20the%20aplot-1.png)

And we can recombine all again :

``` r
p1$plotlist[[1]] <- p2
p1
```

![](p2_multiple_alignments_files/figure-html/reassembling%20align%20and%20tree-1.png)

Pay attention to the arguments provided to the plotting function. In the
various calls to `ggplotTalesMsa`, arguments are gradually removed to
exemplify the type of alignments that are thus displayed.

``` r
p <- ggplotTalesMsa(repeatAlign = repeatAlign,
               talsim = distalrOut$tal.similarity,
               repeatSim = distalrOut$repeat.similarity,
               rvdAlign = rvdAlign,
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatClust" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-1.png)

``` r

p <- ggplotTalesMsa(repeatAlign = repeatAlign,
               talsim = distalrOut$tal.similarity,
               repeatSim = distalrOut$repeat.similarity,
               rvdAlign = rvdAlign,
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatSim" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-2.png)

``` r

p <- ggplotTalesMsa(repeatAlign = repeatAlign,
               talsim = distalrOut$tal.similarity,
               repeatSim = distalrOut$repeat.similarity,
               rvdAlign = NULL,
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatClust" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-3.png)

``` r

# fillType has no effect
p <- ggplotTalesMsa(repeatAlign = repeatAlign,
               talsim = distalrOut$tal.similarity,
               repeatSim = NULL,
               rvdAlign = rvdAlign,
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatSim" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-4.png)

``` r


p <- ggplotTalesMsa(repeatAlign = NULL,
               talsim = distalrOut$tal.similarity,
               repeatSim = NULL,
               rvdAlign = rvdAlign,
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatSim" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-5.png)

``` r


p <- ggplotTalesMsa(repeatAlign = repeatAlign,
               talsim = distalrOut$tal.similarity,
               repeatSim = NULL,
               rvdAlign = NULL,
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatSim" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-6.png)

``` r

p <- ggplotTalesMsa(repeatAlign = NULL,
               talsim = NULL,
               repeatSim = NULL,
               rvdAlign = rvdAlign,
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatSim" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-7.png)

``` r


p <- ggplotTalesMsa(repeatAlign = repeatAlign,
               talsim = NULL,
               repeatSim = NULL,
               rvdAlign = NULL,
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatSim" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-8.png)

``` r




# Single sequence align
p <- ggplotTalesMsa(repeatAlign = repeatAlign[3, , drop = FALSE],
               talsim = distalrOut$tal.similarity,
               repeatSim = distalrOut$repeat.similarity,
               rvdAlign = rvdAlign[3, , drop = FALSE],
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatSim" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-9.png)

``` r

p <- ggplotTalesMsa(repeatAlign = repeatAlign[3, , drop = FALSE],
               talsim = NULL,
               repeatSim = distalrOut$repeat.similarity,
               rvdAlign = rvdAlign[3, , drop = FALSE],
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatClust" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-10.png)

``` r


p <- ggplotTalesMsa(repeatAlign = repeatAlign[3, , drop = FALSE],
               talsim = NULL,
               repeatSim = distalrOut$repeat.similarity,
               rvdAlign = NULL,
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatSim" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-11.png)

``` r

p <- ggplotTalesMsa(repeatAlign = repeatAlign[3, , drop = FALSE],
               talsim = NULL,
               repeatSim = NULL,
               rvdAlign = NULL,
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatSim" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-12.png)

``` r

p <- ggplotTalesMsa(repeatAlign = NULL,
               talsim = NULL,
               repeatSim = NULL,
               rvdAlign = rvdAlign[3, , drop = FALSE],
               repeat.clust.h.cut = 90,
               refgrep = NULL,
               consensusSeq = FALSE,
               fillType = "repeatSim" # "repeatClust"
)
```

![](p2_multiple_alignments_files/figure-html/detailled%20illustration%20of%20ggplotTalesMsa%20others-13.png)

## A couple of utility functions for TALE msa

Compute the consensus of the columns in a TALE alignment:

``` r
taleAlignConsensus(repeatAlign)
#>>  [1] "186" "62"  "48"  "133" "127" "155" "73"  "149" "157" "152" "51"  "118"
#>> [13] "95"  "178" "94"  "57"  "64"  "177" "115" "26"  "26"  "37"  "50"  "64" 
#>> [25] "94"  "26"  "49"  "245"
taleAlignConsensus(rvdAlign)
#>>  [1] "NTERM" "NN"    "ND"    "NN"    "NI"    "NK"    "NN"    "HD"    "NN"   
#>> [10] "NG"    "NG"    "N*"    "HD"    "N*"    "HD"    "NI"    "NN"    "HD"   
#>> [19] "NG"    "HD"    "HD"    "HD"    "NG"    "NN"    "HD"    "HD"    "NG"   
#>> [28] "CTERM"
```

Figure out what elements in the TALE alignment corresponds to the
consensus:

``` r
matchConsensus(repeatAlign)
#>> # A tibble: 84 × 3
#>>    arrayID            positionInArray matchConsensus
#>>    <fct>                        <int> <chr>         
#>>  1 BAI3-1-1_ROI_00003               1 TRUE          
#>>  2 BAI3_ROI_00004                   1 FALSE         
#>>  3 MAI1_ROI_00004                   1 FALSE         
#>>  4 BAI3-1-1_ROI_00003               2 TRUE          
#>>  5 BAI3_ROI_00004                   2 TRUE          
#>>  6 MAI1_ROI_00004                   2 FALSE         
#>>  7 BAI3-1-1_ROI_00003               3 TRUE          
#>>  8 BAI3_ROI_00004                   3 TRUE          
#>>  9 MAI1_ROI_00004                   3 TRUE          
#>> 10 BAI3-1-1_ROI_00003               4 TRUE          
#>> # ℹ 74 more rows
matchConsensus(repeatAlign, returnLong = FALSE)
#>>                    1       2       3      4       5       6       7     
#>> BAI3-1-1_ROI_00003 "TRUE"  "TRUE"  "TRUE" "TRUE"  "TRUE"  "TRUE"  "TRUE"
#>> BAI3_ROI_00004     "FALSE" "TRUE"  "TRUE" "TRUE"  "TRUE"  "TRUE"  "TRUE"
#>> MAI1_ROI_00004     "FALSE" "FALSE" "TRUE" "FALSE" "FALSE" "FALSE" "TRUE"
#>>                    8       9       10      11     12      13      14     
#>> BAI3-1-1_ROI_00003 "TRUE"  "TRUE"  "TRUE"  "TRUE" "TRUE"  "TRUE"  "TRUE" 
#>> BAI3_ROI_00004     "TRUE"  "TRUE"  "TRUE"  "TRUE" "TRUE"  "TRUE"  "TRUE" 
#>> MAI1_ROI_00004     "FALSE" "FALSE" "FALSE" "TRUE" "FALSE" "FALSE" "FALSE"
#>>                    15      16      17     18     19     20     21     22     
#>> BAI3-1-1_ROI_00003 "TRUE"  "TRUE"  "TRUE" "TRUE" "TRUE" "TRUE" "TRUE" "TRUE" 
#>> BAI3_ROI_00004     "TRUE"  "TRUE"  "TRUE" "TRUE" "TRUE" "TRUE" "TRUE" "TRUE" 
#>> MAI1_ROI_00004     "FALSE" "FALSE" "TRUE" "TRUE" "TRUE" "TRUE" "TRUE" "FALSE"
#>>                    23      24     25      26     27     28     
#>> BAI3-1-1_ROI_00003 "TRUE"  "TRUE" "TRUE"  "TRUE" "TRUE" "TRUE" 
#>> BAI3_ROI_00004     "TRUE"  "TRUE" "TRUE"  "TRUE" "TRUE" "FALSE"
#>> MAI1_ROI_00004     "FALSE" "TRUE" "FALSE" "TRUE" "TRUE" "FALSE"
```
