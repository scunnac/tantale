# A comparison between \`runDistal\` and \`distalr\` implementations

e a small test set of tale cds sequences (5 = 2 groups of 2 + outlier)

``` r
library(tantale)
```

    ## 

    ## 

``` r
tellTaleSubjSeqsFile <- system.file("extdata", "bai3_sample_tal_genomic_regions.fasta", package = "tantale", mustWork = T)
tellTaleOutDir <- system.file("extdata", "tellTaleExampleOutput", package = "tantale", mustWork = T)
# 
# tantale::tellTale(subjectFile = tellTaleSubjSeqsFile,
#                           outputDir = tellTaleOutDir
# )
```

## Parse annotale output (rvd, aa, dna, domain type, array id) in a tellTale output dir

``` r
taleParts <- getTaleParts(system.file("extdata", "tellTaleExampleOutput", package = "tantale", mustWork = T))
```

### Run `distalr` to investigate TALE protein domains composition relatedness

``` r
dres <- distalr(taleParts = taleParts, repeats.cluster.h.cut = 10, ncores = 6,
                pairwiseAlnMethod = "Biostrings")

# dres2 <- distalr(taleParts = taleParts, repeats.cluster.h.cut = 10, ncores = 6,
#                 pairwiseAlnMethod = "mmseq2", condaBinPath = "/home/cunnac/bin/miniconda3/condabin/conda")

dres3 <- distalr(taleParts = taleParts, repeats.cluster.h.cut = 10, ncores = 6,
                pairwiseAlnMethod = "DECIPHER")
```
