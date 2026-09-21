# Emulate DisTal in R

This is meant to approximate the results of DisTal in R and is very
similar to
[`runDistal`](https://scunnac.github.io/tantale/reference/runDistal.md).
It still uses the Arlem binary just like DisTal but performs the rest of
the operations with R support and parallelization. Depending on the
`pairwiseAlnMethod` parameter value it is much faster than the orginial
Perl code and returns similar results. Please take a look at the
vignette for tips on how to use it properly.

## Usage

``` r
distalr(
  taleParts,
  repeats.cluster.h.cut = 10,
  ncores = 1,
  pairwiseAlnMethod = "DECIPHER",
  condaBinPath = "auto"
)
```

## Arguments

- taleParts:

  a table of TALE parts as returned by the
  [`getTaleParts`](https://scunnac.github.io/tantale/reference/getTaleParts.md)
  function.

- repeats.cluster.h.cut:

  numeric value to cut the hierarchical clustering tree of the repeat.

- pairwiseAlnMethod:

  Specify the underlying approach for computing pairwise similarities
  between TALE parts amino acid sequences. Must be "Biostrings",
  "mmseq2" or "DECIPHER"

- condaBinPath:

  Path to your Conda binary file if you need to specify a path different
  from the one that is automatically searched by the reticulate package
  functions.

## Value

A list with DisTal output components:

- taleParts: the original input tibble with a 'domCode' column
  corresponding to the unique distal 'code' or label attached to a
  unique domain sequence. Thus all parts with this sequence will have
  the same 'domCode'.

- repeats.code: a data frame of the unique repeat AA sequences and their
  numeric codes

- coded.repeats.str: a list of repeat-coded TALE strings

- repeat.similarity: a long, three columns data frame with pairwise
  similarity scores between repeats

- tal.similarity: a three columns Tals similarity table with pairwise
  similarity scores between TALEs

- repeats.cluster: a data frame containing repeat code and repeat
  clusters.
