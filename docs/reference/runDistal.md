# Run the Distal tool

Run the Distal tool of the
[QueTal](https://doi.org/10.3389/fpls.2015.00545) suite to classify and
compare TAL effectors functionally and phylogenetically.

## Usage

``` r
runDistal(
  fasta.file,
  outdir = NULL,
  treetype = "p",
  repeats.cluster.h.cut = 10,
  overwrite = F
)
```

## Arguments

- fasta.file:

  fasta file containing TALE DNA/AA sequences.

- outdir:

  directory to store disTal output.

- repeats.cluster.h.cut:

  numeric value to cut the hierachycal clustering tree of the repeat.

- overwrite:

  logical indicating whether to rerun disTal or only load the existing
  results.

## Value

A list with DisTal output components:

- repeats.code: a data frame of the unique repeat AA sequences and there
  numeric codes

- coded.repeats.str: a list of repeat-coded TALE strings

- repeat.similarity: a long, three columns data frame with pairwise
  similarity scores between repeats

- tal.similarity: a three columns Tals similarity table with pairwise
  similarity scores between TALEs

- tree: Newick format neighbor-joining tree of TALs constructed based on
  TALEs similarity

- repeats.cluster: a data frame containing repeat code and repeat
  clusters.
