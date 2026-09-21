# Perform multiple alignment of TALE repeat or RVD sequences

Perform multiple alignment of TALE repeat or RVD sequences using the
`--text` mode of the [MAFFT](https://mafft.cbrc.jp/alignment/software/)
Multiple alignment program.

By default, uses the simple scoring matrix defined in [the text mode of
MAFFT](https://mafft.cbrc.jp/alignment/software/textcomparison.html).
Users can optionally provide a custom scoring matrix.

## Usage

``` r
buildRepeatMsa(
  inputSeqs,
  sep = " ",
  distalRepeatSims = NULL,
  mafftOpts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
  mafftPath = system.file("tools", "mafft-linux64", package = "tantale", mustWork = TRUE),
  gapSymbol = NA
)
```

## Arguments

- inputSeqs:

  Any object accepted as input by the
  [`toListOfSplitedStr`](https://scunnac.github.io/tantale/reference/toListOfSplitedStr.md)
  function, such as the path to a fasta file containing the TALE
  sequences to be aligned or the `coded.repeats.str` slot of the object
  returned by the
  [`runDistal`](https://scunnac.github.io/tantale/reference/runDistal.md)
  function. Can also be the return value of the
  [`taleParts2RvdStringSet`](https://scunnac.github.io/tantale/reference/taleParts2RvdStringSet.md)
  function if one wants to align RVD sequences.

- sep:

  Passed to
  [`toListOfSplitedStr()`](https://scunnac.github.io/tantale/reference/toListOfSplitedStr.md)
  to split the TALEs strings in input.

- distalRepeatSims:

  A long, three columns data frame with pairwise similarity scores
  between repeats as available in the `repeat.similarity slot` of the
  object returned by the
  [`runDistal`](https://scunnac.github.io/tantale/reference/runDistal.md)
  function.

- mafftOpts:

  A character string containing additional options for the MAFFT
  command. This is notably useful to tweak the Gap opening and gap
  extension penalties.

- mafftPath:

  Path to a MAFFT installation directory. By default uses the MAFFT
  version included in tantale.

- gapSymbol:

  Specify a alternative symbol for gaps in the alignments.

## Value

A character matrix representing the multiple alignment.
