# 'Nice' plotting a multiple alignment of TALE sequences

Plot TALEs msa in the ggplot2 framework.

## Usage

``` r
ggplotTalesMsa(
  repeatAlign,
  talsim = NULL,
  rvdAlign = NULL,
  repeatSim = NULL,
  repeat.clust.h.cut = 90,
  refgrep = NULL,
  consensusSeq = FALSE,
  fillType = "repeatClust"
)
```

## Arguments

- repeatAlign:

  A multiple Tal repeat sequences alignment in the form of a matrix as
  returned by
  [`buildRepeatMsa`](https://scunnac.github.io/tantale/reference/buildRepeatMsa.md)
  or as one of the elements of the `SeqOfRepsAlignments` slot in the
  return object of the
  [`buildDisTalGroups`](https://scunnac.github.io/tantale/reference/buildDisTalGroups.md)
  function.

- talsim:

  a *three columns Tals similarity table* as obtained with
  [`runDistal`](https://scunnac.github.io/tantale/reference/runDistal.md)
  in the 'tal.similarity' slot of the returned object.

- rvdAlign:

  A multiple Tal RVD sequences alignment in the form of a matrix as
  returned by
  [`convertRepeat2RvdAlign`](https://scunnac.github.io/tantale/reference/convertRepeat2RvdAlign.md),
  [`buildRepeatMsa`](https://scunnac.github.io/tantale/reference/buildRepeatMsa.md)
  or as one of the elements of the `SeqOfRepsAlignments` slot in the
  return object of the
  [`buildDisTalGroups`](https://scunnac.github.io/tantale/reference/buildDisTalGroups.md)
  function.

- repeatSim:

  A long, three columns data frame with pairwise similarity scores
  between repeats as available in the `repeat.similarity slot` of the
  object returned by the
  [`runDistal`](https://scunnac.github.io/tantale/reference/runDistal.md)
  or the the
  [`distalr`](https://scunnac.github.io/tantale/reference/distalr.md)
  functions.

- repeat.clust.h.cut:

  height for tree cutting when defining domain/repeat clusters.

- refgrep:

  Regular expression pattern that will be used to search TALE names to
  select the reference in the alignment.

- consensusSeq:

  (logical) Whether to display the consensus sequence \*\*NOT
  IMPLEMENTED YET\*\*

- fillType:

  Either "repeatClust" or "repeatSim". If both options are possible
  because the necessary information is there (at least a `repeatSim`
  value), this argument will decide what type of 'box color filling' is
  employed and it is either based on the cluster where the repeat falls
  after clustering all the repeat in the alignment or it is based on the
  amino acid similarity between a repeat at a position and the repeat of
  the 'reference' TALE at this position.

## Value

An [`aplot`](https://rdrr.io/pkg/aplot/man/plot-insertion.html) object.

## Details

This function as a similar purpose as
[`heatmap_msa`](https://scunnac.github.io/tantale/reference/heatmap_msa.md)
but has been implemented with
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html). It is
more versatile (takes single row matrices of alignment) but a bit
slower.

The type of plot that you will get will depend on the provided
information in the form of parameter values See the tantale website for
detailed usage cases.

The only mandatory argument is either `repeatAlign` **or** `rvdAlign`.

The plot is printed and returned for further modifications is necessary.
