# Plotting a multiple alignment of TALE sequences

Plot in frame of
[`heatmap.2`](https://rdrr.io/pkg/gplots/man/heatmap.2.html) for Tals
alignment.

## Usage

``` r
heatmap_msa(
  talsim,
  repeatAlign,
  rvdAlign = NULL,
  repeatSim,
  repeat.clust.h.cut = 90,
  refgrep = NULL,
  consensusSeq = FALSE,
  noteColSet = NULL,
  plot.type,
  save.path,
  ...
)
```

## Arguments

- talsim:

  a *three columns Tals similarity table* as obtained with
  [`runDistal`](https://scunnac.github.io/tantale/reference/runDistal.md)
  in the 'tal.similarity' slot of the returned object.

- repeatAlign:

  a multiple Tal repeat sequences alignment in the form of a matrix as
  returned by
  [`buildRepeatMsa`](https://scunnac.github.io/tantale/reference/buildRepeatMsa.md)
  or as one of the elements of the `SeqOfRepsAlignments` slot in the
  return object of the
  [`buildDisTalGroups`](https://scunnac.github.io/tantale/reference/buildDisTalGroups.md)
  function.

- rvdAlign:

  (optional) when the rvds need to be labeled in the plot (plot.type =
  "repeat.similarity" or "repeat.clusters.with.rvd", a multiple Tal
  repeat sequences alignment in the form of a matrix as returned by
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
  function. **(CORRECT???!!!)**

- repeat.clust.h.cut:

  height for tree cutting when plot type in "repeat.clusters".

- refgrep:

  regular expression pattern that will be used to search Tal names to
  select the reference in the alignment.

- consensusSeq:

  (logical) whether to display the consensus sequence when the plot type
  is "repeat.clusters.with.rvd".

- noteColSet:

  In case rvdSim = NULL, vector of 2 colors for rvd alignment, the first
  color is for matched rvds, and the second color is for mismatched
  ones. In the other case, more colors should be supplied.

- plot.type:

  Either `"repeat.similarity"`, `"repeat.clusters"` ,
  `"repeat.clusters.with.rvd"`. Defines the type of plot that will be
  produced by the function. See below for details.

- save.path:

  file path to save the plot. If save.path is NULL, the heatmap will be
  printed. If save.path is specified, the image file will be created
  with the format based on file extension.

- ...:

  any other arguments of
  [`heatmap.2`](https://rdrr.io/pkg/gplots/man/heatmap.2.html)

## Value

the return value of
[`heatmap.2`](https://rdrr.io/pkg/gplots/man/heatmap.2.html)

## Details

"repeat.similarity" plot shows RVD alignment of Tals, a hierarchical
dendrogram reflecting overall similarities between TALEs, similarity
between repeats alignment by the color of cells, and (if rvdSim is
provided) similarity between RVDs alignment in the color of the RVD
labels.

"repeat.clusters" plot shows repeat alignment of Tals with cells filled
with colors representing the repeat clustering group and a hierarchical
dendrogram reflecting overall similarities between TALEs.

"repeat.clusters.with.rvd" plots repeat alignment of Tals with rvd
labeled. If plotting from the outputs of `buildDistalGroups`, you supply
repeatClustID/Similarity alignment to `forMatrix`, with `talsim`,
`forCellNote` - repeat/rvd alignment, and refgrep optionally.

But if you don't have these alignments, you can provide **repeat
alignment** to `forMatrix` with `repeatSim`, the repeat similarity data
frame, the function will convert it into repeatClustID/Similarity
alignment depending on the plot type. In case of *repeat.clusters*, you
may want to adjust the param `repeat.clust.h.cut` to decrease/increase
the number of repeat clusters.
