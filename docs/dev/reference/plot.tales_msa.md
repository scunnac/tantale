# Plot a multiple alignment of TALEs

Draws the alignment as a heatmap: one row per array, one column per
alignment position, with cells coloured by one layer and optionally
labelled with another.

## Usage

``` r
# S3 method for class 'tales_msa'
plot(
  x,
  fill = NULL,
  label = NULL,
  tal_sim = NULL,
  domain_sim = NULL,
  h_cut = 10,
  ref_pattern = NULL,
  consensus = FALSE,
  fill_type = "repeat_clust",
  ...
)
```

## Arguments

- x:

  A
  [`tales_msa`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md)
  object.

- fill:

  Layer whose values colour the cells. Defaults to `"dom_code"` when
  present, otherwise the first available residue layer.

- label:

  Layer whose values are written in the cells. Left unset it defaults to
  `"rvd"` when that is not already the `fill`; pass `NULL` explicitly
  for an unlabelled heatmap.

- tal_sim:

  Pairwise distances between whole TALEs, as the `tale_distances`
  element of a
  [`tales_compare`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
  result. Used to build the tree panel that orders the alignment rows.

- domain_sim:

  Pairwise distances between repeat units, as the `domain_distances`
  element of a
  [`tales_compare`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
  result. Used to group repeats into clusters, and to score each repeat
  against the reference TALE\\s repeat at the same alignment column.

- h_cut:

  Height at which the repeat tree is cut to define clusters. Interpreted
  on a distance scale, so 0 means identical.

- ref_pattern:

  Regular expression matched against the array names to choose the
  reference TALE. Must identify exactly one.

- consensus:

  Whether to add a consensus panel above the alignment. The consensus is
  the most frequent element in each column of the labelled layer, so it
  always matches what the cells say.

- fill_type:

  One of `"repeat_clust"`, `"repeat_sim"` or `"rvd_sim"`. The first two
  colour cells by repeat cluster or by protein-sequence similarity to
  the reference. `"rvd_sim"` colours them instead by how alike each
  RVD's *DNA-binding preference* is to the reference TALE's RVD at that
  position, on a diverging scale over `[-1, 1]`. The repeat- and
  RVD-level views genuinely differ: `HD` and `ND` are distinct repeats
  with identical specificity, while repeats differing only at positions
  12-13 are near-identical proteins targeting different bases.

- ...:

  Unused, present for compatibility with the `plot` generic.

## Value

An [`aplot`](https://rdrr.io/pkg/aplot/man/plot-insertion.html) object.

## Details

A `tales_msa` carries every layer at once – `rvd`, `dom_code` and
whatever else the object holds – so `fill` and `label` name two of them
rather than being passed as separate matrices.

Three things are decided independently, and it helps to read the figure
that way: what each cell *says*, what colour that text is, and what
colour the block behind it is.

**Cell text** is whatever `label` names, or nothing when `label = NULL`.
Termini are relabelled `N-` and `-C`; an unidentified terminus keeps its
`XXXXX` code, which is deliberately not mistakable for an RVD. Repeat
codes are padded to three characters so columns line up.

**Text colour** always answers one question: does this element match the
consensus of its column? Cyan for yes, pink for no. The consensus is the
most frequent element in the column
([`tales_consensus`](https://scunnac.github.io/tantale/dev/reference/tales_consensus.md)),
taken over the labelled layer – so the text colour and the text itself
always describe the same thing.

**Block fill** is what `fill_type` selects, and it is the only part that
can be unavailable:

|                  |                                                                           |                 |
|------------------|---------------------------------------------------------------------------|-----------------|
| **fill_type**    | **shows**                                                                 | **needs**       |
| `"repeat_clust"` | which cluster the repeat falls in, cut at `h_cut`                         | `domain_sim`    |
| `"repeat_sim"`   | protein-sequence similarity to the reference, 0-100                       | `domain_sim`    |
| `"rvd_sim"`      | how alike the RVD's DNA-binding preference is to the reference's, -1 to 1 | a `label` layer |

With no `domain_sim` and no `label`, every block is flat grey: the text
still carries the consensus comparison, but there is nothing to colour
blocks by.

A cell with no value for the chosen layer keeps its text and loses its
colour. In `"rvd_sim"` that is the termini, which have no DNA-binding
preference and so no position on a specificity scale.

**The reference** matters for both similarity fills. `ref_pattern` is
matched against the array names and must identify exactly one, otherwise
the default is used with a warning; by default it is the array with the
most non-gap parts, ties broken alphabetically. The reference row is
marked with a trailing `_#`.

**Two panels may be attached.** Supplying `tal_sim` with more than one
array adds a dendrogram panel on the left; `consensus = TRUE` adds a
consensus panel on top. When either is present the return value is an
`aplot` composition rather than a single ggplot, so modify the alignment
through its `plotlist` element rather than adding layers to the result
directly.

## See also

Other TALE plots:
[`plot.tales()`](https://scunnac.github.io/tantale/dev/reference/plot.tales.md),
[`plot_target_preds()`](https://scunnac.github.io/tantale/dev/reference/plot_target_preds.md),
[`talomes_heatmap()`](https://scunnac.github.io/tantale/dev/reference/talomes_heatmap.md)

## Examples

``` r
aligned <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L),
  alignment_position = c(1L, 2L, 3L, 1L, 3L),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "CTERM")
)
msa <- tales_msa(aligned)
plot(msa)
```
