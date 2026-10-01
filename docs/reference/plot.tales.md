# Plot the domain composition of a set of TALE arrays

A compact, information-rich view of the arrays in a `tales` object: one
point per part, positioned by its place in the array, outlined by domain
type and filled by amino-acid length, with the RVD printed on each
repeat.

A
[`tales_msa`](https://scunnac.github.io/tantale/reference/tales_msa.md)
dispatches to
[`plot.tales_msa`](https://scunnac.github.io/tantale/reference/plot.tales_msa.md)
instead, being the more specific class.

## Usage

``` r
# S3 method for class 'tales'
plot(x, position = c("array", "alignment"), facet_by = "seqnames", ...)
```

## Arguments

- x:

  A [`tales`](https://scunnac.github.io/tantale/reference/tales.md)
  object, as returned by
  [`tales_from_telltales`](https://scunnac.github.io/tantale/reference/tales_from_telltales.md)
  or in the `tales` element of
  [`tales_compare_distal`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md)'s
  output. A data frame of parts is accepted and converted with
  [`tales`](https://scunnac.github.io/tantale/reference/tales.md).

- position:

  Which coordinate to lay the parts out on. `"array"` (default) uses
  `position_in_array`, so each array starts at 1 and runs contiguously.
  `"alignment"` uses `alignment_position`, which requires an aligned
  object (or one demoted from a
  [`tales_msa`](https://scunnac.github.io/tantale/reference/tales_msa.md),
  which keeps the column): gaps then appear as empty columns and shared
  features line up. Aberrant repeats, for instance, are visible as a
  column in the aligned layout and scattered in the unaligned one.

- facet_by:

  Names of the columns whose values split the plot into panels, stacked
  in rows whose height follows the number of arrays. Each column must
  hold one value per array: `"seqnames"` (the default, one panel per
  source sequence) or `"strain"` for a set of genomes, for instance, and
  `c("strain", "seqnames")` for both. `NULL` draws a single panel, as
  does the default when `x` has no `seqnames` column.

- ...:

  Unused, present for compatibility with the `plot` generic.

## Value

The ggplot object, returned invisibly after being printed as a side
effect.

## See also

Other TALE plots:
[`plot.tales_msa()`](https://scunnac.github.io/tantale/reference/plot.tales_msa.md),
[`plot_target_preds()`](https://scunnac.github.io/tantale/reference/plot_target_preds.md),
[`talomes_heatmap()`](https://scunnac.github.io/tantale/reference/talomes_heatmap.md)

## Examples

``` r
x <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
                                      package = "tantale"))
plot(x)

plot(x, facet_by = NULL)
```
