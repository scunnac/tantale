# Plot the domain composition of a set of TALE arrays

A compact, information-rich view of the arrays in a `tales` object: one
point per part, positioned by its place in the array, coloured by domain
type and filled by amino-acid length, with the RVD printed on each
repeat.

A
[`tales_msa`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md)
dispatches to
[`plot.tales_msa`](https://scunnac.github.io/tantale/dev/reference/plot.tales_msa.md)
instead, being the more specific class.

## Usage

``` r
# S3 method for class 'tales'
plot(x, position = c("array", "alignment"), ...)
```

## Arguments

- x:

  A [`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
  object, as returned by
  [`tales_from_telltale`](https://scunnac.github.io/tantale/dev/reference/tales_from_telltale.md)
  or in the `tales` element of
  [`tales_compare`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)'s
  output. A legacy `tale_parts` data frame is accepted and converted.

- position:

  Which coordinate to lay the parts out on. `"array"` (default) uses
  `position_in_array`, so each array starts at 1 and runs contiguously.
  `"alignment"` uses `alignment_position`, which requires an aligned
  object (or one demoted from a
  [`tales_msa`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md),
  which keeps the column): gaps then appear as empty columns and shared
  features line up. Aberrant repeats, for instance, are visible as a
  column in the aligned layout and scattered in the unaligned one.

- ...:

  Unused, present for compatibility with the `plot` generic.

## Value

The ggplot object, invisibly printed as a side effect.

## See also

Other TALE plots:
[`plot.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/plot.tales_msa.md),
[`plot_target_preds()`](https://scunnac.github.io/tantale/dev/reference/plot_target_preds.md),
[`talomes_heatmap()`](https://scunnac.github.io/tantale/dev/reference/talomes_heatmap.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
plot(x)
```
