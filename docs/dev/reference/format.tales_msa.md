# Render a tales_msa object as lines of text

As
[`format.tales`](https://scunnac.github.io/tantale/dev/reference/format.tales.md),
plus the alignment width, and a preview that shows aligned rows rather
than bare sequences: gaps are drawn and every cell padded to a common
width, so the columns line up down the page.

## Usage

``` r
# S3 method for class 'tales_msa'
format(x, n = 2L, gap = "-", ...)
```

## Arguments

- x:

  A
  [`tales_msa`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md)
  object.

- n:

  Number of arrays to show at each end.

- gap:

  What to draw where an array has no residue at a position.

- ...:

  Unused, present for compatibility with the `format` generic.

## Value

A character vector, one element per line.

## See also

Other tales objects:
[`as_tales()`](https://scunnac.github.io/tantale/dev/reference/as_tales.md),
[`format.tales()`](https://scunnac.github.io/tantale/dev/reference/format.tales.md),
[`is_tales()`](https://scunnac.github.io/tantale/dev/reference/is_tales.md),
[`new_tales()`](https://scunnac.github.io/tantale/dev/reference/new_tales.md),
[`print.tales()`](https://scunnac.github.io/tantale/dev/reference/print.tales.md),
[`print.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/print.tales_msa.md),
[`summary.tales()`](https://scunnac.github.io/tantale/dev/reference/summary.tales.md),
[`summary.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/summary.tales_msa.md),
[`tales()`](https://scunnac.github.io/tantale/dev/reference/tales.md),
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_anchor_codes.md),
[`tales_anomalies()`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/dev/reference/tales_assert_complete.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/dev/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/dev/reference/validate_tales.md)
