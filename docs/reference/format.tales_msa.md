# Render a tales_msa object as lines of text

As
[`format.tales`](https://scunnac.github.io/tantale/reference/format.tales.md),
plus the alignment width, and a preview of aligned rows: gaps are drawn
and every cell padded to a common width, so the columns line up down the
page.

## Usage

``` r
# S3 method for class 'tales_msa'
format(x, n = 2L, gap = "-", ...)
```

## Arguments

- x:

  A
  [`tales_msa`](https://scunnac.github.io/tantale/reference/tales_msa.md)
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
[`[.tales()`](https://scunnac.github.io/tantale/reference/sub-.tales.md),
[`as_tales()`](https://scunnac.github.io/tantale/reference/as_tales.md),
[`format.tales()`](https://scunnac.github.io/tantale/reference/format.tales.md),
[`is_tales()`](https://scunnac.github.io/tantale/reference/is_tales.md),
[`new_tales()`](https://scunnac.github.io/tantale/reference/new_tales.md),
[`print.tales()`](https://scunnac.github.io/tantale/reference/print.tales.md),
[`print.tales_msa()`](https://scunnac.github.io/tantale/reference/print.tales_msa.md),
[`summary.tales()`](https://scunnac.github.io/tantale/reference/summary.tales.md),
[`summary.tales_msa()`](https://scunnac.github.io/tantale/reference/summary.tales_msa.md),
[`tales()`](https://scunnac.github.io/tantale/reference/tales.md),
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md),
[`tales_anomalies()`](https://scunnac.github.io/tantale/reference/tales_anomalies.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/reference/tales_assert_complete.md),
[`tales_bind()`](https://scunnac.github.io/tantale/reference/tales_bind.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/reference/validate_tales.md)

## Examples

``` r
aligned <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L),
  alignment_position = c(1L, 2L, 3L, 1L, 3L),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "CTERM")
)
msa <- tales_msa(aligned)
format(msa)
#> [1] "\033[34m<tales_msa>\033[39m 2 arrays, 3 alignment positions"
#> [2] "  layers: rvd"                                              
#> [3] "      rvd"                                                  
#> [4] "  A1  NTERM    HD CTERM"                                    
#> [5] "  A2  NTERM     - CTERM"                                    
```
