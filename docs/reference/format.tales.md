# Render a tales object as lines of text

The representation
[`print.tales`](https://scunnac.github.io/tantale/reference/print.tales.md)
emits, returned as a character vector, so an object's description can go
somewhere other than the console: a log, an error message, a report.

It shows what the class knows and a tibble would not: the number of
arrays as distinct from the number of parts, which residue layers the
object carries, the `dom_code` namespace when it is stamped, and a
preview of the first and last arrays as sequences.

The preview uses `dom_code` when present and `rvd` otherwise. Domain
codes are the more discriminating of the two: two arrays can share an
RVD sequence while being built from different repeats.

This replaces the [`format()`](https://rdrr.io/r/base/format.html) a
`tales` would otherwise inherit from tibble. Call
`format(tibble::as_tibble(x))` for that one.

## Usage

``` r
# S3 method for class 'tales'
format(x, n = 2L, ...)
```

## Arguments

- x:

  A [`tales`](https://scunnac.github.io/tantale/reference/tales.md)
  object.

- n:

  Number of arrays to show at each end. Everything is shown when the
  object holds `2 * n` arrays or fewer.

- ...:

  Unused, present for compatibility with the `format` generic.

## Value

A character vector, one element per line.

## See also

Other tales objects:
[`[.tales()`](https://scunnac.github.io/tantale/reference/sub-.tales.md),
[`as_tales()`](https://scunnac.github.io/tantale/reference/as_tales.md),
[`format.tales_msa()`](https://scunnac.github.io/tantale/reference/format.tales_msa.md),
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
[`tales_names()`](https://scunnac.github.io/tantale/reference/tales_names.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/reference/validate_tales.md)

## Examples

``` r
parts <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L, 3L),
  domain_type = c("N-terminus", "repeat", "C-terminus",
                  "N-terminus", "repeat", "C-terminus"),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "NI", "CTERM")
)
x <- tales(parts)
format(x)
#> [1] "\033[34m<tales>\033[39m 2 arrays, 6 parts"
#> [2] "  layers: rvd   |   1 other column"       
#> [3] "      rvd"                                
#> [4] "  A1  NTERM HD CTERM"                     
#> [5] "  A2  NTERM NI CTERM"                     
```
