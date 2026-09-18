# Render a tales object as lines of text

The representation
[`print.tales`](https://scunnac.github.io/tantale/dev/reference/print.tales.md)
emits, returned as a character vector rather than written out – so an
object's description can go somewhere other than the console: a log, an
error message, a report.

It shows what the class knows and a tibble would not: the number of
arrays as distinct from the number of parts, which residue layers the
object carries, the `dom_code` namespace when it is stamped, and a
preview of the first and last arrays as sequences.

The preview uses `dom_code` when present and `rvd` otherwise. Repeat
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

  A [`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
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
[`as_tales()`](https://scunnac.github.io/tantale/dev/reference/as_tales.md),
[`format.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/format.tales_msa.md),
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
