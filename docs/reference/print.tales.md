# Print a tales object

Emits what
[`format.tales`](https://scunnac.github.io/tantale/reference/format.tales.md)
builds; see there for what is shown and why.

## Usage

``` r
# S3 method for class 'tales'
print(x, n = 2L, ...)
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

`x`, invisibly.

## See also

Other tales objects:
[`[.tales()`](https://scunnac.github.io/tantale/reference/sub-.tales.md),
[`as_tales()`](https://scunnac.github.io/tantale/reference/as_tales.md),
[`format.tales()`](https://scunnac.github.io/tantale/reference/format.tales.md),
[`format.tales_msa()`](https://scunnac.github.io/tantale/reference/format.tales_msa.md),
[`is_tales()`](https://scunnac.github.io/tantale/reference/is_tales.md),
[`new_tales()`](https://scunnac.github.io/tantale/reference/new_tales.md),
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
parts <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L, 3L),
  domain_type = c("N-terminus", "repeat", "C-terminus",
                  "N-terminus", "repeat", "C-terminus"),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "NI", "CTERM")
)
x <- tales(parts)
print(x)
#> <tales> 2 arrays, 6 parts
#>   layers: rvd   |   1 other column
#>       rvd
#>   A1  NTERM HD CTERM
#>   A2  NTERM NI CTERM
```
