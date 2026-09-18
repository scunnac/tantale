# Create a tales object

Builds a `tales` object from a data frame of TALE array parts: one row
per (array, slot). Column names still using the legacy camelCase
spelling (`arrayID`, `positionInArray`, ...) are renamed on the way in.

## Usage

``` r
tales(x, dom_code_namespace = NULL, sanitize = FALSE)
```

## Arguments

- x:

  A data frame with at least `array_id` and `position_in_array` columns,
  plus at least one of `rvd` or `dom_code`. Other recognised columns
  (`domain_type`, `position_in_crd`, `aa_seq`, `dna_seq`, `seqnames`,
  `source_directory`) are validated if present. Any further column is
  preserved untouched.

- dom_code_namespace:

  Optional scalar string, see
  [`tales_namespace`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md).

- sanitize:

  If `TRUE`, arrays carrying biological anomalies – missing sequences,
  impossible terminus arrangements, coordinate disagreements – are
  removed, with a warning naming them and why. If `FALSE` (default) they
  are kept and merely warned about, so odd predictions can still be
  loaded and inspected. Structural corruption is an error either way.
  See
  [`tales_anomalies`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md).

## Value

A validated `tales` object.

## See also

Other tales objects:
[`as_tales()`](https://scunnac.github.io/tantale/dev/reference/as_tales.md),
[`format.tales()`](https://scunnac.github.io/tantale/dev/reference/format.tales.md),
[`format.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/format.tales_msa.md),
[`is_tales()`](https://scunnac.github.io/tantale/dev/reference/is_tales.md),
[`new_tales()`](https://scunnac.github.io/tantale/dev/reference/new_tales.md),
[`print.tales()`](https://scunnac.github.io/tantale/dev/reference/print.tales.md),
[`print.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/print.tales_msa.md),
[`summary.tales()`](https://scunnac.github.io/tantale/dev/reference/summary.tales.md),
[`summary.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/summary.tales_msa.md),
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_anchor_codes.md),
[`tales_anomalies()`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/dev/reference/tales_assert_complete.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/dev/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/dev/reference/validate_tales.md)

## Examples

``` r
parts <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L, 3L),
  domain_type = c("N-terminus", "repeat", "C-terminus",
                  "N-terminus", "repeat", "C-terminus"),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "NI", "CTERM")
)
tales(parts)
#> <tales> 2 arrays, 6 parts
#>   layers: rvd   |   1 other column
#>       rvd
#>   A1  NTERM HD CTERM
#>   A2  NTERM NI CTERM
```
