# Codes marking a TALE array terminus

The values a `rvd` column takes on non-repeat parts. `"NTERM"` and
`"CTERM"` mark identified termini; `"XXXXX"` marks a terminus whose CDS
was detected but for which no HMMer hit was found, so its identity is
unknown (see
[`tell_tales`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)).

## Usage

``` r
tales_anchor_codes()
```

## Value

A character vector of anchor codes.

## Details

These share the `rvd` column with real RVDs, so code that distinguishes
repeats from termini by value should use this function rather than
spelling the codes out.

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
[`tales()`](https://scunnac.github.io/tantale/dev/reference/tales.md),
[`tales_anomalies()`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/dev/reference/tales_assert_complete.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/dev/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/dev/reference/validate_tales.md)

## Examples

``` r
tales_anchor_codes()
#> [1] "NTERM" "CTERM" "XXXXX"
```
