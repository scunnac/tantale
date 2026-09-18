# Report the biological anomalies in a tales object

Lists the arrays that are *odd* rather than *unreadable*: missing
sequence data, impossible domain-type arrangements, coordinate
disagreements, an amino acid sequence paired with more than one RVD, or
an attribute that varies within an array when it should not.

Such arrays are accepted by
[`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md) –
real TALE predictions are messy, and refusing to load them would force
cleaning outside the package and destroy the diagnostic signal.
Construction warns about them; this function tells you which and why;
`tales(x, sanitize = TRUE)` removes them.

## Usage

``` r
tales_anomalies(x)
```

## Arguments

- x:

  A [`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
  object.

## Value

A tibble of `array_id`, `check` and `detail`, one row per anomaly. Zero
rows if the object is clean.

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
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_anchor_codes.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/dev/reference/tales_assert_complete.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/dev/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/dev/reference/validate_tales.md)

## Examples

``` r
# A2 carries two N-termini -- an impossible arrangement.
odd <- data.frame(
  array_id = c("A1", "A1", "A2", "A2", "A2"),
  position_in_array = c(1L, 2L, 1L, 2L, 3L),
  domain_type = c("N-terminus", "repeat",
                  "N-terminus", "N-terminus", "repeat"),
  rvd = c("NTERM", "HD", "NTERM", "NTERM", "NI")
)
x <- suppressWarnings(tales(odd))
tales_anomalies(x)
#> # A tibble: 2 × 3
#>   array_id check               detail                      
#>   <chr>    <chr>               <chr>                       
#> 1 A2       terminus_duplicated more than one N-terminus    
#> 2 A2       terminus_misplaced  N-terminus not at position 1
tales(odd, sanitize = TRUE) # drops A2 instead of merely warning
#> Warning: Dropped 1 array with biological anomalies.
#> ✖ Array: "A2"
#> ℹ Reasons: terminus_duplicated and terminus_misplaced
#> <tales> 1 array, 2 parts
#>   layers: rvd   |   1 other column
#>       rvd
#>   A1  NTERM HD
```
