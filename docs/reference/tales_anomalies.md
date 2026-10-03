# Report the biological anomalies in a tales object

Lists the arrays whose content is biologically odd: missing sequence
data, a structure other than that of a standard TALE (an N-terminus, one
or more repeats and a C-terminus, both termini matched by the profile of
their TALE domain), impossible domain-type arrangements, coordinate
disagreements, an amino acid sequence paired with more than one RVD, or
an attribute that varies within an array when it should not.
Structurally broken input (a duplicated key, a missing required column)
is an error in
[`tales`](https://scunnac.github.io/tantale/reference/tales.md) instead.

Such arrays are accepted by
[`tales`](https://scunnac.github.io/tantale/reference/tales.md) – real
TALE predictions are messy, and refusing to load them would force
cleaning outside the package and destroy the diagnostic signal.
Construction warns about them, this function tells you which and why,
and `tales(x, sanitize = TRUE)` removes them.

The structure checks need a `domain_type` column, and the terminus
profile check an `rvd` column; they are skipped when it is absent. The
structure checks report `terminus_absent` (no N- or no C-terminus part),
`no_repeat` (no repeat part) and `terminus_unmatched` (a terminus coded
`XXXXX`, see
[`tales_anchor_codes`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md)).
They apply to whole arrays: a subset keeping only the repeats is
reported as lacking its termini.

## Usage

``` r
tales_anomalies(x)
```

## Arguments

- x:

  A [`tales`](https://scunnac.github.io/tantale/reference/tales.md)
  object.

## Value

A tibble of `array_id`, `check` and `detail`, one row per anomaly. Zero
rows if the object is clean.

## See also

Other tales objects:
[`[.tales()`](https://scunnac.github.io/tantale/reference/sub-.tales.md),
[`as_tales()`](https://scunnac.github.io/tantale/reference/as_tales.md),
[`format.tales()`](https://scunnac.github.io/tantale/reference/format.tales.md),
[`format.tales_msa()`](https://scunnac.github.io/tantale/reference/format.tales_msa.md),
[`is_tales()`](https://scunnac.github.io/tantale/reference/is_tales.md),
[`new_tales()`](https://scunnac.github.io/tantale/reference/new_tales.md),
[`print.tales()`](https://scunnac.github.io/tantale/reference/print.tales.md),
[`print.tales_msa()`](https://scunnac.github.io/tantale/reference/print.tales_msa.md),
[`summary.tales()`](https://scunnac.github.io/tantale/reference/summary.tales.md),
[`summary.tales_msa()`](https://scunnac.github.io/tantale/reference/summary.tales_msa.md),
[`tales()`](https://scunnac.github.io/tantale/reference/tales.md),
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/reference/tales_assert_complete.md),
[`tales_bind()`](https://scunnac.github.io/tantale/reference/tales_bind.md),
[`tales_names()`](https://scunnac.github.io/tantale/reference/tales_names.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/reference/validate_tales.md)

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
#> # A tibble: 4 × 3
#>   array_id check               detail                      
#>   <chr>    <chr>               <chr>                       
#> 1 A2       terminus_duplicated more than one N-terminus    
#> 2 A1       terminus_absent     no C-terminus               
#> 3 A2       terminus_absent     no C-terminus               
#> 4 A2       terminus_misplaced  N-terminus not at position 1
tales(odd, sanitize = TRUE) # drops A2 instead of merely warning
#> Warning: Dropped 2 arrays with biological anomalies.
#> ✖ Arrays: "A2" and "A1"
#> ℹ Reasons: terminus_duplicated, terminus_absent, and terminus_misplaced
#> <tales> 0 arrays, 0 parts
#>   layers: rvd   |   1 other column
```
