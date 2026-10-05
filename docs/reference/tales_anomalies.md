# Report the biological anomalies in a tales object

Lists the arrays whose content is biologically odd: missing sequence
data, a structure other than that of a standard TALE (an N-terminus, one
or more repeats and a C-terminus, both termini canonical TALE terminal
domains), impossible domain-type arrangements, coordinate disagreements,
an amino acid sequence paired with more than one RVD, or an attribute
that varies within an array when it should not. Structurally broken
input (a duplicated key, a missing required column) is an error in
[`tales`](https://scunnac.github.io/tantale/reference/tales.md) instead.

Such arrays are accepted by
[`tales`](https://scunnac.github.io/tantale/reference/tales.md) – real
TALE predictions are messy, and refusing to load them would force
cleaning outside the package and destroy the diagnostic signal.
Construction warns about them and this function tells you which and why.

Each anomaly is of one of two kinds. `"noncanonical"` marks a TALE whose
structure departs from the standard one: a terminus that is not a
canonical TALE terminal domain (`terminus_noncanonical`), such as the
C-terminus of a truncTALE, or no terminus part on one side
(`terminus_absent`). Such TALEs can be real and interesting.
`"integrity"` marks every other anomaly: data that are inconsistent or
incomplete, or an array without repeats. `tales(x, sanitize = TRUE)`
removes the arrays with an anomaly of kind `"integrity"` and keeps the
others; `tales(x, sanitize = "canonical")` keeps canonical TALEs only.

The structure checks need a `domain_type` column, and the terminus
profile check an `rvd` column; they are skipped when it is absent. The
structure checks report `terminus_absent` (no N- or no C-terminus part),
`no_repeat` (no repeat part) and `terminus_noncanonical` (a terminus
coded `XXXXX`, see
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

A tibble of `array_id`, `check`, `kind` (`"integrity"` or
`"noncanonical"`) and `detail`, one row per anomaly. Zero rows if the
object is clean.

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
#> # A tibble: 4 × 4
#>   array_id check               kind         detail                      
#>   <chr>    <chr>               <chr>        <chr>                       
#> 1 A1       terminus_absent     noncanonical no C-terminus               
#> 2 A2       terminus_absent     noncanonical no C-terminus               
#> 3 A2       terminus_duplicated integrity    more than one N-terminus    
#> 4 A2       terminus_misplaced  integrity    N-terminus not at position 1
tales(odd, sanitize = TRUE) # drops A2 instead of merely warning
#> Warning: Dropped 1 array with inconsistent or incomplete data.
#> ✖ Array: "A2"
#> ℹ Reasons: terminus_duplicated and terminus_misplaced
#> Kept 1 non-canonical array.
#> ℹ Array: "A1"
#> ℹ `tales_anomalies()` lists it, with the reason.
#> <tales> 1 array, 2 parts
#>   layers: rvd   |   1 other column
#>       rvd
#>   A1  NTERM HD

# B2 is a truncTALE-like array: its C-terminus is not canonical.
parts <- data.frame(
  array_id = rep(c("B1", "B2"), each = 3),
  position_in_array = rep(1:3, 2),
  domain_type = rep(c("N-terminus", "repeat", "C-terminus"), 2),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "NI", "XXXXX")
)
y <- tales(parts, sanitize = TRUE)          # keeps B2, with a message
#> Kept 1 non-canonical array.
#> ℹ Array: "B2"
#> ℹ `tales_anomalies()` lists it, with the reason.
tales_anomalies(y)
#> # A tibble: 1 × 4
#>   array_id check                 kind         detail                            
#>   <chr>    <chr>                 <chr>        <chr>                             
#> 1 B2       terminus_noncanonical noncanonical C-terminus is not a canonical TAL…
z <- tales(parts, sanitize = "canonical")   # drops B2
#> Warning: Dropped 1 array that is not a canonical TALE.
#> ✖ Array: "B2"
#> ℹ Reason: terminus_noncanonical
unique(z$array_id)
#> [1] "B1"
```
