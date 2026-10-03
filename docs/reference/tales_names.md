# The names of the TALEs in a tales object

Each array in a `tales` object is one TALE, identified by its
`array_id`. This returns those identifiers once each, in the order the
arrays appear in `x`. They are also the names of the vector that
[`tales_rvd_strings`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md)
returns. A `tales_msa` has the same arrays as the `tales` object it was
aligned from.

## Usage

``` r
tales_names(x)
```

## Arguments

- x:

  A `tales` object.

## Value

A character vector of array identifiers.

## Details

[`names()`](https://rdrr.io/r/base/names.html) is left to its usual
meaning: a `tales` object is a data frame, so `names(x)` returns its
column names.

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
[`tales_anomalies()`](https://scunnac.github.io/tantale/reference/tales_anomalies.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/reference/tales_assert_complete.md),
[`tales_bind()`](https://scunnac.github.io/tantale/reference/tales_bind.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/reference/validate_tales.md)

## Examples

``` r
rvd_fasta <- system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",
                         package = "tantale")
x <- as_tales(rvd_fasta, sep = "-")
tales_names(x)
#>  [1] "TalA_BAI3"     "TalA_CFBP1947" "TalA_MAI1"     "TalA_MAI106"  
#>  [5] "TalA_MAI129"   "TalA_MAI134"   "TalA_MAI145"   "TalA_MAI68"   
#>  [9] "TalA_MAI73"    "TalA_MAI95"    "TalA_MAI99"   
length(tales_names(x)) # the number of TALEs
#> [1] 11
```
