# Codes marking a TALE array terminus

The values a `rvd` column takes on non-repeat parts.

## Usage

``` r
tales_anchor_codes()
```

## Value

A named character vector of the three codes.

## Details

AnnoTALE reports as N-terminus whatever the ORF encodes upstream of the
first repeat, and as C-terminus whatever it encodes downstream of the
last one.
[`tell_tales`](https://scunnac.github.io/tantale/reference/tell_tales.md)
searches each of these segments with the TALE N- or C-terminal protein
profile (`hmmsearch`, E-value at most `terminus_max_evalue`). `"NTERM"`
and `"CTERM"` mark a segment that matches its profile, a canonical TALE
terminal domain, complete or truncated. `"XXXXX"` marks a segment that
does not match, typically unrelated sequence where the ORF starts or
ends inside a frameshifted region. An array for which AnnoTALE reported
no segment on one side has no terminus part on that side.

These share the `rvd` column with real RVDs, so code that distinguishes
repeats from termini by value should use this function rather than
spelling the codes out.

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
[`tales_anomalies()`](https://scunnac.github.io/tantale/reference/tales_anomalies.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/reference/tales_assert_complete.md),
[`tales_bind()`](https://scunnac.github.io/tantale/reference/tales_bind.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/reference/validate_tales.md)

## Examples

``` r
tales_anchor_codes()
#>      N-      -C      ?? 
#> "NTERM" "CTERM" "XXXXX" 
```
