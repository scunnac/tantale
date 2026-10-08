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
profile (`hmmsearch`). `"NTERM"` and `"CTERM"` mark a canonical TALE
terminal domain, one that can be expected to do its usual job: the match
has an E-value of at most `terminus_max_evalue`, covers at least
`terminus_min_cover` of the profile (0.9 by default) and reaches the end
of the profile that adjoins the repeats. `"XXXXX"` marks any other
segment. An array for which AnnoTALE reported no segment on one side has
no terminus part at the protein level on that side.

A terminus can fail to be canonical in several ways, which
`array_report.tsv` tells apart (its `*_aa_*` and `*_dna_*` columns, see
[`tell_tales`](https://scunnac.github.io/tantale/reference/tell_tales.md)):
the segment may be unrelated sequence where the ORF starts or ends
inside a frameshifted region, its repeat-side part may be in another
reading frame, or the terminus may have lost a large part of its domain.
The truncTALEs of *Xanthomonas oryzae* have lost the activation domain,
and their C-termini are coded `"XXXXX"`. A smaller internal deletion,
such as the one in the N-terminus of TalC, a major TALE of African *X.
oryzae* pv. *oryzae*, leaves a terminus canonical.

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
[`tales_names()`](https://scunnac.github.io/tantale/reference/tales_names.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/reference/validate_tales.md)

## Examples

``` r
tales_anchor_codes()
#>      N-      -C      ?? 
#> "NTERM" "CTERM" "XXXXX" 
```
