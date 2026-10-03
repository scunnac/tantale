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
profile (`hmmsearch`, E-value at most `terminus_max_evalue`, the match
reaching the end of the profile that adjoins the repeats). `"NTERM"` and
`"CTERM"` mark a segment that matches its profile, a canonical TALE
terminal domain, complete or truncated at its far end. `"XXXXX"` marks a
segment that does not match, typically unrelated sequence where the ORF
starts or ends inside a frameshifted region, or a terminus whose
repeat-side part a frameshift has put in another reading frame. An array
for which AnnoTALE reported no segment on one side has no terminus part
on that side.

The codes record sequence relatedness only. A terminus shorter than the
canonical one is coded `"NTERM"` or `"CTERM"` as long as it matches its
profile, which it can do over its whole length: an internal deletion, or
a C-terminus that stops early, still aligns with the part of the profile
it keeps. Such a terminus has probably lost functional regions. The
N-terminal region carries the type III secretion signal and, next to the
repeats, the degenerate repeats that bind the thymine preceding the
target; the C-terminal region carries the nuclear localisation signals
and, at its far end, the transcription activation domain. The truncTALEs
of *Xanthomonas oryzae* have lost the activation domain, and their
C-termini are coded `"CTERM"`. The length of the terminus parts,
`nchar(aa_seq)`, is the quickest way to spot such arrays.

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
