# Domain-coded strings, one per TALE array

Renders each array as a separated string of its `dom_code`s in part
order, where each distinct domain sequence is one "residue".
[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md)
and
[`tales_tale_distances()`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md)
work on arrays in this form.

## Usage

``` r
tales_coded_strings(x, sep = " ", repeats_only = FALSE)
```

## Arguments

- x:

  A [tales](https://scunnac.github.io/tantale/reference/tales.md) object
  carrying a `dom_code` column.

- sep:

  Separator between codes. Defaults to `" "`; see Details.

- repeats_only:

  Drop the terminus parts, keeping only repeats. `FALSE` by default; see
  Details. Needs a `domain_type` column.

## Value

A
[Biostrings::BStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html),
named by `array_id`.

## Details

This is the sibling of
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md)
and takes the same two arguments, but **both defaults differ**, because
the two projections feed different consumers:

|        |                                                                                           |                                      |
|--------|-------------------------------------------------------------------------------------------|--------------------------------------|
|        | [`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md) | `tales_coded_strings()`              |
| `sep`  | `"-"`, the AnnoTALE convention                                                            | `" "`, what MAFFT `--text` splits on |
| filter | `repeats_only = TRUE`                                                                     | `repeats_only = FALSE`               |

The separator is free to choose here in a way it is not for RVDs: a
`dom_code` is a bare integer rendered as text, so `"1 2 3"` and
`"1-2-3"` are equally unambiguous. It defaults to a space because that
is what MAFFT's text mode expects.

The termini are kept by default, where
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md)
drops them. Target prediction concerns the repeat domain only, so
dropping them there is right; alignment is the consumer here, and the
two termini are the most reliable anchors an alignment of TALE arrays
has. Set `repeats_only = TRUE` to compare bare repeat arrays.

## See also

[`tales_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_domain_codes.md)
for the code-to-sequence lookup,
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md)
for the sibling projection.

Other tales projections:
[`tales_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_domain_codes.md),
[`tales_get_dna_seq()`](https://scunnac.github.io/tantale/reference/tales_get_dna_seq.md),
[`tales_get_protein_seq()`](https://scunnac.github.io/tantale/reference/tales_get_protein_seq.md),
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md),
[`tales_to_universalmotif()`](https://scunnac.github.io/tantale/reference/tales_to_universalmotif.md)

## Examples

``` r
x <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
                                      package = "tantale"))
xa <- tales_assign_domain_codes(x)
tales_coded_strings(xa)[1]
#> BStringSet object of length 1:
#>     width seq                                               names               
#> [1]    77 43 19 12 19 24 32 2 14 ...33 8 16 3 12 7 29 11 47 ROI_00001
tales_coded_strings(xa, sep = "-", repeats_only = TRUE)[1]
#> BStringSet object of length 1:
#>     width seq                                               names               
#> [1]    71 19-12-19-24-32-2-14-1-2...19-33-8-16-3-12-7-29-11 ROI_00001
```
