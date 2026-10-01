# Whole-array DNA sequence, one per TALE array

The DNA sibling of
[`tales_get_protein_seq()`](https://scunnac.github.io/tantale/reference/tales_get_protein_seq.md):
reassembles each array's parts, in order, into one full-length
nucleotide sequence, by pasting `dna_seq` together directly. See
[`tales_get_protein_seq()`](https://scunnac.github.io/tantale/reference/tales_get_protein_seq.md)'s
Details for why `position_in_array` is the ordering key even for a
`tales_msa`.

## Usage

``` r
tales_get_dna_seq(x)
```

## Arguments

- x:

  A [tales](https://scunnac.github.io/tantale/reference/tales.md) object
  carrying a `dna_seq` column.

## Value

A
[Biostrings::DNAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html),
named by `array_id`.

## See also

[`tales_get_protein_seq()`](https://scunnac.github.io/tantale/reference/tales_get_protein_seq.md)
for the protein sequence sibling.

Other tales projections:
[`tales_coded_strings()`](https://scunnac.github.io/tantale/reference/tales_coded_strings.md),
[`tales_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_domain_codes.md),
[`tales_get_protein_seq()`](https://scunnac.github.io/tantale/reference/tales_get_protein_seq.md),
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md),
[`tales_to_universalmotif()`](https://scunnac.github.io/tantale/reference/tales_to_universalmotif.md)

## Examples

``` r
x <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
                                      package = "tantale"))
tales_get_dna_seq(x)[1]
#> DNAStringSet object of length 1:
#>     width seq                                               names               
#> [1]  4305 ATGGATCCCATTCGTCCGCGCAG...TGAGGGAGCTATTGCCTCAGTGA ROI_00001
```
