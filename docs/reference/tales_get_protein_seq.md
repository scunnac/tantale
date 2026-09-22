# Whole-array protein sequence, one per TALE array

Reassembles each array's parts – N-terminus, repeats and C-terminus, in
order – into one full-length amino acid sequence, by pasting `aa_seq`
together directly (no separator, unlike
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md)'s
hyphen-joined RVDs: this is a real protein sequence, not a token
string).

## Usage

``` r
tales_get_protein_seq(x)
```

## Arguments

- x:

  A [tales](https://scunnac.github.io/tantale/reference/tales.md) object
  carrying an `aa_seq` column.

## Value

An
[Biostrings::AAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html),
named by `array_id`.

## Details

Ordered by `position_in_array`, not `alignment_position`: an alignment
never reorders an array's parts, only inserts gaps between them (see
`dev/class-design.md`), and a `tales_msa`'s gaps are never rows to begin
with – a gap is a column with no row for that array. So
`position_in_array` gives the same part order `alignment_position`
would, without needing to skip anything, and works unchanged whether `x`
is a bare `tales` or a `tales_msa`.

## See also

[`tales_get_dna_seq()`](https://scunnac.github.io/tantale/reference/tales_get_dna_seq.md)
for the DNA sequence sibling.

Other tales projections:
[`repeat_to_rvd_map()`](https://scunnac.github.io/tantale/reference/repeat_to_rvd_map.md),
[`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/reference/repeat_to_rvd_map_distalr.md),
[`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/reference/tale_parts_to_rvd.md),
[`tales_coded_strings()`](https://scunnac.github.io/tantale/reference/tales_coded_strings.md),
[`tales_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_domain_codes.md),
[`tales_get_dna_seq()`](https://scunnac.github.io/tantale/reference/tales_get_dna_seq.md),
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md),
[`tales_to_universalmotif()`](https://scunnac.github.io/tantale/reference/tales_to_universalmotif.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
tales_get_protein_seq(x)[1]
#> AAStringSet object of length 1:
#>     width seq                                               names               
#> [1]  1434 MDPIRPRRPSPAREILPGPQPDR...GAADDFPAFNEEELAWLRELLPQ ROI_00001
```
