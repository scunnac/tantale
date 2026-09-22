# Generate a mapping between Distal repeat IDs and their cognate RVD.

Uses Distal repeat sequences and RVD sequences from a set of TALEs to
return the association between repeat ID and RVD.

## Usage

``` r
repeat_to_rvd_map(repeat_vecs, rvd_vecs)
```

## Arguments

- repeat_vecs:

  Expects a list of Distal repeat IDs character vectors. Each **named**
  element corresponding to a TALE.

- rvd_vecs:

  A named list of RVD vectors, one per array, parallel to `repeat_vecs`.

## Value

A two columns repeatID - RVD data frame.

## Details

Care must be taken that TALEs in the two sets of sequences have the same
name. In addition, the function tries hard to make sure that the two
sets of sequences are identical in every ways but the individual
'values' they contain. It is therefore notably important to make sure
that the sequences are consistent in whether they include N-term and
C-term domains IDs/Tags or not.

## See also

Other tales projections:
[`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/reference/repeat_to_rvd_map_distalr.md),
[`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/reference/tale_parts_to_rvd.md),
[`tales_coded_strings()`](https://scunnac.github.io/tantale/reference/tales_coded_strings.md),
[`tales_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_domain_codes.md),
[`tales_get_dna_seq()`](https://scunnac.github.io/tantale/reference/tales_get_dna_seq.md),
[`tales_get_protein_seq()`](https://scunnac.github.io/tantale/reference/tales_get_protein_seq.md),
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md),
[`tales_to_universalmotif()`](https://scunnac.github.io/tantale/reference/tales_to_universalmotif.md)

## Examples

``` r
repeat_vecs <- list(A1 = c("12", "45", "12"), A2 = c("45", "78"))
rvd_vecs <- list(A1 = c("HD", "NI", "HD"), A2 = c("NI", "NG"))
repeat_to_rvd_map(repeat_vecs, rvd_vecs)
#> # A tibble: 3 × 2
#>   repeatID RVD  
#>   <chr>    <chr>
#> 1 12       HD   
#> 2 45       NI   
#> 3 78       NG   
```
