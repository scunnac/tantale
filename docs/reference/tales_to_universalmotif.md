# Convert a tales object to a list of predicted DNA-binding-specificity PWMs

One `universalmotif` position weight matrix per array, built by looking
up each repeat's RVD in
[`rvd_dna_specificity`](https://scunnac.github.io/tantale/reference/rvd_dna_specificity.md)
and stacking the rows in repeat order. The conversion
[`tales_compare_functal`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md)
is built on, exposed on its own so any `universalmotif` function –
`compare_motifs()`, `motif_tree()`, `view_motifs()`, `scan_sequences()`,
`merge_motifs()`, ... – can be run directly on real TALE binding models,
not just the one comparison
[`tales_compare_functal()`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md)
does with them. Ledger §12b lists several of these as concrete
follow-ups.

## Usage

``` r
tales_to_universalmotif(x)
```

## Arguments

- x:

  A [`tales`](https://scunnac.github.io/tantale/reference/tales.md)
  object carrying an `rvd` column.

## Value

A named list of `universalmotif` objects, one per array, named by
`array_id` – a plain list, since that is what
`compare_motifs()`/`motif_tree()` themselves accept.

## Details

Only `rvd`, in repeat order, is used (via
[`tales_rvd_strings`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md),
which drops the two termini by default – DNA-binding specificity is a
property of the repeat region; QueTAL FuncTAL, the tool this table comes
from, never scored termini either, since its own RVD extraction only
ever found repeats). An RVD absent from
[`rvd_dna_specificity`](https://scunnac.github.io/tantale/reference/rvd_dna_specificity.md)
is scored with a flat, uninformative row rather than dropped, so an
unusual RVD costs a comparison specificity rather than an error.

## See also

[`tales_compare_functal()`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md),
the one comparison built on this;
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md),
the sibling projection at the string layer.

Other tales projections:
[`repeat_to_rvd_map()`](https://scunnac.github.io/tantale/reference/repeat_to_rvd_map.md),
[`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/reference/repeat_to_rvd_map_distalr.md),
[`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/reference/tale_parts_to_rvd.md),
[`tales_coded_strings()`](https://scunnac.github.io/tantale/reference/tales_coded_strings.md),
[`tales_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_domain_codes.md),
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
motifs <- tales_to_universalmotif(x)
motifs[[1]]
#> 
#>        Motif name:   ROI_00001
#>          Alphabet:   ACGT
#>              Type:   PCM
#>          Total IC:   15.98
#>       Pseudocount:   0
#>      Target sites:   72
#> 
#>   [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#> A   30    5   30   15   15   50    1    5   15    50     5    30    15    50
#> C   10   10   10   50   50   10    4   10   50    10    10    10    50    10
#> G   30    5   30    5    5    5    1    5    5     5     5    30     5     5
#> T    2   50    2    5    5    5    4   50    5     5    50     2     5     5
#>   [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
#> A     5    50     5    30     5    15    50    50     5    15    30     5
#> C    10    10    10    10    10    50    10    10    10    50    10    10
#> G     5     5     5    30     5     5     5     5     5     5    30     5
#> T    50     5    50     2    50     5     5     5    50     5     2    50
```
