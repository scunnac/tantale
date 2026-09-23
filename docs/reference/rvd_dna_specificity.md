# RVD-to-DNA-binding-specificity weights

One row per RVD, giving its relative binding weight for each of the four
bases. This is the table
[`tales_compare_functal`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md)
stacks, one row per repeat, to build an array's position weight matrix.

## Usage

``` r
rvd_dna_specificity
```

## Format

A tibble with 404 rows and 5 columns:

- rvd:

  The two-letter RVD code. `"N*"` and `"H*"` are RVDs whose residue 13
  is missing; `"OO"` is position 0, the base just before the first
  repeat's target, where TALEs prefer a T; `"XX"` is the flat fallback
  row used for an unrecognised RVD.

- A, C, G, T:

  Relative binding weight for that base. The weights do not sum to one:
  [`tales_compare_functal`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md)
  hands them to
  [`create_motif`](https://rdrr.io/pkg/universalmotif/man/create_motif.html)
  as raw counts (`type = "PCM"`), which normalises internally.

## Details

Taken verbatim from QueTAL FuncTAL's own table (shipped as
`legacy/QueTAL_v1.1/FuncTAL/Info/2014mat18` in the installed package):
the values are unchanged, only a header and column names added.
Conceptually related to but distinct from the internal `rvdSimDf` used
by
[`tales_align`](https://scunnac.github.io/tantale/reference/tales_align.md)'s
RVD scoring: that one is a *derived* RVD-vs-RVD similarity (a
correlation between two RVDs' base-preference profiles, itself built
from TALVEZ's much smaller 17-RVD `mat1`), used to score repeat
*substitutions* during sequence alignment. This one is the *raw*,
per-RVD base preference itself, one level upstream, used to build a
whole array's binding-specificity model for comparison against another
array's. The two tables serve different consumers and cannot replace
each other.

## References

Pérez-Quintero A.L. et al. (2015). QueTAL: a suite of tools to classify
and compare TAL effectors functionally and phylogenetically. *Frontiers
in Plant Science* **6**, 545.
[doi:10.3389/fpls.2015.00545](https://doi.org/10.3389/fpls.2015.00545)

## See also

[`tales_compare_functal`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md),
its only consumer.

Other pairwise distances:
[`[.pairwise_distances()`](https://scunnac.github.io/tantale/reference/sub-.pairwise_distances.md),
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/reference/is_pairwise_distances.md),
[`pairwise_distances()`](https://scunnac.github.io/tantale/reference/pairwise_distances.md),
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_assign_domain_codes.md),
[`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md),
[`tales_compare_functal()`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md),
[`tales_domain_distances()`](https://scunnac.github.io/tantale/reference/tales_domain_distances.md),
[`tales_group_hclust()`](https://scunnac.github.io/tantale/reference/tales_group_hclust.md),
[`tales_group_kmedoids()`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/reference/validate_pairwise_distances.md)

## Examples

``` r
rvd_dna_specificity[rvd_dna_specificity$rvd %in% c("HD", "NI", "NG", "NN"), ]
#> # A tibble: 4 × 5
#>   rvd       A     C     G     T
#>   <chr> <dbl> <dbl> <dbl> <dbl>
#> 1 HD       15    50     5     5
#> 2 NG        5    10     5    50
#> 3 NI       50    10     5     5
#> 4 NN       30    10    30     2
```
