# Compute TALE and domain relatedness by domain-sequence alignment (DisTAL)

Quantifies how TALE arrays, and the individual domains they are built
from – repeats and the two termini alike – relate to one another by
aligning their domain sequences. An R re-implementation of the original
DisTAL Perl program, including the array alignment step (the ARLEM
model, see
[`tales_tale_distances`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md)),
with parallel domain alignment, which makes it much faster (the exact
speedup depends on `aln_method`).

## Usage

``` r
tales_compare_distal(
  x,
  ncores = 1,
  aln_method = "DECIPHER",
  conda_bin = "auto"
)
```

## Arguments

- x:

  A [`tales`](https://scunnac.github.io/tantale/reference/tales.md)
  object whose parts carry amino acid sequences.

- ncores:

  Number of cores for the pairwise alignment step.

- aln_method:

  Approach for pairwise distances between part amino acid sequences:
  `"DECIPHER"` (default), `"Biostrings"` or `"mmseq2"`. See
  [`tales_domain_distances`](https://scunnac.github.io/tantale/reference/tales_domain_distances.md).

- conda_bin:

  Path to a Conda binary, if `reticulate` cannot find it.

## Value

A list of three objects, all describing the same run:

- `tales`: the input, with a `dom_code` column added.

- `domain_distances`: a
  [`domain_distances`](https://scunnac.github.io/tantale/reference/pairwise_distances.md)
  between distinct domains, keyed by `dom_code`.

- `tale_distances`: a
  [`tale_distances`](https://scunnac.github.io/tantale/reference/pairwise_distances.md)
  between whole arrays, keyed by `array_id`.

## Details

The name says which comparison this is:
[`tales_compare_functal`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md)
compares TALEs by their predicted DNA-binding specificity instead.

It returns the two expensive products (the pairwise protein alignment
between distinct domains, and the alignment of the coded arrays) with
the coded input; anything else is cheap to derive from these.

This is where `dom_code` is minted, over the whole set of parts
supplied, and where `tales`/`domain_distances` are stamped with a
namespace identifying that set — see
[`tales_namespace`](https://scunnac.github.io/tantale/reference/tales_namespace.md).
Passing a subset later is safe; re-running on a different part set mints
different codes, and the differing namespace is what stops the two being
joined by mistake. `tale_distances` carries no such stamp: it is keyed
by `array_id`, a name that keeps its meaning across runs.

## References

Pérez-Quintero A.L. et al. (2015). QueTAL: a suite of tools to classify
and compare TAL effectors functionally and phylogenetically. *Frontiers
in Plant Science* **6**, 545.
[doi:10.3389/fpls.2015.00545](https://doi.org/10.3389/fpls.2015.00545)

Abouelhoda M.I., Giegerich R., Behzadi B., Steyaert J.-M. (2009).
Alignment of minisatellite maps based on run-length encoding scheme.
*Journal of Bioinformatics and Computational Biology* **7**(2), 287–308.
[doi:10.1142/S0219720009004060](https://doi.org/10.1142/S0219720009004060)

## See also

[`tales_group_hclust`](https://scunnac.github.io/tantale/reference/tales_group_hclust.md)/[`tales_group_kmedoids`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md)
to cluster arrays from the returned `tale_distances`;
[`tales_compare_functal`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md),
comparing TALEs by predicted DNA-binding specificity instead of
domain-sequence relatedness.

Other pairwise distances:
[`[.pairwise_distances()`](https://scunnac.github.io/tantale/reference/sub-.pairwise_distances.md),
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/reference/is_pairwise_distances.md),
[`pairwise_distances()`](https://scunnac.github.io/tantale/reference/pairwise_distances.md),
[`rvd_dna_specificity`](https://scunnac.github.io/tantale/reference/rvd_dna_specificity.md),
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_assign_domain_codes.md),
[`tales_compare_functal()`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md),
[`tales_domain_distances()`](https://scunnac.github.io/tantale/reference/tales_domain_distances.md),
[`tales_group_hclust()`](https://scunnac.github.io/tantale/reference/tales_group_hclust.md),
[`tales_group_kmedoids()`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/reference/validate_pairwise_distances.md)

## Examples

``` r
x <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
                                      package = "tantale"))
cmp <- tales_compare_distal(x)
#> Computing a distance matrix between TALE parts amino acid sequences using:
#> DECIPHER
#> Deriving domain substitution costs that meet the triangle inequality (Minkowski
#> distance between domain distance profiles).
#> Aligning 4 TALE arrays pairwise (6 pairs).
#> Finished computing TALE and repeat relatedness.
names(cmp)
#> [1] "tales"            "domain_distances" "tale_distances"  
cmp$tale_distances
#> # A tibble: 16 × 5
#>    id1       id2       dissim arlem_score max_length
#>    <chr>     <chr>      <dbl>       <dbl>      <int>
#>  1 ROI_00001 ROI_00001   0              0         28
#>  2 ROI_00002 ROI_00001   5.64         158         28
#>  3 ROI_00003 ROI_00001   4.96         139         28
#>  4 ROI_00004 ROI_00001   5.29         148         28
#>  5 ROI_00001 ROI_00002   5.64         158         28
#>  6 ROI_00002 ROI_00002   0              0         16
#>  7 ROI_00003 ROI_00002   6            168         28
#>  8 ROI_00004 ROI_00002   5.62         135         24
#>  9 ROI_00001 ROI_00003   4.96         139         28
#> 10 ROI_00002 ROI_00003   6            168         28
#> 11 ROI_00003 ROI_00003   0              0         28
#> 12 ROI_00004 ROI_00003   5.43         152         28
#> 13 ROI_00001 ROI_00004   5.29         148         28
#> 14 ROI_00002 ROI_00004   5.62         135         24
#> 15 ROI_00003 ROI_00004   5.43         152         28
#> 16 ROI_00004 ROI_00004   0              0         24
```
