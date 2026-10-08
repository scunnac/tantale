# Group TALEs by hierarchical clustering of their pairwise distance

Clusters TALE arrays with
[`hclust`](https://rdrr.io/r/stats/hclust.html) on their pairwise
distance, cuts the tree into exactly `k` groups, and returns the
[`tales`](https://scunnac.github.io/tantale/reference/tales.md) object
the distances were computed from with the result attached as its `group`
column.

## Usage

``` r
tales_group_hclust(x, tale_distances, k = NULL, plot_tree = FALSE)
```

## Arguments

- x:

  A [tales](https://scunnac.github.io/tantale/reference/tales.md) object
  – the one whose comparison produced `tale_distances`.

- tale_distances:

  A
  [tale_distances](https://scunnac.github.io/tantale/reference/pairwise_distances.md)
  object, as returned by
  [`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md).
  A plain data frame with the same columns is also accepted and coerced.

- k:

  Integer, the number of groups to cut the tree into.

- plot_tree:

  Logical, whether to draw the `ggtree` dendrogram: colored by group,
  with a dashed line at the cut height. `FALSE` by default.

## Value

`x` with an added (or replaced) `group` column.

## Details

The clustering is computed from `tale_distances`, but the result belongs
on the `tales` object the distances were computed from, so that is what
comes back. `group` is a recognised `tales` column, validated as
constant within an array – it is an array-level property, like
`seqnames`.

Taking `x` rather than returning a bare lookup table is what makes the
correspondence checkable: the array names in `tale_distances` must be
the array names in `x`, and this is the only place that can be verified.
A mismatch is an error rather than a silent `NA` group, because a
partly-grouped object is the kind of thing that fails much later and
confusingly.

The bare mapping is still one line away if you want it:
`unique(out[c("array_id", "group")])`.

The tree is built directly on the distances
(`stats::hclust(stats::as.dist(distMat))`), matching the original DisTAL
clustering this package reimplements. The tree is cut with
`stats::cutree(tree, k = k)`, which always succeeds, including when a
tie in merge heights would make a height-based cut ambiguous.

## See also

[`tales_group_kmedoids()`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md),
the alternative method;
[`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md),
which produces both inputs.

Other pairwise distances:
[`[.pairwise_distances()`](https://scunnac.github.io/tantale/reference/sub-.pairwise_distances.md),
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/reference/is_pairwise_distances.md),
[`pairwise_distances()`](https://scunnac.github.io/tantale/reference/pairwise_distances.md),
[`rvd_dna_specificity`](https://scunnac.github.io/tantale/reference/rvd_dna_specificity.md),
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_assign_domain_codes.md),
[`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md),
[`tales_compare_functal()`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md),
[`tales_domain_distances()`](https://scunnac.github.io/tantale/reference/tales_domain_distances.md),
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
grouped <- tales_group_hclust(cmp$tales, cmp$tale_distances, k = 2)
unique(grouped[c("array_id", "group")])
#> # A tibble: 4 × 2
#>   array_id  group
#>   <chr>     <int>
#> 1 ROI_00001     1
#> 2 ROI_00002     2
#> 3 ROI_00003     1
#> 4 ROI_00004     1
```
