# Group TALEs by similarity

Classifies TALE arrays into groups by hierarchical or k-medoids
clustering of their pairwise similarity, and returns the
[tales](https://scunnac.github.io/tantale/dev/reference/tales.md) object
with the result attached as its `group` column.

## Usage

``` r
tales_group(
  x,
  tal_sim,
  plot_tree = FALSE,
  k = NULL,
  k_range = NULL,
  method = "k-medoids"
)
```

## Arguments

- x:

  A [tales](https://scunnac.github.io/tantale/dev/reference/tales.md)
  object – the one whose comparison produced `tal_sim`.

- tal_sim:

  A
  [tale_distances](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md)
  object, as returned by
  [`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md).
  A plain data frame using the legacy `TAL1`/`TAL2`/`Sim` column names
  is also accepted and coerced.

- plot_tree:

  Logical, whether to plot the hclust tree. With `method = "k-medoids"`
  no tree is drawn; a silhouette-value plot is produced instead.

- k:

  Integer, the number of groups to classify arrays into. With
  `method = "k-medoids"`, `k = "auto"` picks the optimum automatically
  and `k = NULL` prompts for it interactively. The automatic pick is
  worth checking rather than trusting.

- k_range:

  Integer vector of length 2 giving the range of `k` to test. Only used
  when `method = "k-medoids"`; the minimum is 2.

- method:

  One of `"hclust"` (see
  [`cutree`](https://rdrr.io/r/stats/cutree.html)) or `"k-medoids"` (see
  [`pam`](https://rdrr.io/pkg/cluster/man/pam.html)).

## Value

`x` with an added (or replaced) `group` column.

## Details

The clustering is computed from `tal_sim`, but the result belongs on the
`tales` object the distances were computed from, so that is what comes
back. `group` is a recognised `tales` column, validated as constant
within an array – it is an array-level property, like `seqnames`.

Taking `x` rather than returning a bare lookup table is what makes the
correspondence checkable: the array names in `tal_sim` must be the array
names in `x`, and this is the only place that can be verified. A
mismatch is an error rather than a silent `NA` group, because a
partly-grouped object is the kind of thing that fails much later and
confusingly.

The bare mapping is still one line away if you want it:
`unique(out[c("array_id", "group")])`.

## See also

[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md),
which produces both inputs.

Other pairwise distances:
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/dev/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/dev/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/is_pairwise_distances.md),
[`pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md),
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md),
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md),
[`tales_domain_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_distances.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_tale_distances.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/validate_pairwise_distances.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
cmp <- tales_compare(x)
#> Computing a distance matrix between TALE parts amino acid sequences using:
#> DECIPHER
#> Generate an ARLEM cost matrix which meets triangle inequality criteria by
#> computing the minkowski distance between pairwise distance vectors.
#> Running ARLEM version 1.0 :
#> Copyright by Mohamed I. Abouelhoda
#> Plz. cite Abouelhoda, Giegerich, Behzadi, and Steyaert
#> Finished computing TALE and repeat relatedness.
grouped <- tales_group(cmp$tales, cmp$tale_distances, method = "hclust", k = 2)
unique(grouped[c("array_id", "group")])
#> # A tibble: 4 × 2
#>   array_id  group
#>   <chr>     <int>
#> 1 ROI_00001     1
#> 2 ROI_00002     2
#> 3 ROI_00003     2
#> 4 ROI_00004     1
```
