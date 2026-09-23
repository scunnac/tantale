# Group TALEs by k-medoids clustering of their pairwise distance

Clusters TALE arrays with
[`pam`](https://rdrr.io/pkg/cluster/man/pam.html) for every candidate
`k` in `k_range`, then returns the
[`tales`](https://scunnac.github.io/tantale/reference/tales.md) object
the distances were computed from with the chosen clustering's result
attached as its `group` column.

## Usage

``` r
tales_group_kmedoids(
  x,
  tale_distances,
  k_range = NULL,
  k = NULL,
  seed = 7,
  plot_silhouette = TRUE
)
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
  A plain data frame using the legacy `TAL1`/`TAL2`/`Sim` column names
  is also accepted and coerced.

- k_range:

  Integer vector of candidate values of `k` to evaluate. Each must be
  between 1 and one less than the number of arrays being grouped
  ([`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html)'s own
  requirement).

- k:

  See Details.

- seed:

  Passed to [`set.seed()`](https://rdrr.io/r/base/Random.html) before
  every [`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html)
  call, so the same candidate always clusters the same way. Previously a
  hardcoded `7`; now a documented, overridable default.

- plot_silhouette:

  Logical, whether to draw the silhouette-vs-k plot.

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

Unlike a hierarchical tree, PAM has no single structure that can be cut
at an arbitrary `k` after the fact – `k` is a parameter to the
clustering itself. So one clustering is computed per candidate in
`k_range`, and `k` picks which of those to keep:

- `k` a single number uses that candidate directly.

- `k = "auto"` picks the elbow of the silhouette-vs-k curve (a partial,
  first-step application of the Kneedle algorithm – see
  `.tales_group_kmedoids_elbow()` – good enough in practice to be worth
  keeping, not a validated implementation of the full method).

- `k = NULL` (the default) shows the silhouette plot and asks for a
  number at the console – but only when
  [`interactive()`](https://rdrr.io/r/base/interactive.html) is `TRUE`.
  In a script, a test or a vignette render, `k = NULL` errors instead of
  blocking on input that will never arrive.

## See also

[`tales_group_hclust()`](https://scunnac.github.io/tantale/reference/tales_group_hclust.md),
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
[`tales_group_hclust()`](https://scunnac.github.io/tantale/reference/tales_group_hclust.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/reference/validate_pairwise_distances.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
cmp <- tales_compare_distal(x)
#> Computing a distance matrix between TALE parts amino acid sequences using:
#> DECIPHER
#> Deriving domain substitution costs that meet the triangle inequality (Minkowski
#> distance between domain distance profiles).
#> Aligning 4 TALE arrays pairwise (6 pairs).
#> Finished computing TALE and repeat relatedness.
grouped <- tales_group_kmedoids(cmp$tales, cmp$tale_distances,
                                k_range = 2:3, k = 2)

#> Number of groups is decided based on the provided value of k: 2
unique(grouped[c("array_id", "group")])
#> # A tibble: 4 × 2
#>   array_id  group
#>   <chr>     <int>
#> 1 ROI_00001     1
#> 2 ROI_00002     1
#> 3 ROI_00003     2
#> 4 ROI_00004     1
```
