# Subset a pairwise_distances object

Subsets rows and columns like an ordinary tibble, but the class travels
only while the result still satisfies the contract: `id1`/`id2` (both
character) and `dissim` (numeric) all present. Dropping any of them
leaves something that can no longer be described as a pairwise
similarity table, and the class quietly steps out of the way rather than
continuing to claim invariants it can no longer keep – the result is a
plain tibble, not an error.

## Usage

``` r
# S3 method for class 'pairwise_distances'
x[...]
```

## Arguments

- x:

  A `pairwise_distances` object.

- ...:

  Passed on to the tibble/data frame method.

## Value

The subset: still a `pairwise_distances` (or its
`tale_distances`/`domain_distances` subclass) if the contract holds,
otherwise a plain tibble.

## Details

Squareness is never checked here, on purpose (see
[`distances_assert_square`](https://scunnac.github.io/tantale/reference/distances_assert_square.md)):
subsetting rows is a normal, everyday way to end up with a non-square
table, and re-checking that on every `[` call would outlaw it.

## See also

Other pairwise distances:
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
[`tales_group_kmedoids()`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/reference/validate_pairwise_distances.md)

## Examples

``` r
d <- pairwise_distances(data.frame(
  id1 = c("A1", "A1", "A2", "A2"),
  id2 = c("A1", "A2", "A1", "A2"),
  dissim = c(0, 35, 35, 0)
))
is_pairwise_distances(d[1:2, ]) # row subsetting never breaks the contract
#> [1] TRUE
class(d[, "id1"])               # dropping dissim/id2: the class steps aside
#> [1] "tbl_df"     "tbl"        "data.frame"
```
