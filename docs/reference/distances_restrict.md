# Restrict a distance table to a set of entities

Filters **both** id columns, which is what keeps the result square.

## Usage

``` r
distances_restrict(x, ids)
```

## Arguments

- x:

  A `pairwise_distances` object.

- ids:

  A character vector of entity ids to keep.

## Value

A `pairwise_distances` over `ids` only.

## See also

Other pairwise distances:
[`[.pairwise_distances()`](https://scunnac.github.io/tantale/reference/sub-.pairwise_distances.md),
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/reference/distances_assert_square.md),
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
  id1 = c("A1", "A1", "A1", "A2", "A2", "A2", "A3", "A3", "A3"),
  id2 = c("A1", "A2", "A3", "A1", "A2", "A3", "A1", "A2", "A3"),
  dissim = c(0, 35, 60, 35, 0, 40, 60, 40, 0)
))
distances_restrict(d, c("A1", "A2"))
#> # A tibble: 4 × 3
#>   id1   id2   dissim
#>   <chr> <chr>  <dbl>
#> 1 A1    A1         0
#> 2 A1    A2        35
#> 3 A2    A1        35
#> 4 A2    A2         0
```
