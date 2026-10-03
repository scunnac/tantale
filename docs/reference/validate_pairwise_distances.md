# Validate a pairwise distance table

Checks the column contract, and nothing else. Squareness, the diagonal
and symmetry are preconditions of the methods that need them, checked by
[`distances_assert_square`](https://scunnac.github.io/tantale/reference/distances_assert_square.md):
a table filtered on one id column is a legitimate intermediate step.

## Usage

``` r
validate_pairwise_distances(x)
```

## Arguments

- x:

  A `pairwise_distances` object.

## Value

`x`, invisibly, if valid; otherwise an error.

## See also

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
[`tales_group_kmedoids()`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md)

## Examples

``` r
d <- data.frame(
  id1 = c("A1", "A1", "A2", "A2"),
  id2 = c("A1", "A2", "A1", "A2"),
  dissim = c(0, 35, 35, 0)
)
validate_pairwise_distances(pairwise_distances(d))

# Without its dissim column, a table is not a pairwise_distances
try(validate_pairwise_distances(d[, c("id1", "id2")]))
#> Error in validate_pairwise_distances(d[, c("id1", "id2")]) : 
#>   A <pairwise_distances> object requires the column dissim.
```
