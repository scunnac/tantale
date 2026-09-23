# Assert that a distance table is complete and square

Checks that the table holds every ordered pair of the ids it contains —
`n^2` rows for `n` ids. Phrased over the ids *present*, so that a
symmetric subset stays square: filtering both id columns to the same set
of entities preserves this, while filtering one of them does not.

## Usage

``` r
distances_assert_square(x, arg = "x")
```

## Arguments

- x:

  A `pairwise_distances` object.

- arg:

  Name of the argument being checked, for the error message.

## Value

`x`, invisibly.

## Details

Squareness is a **precondition** of the functions that need it. The
class itself does not require it, because filtering the two id columns
one after the other passes through a legitimate non-square intermediate.

## See also

Other pairwise distances:
[`[.pairwise_distances()`](https://scunnac.github.io/tantale/reference/sub-.pairwise_distances.md),
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/reference/as.matrix.pairwise_distances.md),
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
distances_assert_square(d)

try(distances_assert_square(d[1:3, ])) # missing the A2-A2 pair -- errors
#> Error in distances_assert_square(d[1:3, ]) : 
#>   `x` must hold every pair of the 2 ids it contains (4 rows), but has 3.
#> ℹ Filtering both id columns to the same entities keeps a table square;
#>   filtering one does not.
```
