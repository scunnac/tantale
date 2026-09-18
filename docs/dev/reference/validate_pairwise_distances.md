# Validate a pairwise similarity table

Checks the column contract, and nothing else. Squareness, the diagonal
and symmetry are deliberately *not* checked here: the package filters
these tables asymmetrically on purpose, so those properties are
preconditions of the methods that need them — see
[`distances_assert_square`](https://scunnac.github.io/tantale/dev/reference/distances_assert_square.md).

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
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/dev/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/dev/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/is_pairwise_distances.md),
[`pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md),
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md),
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md),
[`tales_domain_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_distances.md),
[`tales_group()`](https://scunnac.github.io/tantale/dev/reference/tales_group.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_tale_distances.md)
