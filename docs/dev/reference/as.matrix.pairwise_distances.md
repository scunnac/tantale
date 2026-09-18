# Render a similarity table as a square matrix

Materialises the wide form, with ids as both row and column names,
sorted. This replaces the hand-written
`acast(x, id1 ~ id2, value.var = "sim")` that appears at four call sites
in the package.

## Usage

``` r
# S3 method for class 'pairwise_distances'
as.matrix(x, value = PAIRWISE_DISTANCES_VALUE_COL, ...)
```

## Arguments

- x:

  A `pairwise_distances` object.

- value:

  Name of the column to fill cells with. Defaults to `dissim`.

- ...:

  Ignored.

## Value

A numeric matrix, square, with sorted ids as dimnames.

## See also

Other pairwise distances:
[`distances_assert_square()`](https://scunnac.github.io/tantale/dev/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/dev/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/is_pairwise_distances.md),
[`pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md),
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md),
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md),
[`tales_domain_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_distances.md),
[`tales_group()`](https://scunnac.github.io/tantale/dev/reference/tales_group.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_tale_distances.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/validate_pairwise_distances.md)

## Examples

``` r
d <- pairwise_distances(data.frame(
  id1 = c("A1", "A1", "A2", "A2"),
  id2 = c("A1", "A2", "A1", "A2"),
  dissim = c(0, 35, 35, 0)
))
as.matrix(d)
#>    A1 A2
#> A1  0 35
#> A2 35  0
```
