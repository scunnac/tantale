# Create a pairwise similarity table

The long form of a square pairwise similarity over one entity set: one
row per ordered pair of entities, with a `sim` score.

## Usage

``` r
pairwise_distances(x, dom_code_namespace = NULL)

tale_distances(x, dom_code_namespace = NULL)

domain_distances(x, dom_code_namespace = NULL)
```

## Arguments

- x:

  A data frame with `id1`, `id2` and `dissim` columns. Legacy spellings
  and the similarity vocabulary are accepted and converted: a table
  carrying `Sim` or `normArlemScore` instead of a distance is folded
  into `dissim`, and those restatements are then dropped so only one
  copy of the quantity is stored. Further columns (`arlem_score`,
  `max_length`, or anything else) are preserved.

- dom_code_namespace:

  Optional scalar string identifying the run whose `dom_code` values
  this table is keyed by, see
  [`tales_namespace`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md).
  Relevant for `domain_distances()`, whose ids *are* `dom_code`s.

## Value

A validated `pairwise_distances` object.

## Details

`tale_distances()` and `domain_distances()` are the entity-specific
flavours: similarity between whole TALE arrays, and between individual
repeat units. They add no structure, only semantics — every method is
written once on the parent.

Legacy column spellings (`TAL1`/`TAL2`, `RepU1`/`RepU2`, `Sim`,
`Dissim`, `arlemScore`, ...) are renamed on the way in.

## See also

Other pairwise distances:
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/dev/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/dev/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/is_pairwise_distances.md),
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md),
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md),
[`tales_domain_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_distances.md),
[`tales_group()`](https://scunnac.github.io/tantale/dev/reference/tales_group.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_tale_distances.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/validate_pairwise_distances.md)

## Examples

``` r
d <- data.frame(
  id1 = c("A1", "A1", "A2", "A2"),
  id2 = c("A1", "A2", "A1", "A2"),
  dissim = c(0, 35, 35, 0)
)
pairwise_distances(d)

# Legacy spellings are recognised and folded in.
legacy <- data.frame(TAL1 = "A1", TAL2 = "A2", Sim = 65)
tale_distances(legacy)
```
