# Create a pairwise distance table

The long form of a square pairwise distance over one entity set: one row
per ordered pair of entities, with a `dissim` value (0 for identical
entities, larger for more different ones).

## Usage

``` r
pairwise_distances(x, dom_code_namespace = NULL)

tale_distances(x, dom_code_namespace = NULL)

domain_distances(x, dom_code_namespace = NULL)
```

## Arguments

- x:

  A data frame with `id1`, `id2` and `dissim` columns. A similarity is
  accepted and converted: a table carrying `sim` or `norm_arlem_score`
  instead of a distance is folded into `dissim`, and those restatements
  are then dropped so only one copy of the quantity is stored. Further
  columns (`arlem_score`, `max_length`, or anything else) are preserved.

- dom_code_namespace:

  Optional scalar string identifying the run whose `dom_code` values
  this table is keyed by, see
  [`tales_namespace`](https://scunnac.github.io/tantale/reference/tales_namespace.md).
  Relevant for `domain_distances()`, whose ids *are* `dom_code`s.

## Value

A validated `pairwise_distances` object.

## Details

`tale_distances()` and `domain_distances()` are the entity-specific
flavours: distances between whole TALE arrays, and between distinct
domains, repeats and the two termini alike. For `domain_distances`,
`dissim` is the percentage of amino acids that differ between two
domains (see
[`tales_domain_distances`](https://scunnac.github.io/tantale/reference/tales_domain_distances.md));
for `tale_distances`, it is the array alignment cost between two TALEs
(see
[`tales_tale_distances`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md)).
The two flavours add no structure, only meaning: every method is written
once on the parent.

## See also

Other pairwise distances:
[`[.pairwise_distances()`](https://scunnac.github.io/tantale/reference/sub-.pairwise_distances.md),
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/reference/is_pairwise_distances.md),
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
d <- data.frame(
  id1 = c("A1", "A1", "A2", "A2"),
  id2 = c("A1", "A2", "A1", "A2"),
  dissim = c(0, 35, 35, 0)
)
pairwise_distances(d)
domain_distances(d, dom_code_namespace = "example")

# A similarity is folded into the distance.
similar <- data.frame(id1 = "A1", id2 = "A2", sim = 65)
tale_distances(similar)
```
