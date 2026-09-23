# Combine tales objects

Row-binds one or more
[`tales`](https://scunnac.github.io/tantale/reference/tales.md) objects
into one, reconciling the invariants that a plain
[`bind_rows`](https://dplyr.tidyverse.org/reference/bind_rows.html)
would not check: `array_id` uniqueness across inputs, the `dom_code`
namespace, and the `group` column. A single input is accepted: it goes
through the same reconciliation and re-validation, harmlessly.

## Usage

``` r
tales_bind(..., on_namespace_mismatch = c("recode", "error"), sanitize = FALSE)
```

## Arguments

- ...:

  One or more
  [`tales`](https://scunnac.github.io/tantale/reference/tales.md)
  objects (not `tales_msa`).

- on_namespace_mismatch:

  What to do when the inputs' `dom_code` namespaces (see
  [`tales_namespace`](https://scunnac.github.io/tantale/reference/tales_namespace.md))
  disagree. `"recode"` (default) drops the stale `dom_code` column and
  calls
  [`tales_assign_domain_codes`](https://scunnac.github.io/tantale/reference/tales_assign_domain_codes.md)
  on the bound result, over the union of `aa_seq`; `"error"` aborts
  instead, for callers who want that stricter behaviour.

- sanitize:

  Passed to
  [`tales`](https://scunnac.github.io/tantale/reference/tales.md) for
  the final construction.

## Value

A validated
[`tales`](https://scunnac.github.io/tantale/reference/tales.md) object.

## Details

There is no [`c()`](https://rdrr.io/r/base/c.html) method: base
[`c()`](https://rdrr.io/r/base/c.html) can silently drop attributes when
a `tales` is mixed with another object, and its signature has no room
for `on_namespace_mismatch`.

`tales_msa` inputs are refused: `alignment_width` and
`alignment_position` are a coordinate system specific to one alignment
run, and binding two runs' matrices would produce an object that
misdescribes what a given column means in each row. Demote explicitly
with
[`as_tales`](https://scunnac.github.io/tantale/reference/as_tales.md)
first if that is really what you want.

`group` is dropped from the result whenever present on any input: it is
a clustering result over one specific distance matrix and one specific
set of arrays, so two "group 1"s from separate
[`tales_group_hclust`](https://scunnac.github.io/tantale/reference/tales_group_hclust.md)/[`tales_group_kmedoids`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md)
calls are not comparable. Recompute it with either on the bound result,
after re-running
[`tales_compare_distal`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md).

`tale_distances`/`domain_distances` are untouched: this function only
binds `tales` data. A companion distance table from either input does
not describe the bound object; re-run
[`tales_compare_distal`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md)
if you need one.

## See also

Other tales objects:
[`[.tales()`](https://scunnac.github.io/tantale/reference/sub-.tales.md),
[`as_tales()`](https://scunnac.github.io/tantale/reference/as_tales.md),
[`format.tales()`](https://scunnac.github.io/tantale/reference/format.tales.md),
[`format.tales_msa()`](https://scunnac.github.io/tantale/reference/format.tales_msa.md),
[`is_tales()`](https://scunnac.github.io/tantale/reference/is_tales.md),
[`new_tales()`](https://scunnac.github.io/tantale/reference/new_tales.md),
[`print.tales()`](https://scunnac.github.io/tantale/reference/print.tales.md),
[`print.tales_msa()`](https://scunnac.github.io/tantale/reference/print.tales_msa.md),
[`summary.tales()`](https://scunnac.github.io/tantale/reference/summary.tales.md),
[`summary.tales_msa()`](https://scunnac.github.io/tantale/reference/summary.tales_msa.md),
[`tales()`](https://scunnac.github.io/tantale/reference/tales.md),
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md),
[`tales_anomalies()`](https://scunnac.github.io/tantale/reference/tales_anomalies.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/reference/tales_assert_complete.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/reference/validate_tales.md)

## Examples

``` r
a <- tales(data.frame(
  array_id = "A1", position_in_array = 1:2,
  rvd = c("NTERM", "HD")
))
b <- tales(data.frame(
  array_id = "A2", position_in_array = 1:2,
  rvd = c("NTERM", "NI")
))
tales_bind(a, b)
#> <tales> 2 arrays, 4 parts
#>   layers: rvd
#>       rvd
#>   A1  NTERM HD
#>   A2  NTERM NI
```
