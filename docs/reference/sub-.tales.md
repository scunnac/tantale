# Subset a tales object

Subsets rows and columns like an ordinary tibble, but the class travels
only while the result still satisfies the `tales` contract.

## Usage

``` r
# S3 method for class 'tales'
x[...]
```

## Arguments

- x:

  A `tales` object.

- ...:

  Passed on to the tibble/data frame method.

## Value

The subset: still a `tales` (or `tales_msa`) if the contract holds,
otherwise a plain tibble.

## Details

Row subsetting never breaks anything: every invariant a `tales` checks
is closed under keeping a subset of rows, so filtering to one array, or
to its repeats only, is still a valid `tales`. Column subsetting is
where it degrades: dropping `array_id`, or both residue columns (`rvd`,
`dom_code`) at once, leaves something that can no longer be described as
a `tales`, and the class quietly steps out of the way rather than
continuing to claim invariants it can no longer keep – the result is a
plain tibble, not an error. Dropping an optional column (`seqnames`,
`aa_seq`, ...) has no such consequence. A `tales_msa` degrades one step
at a time: losing `alignment_position` alone steps it back to a plain
`tales`, not all the way to a tibble.

Attributes travel with a valid subset. `dom_code_namespace` describes
the run the codes came from, not which rows happen to be kept right now,
so cutting an object down with `[` keeps its namespace unchanged.

## See also

Other tales objects:
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
[`tales_bind()`](https://scunnac.github.io/tantale/reference/tales_bind.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/reference/validate_tales.md)

## Examples

``` r
x <- tales(data.frame(
  array_id = c("A1", "A1"), position_in_array = c(1L, 2L),
  rvd = c("NTERM", "HD")
))
is_tales(x[1, ]) # row subsetting never breaks the contract
#> [1] TRUE
class(x[, -1])   # dropping array_id: the class steps out of the way
#> [1] "tbl_df"     "tbl"        "data.frame"
```
