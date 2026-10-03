# Summarise a tales_msa object

What an alignment looks like, in numbers: how gappy it is, and how much
the arrays actually agree.

## Usage

``` r
# S3 method for class 'tales_msa'
summary(object, ...)

# S3 method for class 'summary.tales_msa'
print(x, ...)
```

## Arguments

- object:

  A
  [`tales_msa`](https://scunnac.github.io/tantale/reference/tales_msa.md)
  object.

- ...:

  Unused.

- x:

  A `summary.tales_msa` object.

## Value

An object of class `summary.tales_msa`.

## Details

**Columns with no consensus** is the measure worth having. A column has
no consensus when no element is strictly more common than the rest (see
[`tales_consensus`](https://scunnac.github.io/tantale/reference/tales_consensus.md)),
so the count says how much of the alignment the arrays agree on.

It is reported per layer, and the two usually differ in a way that means
something. Repeats that are distinct proteins can carry the same RVD, so
a column can lack a `dom_code` consensus while having a clear `rvd` one:
the arrays disagree about the repeat and agree about the base it binds.
The reverse is rarer and more interesting.

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
[`tales()`](https://scunnac.github.io/tantale/reference/tales.md),
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md),
[`tales_anomalies()`](https://scunnac.github.io/tantale/reference/tales_anomalies.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/reference/tales_assert_complete.md),
[`tales_bind()`](https://scunnac.github.io/tantale/reference/tales_bind.md),
[`tales_names()`](https://scunnac.github.io/tantale/reference/tales_names.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/reference/validate_tales.md)

## Examples

``` r
aligned <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L),
  alignment_position = c(1L, 2L, 3L, 1L, 3L),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "CTERM")
)
msa <- tales_msa(aligned)
summary(msa)
#> <tales_msa> summary
#>   arrays / width            2 / 3
#>   gaps                      1 of 6 cells (17%)
#>   columns with no gap       2 of 3
#>   no consensus (rvd)        1 of 3 columns
#>   matching consensus (rvd)  4 of 4 scorable cells
#>   anomalies                 none
```
