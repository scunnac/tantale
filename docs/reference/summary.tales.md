# Summarise a tales object

The measures worth knowing about a set of TALE arrays before doing
anything with them, and the ones too expensive for
[`print.tales`](https://scunnac.github.io/tantale/reference/print.tales.md).

## Usage

``` r
# S3 method for class 'tales'
summary(object, ...)

# S3 method for class 'summary.tales'
print(x, ...)
```

## Arguments

- object:

  A [`tales`](https://scunnac.github.io/tantale/reference/tales.md)
  object.

- ...:

  Unused.

- x:

  A `summary.tales` object.

## Value

An object of class `summary.tales`: a list of the measures, so they can
be used as well as read.

## Details

**Distinct domains** counts different `dom_code`s against the number of
parts, broken down by domain type. These are *domains*, repeats and
termini alike: a `dom_code` identifies any distinct part sequence, and
the two termini are parts like the repeats are. In the example output
shipped with the package (see
[`tales_from_telltales`](https://scunnac.github.io/tantale/reference/tales_from_telltales.md)),
47 distinct codes cover 39 repeats and 8 termini, so reading the total
as a count of distinct repeats would overstate repeat diversity.

The ratio to parts is a property of the biology: TALEs reuse repeats
heavily, within an array and between arrays. It is also what decides the
cost of
[`tales_compare_distal`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md),
whose pairwise comparison runs over distinct domains.

**Repeats per array** counts repeats only, excluding the termini, so it
is the number that determines how long a target box each TALE
recognises. Note that
[`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md)'s
log reports array length as the number of domain hits, which includes
the two termini and is therefore larger by up to two.

**Complete arrays** are those where both an N- and a C-terminus were
identified. An incomplete one is not necessarily wrong (the array may
sit at the end of a contig, or a terminus may simply not have been
detected), but completeness is a precondition several downstream
functions depend on.

**Anomalies** come from
[`tales_anomalies`](https://scunnac.github.io/tantale/reference/tales_anomalies.md).
They are worth computing but too slow to compute every time an object is
printed, which is why they appear here.

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
[`summary.tales_msa()`](https://scunnac.github.io/tantale/reference/summary.tales_msa.md),
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
parts <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L, 3L),
  domain_type = c("N-terminus", "repeat", "C-terminus",
                  "N-terminus", "repeat", "C-terminus"),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "NI", "CTERM")
)
x <- tales(parts)
summary(x)
#> <tales> summary
#>   arrays / parts            2 / 6
#>   distinct RVDs             2
#>   repeats per array         min 1   median 1   max 1
#>   arrays with both termini  2 of 2
#>   anomalies                 none
```
