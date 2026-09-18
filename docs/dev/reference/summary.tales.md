# Summarise a tales object

The measures worth knowing about a set of TALE arrays before doing
anything with them, and the ones too expensive for
[`print.tales`](https://scunnac.github.io/tantale/dev/reference/print.tales.md).

## Usage

``` r
# S3 method for class 'tales'
summary(object, ...)

# S3 method for class 'summary.tales'
print(x, ...)
```

## Arguments

- object:

  A [`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
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
parts, broken down by domain type. *Domains*, not repeats: a `dom_code`
identifies any distinct part sequence, and the two termini are parts
like the repeats are. On the reference fixture 251 distinct codes cover
180 repeats and 71 termini, so calling the total a repeat count would
overstate it by nearly a third.

The ratio to parts is a property of the biology: TALEs reuse repeats
heavily, within an array and between arrays. It is also what decides the
cost of
[`tales_compare`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md),
whose pairwise comparison runs over distinct domains rather than over
parts.

**Repeats per array** counts repeats only, excluding the termini, so it
is the number that determines how long a target box each TALE
recognises. Note that
[`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)'s
log reports array length as the number of domain hits, which includes
the two termini and is therefore larger by up to two.

**Complete arrays** are those where both an N- and a C-terminus were
identified. An incomplete one is not necessarily wrong – the array may
sit at the end of a contig, or a terminus may simply not have been
detected – but it is the precondition several downstream functions
depend on.

**Anomalies** come from
[`tales_anomalies`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md),
and are the reason this is a
[`summary()`](https://rdrr.io/r/base/summary.html) rather than part of
printing: they are worth computing, and too slow to compute every time
an object is echoed.

## See also

Other tales objects:
[`as_tales()`](https://scunnac.github.io/tantale/dev/reference/as_tales.md),
[`format.tales()`](https://scunnac.github.io/tantale/dev/reference/format.tales.md),
[`format.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/format.tales_msa.md),
[`is_tales()`](https://scunnac.github.io/tantale/dev/reference/is_tales.md),
[`new_tales()`](https://scunnac.github.io/tantale/dev/reference/new_tales.md),
[`print.tales()`](https://scunnac.github.io/tantale/dev/reference/print.tales.md),
[`print.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/print.tales_msa.md),
[`summary.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/summary.tales_msa.md),
[`tales()`](https://scunnac.github.io/tantale/dev/reference/tales.md),
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_anchor_codes.md),
[`tales_anomalies()`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/dev/reference/tales_assert_complete.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md),
[`tales_requirements()`](https://scunnac.github.io/tantale/dev/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/dev/reference/validate_tales.md)
