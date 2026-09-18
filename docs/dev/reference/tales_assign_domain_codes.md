# Assign a domain code to every distinct part sequence

Step 1 of
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md).
Gives each distinct `aa_seq` in `x` an integer code, recorded in a
`dom_code` column, so that the parts can be compared once each rather
than once per occurrence.

A code names a distinct **domain** sequence, not a repeat: the N- and
C-termini are parts like the repeats are, and they get codes too.

## Usage

``` r
tales_assign_domain_codes(x)
```

## Arguments

- x:

  A [tales](https://scunnac.github.io/tantale/dev/reference/tales.md)
  object carrying an `aa_seq` column.

## Value

`x` with a `dom_code` column and a `dom_code_namespace` stamp.

## Details

TALEs reuse domains heavily, within an array and between arrays, so the
number of distinct sequences is far smaller than the number of parts.
That ratio is what makes
[`tales_domain_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_distances.md)
affordable – it compares distinct domains, not parts.

## The codes are only meaningful within one call

Codes are assigned with
[`dplyr::cur_group_id()`](https://dplyr.tidyverse.org/reference/context.html)
over the distinct `aa_seq` values **present in `x`**. Add an array,
remove one, or reorder the sequences, and the same protein can get a
different number. They are positions in this table's own vocabulary, not
identifiers of anything.

So **comparing codes between two calls is an error**, and a similarity
table keyed by one call's codes must never be used with another call's.
This is not a caution to remember: it is enforced. Every object minted
here is stamped with a `dom_code_namespace` derived from the sequences
that produced it, and the classes refuse to join across namespaces. Read
the stamp with
[`tales_namespace()`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md).

## See also

[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md),
which runs all three steps;
[`tales_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_codes.md),
which reads the code-to-sequence table back out of an object that
already has codes.

Other pairwise distances:
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/dev/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/dev/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/is_pairwise_distances.md),
[`pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md),
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md),
[`tales_domain_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_distances.md),
[`tales_group()`](https://scunnac.github.io/tantale/dev/reference/tales_group.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_tale_distances.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/validate_pairwise_distances.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
xa <- tales_assign_domain_codes(x)
tales_namespace(xa)
#> [1] "7377e3f80aa36898"
```
