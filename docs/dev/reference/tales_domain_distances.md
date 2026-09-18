# Pairwise distances between distinct TALE domains

Step 2 of
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md).
Aligns every distinct domain sequence against every other and returns
their pairwise dissimilarity.

## Usage

``` r
tales_domain_distances(
  x,
  aln_method = "DECIPHER",
  ncores = 1,
  conda_bin = "auto"
)
```

## Arguments

- x:

  A [tales](https://scunnac.github.io/tantale/dev/reference/tales.md)
  object carrying `dom_code`, as returned by
  [`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md).

- aln_method:

  One of `"DECIPHER"` (the default), `"Biostrings"` or `"mmseq2"`.
  `"mmseq2"` runs in the conda environment.

- ncores:

  Number of cores for the pairwise alignment.

- conda_bin:

  Passed to `reticulate`, for `aln_method = "mmseq2"`.

## Value

A
[domain_distances](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md)
object, stamped with `x`'s namespace.

## Details

This is the expensive step, and the one worth having on its own: the
repeat-level distances answer questions about repeat diversity that need
no TALE-level alignment at all, and computing them does not require
running ARLEM.

Distances are between **distinct domains**, keyed by `dom_code`, so the
cost goes with the number of distinct sequences rather than the number
of parts. See
[summary()](https://scunnac.github.io/tantale/dev/reference/summary.tales.md)
for that ratio on a given object.

## See also

[`tales_tale_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_tale_distances.md),
which consumes this.

Other pairwise distances:
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/dev/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/dev/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/is_pairwise_distances.md),
[`pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md),
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md),
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md),
[`tales_group()`](https://scunnac.github.io/tantale/dev/reference/tales_group.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_tale_distances.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/validate_pairwise_distances.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
xa <- tales_assign_domain_codes(x)
tales_domain_distances(xa)
#> Computing a distance matrix between TALE parts amino acid sequences using:
#> DECIPHER
```
