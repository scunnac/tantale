# Pairwise distances between distinct TALE domains

Step 2 of
[`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md).
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

  A [tales](https://scunnac.github.io/tantale/reference/tales.md) object
  carrying `dom_code`, as returned by
  [`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_assign_domain_codes.md).

- aln_method:

  One of `"DECIPHER"` (the default), `"Biostrings"` or `"mmseq2"`.
  `"mmseq2"` runs in the conda environment.

- ncores:

  Number of cores for the pairwise alignment.

- conda_bin:

  Passed to `reticulate`, for `aln_method = "mmseq2"`.

## Value

A
[domain_distances](https://scunnac.github.io/tantale/reference/pairwise_distances.md)
object, stamped with `x`'s namespace.

## Details

This is the expensive step, and the one worth having on its own: the
domain-level distances answer questions about domain diversity without
any array-level alignment.

Distances are between **distinct domains**, keyed by `dom_code`, so the
cost goes with the number of distinct sequences rather than the number
of parts. See
[summary()](https://scunnac.github.io/tantale/reference/summary.tales.md)
for that ratio on a given object.

`dissim` follows DisTAL's definition: the percentage of amino acids that
change between two domains, normalised by the longer one, with residues
one domain lacks counted as changes. A 20-residue half-repeat that is an
exact prefix of a 34-residue repeat is therefore 14/34, about 41
percent, from it. The `"DECIPHER"` and `"mmseq2"` backends compute this;
`"Biostrings"` uses a global alignment with gap penalties and scores
length differences more severely.

## See also

[`tales_tale_distances()`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md),
which consumes this.

Other pairwise distances:
[`[.pairwise_distances()`](https://scunnac.github.io/tantale/reference/sub-.pairwise_distances.md),
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/reference/is_pairwise_distances.md),
[`pairwise_distances()`](https://scunnac.github.io/tantale/reference/pairwise_distances.md),
[`rvd_dna_specificity`](https://scunnac.github.io/tantale/reference/rvd_dna_specificity.md),
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_assign_domain_codes.md),
[`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md),
[`tales_compare_functal()`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md),
[`tales_group_hclust()`](https://scunnac.github.io/tantale/reference/tales_group_hclust.md),
[`tales_group_kmedoids()`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md),
[`tales_tale_distances()`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/reference/validate_pairwise_distances.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
xa <- tales_assign_domain_codes(x)
tales_domain_distances(xa)
#> Computing a distance matrix between TALE parts amino acid sequences using:
#> DECIPHER
```
