# Pairwise distances between TALE arrays

Step 3 of
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md).
Aligns the arrays against each other as strings of domain codes, using
the domain distances from step 2 as the cost of substituting one domain
for another.

## Usage

``` r
tales_tale_distances(x, domain_distances)
```

## Arguments

- x:

  A [tales](https://scunnac.github.io/tantale/dev/reference/tales.md)
  object carrying `dom_code`.

- domain_distances:

  A
  [domain_distances](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md)
  object for the same `x`, as returned by
  [`tales_domain_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_distances.md).

## Value

A
[tale_distances](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md)
object, keyed by `array_id`.

## Details

The two arguments are not independent, and that is the point. ARLEM
aligns each array's sequence of `dom_code`s; what it costs to align one
domain against a different one is taken from `domain_distances`, so the
TALE-level comparison is built on the domain-level one rather than
computed beside it.

The domain distances are first passed through a Minkowski distance
(`p = 3.5`) between their rows and rescaled to 0-100. That step is not
cosmetic: ARLEM needs a cost matrix satisfying the triangle inequality,
and raw pairwise alignment dissimilarities do not.

## Both arguments must come from the same call

`domain_distances` is keyed by `dom_code`, and those codes mean what
they mean only within the call that minted them
([`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md)).
Passing distances from one run with codes from another silently compares
the wrong domains, which is why both objects carry a namespace stamp and
this function refuses when they disagree.

## See also

[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md),
which runs all three steps.

Other pairwise distances:
[`as.matrix.pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.pairwise_distances.md),
[`distances_assert_square()`](https://scunnac.github.io/tantale/dev/reference/distances_assert_square.md),
[`distances_restrict()`](https://scunnac.github.io/tantale/dev/reference/distances_restrict.md),
[`is_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/is_pairwise_distances.md),
[`pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md),
[`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md),
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md),
[`tales_domain_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_distances.md),
[`tales_group()`](https://scunnac.github.io/tantale/dev/reference/tales_group.md),
[`validate_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/validate_pairwise_distances.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
xa <- tales_assign_domain_codes(x)
dd <- tales_domain_distances(xa)
#> Computing a distance matrix between TALE parts amino acid sequences using:
#> DECIPHER
tales_tale_distances(xa, dd)
#> Generate an ARLEM cost matrix which meets triangle inequality criteria by
#> computing the minkowski distance between pairwise distance vectors.
#> Running ARLEM version 1.0 :
#> Copyright by Mohamed I. Abouelhoda
#> Plz. cite Abouelhoda, Giegerich, Behzadi, and Steyaert
```
