# Normalise a pairwise table to what MAFFT's –textmatrix wants

MAFFT scores matches, so it wants a *similarity*: higher means more
alike. Takes the
[`pairwise_distances`](https://scunnac.github.io/tantale/reference/pairwise_distances.md)
columns (`id1`, `id2`, `dissim`, or `sim`), inverting the distance on
the way in.

## Usage

``` r
.as_mafft_score_table(x)
```

## Arguments

- x:

  A data frame of pairwise scores.

## Value

A three-column data frame named `id1`, `id2`, `sim`.
