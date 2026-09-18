# Align TALE sequences with MAFFT text mode

Implementation behind
[`tales_align`](https://scunnac.github.io/tantale/dev/reference/tales_align.md)
and the deprecated
[`tales_align`](https://scunnac.github.io/tantale/dev/reference/tales_align.md).
Internal so that package code can call it without tripping the
deprecation warning. Normalise a pairwise table to what MAFFT's
–textmatrix wants

## Usage

``` r
.as_mafft_score_table(x)
```

## Arguments

- x:

  A data frame of pairwise scores.

## Value

A three-column data frame named `id1`, `id2`, `sim`.

## Details

MAFFT scores matches, so it wants a *similarity*: higher means more
alike. Accepts either the canonical
[`pairwise_distances`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md)
vocabulary (`id1`, `id2`, `dissim`) or the legacy one (`RepU1`, `RepU2`,
`Sim`), inverting the distance on the way in.
