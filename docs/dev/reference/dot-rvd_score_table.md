# RVD similarity matrix for aligning RVD sequences

The repeat-level similarity matrix is keyed by `dom_code`, which is
meaningless for an RVD sequence, so RVD alignments have had no scoring
matrix at all – MAFFT treated `NI` and `NN` as no more alike than `NI`
and `HD`. This supplies the missing one, from the Spearman correlation
of each RVD's A/C/G/T preference profile (`rvdSimDf`, derived from
TALVEZ's `mat1`).

## Usage

``` r
.rvd_score_table(residues)
```

## Arguments

- residues:

  Character vector of the RVDs present in the alignment.

## Value

A data frame of `id1`, `id2`, `sim` covering every ordered pair of
`residues`.

## Details

`XX` – a terminus detected but not identified – has a uniform base
profile, so its correlation with everything is undefined. Those cells
are filled as **neutral** (0): we know nothing about it, so it should
neither attract nor repel.

The one exception is `XX` against itself, which is set to the maximum.
That is forced, not chosen: MAFFT produces unusable output when the
diagonal is not high, and a low diagonal would anyway assert that an
`XX` must *not* align with an `XX`, which is a stronger claim than
ignorance.

Scale is irrelevant – MAFFT normalises the matrix, so a linear rescale
or offset leaves the alignment unchanged. Only relative structure
matters, so the correlations are used as they are.
