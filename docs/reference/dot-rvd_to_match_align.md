# Recode an RVD alignment as similarity to a reference row

Substitutes each RVD with a score expressing how similar its DNA-binding
preference is to the RVD of a reference TALE, column by column. This is
the RVD-level counterpart of `.domain_to_sim_align()`, which works on
protein sequence similarity instead: the two come apart, since domains
can be sequence-divergent yet share an RVD, or near-identical yet differ
at positions 12-13.

## Usage

``` r
.rvd_to_match_align(rvd_align, rvd_sims = rvdSimDf, ref_tag = NULL)
```

## Arguments

- rvd_align:

  A character matrix of aligned RVDs.

- rvd_sims:

  A data frame of pairwise RVD similarity with columns `rvd1`, `rvd2`
  and `Cor`. Defaults to the package's internal `rvdSimDf`.

- ref_tag:

  Pattern selecting the reference row; see `.pick_ref_name()`.

## Value

A numeric matrix with the dimensions and dimnames of `rvd_align`.

## Details

Wired through `fill_type = "rvd_sim"` in
[`plot.tales_msa`](https://scunnac.github.io/tantale/reference/plot.tales_msa.md)
– the RVD-level counterpart of `"domain_sim"`. It is the only consumer
of the internal `rvdSimDf` dataset.
