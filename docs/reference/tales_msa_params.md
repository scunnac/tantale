# Settings a TALE alignment was made with

What
[`tales_align`](https://scunnac.github.io/tantale/reference/tales_align.md)
aligned on and how, kept with the alignment so that it can be reproduced
or interpreted later, for example to tell whether a `tales_msa` holding
both `rvd` and `dom_code` columns was aligned on its RVDs or on its
domain codes. Subsetting keeps these settings; returning to a plain
`tales`
([`as_tales`](https://scunnac.github.io/tantale/reference/as_tales.md))
drops them.

## Usage

``` r
tales_msa_params(x)
```

## Arguments

- x:

  A `tales_msa` object.

## Value

A list with `residue_col`, the column aligned on; `domain_distances`,
the scoring given to
[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md)
(`NULL`, `"rvd"`, or the matrix, table or file supplied); and
`mafft_opts`, the options MAFFT ran with. `NULL` for a `tales_msa` built
by hand with
[`tales_msa`](https://scunnac.github.io/tantale/reference/tales_msa.md).

## See also

Other TALE alignment:
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/reference/as.matrix.tales_msa.md),
[`is_tales_msa()`](https://scunnac.github.io/tantale/reference/is_tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md),
[`tales_consensus()`](https://scunnac.github.io/tantale/reference/tales_consensus.md),
[`tales_consensus_match()`](https://scunnac.github.io/tantale/reference/tales_consensus_match.md),
[`tales_msa()`](https://scunnac.github.io/tantale/reference/tales_msa.md),
[`tales_msa_width()`](https://scunnac.github.io/tantale/reference/tales_msa_width.md),
[`validate_tales_msa()`](https://scunnac.github.io/tantale/reference/validate_tales_msa.md)

## Examples

``` r
# \donttest{
# Needs MAFFT, resolved from the tantale conda environment on first use.
rvd_fasta <- system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",
                         package = "tantale")
msa <- tales_align(as_tales(rvd_fasta, sep = "-"), residue_col = "rvd")
#> Now running MAFFT (Copyright 2002-2007 Kazutaka Katoh) on TALE array sequences.
tales_msa_params(msa)$residue_col
#> [1] "rvd"
# }
```
