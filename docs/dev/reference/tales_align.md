# Align the repeat arrays of a tales object

Aligns TALE arrays on one of their residue layers with MAFFT's text
mode, returning the alignment as a
[`tales_msa`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md)
— the input object with an `alignment_position` column added, so every
other layer (`rvd`, `dom_code`, `aa_seq`, ...) remains available.

## Usage

``` r
tales_align(
  x,
  residue_col = c("rvd", "dom_code"),
  repeat_sims = NULL,
  mafft_opts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
  mafft_path = NULL,
  mafft_verbose = FALSE,
  ...
)
```

## Arguments

- x:

  A [`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
  object holding complete arrays.

- residue_col:

  Which layer to align on: `"rvd"` (default) or `"dom_code"`. Given
  explicitly rather than guessed from the values.

- repeat_sims:

  Scoring matrix for the residues being aligned. `NULL` (default) means
  none. Pass `"rvd"` to opt in to the built-in RVD similarity matrix
  when aligning RVDs, or a
  [`domain_distances`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md)
  object when aligning repeat codes. Optional similarity table passed to
  MAFFT as a scoring matrix, as accepted by `tales_align`.

- mafft_opts:

  Command-line options handed to MAFFT. The default,
  `"--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1"`,
  asks for an accurate local alignment and sets the gap penalties.

  `--op` is the cost of opening a gap and `--ep` the cost of extending
  one. They are set unusually here – free to open, expensive to extend –
  because of what is being aligned. A TALE gains or loses whole repeats,
  so a gap should start wherever it needs to; what should be discouraged
  is one long gap swallowing a stretch of repeats that really do
  correspond. Raise `--op` if the alignment fragments into too many
  small gaps.

- mafft_path:

  Where to find MAFFT. `NULL`, the default, uses the `tantale` conda
  environment, creating it on first use. Give the root of a standalone
  MAFFT directory instead (one holding `mafft.bat` with the helpers
  under `mafftdir/libexec`) to use your own copy – but note that the
  version matters: MAFFT changed how it aligns text-mode sequences after
  7.4x, and later releases leave the termini of a TALE alignment
  unanchored.

- mafft_verbose:

  Whether to let MAFFT write to the console. It reports its banner, the
  strategy it chose and its progress through the sequences, which is
  dozens of lines per alignment and rarely what you want. Left `FALSE`
  that output is captured rather than discarded, and replayed if the
  alignment fails – so silence costs nothing diagnostically.

- ...:

  Further arguments to the MAFFT runner, chiefly `gap_symbol`, the value
  gaps take in the returned matrix (`NA` by default).

## Value

A `tales_msa` object.

## Details

The mapping back from MAFFT's output is **positional**: the k-th non-gap
cell of an aligned row is the k-th part fed in. That is well defined
only because this function builds MAFFT's input from `x` itself, which
is why
[`tales_assert_complete`](https://scunnac.github.io/tantale/dev/reference/tales_assert_complete.md)
is enforced first.

## See also

Other TALE alignment:
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.tales_msa.md),
[`is_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/is_tales_msa.md),
[`tales_consensus()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus.md),
[`tales_consensus_match()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus_match.md),
[`tales_msa()`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md),
[`tales_width()`](https://scunnac.github.io/tantale/dev/reference/tales_width.md),
[`validate_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/validate_tales_msa.md)

## Examples

``` r
# \donttest{
# Needs MAFFT, resolved from the tantale conda environment on first use.
rvd_fasta <- system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",
                         package = "tantale")
x <- as_tales(rvd_fasta, sep = "-")
msa <- tales_align(x, residue_col = "rvd")
#> Now running MAFFT (Copyright 2002-2007 Kazutaka Katoh) on TALE array sequences.
as.matrix(msa)
#>               1    2    3    4    5    6    7    8    9    10   11   12   13  
#> TalA_BAI3     "NN" "NG" "NN" "HD" "HD" "NI" "N*" "NG" "HD" "NI" "NG" "NN" "HD"
#> TalA_CFBP1947 "NN" "NG" "NN" "HD" "HD" "NI" "N*" "NG" "HD" "NI" "NG" "NN" "HD"
#> TalA_MAI1     "NN" "N*" "NN" "HD" "HD" "NI" "N*" "NG" "HD" "NI" "NG" "NN" "HD"
#> TalA_MAI106   "NN" "N*" "NN" "HD" NA   NA   NA   NA   "HD" "NI" "NG" "NN" "HD"
#> TalA_MAI129   "NN" "N*" "NN" "HD" NA   NA   NA   NA   "HD" "NI" "NG" "NN" "HD"
#> TalA_MAI134   "NN" "N*" "NN" "HD" "HD" "NI" "N*" "NG" "HD" "NI" "NG" "NN" "HD"
#> TalA_MAI145   "NN" "N*" "NN" "HD" NA   NA   NA   NA   "HD" "NI" "NG" "NN" "HD"
#> TalA_MAI68    "NN" "N*" "NN" "HD" NA   NA   NA   NA   "HD" "NI" "NG" "NN" "HD"
#> TalA_MAI73    "NN" "N*" "NN" "HD" NA   NA   NA   NA   "HD" "NI" "NG" "NN" "HD"
#> TalA_MAI95    "NN" "N*" "NN" "HD" NA   NA   NA   NA   "HD" "NI" "NG" "NN" "HD"
#> TalA_MAI99    "NN" "N*" "NN" "HD" NA   NA   NA   NA   "HD" "NI" "NG" "NN" "HD"
#>               14   15   16   17   18   19   20   21   22   23   24   25   26  
#> TalA_BAI3     "NI" "NG" "NI" "NG" "NN" "NG" "HD" "NI" "NI" "NG" "HD" "NN" "NG"
#> TalA_CFBP1947 "NI" "NG" "NI" "NG" "NN" "NG" "HD" "NI" "NI" "NG" "HD" "NN" "NG"
#> TalA_MAI1     "NS" "NG" "NI" "NG" "NN" "NG" "HD" "NI" "NI" "NG" "HD" "NN" "NG"
#> TalA_MAI106   "NS" "NG" "NI" "N*" "NN" "NG" "HD" "NI" "NI" "NG" "HD" "NN" "NG"
#> TalA_MAI129   "NS" "NG" "NI" "N*" "NN" "NG" "HD" "NI" "NI" "NG" "HD" "NN" "NG"
#> TalA_MAI134   "NS" "NG" "NI" "NG" "NN" "NG" "HD" "NI" "NI" "NG" "HD" "NN" "NG"
#> TalA_MAI145   "NS" "NG" "NI" "N*" "NN" "NG" "HD" "NI" "NI" "NG" "HD" "NN" "NG"
#> TalA_MAI68    "NS" "NG" "NI" "N*" "NN" "NG" "HD" "NI" "NI" "NG" "HD" "NN" "NG"
#> TalA_MAI73    "NS" "NG" "NI" "N*" "NN" "NG" "HD" "NI" "NI" "NG" "HD" "NN" "NG"
#> TalA_MAI95    "NS" "NG" "NI" "N*" "NN" "NG" "HD" "NI" "NI" "NG" "HD" "NN" "NG"
#> TalA_MAI99    "NS" "NG" "NI" "N*" "NN" "NG" "HD" "NI" "NI" "NG" "HD" "NN" "NG"
# }
```
