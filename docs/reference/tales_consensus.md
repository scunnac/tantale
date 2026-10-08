# Compute a consensus from a TALE msa

Pick the most frequent element in each column of the alignment matrix.

## Usage

``` r
tales_consensus(align)
```

## Arguments

- align:

  A TALE alignment as a character matrix, one row per array and one
  column per alignment position, with `NA` for gaps.
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) on a
  [`tales_msa`](https://scunnac.github.io/tantale/reference/tales_msa.md)
  returns one (see
  [`as.matrix.tales_msa`](https://scunnac.github.io/tantale/reference/as.matrix.tales_msa.md)),
  filled with RVDs or domain codes.

## Value

A vector of consensus elements in each column of `align`.

## Details

A column has a consensus only when one element is strictly more common
than every other. Where two or more are tied for most frequent – as
happens whenever each array carries a different repeat at that position
– the result is `NA`, because there is no agreement to report. `NA` is
likewise returned when the most common thing at a position is a gap.

## See also

Other TALE alignment:
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/reference/as.matrix.tales_msa.md),
[`is_tales_msa()`](https://scunnac.github.io/tantale/reference/is_tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md),
[`tales_consensus_match()`](https://scunnac.github.io/tantale/reference/tales_consensus_match.md),
[`tales_msa()`](https://scunnac.github.io/tantale/reference/tales_msa.md),
[`tales_msa_params()`](https://scunnac.github.io/tantale/reference/tales_msa_params.md),
[`tales_msa_width()`](https://scunnac.github.io/tantale/reference/tales_msa_width.md),
[`validate_tales_msa()`](https://scunnac.github.io/tantale/reference/validate_tales_msa.md)

## Examples

``` r
# column 1: HD is a clear majority. column 2: a three-way tie, no consensus.
aln <- matrix(c("HD", "HD", "NI",
               "NG", "NI", "HD"),
             nrow = 3, dimnames = list(c("A1", "A2", "A3"), NULL))
aln
#>    [,1] [,2]
#> A1 "HD" "NG"
#> A2 "HD" "NI"
#> A3 "NI" "HD"
tales_consensus(aln)
#> [1] "HD" NA  

# From a tales_msa, such as tales_align() returns
msa <- tales_msa(data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2", "A3", "A3", "A3"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L, 1L, 2L, 3L),
  alignment_position = c(1L, 2L, 3L, 1L, 3L, 1L, 2L, 3L),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "CTERM", "NTERM", "HD", "CTERM")
))
tales_consensus(as.matrix(msa))
#> [1] "NTERM" "HD"    "CTERM"
```
