# Compute a consensus from a TALE msa

Pick the most frequent element in each column of the alignment matrix.

## Usage

``` r
tales_consensus(align)
```

## Arguments

- align:

  A multiple Tal sequences alignment in the form of a matrix.

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
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.tales_msa.md),
[`is_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/is_tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/dev/reference/tales_align.md),
[`tales_consensus_match()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus_match.md),
[`tales_msa()`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md),
[`tales_width()`](https://scunnac.github.io/tantale/dev/reference/tales_width.md),
[`validate_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/validate_tales_msa.md)

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
```
