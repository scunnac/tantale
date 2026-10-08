# Do elements in a TALE msa match the consensus?

Compute a logical matrix corresponding to the input `align` input with
`TRUE` where an element matches the consensus at that position and
`FALSE` where it does not. Columns with no consensus – see
[`tales_consensus`](https://scunnac.github.io/tantale/reference/tales_consensus.md)
– are `NA` throughout, since there is nothing there to match.

## Usage

``` r
tales_consensus_match(align, long = TRUE)
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

- long:

  Set to `TRUE` (default) to return a long tibble, or `FALSE` to return
  a logical matrix with the same shape as `align`.

## Value

A multiple Tal sequences alignment in the form of a matrix filled with
logical values if `long` is `FALSE` and a long tibble representing the
original alignment otherwise (default).

## See also

Other TALE alignment:
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/reference/as.matrix.tales_msa.md),
[`is_tales_msa()`](https://scunnac.github.io/tantale/reference/is_tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md),
[`tales_consensus()`](https://scunnac.github.io/tantale/reference/tales_consensus.md),
[`tales_msa()`](https://scunnac.github.io/tantale/reference/tales_msa.md),
[`tales_msa_params()`](https://scunnac.github.io/tantale/reference/tales_msa_params.md),
[`tales_msa_width()`](https://scunnac.github.io/tantale/reference/tales_msa_width.md),
[`validate_tales_msa()`](https://scunnac.github.io/tantale/reference/validate_tales_msa.md)

## Examples

``` r
aln <- matrix(c("HD", "HD", "NI",
               "NG", "NI", "HD"),
             nrow = 3, dimnames = list(c("A1", "A2", "A3"), NULL))
tales_consensus_match(aln, long = FALSE)
#>     [,1] [,2]
#> A1  TRUE   NA
#> A2  TRUE   NA
#> A3 FALSE   NA
tales_consensus_match(aln)
#> # A tibble: 6 × 3
#>   array_id alignment_position tales_consensus_match
#>   <fct>    <fct>              <lgl>                
#> 1 A1       A                  TRUE                 
#> 2 A2       A                  TRUE                 
#> 3 A3       A                  FALSE                
#> 4 A1       B                  NA                   
#> 5 A2       B                  NA                   
#> 6 A3       B                  NA                   
```
