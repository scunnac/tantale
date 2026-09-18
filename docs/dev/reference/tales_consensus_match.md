# Do elements in a TALE msa match the consensus?

Compute a logical matrix corresponding to the input `align` input with
`TRUE` where an element matches the consensus at that position and
`FALSE` where it does not. Columns with no consensus – see
[`tales_consensus`](https://scunnac.github.io/tantale/dev/reference/tales_consensus.md)
– are `NA` throughout, since there is nothing there to match.

## Usage

``` r
tales_consensus_match(align, long = TRUE)
```

## Arguments

- align:

  A multiple Tal sequences alignment in the form of a matrix.

- long:

  Set to `TRUE` (default) to return a long tibble, or `FALSE` to return
  a logical matrix with the same shape as `align`.

## Value

A multiple Tal sequences alignment in the form of a matrix filled with
logical values if `long` is `FALSE` and a long tibble representing the
original alignment otherwise (default).

## See also

Other TALE alignment:
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.tales_msa.md),
[`is_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/is_tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/dev/reference/tales_align.md),
[`tales_consensus()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus.md),
[`tales_msa()`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md),
[`tales_width()`](https://scunnac.github.io/tantale/dev/reference/tales_width.md),
[`validate_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/validate_tales_msa.md)

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
#>   array_id position_in_array tales_consensus_match
#>   <fct>                <int> <lgl>                
#> 1 A1                       1 TRUE                 
#> 2 A2                       1 TRUE                 
#> 3 A3                       1 FALSE                
#> 4 A1                       2 NA                   
#> 5 A2                       2 NA                   
#> 6 A3                       2 NA                   
```
