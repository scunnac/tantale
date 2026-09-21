# Do elements in a TALE msa match the consensus?

Compute a logical matrix corresponding to the input `align` input with
`TRUE` if an element match the consensus element at that position or
`FALSE` otherwise.

## Usage

``` r
matchConsensus(align, returnLong = TRUE)
```

## Arguments

- align:

  A multiple Tal sequences alignment in the form of a matrix.

## Value

A multiple Tal sequences alignment in the form of a matrix filled with
logical values if is `FALSE` and a long tibble representing the original
alignment otherwise (default).
