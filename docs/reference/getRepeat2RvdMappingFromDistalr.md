# Generate a mapping between Distal repeat IDs and their cognate RVD.

Uses Distal repeat sequences and RVD sequences from a set of TALEs
analyzed with the
[`distalr`](https://scunnac.github.io/tantale/reference/distalr.md)
function to return the association between repeat ID and RVD.

## Usage

``` r
getRepeat2RvdMappingFromDistalr(distalrTaleParts)
```

## Arguments

- distalrTaleParts:

  The taleParts object in a
  [`distalr`](https://scunnac.github.io/tantale/reference/distalr.md)
  output.

## Value

A two columns repeatID - RVD data frame.
