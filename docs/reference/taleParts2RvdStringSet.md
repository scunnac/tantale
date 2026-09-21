# Generates a RVD sequences set from a taleParts object

Uses a taleParts object in a
[`distalr`](https://scunnac.github.io/tantale/reference/distalr.md)
output to return a `BStringSet` of RVD sequences. RVDs are separated by
the character specified in the `sep` parameter.

## Usage

``` r
taleParts2RvdStringSet(taleParts, sep = "-")
```

## Arguments

- sep:

  Used as a RVD separatator

- distalrTaleParts:

  The taleParts object in a
  [`distalr`](https://scunnac.github.io/tantale/reference/distalr.md)
  output.

## Value

A two columns repeatID - RVD data frame.

## Details

Uses Distal repeat sequences and RVD sequences from a set of TALEs
analyzed with the
[`distalr`](https://scunnac.github.io/tantale/reference/distalr.md)
function to return the association between repeat ID and RVD.
