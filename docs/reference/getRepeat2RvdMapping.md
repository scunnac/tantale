# Generate a mapping between Distal repeat IDs and their cognate RVD.

Uses Distal repeat sequences and RVD sequences from a set of TALEs to
return the association between repeat ID and RVD.

## Usage

``` r
getRepeat2RvdMapping(talesRepeatVectors, talesRvdVectors)
```

## Arguments

- talesRepeatVectors:

  Expects a list of Distal RVDs character vectors. Each **named**
  element corresponding to a TALE.

## Value

A two columns repeatID - RVD data frame.

## Details

Care must be taken that TALEs in the two sets of sequences have the same
name. In addition, the function tries hard to make sure that the two
sets of sequences are identical in every ways but the individual
'values' they contain. It is therefore notably important to make sure
that the sequences are consistent in whether they include N-term and
C-term domains IDs/Tags or not.
