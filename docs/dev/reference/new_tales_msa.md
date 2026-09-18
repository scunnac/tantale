# Low-level constructor for a tales_msa object

Attaches the class without validating. Use
[`tales_msa`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md)
unless the invariants are already established.

## Usage

``` r
new_tales_msa(x, alignment_width = NULL)
```

## Arguments

- x:

  A `tales` object carrying an `alignment_position` column.

- alignment_width:

  Integer width of the alignment. Stored as an attribute rather than
  derived, because subsetting arrays can empty the last column and would
  silently shrink a derived value.

## Value

A `tales_msa` object.
