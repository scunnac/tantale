# Width of a TALE alignment

The number of columns in the alignment, including those that are all
gaps in the object at hand. Stored rather than derived: subsetting
arrays can empty the last column, which would silently shrink
`max(alignment_position)`.

## Usage

``` r
tales_width(x)
```

## Arguments

- x:

  A `tales_msa` object.

## Value

An integer scalar, or `NULL` if unset.

## See also

Other TALE alignment:
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.tales_msa.md),
[`is_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/is_tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/dev/reference/tales_align.md),
[`tales_consensus()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus.md),
[`tales_consensus_match()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus_match.md),
[`tales_msa()`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md),
[`validate_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/validate_tales_msa.md)
