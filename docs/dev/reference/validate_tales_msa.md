# Validate a tales_msa object

Checks every
[`validate_tales`](https://scunnac.github.io/tantale/dev/reference/validate_tales.md)
invariant, then those specific to an alignment. As for `tales`, only
properties closed under row subsetting are checked here; grid
completeness is a precondition of the functions that need it.

## Usage

``` r
validate_tales_msa(x)
```

## Arguments

- x:

  A `tales_msa` object.

## Value

`x`, invisibly, if valid; otherwise an error.

## See also

Other TALE alignment:
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.tales_msa.md),
[`is_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/is_tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/dev/reference/tales_align.md),
[`tales_consensus()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus.md),
[`tales_consensus_match()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus_match.md),
[`tales_msa()`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md),
[`tales_width()`](https://scunnac.github.io/tantale/dev/reference/tales_width.md)
