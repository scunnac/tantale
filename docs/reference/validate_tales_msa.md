# Validate a tales_msa object

Checks every
[`validate_tales`](https://scunnac.github.io/tantale/reference/validate_tales.md)
invariant, then those specific to an alignment. As for `tales`, only
properties that still hold on any subset of rows are checked here; grid
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
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/reference/as.matrix.tales_msa.md),
[`is_tales_msa()`](https://scunnac.github.io/tantale/reference/is_tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md),
[`tales_consensus()`](https://scunnac.github.io/tantale/reference/tales_consensus.md),
[`tales_consensus_match()`](https://scunnac.github.io/tantale/reference/tales_consensus_match.md),
[`tales_msa()`](https://scunnac.github.io/tantale/reference/tales_msa.md),
[`tales_msa_width()`](https://scunnac.github.io/tantale/reference/tales_msa_width.md)

## Examples

``` r
aligned <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L),
  alignment_position = c(1L, 2L, 3L, 1L, 3L),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "CTERM")
)
validate_tales_msa(tales_msa(aligned))

# Two parts of A1 placed in the same alignment column
aligned$alignment_position[2] <- 1L
try(tales_msa(aligned))
#> Error in validate_tales_msa(new_tales_msa(x, alignment_width = alignment_width)) : 
#>   alignment_position must be unique within an array.
```
