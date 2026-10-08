# Width of a TALE alignment

The number of columns in the alignment, including those that are all
gaps in the object at hand. It is stored with the object, because
subsetting arrays can empty the last column, which would silently shrink
`max(alignment_position)`. A plain `tales` has no width.

## Usage

``` r
tales_msa_width(x)
```

## Arguments

- x:

  A `tales_msa` object.

## Value

An integer scalar, or `NULL` if unset.

## See also

Other TALE alignment:
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/reference/as.matrix.tales_msa.md),
[`is_tales_msa()`](https://scunnac.github.io/tantale/reference/is_tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md),
[`tales_consensus()`](https://scunnac.github.io/tantale/reference/tales_consensus.md),
[`tales_consensus_match()`](https://scunnac.github.io/tantale/reference/tales_consensus_match.md),
[`tales_msa()`](https://scunnac.github.io/tantale/reference/tales_msa.md),
[`tales_msa_params()`](https://scunnac.github.io/tantale/reference/tales_msa_params.md),
[`validate_tales_msa()`](https://scunnac.github.io/tantale/reference/validate_tales_msa.md)

## Examples

``` r
# A2 has no C-terminus, so the last column holds A1's alone.
aligned <- data.frame(
  array_id = c("A1", "A1", "A1", "A1", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 4L, 1L, 2L),
  alignment_position = c(1L, 2L, 3L, 4L, 1L, 2L),
  rvd = c("NTERM", "HD", "NI", "CTERM", "NTERM", "HD")
)
msa <- tales_msa(aligned)
tales_msa_width(msa)
#> [1] 4
# Keeping A2 alone empties columns 3 and 4; the width is still 4.
a2 <- msa[msa$array_id == "A2", ]
max(a2$alignment_position)
#> [1] 2
tales_msa_width(a2)
#> [1] 4
```
