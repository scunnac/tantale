# Is this a tales_msa object?

Is this a tales_msa object?

## Usage

``` r
is_tales_msa(x)
```

## Arguments

- x:

  An object.

## Value

A logical scalar.

## See also

Other TALE alignment:
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/reference/as.matrix.tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md),
[`tales_consensus()`](https://scunnac.github.io/tantale/reference/tales_consensus.md),
[`tales_consensus_match()`](https://scunnac.github.io/tantale/reference/tales_consensus_match.md),
[`tales_msa()`](https://scunnac.github.io/tantale/reference/tales_msa.md),
[`tales_msa_width()`](https://scunnac.github.io/tantale/reference/tales_msa_width.md),
[`validate_tales_msa()`](https://scunnac.github.io/tantale/reference/validate_tales_msa.md)

## Examples

``` r
aligned <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L),
  alignment_position = c(1L, 2L, 3L, 1L, 3L),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "CTERM")
)
msa <- tales_msa(aligned)
is_tales_msa(msa)
#> [1] TRUE
is_tales_msa(as_tales(msa)) # back to a plain tales
#> [1] FALSE
```
