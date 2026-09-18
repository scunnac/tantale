# Create a tales_msa object

A [`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
object plus an `alignment_position` column. Gaps are implicit: a gap is
simply the absence of a row at that (`array_id`, `alignment_position`).

## Usage

``` r
tales_msa(x, alignment_width = NULL, dom_code_namespace = NULL)
```

## Arguments

- x:

  A data frame with the
  [`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
  columns plus `alignment_position`.

- alignment_width:

  Integer alignment width; defaults to `max(alignment_position)`.

- dom_code_namespace:

  Optional scalar string, see
  [`tales_namespace`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md).

## Value

A validated `tales_msa` object.

## Details

Normally produced by
[`tales_align`](https://scunnac.github.io/tantale/dev/reference/tales_align.md)
rather than called directly.

## See also

Other TALE alignment:
[`as.matrix.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/as.matrix.tales_msa.md),
[`is_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/is_tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/dev/reference/tales_align.md),
[`tales_consensus()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus.md),
[`tales_consensus_match()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus_match.md),
[`tales_width()`](https://scunnac.github.io/tantale/dev/reference/tales_width.md),
[`validate_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/validate_tales_msa.md)

## Examples

``` r
# A hand-built alignment: A2 has a gap at position 2 relative to A1.
aligned <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L),
  alignment_position = c(1L, 2L, 3L, 1L, 3L),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "CTERM")
)
msa <- tales_msa(aligned)
tales_width(msa)
#> [1] 3
as.matrix(msa)
#>    1       2    3      
#> A1 "NTERM" "HD" "CTERM"
#> A2 "NTERM" NA   "CTERM"
```
