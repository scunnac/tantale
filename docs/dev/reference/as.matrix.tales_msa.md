# Render a TALE alignment as a matrix

Materialises the rectangular form: one row per array, one column per
alignment position, cells holding the requested layer. This is the only
place the gapped matrix is built — the long object is the canonical
storage.

## Usage

``` r
# S3 method for class 'tales_msa'
as.matrix(x, value = NULL, gap = NA, ...)
```

## Arguments

- x:

  A `tales_msa` object.

- value:

  Name of the column to fill cells with. Defaults to the first available
  of `rvd`, `dom_code`.

- gap:

  Value to use for gaps. Defaults to `NA`.

- ...:

  Ignored.

## Value

A character matrix with arrays as rows and alignment positions as
columns.

## See also

Other TALE alignment:
[`is_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/is_tales_msa.md),
[`tales_align()`](https://scunnac.github.io/tantale/dev/reference/tales_align.md),
[`tales_consensus()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus.md),
[`tales_consensus_match()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus_match.md),
[`tales_msa()`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md),
[`tales_width()`](https://scunnac.github.io/tantale/dev/reference/tales_width.md),
[`validate_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/validate_tales_msa.md)

## Examples

``` r
aligned <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2"),
  position_in_array = c(1L, 2L, 3L, 1L, 2L),
  alignment_position = c(1L, 2L, 3L, 1L, 3L),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "CTERM")
)
msa <- tales_msa(aligned)
as.matrix(msa) # gaps as NA
#>    1       2    3      
#> A1 "NTERM" "HD" "CTERM"
#> A2 "NTERM" NA   "CTERM"
as.matrix(msa, gap = "-") # gaps as "-"
#>    1       2    3      
#> A1 "NTERM" "HD" "CTERM"
#> A2 "NTERM" "-"  "CTERM"
```
