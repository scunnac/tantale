# Generate a mapping between Distal repeat IDs and their cognate RVD

Returns the association between each domain code (`dom_code`) and its
RVD, one row per distinct pair.

## Usage

``` r
repeat_to_rvd_map_distalr(tale_parts)
```

## Arguments

- tale_parts:

  A [`tales`](https://scunnac.github.io/tantale/reference/tales.md)
  object or data frame with `dom_code` and `rvd` columns, such as the
  `tales` element of
  [`tales_compare_distal`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md)'s
  result.

## Value

A two columns repeatID - RVD data frame.

## See also

Other tales projections:
[`repeat_to_rvd_map()`](https://scunnac.github.io/tantale/reference/repeat_to_rvd_map.md),
[`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/reference/tale_parts_to_rvd.md),
[`tales_coded_strings()`](https://scunnac.github.io/tantale/reference/tales_coded_strings.md),
[`tales_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_domain_codes.md),
[`tales_get_dna_seq()`](https://scunnac.github.io/tantale/reference/tales_get_dna_seq.md),
[`tales_get_protein_seq()`](https://scunnac.github.io/tantale/reference/tales_get_protein_seq.md),
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md),
[`tales_to_universalmotif()`](https://scunnac.github.io/tantale/reference/tales_to_universalmotif.md)

## Examples

``` r
parts <- data.frame(dom_code = c(1, 2, 1, 3), rvd = c("HD", "NI", "HD", "NG"))
repeat_to_rvd_map_distalr(parts)
#>   repeatID RVD
#> 1        1  HD
#> 2        2  NI
#> 3        3  NG
```
