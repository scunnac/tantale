# Generates a RVD sequences set from a tale_parts object

Uses a tale_parts object in a
[`tales_compare_distal`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md)
output to return a
[`BStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
of RVD sequences, one per array, ordered by `position_in_array` and
joined with `sep`.

## Usage

``` r
tale_parts_to_rvd(tale_parts, sep = "-", rvd_only = FALSE)
```

## Arguments

- tale_parts:

  The tale_parts object in a
  [`tales_compare_distal`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md)
  output.

- sep:

  Separator joining consecutive RVDs within an array's string.

- rvd_only:

  Return only RVDs and omit N- and C-terminal domains.

## Value

A
[`BStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html),
one element per array, named by `array_id`.

## See also

Other tales projections:
[`repeat_to_rvd_map()`](https://scunnac.github.io/tantale/reference/repeat_to_rvd_map.md),
[`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/reference/repeat_to_rvd_map_distalr.md),
[`tales_coded_strings()`](https://scunnac.github.io/tantale/reference/tales_coded_strings.md),
[`tales_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_domain_codes.md),
[`tales_get_dna_seq()`](https://scunnac.github.io/tantale/reference/tales_get_dna_seq.md),
[`tales_get_protein_seq()`](https://scunnac.github.io/tantale/reference/tales_get_protein_seq.md),
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md),
[`tales_to_universalmotif()`](https://scunnac.github.io/tantale/reference/tales_to_universalmotif.md)

## Examples

``` r
parts <- data.frame(
  array_id = c("A1", "A1", "A1", "A2", "A2"),
  position_in_array = c(1, 2, 3, 1, 2),
  rvd = c("NTERM", "HD", "CTERM", "NTERM", "NI")
)
tale_parts_to_rvd(parts)
#> BStringSet object of length 2:
#>     width seq                                               names               
#> [1]    14 NTERM-HD-CTERM                                    A1
#> [2]     8 NTERM-NI                                          A2
```
