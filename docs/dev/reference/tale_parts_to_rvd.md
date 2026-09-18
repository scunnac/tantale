# Generates a RVD sequences set from a tale_parts object

Uses a tale_parts object in a
[`tales_compare`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
output to return a
[`BStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
of RVD sequences. RVDs are separated by the character specified in the
`sep` parameter.

## Usage

``` r
tale_parts_to_rvd(tale_parts, sep = "-", rvd_only = FALSE)
```

## Arguments

- tale_parts:

  The tale_parts object in a
  [`tales_compare`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
  output.

- sep:

  Used as a RVD separatator

- rvd_only:

  Retrun only RVDs and ommit N- and C- terminal domains

## Value

A two columns repeatID - RVD data frame.

## Details

Uses Distal repeat sequences and RVD sequences from a set of TALEs
analyzed with the
[`tales_compare`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
function to return the association between repeat ID and RVD.

## See also

Other tales projections:
[`repeat_to_rvd_map()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map.md),
[`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map_distalr.md),
[`tales_coded_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_coded_strings.md),
[`tales_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_codes.md),
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_rvd_strings.md)
