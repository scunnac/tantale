# Generate a mapping between Distal repeat IDs and their cognate RVD.

Uses Distal repeat sequences and RVD sequences from a set of TALEs
analyzed with the
[`tales_compare`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
function to return the association between repeat ID and RVD.

## Usage

``` r
repeat_to_rvd_map_distalr(tale_parts)
```

## Arguments

- tale_parts:

  The tale_parts object in a
  [`tales_compare`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
  output.

## Value

A two columns repeatID - RVD data frame.

## See also

Other tales projections:
[`repeat_to_rvd_map()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map.md),
[`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/dev/reference/tale_parts_to_rvd.md),
[`tales_coded_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_coded_strings.md),
[`tales_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_codes.md),
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_rvd_strings.md)
