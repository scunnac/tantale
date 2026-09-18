# Generate a mapping between Distal repeat IDs and their cognate RVD.

Uses Distal repeat sequences and RVD sequences from a set of TALEs to
return the association between repeat ID and RVD.

## Usage

``` r
repeat_to_rvd_map(repeat_vecs, rvd_vecs)
```

## Arguments

- repeat_vecs:

  Expects a list of Distal RVDs character vectors. Each **named**
  element corresponding to a TALE.

- rvd_vecs:

  A named list of RVD vectors, one per array, parallel to `repeat_vecs`.

## Value

A two columns repeatID - RVD data frame.

## Details

Care must be taken that TALEs in the two sets of sequences have the same
name. In addition, the function tries hard to make sure that the two
sets of sequences are identical in every ways but the individual
'values' they contain. It is therefore notably important to make sure
that the sequences are consistent in whether they include N-term and
C-term domains IDs/Tags or not.

## See also

Other tales projections:
[`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map_distalr.md),
[`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/dev/reference/tale_parts_to_rvd.md),
[`tales_coded_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_coded_strings.md),
[`tales_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_codes.md),
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_rvd_strings.md)
