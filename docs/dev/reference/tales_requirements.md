# Column requirements of the tales consumers

What each function needs of a
[`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
object. This is distinct from the class's own column contract, which
says what makes an object *valid*: a valid `tales` may still lack what a
given operation needs.

## Usage

``` r
tales_requirements()
```

## Value

A tibble of function, requirement kind, and columns.

## Details

- [`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md):

  `aa_seq` or `dna_seq` (translated)

- [`tales_align()`](https://scunnac.github.io/tantale/dev/reference/tales_align.md):

  whichever of `rvd` / `dom_code` `residue_col` names

- [`tales_rvd_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_rvd_strings.md):

  `rvd`

- [`tales_coded_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_coded_strings.md):

  `dom_code`

- [`tales_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_codes.md):

  `dom_code` and `aa_seq`; `rvd` included when present

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a tales:

  `rvd`, `aa_seq`, `domain_type`; `seqnames` adds a facet,
  `alignment_position` enables the aligned layout

- [`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/dev/reference/tale_parts_to_rvd.md):

  `rvd`

- [`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map_distalr.md):

  `dom_code`, `rvd`

## See also

Other tales objects:
[`as_tales()`](https://scunnac.github.io/tantale/dev/reference/as_tales.md),
[`format.tales()`](https://scunnac.github.io/tantale/dev/reference/format.tales.md),
[`format.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/format.tales_msa.md),
[`is_tales()`](https://scunnac.github.io/tantale/dev/reference/is_tales.md),
[`new_tales()`](https://scunnac.github.io/tantale/dev/reference/new_tales.md),
[`print.tales()`](https://scunnac.github.io/tantale/dev/reference/print.tales.md),
[`print.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/print.tales_msa.md),
[`summary.tales()`](https://scunnac.github.io/tantale/dev/reference/summary.tales.md),
[`summary.tales_msa()`](https://scunnac.github.io/tantale/dev/reference/summary.tales_msa.md),
[`tales()`](https://scunnac.github.io/tantale/dev/reference/tales.md),
[`tales_anchor_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_anchor_codes.md),
[`tales_anomalies()`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md),
[`tales_assert_complete()`](https://scunnac.github.io/tantale/dev/reference/tales_assert_complete.md),
[`tales_namespace()`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md),
[`validate_tales()`](https://scunnac.github.io/tantale/dev/reference/validate_tales.md)
