# RVD strings, one per TALE array

Renders each array as a separated string of its RVDs in part order — the
form target-prediction tools consume.

## Usage

``` r
tales_rvd_strings(x, sep = "-", rvd_only = TRUE)
```

## Arguments

- x:

  A [`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
  object carrying an `rvd` column.

- sep:

  Separator between RVDs. Defaults to `"-"`, the convention used by
  AnnoTALE and by this package's own sample files.

- rvd_only:

  Drop the terminus parts, whose `rvd` holds an anchor code rather than
  a real RVD (see
  [`tales_anchor_codes`](https://scunnac.github.io/tantale/dev/reference/tales_anchor_codes.md)).
  `TRUE` by default, since target prediction concerns the central repeat
  domain only.

## Value

A
[`BStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html),
named by `array_id`.

## See also

[`tales_coded_strings`](https://scunnac.github.io/tantale/dev/reference/tales_coded_strings.md)

Other tales projections:
[`repeat_to_rvd_map()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map.md),
[`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map_distalr.md),
[`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/dev/reference/tale_parts_to_rvd.md),
[`tales_coded_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_coded_strings.md),
[`tales_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_codes.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
tales_rvd_strings(x)[1]
#> BStringSet object of length 1:
#>     width seq                                               names               
#> [1]    77 NN-NG-NN-HD-HD-NI-N*-NG...NG-HD-NI-NI-NG-HD-NN-NG ROI_00001
tales_rvd_strings(x, rvd_only = FALSE)[1] # keeps NTERM/CTERM markers
#> BStringSet object of length 1:
#>     width seq                                               names               
#> [1]    89 NTERM-NN-NG-NN-HD-HD-NI...NI-NI-NG-HD-NN-NG-CTERM ROI_00001
```
