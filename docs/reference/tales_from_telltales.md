# Build a tales object from a tell_tales run directory

Reads the AnnoTALE/telltale part files of a single
[`tell_tales`](https://scunnac.github.io/tantale/reference/tell_tales.md)
output directory and returns a validated
[`tales`](https://scunnac.github.io/tantale/reference/tales.md) object.

## Usage

``` r
tales_from_telltales(telltale_dir, sanitize = FALSE)
```

## Arguments

- telltale_dir:

  Path to a single
  [`tell_tales`](https://scunnac.github.io/tantale/reference/tell_tales.md)
  output directory.

- sanitize:

  If `TRUE`, arrays carrying biological anomalies are removed with a
  warning naming them and why; if `FALSE` (default) they are kept and
  merely warned about. See
  [`tales_anomalies`](https://scunnac.github.io/tantale/reference/tales_anomalies.md).

## Value

A validated `tales` object.

## Details

One row per part: the N-terminus, each repeat and the C-terminus, as
AnnoTALE split the array's longest ORF. The `rvd` column holds the RVD
of a repeat, or a terminus code (see
[`tales_anchor_codes`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md)):
`NTERM`/`CTERM` when the terminus matches the TALE terminal-domain
protein profile, `XXXXX` when it does not. When AnnoTALE reported no
terminus on one side, the array has no part there, with a warning. An
array whose protein and DNA parts disagree is left out, with a warning.

A candidate array with a TALE terminus DNA hit that AnnoTALE could not
split into parts is absent from the result, with a warning naming it.
[`tales_anomalies`](https://scunnac.github.io/tantale/reference/tales_anomalies.md)
cannot report such an array, since it is not in the object. After
frameshift correction this usually means the array was corrected against
a distant reference; see `max_comparisons` in
[`tell_tales`](https://scunnac.github.io/tantale/reference/tell_tales.md).

The directory must have been written by the current version of
[`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md),
whose `array_report.tsv` holds the terminus check; an older one is an
error.

The result carries no `dom_code`: that surrogate key is minted later, by
the relatedness computation, over the whole set of parts being analysed.

## See also

Other TALE discovery:
[`correct_tales()`](https://scunnac.github.io/tantale/reference/correct_tales.md),
[`tale_annotations`](https://scunnac.github.io/tantale/reference/tale_annotations.md),
[`tales_from_annotale()`](https://scunnac.github.io/tantale/reference/tales_from_annotale.md),
[`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md)

## Examples

``` r
tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
                                 package = "tantale"))
#> <tales> 4 arrays, 96 parts
#>   layers: rvd   |   6 other columns
#>              rvd
#>   ROI_00001  NTERM NN NG NN HD HD NI N* NG HD NI NG NN HD NI NG NI NG NN N ...
#>   ROI_00002  NTERM NN HD NI NN HD NG HD HD NG NG NI NG NI NG CTERM
#>   ROI_00003  NTERM NN ND NN NI NK NN HD NN NG NG N* HD N* HD NI NN HD NG H ...
#>   ROI_00004  NTERM NI HD NN NS NN NG HD NG HD NG NN NG HD NS HD NI NG HD H ...
```
