# Build a tales object from a tell_tales run directory

Reads the AnnoTALE/telltale part files of a single
[`tell_tales`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)
output directory and returns a validated
[`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
object.

## Usage

``` r
tales_from_telltale(telltale_dir, sanitize = FALSE)
```

## Arguments

- telltale_dir:

  Path to a single
  [`tell_tales`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)
  output directory.

- sanitize:

  If `TRUE`, arrays carrying biological anomalies are removed with a
  warning naming them and why; if `FALSE` (default) they are kept and
  merely warned about. See
  [`tales_anomalies`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md).

## Value

A validated `tales` object.

## Details

The result carries no `dom_code`: that surrogate key is minted later, by
the relatedness computation, over the whole set of parts being analysed.

## See also

Other TALE discovery:
[`correct_tales()`](https://scunnac.github.io/tantale/dev/reference/correct_tales.md),
[`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)

## Examples

``` r
tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                package = "tantale"))
#> <tales> 4 arrays, 96 parts
#>   layers: rvd   |   6 other columns
#>              rvd
#>   ROI_00001  NTERM NN NG NN HD HD NI N* NG HD NI NG NN HD NI NG NI NG NN N ...
#>   ROI_00002  NTERM NN HD NI NN HD NG HD HD NG NG NI NG NI NG CTERM
#>   ROI_00003  NTERM NN ND NN NI NK NN HD NN NG NG N* HD N* HD NI NN HD NG H ...
#>   ROI_00004  NTERM NI HD NN NS NN NG HD NG HD NG NN NG HD NS HD NI NG HD H ...
```
