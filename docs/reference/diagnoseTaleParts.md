# Report on potential 'pseudo TALEs' in a taleParts object

NOT TESTED!!!! This displays a compact but information rich view of the
TALEs stored in a taleParts object.

## Usage

``` r
diagnoseTaleParts(taleParts, sanitize = FALSE)
```

## Arguments

- taleParts:

  a table of TALE parts as returned by the
  [`getTaleParts`](https://scunnac.github.io/tantale/reference/getTaleParts.md)
  function or
  [`distalr`](https://scunnac.github.io/tantale/reference/distalr.md)

- sanitize:

  If `FALSE`, will return all the arrays with at least one part with a
  missing sequence. If `TRUE`, will return all the arrays that have no
  part with a missing sequence.

## Value

a taleParts object
