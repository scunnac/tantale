# Structural part of the CRD contract: the coordinate must be a key

Placement (NA on exactly the non-repeat parts) and order agreement are
*biological* checks and live in
[`.tales_anomalies()`](https://scunnac.github.io/tantale/dev/reference/dot-tales_anomalies.md);
only uniqueness is structural, since a repeated coordinate makes the
array unindexable.

## Usage

``` r
.tales_check_crd_unique(x)
```
