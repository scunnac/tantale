# Biological anomalies in a tales object

Collects the array-level anomalies that make an object *odd* rather than
*unreadable*: missing sequence data, impossible domain-type
arrangements, coordinate disagreements, an amino acid sequence paired
with more than one RVD, attributes that should be constant within an
array but are not.

## Usage

``` r
.tales_anomalies(x)
```

## Arguments

- x:

  A data frame with at least the tales key columns.

## Value

A tibble with one row per (array, anomaly): `array_id`, `check` and
`detail`. Zero rows if the object is clean.

## Details

These are deliberately **not** errors. Real TALE predictions are messy,
and a class that refuses to load them forces cleaning outside the
package and destroys exactly the diagnostic signal a user wants. They
are reported as a warning on construction and can be removed with
`sanitize = TRUE`.

Structural violations – a duplicated key, an `NA` `array_id`, a broken
`aa_seq`/`dom_code` bijection – are a different matter and remain hard
errors: the table cannot be interpreted at all, and downstream code
would silently compute wrong answers rather than merely odd ones.
