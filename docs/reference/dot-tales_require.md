# Check a tales against the documented requirements of a consumer

Reads `TALES_REQUIREMENTS` so every function reports a missing column
the same way, rather than each inventing its own guard.

## Usage

``` r
.tales_require(x, fn)
```

## Arguments

- x:

  A `tales` object.

- fn:

  Name of the consumer, as a key of `TALES_REQUIREMENTS`.

## Value

`x`, invisibly. Aborts if a requirement is unmet.
