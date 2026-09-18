# The domain code lookup table

One row per distinct domain sequence: its `dom_code`, the amino acid
sequence it stands for, and the RVD carried by that domain.

## Usage

``` r
tales_domain_codes(x)
```

## Arguments

- x:

  A [`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
  object carrying a `dom_code` column.

## Value

A [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) with
`dom_code`, `aa_seq` and `rvd` columns, one row per code.

## See also

[`tales_coded_strings`](https://scunnac.github.io/tantale/dev/reference/tales_coded_strings.md)

Other tales projections:
[`repeat_to_rvd_map()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map.md),
[`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map_distalr.md),
[`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/dev/reference/tale_parts_to_rvd.md),
[`tales_coded_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_coded_strings.md),
[`tales_rvd_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_rvd_strings.md)

## Examples

``` r
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
xa <- tales_assign_domain_codes(x)
head(tales_domain_codes(xa))
#> # A tibble: 6 × 3
#>   dom_code aa_seq                             rvd  
#>   <chr>    <chr>                              <chr>
#> 1 1        LPPDQVVAIASNGGGKQALETVQRLLPVLCQAHG NG   
#> 2 10       LTPAQVVAIASNDGGKQALETVQRLLPVLCQAHG ND   
#> 3 11       LTPAQVVAIASNGGGKQALE               NG   
#> 4 12       LTPAQVVAIASNGGGKQALETVQRLLPVLCQAHG NG   
#> 5 13       LTPAQVVAIASNGGGKQALETVQRLLPVLCQARG NG   
#> 6 14       LTPAQVVAIASNGGKQALETVQRLLPVLCQAHG  N*   
```
