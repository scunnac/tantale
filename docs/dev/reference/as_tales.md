# Coerce sequences of TALE parts to a tales object

Turns `sep`-separated TALE sequences — RVD strings, or repeat-code
strings — into a
[`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)
object, taking sequence order as `position_in_array` and sequence names
as `array_id`.

## Usage

``` r
as_tales(x, ...)

# S3 method for class 'data.frame'
as_tales(x, ...)

# Default S3 method
as_tales(x, sep = "-", residue_col = c("rvd", "dom_code"), ...)
```

## Arguments

- x:

  A path to a fasta file, a `BStringSet`/`AAStringSet`, a list of
  strings, or a data frame (which is passed to
  [`tales`](https://scunnac.github.io/tantale/dev/reference/tales.md)).

- ...:

  Passed to methods.

- sep:

  Separator between elements of a sequence. Use `"-"` for RVD sequences
  and `" "` for repeat-code strings.

- residue_col:

  Which residue column the parsed elements become: `"rvd"` (default) or
  `"dom_code"`. Given explicitly rather than guessed from the values.

## Value

A validated `tales` object.

## Details

The result is deliberately column-poor: a bare sequence file carries no
domain types, amino acid sequences or source contigs, so only
`array_id`, `position_in_array` and the chosen residue column are
produced. That is a valid `tales` — see `dev/class-design.md` §2.3 for
why `seqnames` and the rest are optional.

## See also

Other tales objects:
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
[`tales_requirements()`](https://scunnac.github.io/tantale/dev/reference/tales_requirements.md),
[`validate_tales()`](https://scunnac.github.io/tantale/dev/reference/validate_tales.md)

## Examples

``` r
rvd_fasta <- system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",
                         package = "tantale")
as_tales(rvd_fasta, sep = "-")
#> <tales> 11 arrays, 258 parts
#>   layers: rvd
#>                  rvd
#>   TalA_BAI3      NN NG NN HD HD NI N* NG HD NI NG NN HD NI NG NI NG NN NG  ...
#>   TalA_CFBP1947  NN NG NN HD HD NI N* NG HD NI NG NN HD NI NG NI NG NN NG  ...
#>   ...            ...
#>   TalA_MAI95     NN N* NN HD HD NI NG NN HD NS NG NI N* NN NG HD NI NI NG  ...
#>   TalA_MAI99     NN N* NN HD HD NI NG NN HD NS NG NI N* NN NG HD NI NI NG  ...
```
