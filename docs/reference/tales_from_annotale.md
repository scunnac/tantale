# Build a tales object from AnnoTALE's own TALE predictions

Reads the output of AnnoTALE's analyze stage, as written by
[`run_annotale_predict`](https://scunnac.github.io/tantale/reference/run_annotale_predict.md),
and returns a validated
[`tales`](https://scunnac.github.io/tantale/reference/tales.md) object
with the same columns as
[`tales_from_telltales`](https://scunnac.github.io/tantale/reference/tales_from_telltales.md).

## Usage

``` r
tales_from_annotale(
  annotale_dir,
  terminus_max_evalue = 1e-05,
  sanitize = FALSE,
  hmm_dir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = TRUE)
)
```

## Arguments

- annotale_dir:

  A directory holding AnnoTALE analyze's `TALE_Protein_parts.fasta`,
  `TALE_DNA_parts.fasta` and `TALE_RVDs.fasta`, in it or in a
  subdirectory: the `output_dir` of
  [`run_annotale_predict()`](https://scunnac.github.io/tantale/reference/run_annotale_predict.md)
  will do.

- terminus_max_evalue:

  Maximum `hmmsearch` E-value for a terminal segment to be coded
  `NTERM`/`CTERM`.

- sanitize:

  If `TRUE`, arrays carrying biological anomalies are removed with a
  warning naming them and why; if `FALSE` (default) they are kept and
  merely warned about. See
  [`tales_anomalies`](https://scunnac.github.io/tantale/reference/tales_anomalies.md).

- hmm_dir:

  Directory holding the TALE terminus protein profiles.

## Value

A validated `tales` object.

## Details

AnnoTALE predict finds TALE genes in a genome; analyze splits each one
into its N-terminal region, repeats and C-terminal region, and reads the
RVD of every repeat.
[`tell_tales`](https://scunnac.github.io/tantale/reference/tell_tales.md)
runs analyze too, on ORFs it delimits itself from nhmmer hits of the
TALE DNA profiles, so the two need not report the same set of TALEs.

As in
[`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md),
each terminal segment is searched with the TALE N- or C-terminal protein
profile (`hmmsearch`, from the tantale environment; see
[`tantale_setup`](https://scunnac.github.io/tantale/reference/tantale_setup.md)).
The `rvd` column holds `NTERM`/`CTERM` for a segment that matches its
profile and `XXXXX` for one that does not (see
[`tales_anchor_codes`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md)).

`array_id` is AnnoTALE's name for the TALE (`MAI1-tempTALE1`), without
the location AnnoTALE appends to it. `seqnames` is filled when predict's
GFF3 file is found under `annotale_dir`.

## See also

Other TALE discovery:
[`correct_tales()`](https://scunnac.github.io/tantale/reference/correct_tales.md),
[`tales_from_telltales()`](https://scunnac.github.io/tantale/reference/tales_from_telltales.md),
[`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md)

## Examples

``` r
# \donttest{
# Needs the tantale environment for hmmsearch (see tantale_setup()).
tales_from_annotale(system.file("extdata", "annotaleExampleOutput",
                                package = "tantale"))
#> <tales> 4 arrays, 96 parts
#>   layers: rvd   |   6 other columns
#>                                 rvd
#>   bai3_sample_tal_genomic_regions-tempTALE1  NTERM NN NG NN HD HD NI N* NG HD NI NG NN  ...
#>   bai3_sample_tal_genomic_regions-tempTALE2  NTERM NN HD NI NN HD NG HD HD NG NG NI NG  ...
#>   bai3_sample_tal_genomic_regions-tempTALE3  NTERM NN ND NN NI NK NN HD NN NG NG N* HD  ...
#>   bai3_sample_tal_genomic_regions-tempTALE4  NTERM NI HD NN NS NN NG HD NG HD NG NN NG  ...
# }
```
