# Runs the "predict" and "analyze" steps of AnnoTALE on a fasta file

An R wrapper around the
[AnnoTALE](https://www.ncbi.nlm.nih.gov/pubmed/26876161) 'AnnoTALE.jar
predict' and 'AnnoTALE.jar analyze' shell calls. The whole AnnoTALE
workflow can be completed by a subsequent call to the
[`run_annotale_build`](https://scunnac.github.io/tantale/reference/run_annotale_build.md)
function.

## Usage

``` r
run_annotale_predict(
  fasta_file,
  output_dir = getwd(),
  prefix = NULL,
  annotale_jar = system.file("tools", "AnnoTALEcli-1.5.jar", package = "tantale",
    mustWork = T)
)
```

## Arguments

- fasta_file:

  Path to a fasta file containing DNA sequences (e.g. a genome assembly)
  to be analyzed for TALE content.

- output_dir:

  Directory where output will be written (created if it does not exist).

- prefix:

  A scalar character vector containing a prefix that will be appended to
  TALE names by AnnoTALE. If not supplied, the function will try to
  guess the prefix from the input file name.

- annotale_jar:

  Path to the AnnoTALE jar file if you want to use another version than
  the one provided with tantale.

## Value

`0`, invisibly; called for the files it writes to `output_dir`. If
either stage exits with a non-zero status, the function stops with an
error of class `tantale_error_annotale_failed`.

## See also

Other external TALE tools:
[`run_annotale_build()`](https://scunnac.github.io/tantale/reference/run_annotale_build.md)

## Examples

``` r
# \donttest{
# Needs a Java runtime.
fasta <- system.file("extdata", "MAI1.fa", package = "tantale")
out <- file.path(tempdir(), "annotale_predict_example")
run_annotale_predict(fasta, output_dir = out)
#> Running AnnoTALE predict for "MAI1"
#>   java -jar
#>   '/home/cunnac/Lab-Related/MyScripts/tantale/inst/tools/AnnoTALEcli-1.5.jar'
#>   predict g='/home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/MAI1.fa'
#>   s='MAI1' outdir='/tmp/Rtmp5o83zd/annotale_predict_example/Predict'
#> Running AnnoTALE analyze for "MAI1"
#>   java -jar
#>   /home/cunnac/Lab-Related/MyScripts/tantale/inst/tools/AnnoTALEcli-1.5.jar
#>   analyze
#>   t='/tmp/Rtmp5o83zd/annotale_predict_example/Predict/TALE_DNA_sequences_(MAI1).fasta'
#>   outdir='/tmp/Rtmp5o83zd/annotale_predict_example/Analyze'
list.files(file.path(out, "Predict"))
#> [1] "GFF__TALE_predictions_(MAI1).gff3"   "Genbank__TALE_predictions_(MAI1).gb"
#> [3] "TALE_DNA_sequences_(MAI1).fasta"     "TALE_protein_sequences_(MAI1).fasta"
#> [5] "protocol_predict.txt"               
# }
```
