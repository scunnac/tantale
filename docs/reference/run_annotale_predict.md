# Runs the "predict" and "analyze" steps of AnnoTALE on a fasta file

An R wrapper around the [AnnoTALE](https://doi.org/10.1038/srep21077)
'AnnoTALE.jar predict' and 'AnnoTALE.jar analyze' shell calls. The whole
AnnoTALE workflow can be completed by a subsequent call to the
[`run_annotale_build`](https://scunnac.github.io/tantale/reference/run_annotale_build.md)
function.

## Usage

``` r
run_annotale_predict(
  fasta_file,
  output_dir = tempfile("annotale_predict_"),
  prefix = NULL,
  opt_param = "Sensitive=false",
  java_args = "",
  annotale_jar = .tantale_tool("annotale")
)
```

## Arguments

- fasta_file:

  Path to a fasta file containing DNA sequences (e.g. a genome assembly)
  to be analyzed for TALE content.

- output_dir:

  Directory where output will be written, created if it does not exist.
  The default is a new directory under
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html), which R deletes
  when the session ends: give a path to keep the results.

- prefix:

  A scalar character vector containing a prefix that will be appended to
  TALE names by AnnoTALE. If not supplied, the function will try to
  guess the prefix from the input file name.

- opt_param:

  A single string of options for the predict stage, as `key=value` pairs
  separated by spaces. AnnoTALE predict has one: `Sensitive`, `false` by
  default; `"Sensitive=true"` runs its sensitive scan. The keys this
  function sets itself (`g`, `s`, `outdir`) are refused: use
  `fasta_file`, `prefix` and `output_dir`. The analyze stage has no
  option of its own beyond a run name.

- java_args:

  A single string of options for the Java virtual machine, placed before
  `-jar` in both stages, such as `"-Xmx8G"` to raise its memory limit.
  The default, `""`, leaves Java's own defaults.

- annotale_jar:

  Path to the AnnoTALE jar file. The default is the copy
  [`tantale_setup()`](https://scunnac.github.io/tantale/reference/tantale_setup.md)
  downloads; give a path to use another version.

## Value

`output_dir`, invisibly, so that the call can be passed straight to
[`tales_from_annotale()`](https://scunnac.github.io/tantale/reference/tales_from_annotale.md).
If either stage exits with a non-zero status, the function stops with an
error of class `tantale_error_annotale_failed`.

## See also

Other external TALE tools:
[`run_annotale_assign()`](https://scunnac.github.io/tantale/reference/run_annotale_assign.md),
[`run_annotale_build()`](https://scunnac.github.io/tantale/reference/run_annotale_build.md),
[`run_annotale_load_classes()`](https://scunnac.github.io/tantale/reference/run_annotale_load_classes.md)

## Examples

``` r
# \donttest{
# Needs a Java runtime, and the tools and genomes tantale_setup() downloads.
fasta <- tantale_genome("MAI1")
out <- file.path(tempdir(), "annotale_predict_example")
run_annotale_predict(fasta, output_dir = out)
#> Running AnnoTALE predict for "MAI1"
#>   java -jar
#>   '/home/cunnac/snap/codium/495/.local/share/R/tantale/tools-1/AnnoTALEcli-1.5.jar'
#>   predict Sensitive=false
#>   g='/home/cunnac/snap/codium/495/.local/share/R/tantale/genomes-1/MAI1.fa'
#>   s='MAI1' outdir='/tmp/RtmpcjOv6d/annotale_predict_example/Predict'
#> Running AnnoTALE analyze for "MAI1"
#>   java -jar
#>   '/home/cunnac/snap/codium/495/.local/share/R/tantale/tools-1/AnnoTALEcli-1.5.jar'
#>   analyze
#>   t='/tmp/RtmpcjOv6d/annotale_predict_example/Predict/TALE_DNA_sequences_(MAI1).fasta'
#>   outdir='/tmp/RtmpcjOv6d/annotale_predict_example/Analyze'
list.files(file.path(out, "Predict"))
#> [1] "GFF__TALE_predictions_(MAI1).gff3"   "Genbank__TALE_predictions_(MAI1).gb"
#> [3] "TALE_DNA_sequences_(MAI1).fasta"     "TALE_protein_sequences_(MAI1).fasta"
#> [5] "protocol_predict.txt"               
# }
```
