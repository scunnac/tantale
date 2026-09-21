# Runs the "predict" and "analyze" steps of AnnoTALE on a fasta file.

A R wrapper around the
[AnnoTALE](https://www.ncbi.nlm.nih.gov/pubmed/26876161) 'AnnoTALE.jar
predict' and 'AnnoTALE.jar analyze' shell calls. The whole AnnoTALE
workflow can be completed by a subsequent call to the
[`buildAnnoTALE`](https://scunnac.github.io/tantale/reference/buildAnnoTALE.md)
function.

## Usage

``` r
analyzeAnnoTALE(
  inputFastaFile,
  outputDir = getwd(),
  prefix = NULL,
  annoTALE = system.file("tools", "AnnoTALEcli-1.5.jar", package = "tantale", mustWork =
    T)
)
```

## Arguments

- inputFastaFile:

  Path to a fasta file containing DNA (?) sequences to be analyzed for
  TALE content.

- outputDir:

  Directory where output will be written (created if does not exist).

- prefix:

  A scalar character vector containing a prefix that will be appended to
  TALE names by AnnoTALE. If not supplied, the function will try to
  guess the prefix from the input file name.

- annoTALE:

  Path to the AnnoTALE jar file if you want to use another version than
  the one provided with tantale.

## Value

Returns invisibly the edit code of the shell call to the last AnnoTALE
step (ie '0' if successful).
