# Run TALE target predictions on DNA sequence(s) using PrediTale.

A R wrapper around the
[PrediTale](https://www.jstacs.de/index.php/PrediTALE) 'PrediTALE.jar
preditale' module. Takes a list of TALE RVD sequences and a fasta file
of DNA sequences and runs PrediTale.

## Usage

``` r
preditale(
  rvd_seqs,
  subj_file,
  opt_param = "",
  output_dir = NULL,
  predictor_path = system.file("tools", "PrediTALE.jar", package = "tantale", mustWork =
    T)
)
```

## Arguments

- rvd_seqs:

  Tale RVD sequences are supplied as either a fasta file (atomic
  character vector) with Tale info (name) in title and sequences of RVD
  as a space or '-' separeted string or as a Biostrings XStringSet with
  sequences of RVD similarly formated. See the
  [PrediTale](https://www.jstacs.de/index.php/PrediTALE) man page for
  how to encode RVDs present on aberrant repeats.

- subj_file:

  Expects a character vector specifying the path to the fasta file
  holding subject DNA sequence(s).

- opt_param:

  An atomic character vector specifying optionnal parameters for
  PrediTALE.jar preditale (eg "Strand=\\forward strand\\").

- output_dir:

  Expects a character vector specifying the path to an output directory.
  If not supplied, output files will be temporary.

- predictor_path:

  If you want to use another version of "PrediTALE.jar" than the one
  supplied with tantale, specify its path here.

## Value

A tibble with the EBE predictions. **Note that column names have been
modified** relative to the column names found in the originale
programs's output in order to homogenize column names across TALE target
prediction programs in tantale

## See also

Other target prediction:
[`plot_target_preds()`](https://scunnac.github.io/tantale/dev/reference/plot_target_preds.md),
[`tales_predict_targets()`](https://scunnac.github.io/tantale/dev/reference/tales_predict_targets.md),
[`talvez()`](https://scunnac.github.io/tantale/dev/reference/talvez.md)

## Examples

``` r
# \donttest{
# Needs a Java runtime.
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
rvds <- tales_rvd_strings(x)
subj <- system.file("extdata", "cladeIII_sweet_promoters.fasta",
                    package = "tantale")
preds <- preditale(rvd_seqs = rvds, subj_file = subj)
head(preds)
#> # A tibble: 1 × 9
#>   subjSeqId            start   end strand score ebeSeq         pval rvds  taleId
#>   <chr>                <dbl> <dbl> <chr>  <dbl> <chr>         <dbl> <chr> <chr> 
#> 1 SWEET11p_93-11_Sense   210   236 +      0.242 TATAAAAATG… 5.86e-5 NN-N… ROI_0…
# }
```
