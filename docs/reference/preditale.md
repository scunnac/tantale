# Run TALE target predictions on DNA sequence(s) using PrediTale

An R wrapper around the
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
  java_args = "-Xms512M -Xmx2G",
  predictor_path = .tantale_tool("preditale")
)
```

## Arguments

- rvd_seqs:

  Tale RVD sequences are supplied as either a fasta file (atomic
  character vector) with Tale info (name) in title and sequences of RVD
  as a space or '-' separated string or as a Biostrings XStringSet with
  sequences of RVD similarly formatted. See the
  [PrediTale](https://www.jstacs.de/index.php/PrediTALE) man page for
  how to encode RVDs present on aberrant repeats.

- subj_file:

  Expects a character vector specifying the path to the fasta file
  holding subject DNA sequence(s).

- opt_param:

  A single string of options for PrediTALE, as `key=value` pairs
  separated by spaces, such as `"Strand=\"forward strand\""`. The
  default, `""`, keeps PrediTALE's own: target sites on both strands
  (`Strand`) with a penalty of 0.01 on the reverse one (`r`), and a
  prediction threshold (`t`) from a significance level of 1e-4 (`sl`),
  estimated on a sub-sample of the subject sequences (`b`). The
  [PrediTALE](https://www.jstacs.de/index.php/PrediTALE) page lists the
  alternatives (a number of expected sites, `n`; dedicated background
  sequences, `bs`). The keys this function sets itself (`TALEs`, `s`,
  `outdir`) are refused: use `rvd_seqs`, `subj_file` and `output_dir`.

- output_dir:

  Expects a character vector specifying the path to an output directory.
  If not supplied, output files will be temporary.

- java_args:

  A single string of options for the Java virtual machine, placed before
  `-jar`. The default starts Java with 512 MB of memory and lets it grow
  to 2 GB; raise `-Xmx` if PrediTALE runs out of memory on a large
  subject.

- predictor_path:

  Path to "PrediTALE.jar". The default is the copy
  [`tantale_setup()`](https://scunnac.github.io/tantale/reference/tantale_setup.md)
  downloads; give a path to use another version.

## Value

A tibble with the EBE predictions. **Note that column names have been
modified** relative to the column names found in the original program's
output in order to homogenize column names across TALE target prediction
programs in tantale

## See also

Other target prediction:
[`plot_target_preds()`](https://scunnac.github.io/tantale/reference/plot_target_preds.md),
[`tales_predict_targets()`](https://scunnac.github.io/tantale/reference/tales_predict_targets.md),
[`talvez()`](https://scunnac.github.io/tantale/reference/talvez.md)

## Examples

``` r
# \donttest{
# Needs a Java runtime.
x <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
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
