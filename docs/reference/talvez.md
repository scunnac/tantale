# Run TALE target predictions on DNA sequence(s) using Talvez

An R wrapper around the
[Talvez](https://doi.org/10.1371/journal.pone.0068464) predictor Perl
script. Takes a list of TALE RVD sequences and a fasta file of DNA
sequences and runs Talvez.

## Usage

``` r
talvez(
  rvd_seqs,
  subj_file,
  opt_param = "-t 0 -l 19",
  output_dir = NULL,
  talvez_dir = system.file("tools", "TALVEZ_3.2", package = "tantale", mustWork = T),
  conda_bin = "auto"
)
```

## Arguments

- rvd_seqs:

  Tale RVD sequences are supplied as either a fasta file (atomic
  character vector) with Tale info (name) in title and sequences of RVD
  as a space or '-' separated string or as a Biostrings XStringSet with
  sequences of RVD similarly formatted.

- subj_file:

  Expects a character vector specifying the path to the fasta file
  holding subject DNA sequence(s). Talvez does not accept sequences
  wrapped at a fixed width. The function unwraps them with Biostrings,
  but this fails for subject sequences longer than 20 kb: unwrap such
  sequences beforehand.

- opt_param:

  An atomic character vector specifying optional parameters for the
  Talvez script (eg "-t 0 -l 19"). **These may not include** the '-e'
  and '-z' options specifying the matrix files.

- output_dir:

  Expects a character vector specifying the path to an output directory.
  If not supplied, output files will be temporary.

- talvez_dir:

  If you want to use another version of Talvez than the one supplied
  with tantale, specify the path of the directory containing the
  necessary files here.

- conda_bin:

  Path to your Conda binary file if you need to specify a path different
  from the one that is automatically searched by the reticulate package
  functions.

## Value

A tibble with one row per predicted EBE (effector binding element), the
columns renamed from the tool's own output so that `talvez` and
[`preditale`](https://scunnac.github.io/tantale/reference/preditale.md)
return the same ones in the same order: `taleId`, `rvds`, `subjSeqId`,
`start`, `end`, `strand`, `ebeSeq` and `score`, followed by the tool's
own column: `rank`, the site's rank among the TALE's predictions by
score.

## Details

This wrapper uses the RVD-nucleotide specificity matrices shipped with
Talvez; custom matrices are not supported. The Talvez script runs in a
conda environment providing its dependencies, created automatically on
first use.

## See also

Other target prediction:
[`plot_target_preds()`](https://scunnac.github.io/tantale/reference/plot_target_preds.md),
[`preditale()`](https://scunnac.github.io/tantale/reference/preditale.md),
[`tales_predict_targets()`](https://scunnac.github.io/tantale/reference/tales_predict_targets.md)

## Examples

``` r
# \donttest{
# Needs the tantale conda environment, built on first use.
x <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
                                      package = "tantale"))
rvds <- tales_rvd_strings(x)
subj <- system.file("extdata", "cladeIII_sweet_promoters.fasta",
                    package = "tantale")
preds <- talvez(rvd_seqs = rvds, subj_file = subj)
#> Invoking Talvez using the following command:
#> '/home/cunnac/mamba/envs/tantale/bin/perl' TALVEZ_3.2.pl -t 0 -l 19 -e mat1 -z
#> mat2 'rvdSeqsTalvez_a465d15b3a8a9.tsv' 'cladeIII_sweet_promoters.fasta'
head(preds)
#> # A tibble: 6 × 9
#>   taleId    rvds                 subjSeqId start   end strand ebeSeq score  rank
#>   <chr>     <chr>                <chr>     <dbl> <dbl> <chr>  <chr>  <dbl> <dbl>
#> 1 ROI_00001 NN-NG-NN-HD-HD-NI-N… SWEET11p…   980  1006 -      TGTAC… 10.0      1
#> 2 ROI_00001 NN-NG-NN-HD-HD-NI-N… SWEET14p…   160   186 -      TGTTT…  6.35     2
#> 3 ROI_00002 NN-HD-NI-NN-HD-NG-H… SWEET11p…  1395  1409 +      TGTAC…  6.96     1
#> 4 ROI_00002 NN-HD-NI-NN-HD-NG-H… SWEET11p…  1489  1503 +      TCCAG…  6.66     2
#> 5 ROI_00002 NN-HD-NI-NN-HD-NG-H… SWEET11p…   135   149 +      TGCAT…  6.53     3
#> 6 ROI_00002 NN-HD-NI-NN-HD-NG-H… SWEET11p…   130   144 +      TGCAT…  6.53     4
# }
```
