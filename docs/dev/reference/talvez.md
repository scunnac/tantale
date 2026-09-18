# Run TALE target predictions on DNA sequence(s) using Talvez

A R wrapper around the
[Talvez](https://doi.org/10.1371/journal.pone.0068464) predictor perl
script. Takes a list of TALE RVD sequences and a fasta file of DNA
sequences and runs Talvez

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
  as a space or '-' separeted string or as a Biostrings XStringSet with
  sequences of RVD similarly formated.

- subj_file:

  Expects a character vector specifying the path to the fasta file
  holding subject DNA sequence(s). Talvez forbids to have sequences in
  the file wrapped at a fixed width. The function uses Biostrings to
  unwrap them but if you use subject sequences longer than 20kb, this
  will fail and you are advised to unwrap your sequences before hand.

- opt_param:

  An atomic character vector specifying optionnal parameters for the
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

A tibble with the EBE predictions. **Note that column names have been
modified** relative to the column names found in the originale
programs's output in order to homogenize column names across TALE target
prediction programs in tantale.

## Details

Note that this talvez wrapper, uses the RVD - Nucleotide specificity
matrices used with talvez, it is not possible to use custom ones. Note
also that the talvez script is run in a conda environment providing the
necessary dependencies. This environment will be created automatically
if necessary.

## See also

Other target prediction:
[`plot_target_preds()`](https://scunnac.github.io/tantale/dev/reference/plot_target_preds.md),
[`preditale()`](https://scunnac.github.io/tantale/dev/reference/preditale.md),
[`tales_predict_targets()`](https://scunnac.github.io/tantale/dev/reference/tales_predict_targets.md)

## Examples

``` r
# \donttest{
# Needs the tantale conda environment, built on first use.
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
rvds <- tales_rvd_strings(x)
subj <- system.file("extdata", "cladeIII_sweet_promoters.fasta",
                    package = "tantale")
preds <- talvez(rvd_seqs = rvds, subj_file = subj)
#> Invoking Talvez using the following command:
#> '/home/cunnac/mamba/envs/tantale/bin/perl' TALVEZ_3.2.pl -t 0 -l 19 -e mat1 -z
#> mat2 'rvdSeqsTalvez_1edd632aab20e9.tsv' 'cladeIII_sweet_promoters.fasta'
head(preds)
#> # A tibble: 6 × 9
#>   taleId    rvds                 subjSeqId score strand start   end ebeSeq  rank
#>   <chr>     <chr>                <chr>     <dbl> <chr>  <dbl> <dbl> <chr>  <dbl>
#> 1 ROI_00001 NN-NG-NN-HD-HD-NI-N… SWEET11p… 10.0  -        980  1006 TGTAC…     1
#> 2 ROI_00001 NN-NG-NN-HD-HD-NI-N… SWEET14p…  6.35 -        160   186 TGTTT…     2
#> 3 ROI_00002 NN-HD-NI-NN-HD-NG-H… SWEET11p…  6.96 +       1395  1409 TGTAC…     1
#> 4 ROI_00002 NN-HD-NI-NN-HD-NG-H… SWEET11p…  6.66 +       1489  1503 TCCAG…     2
#> 5 ROI_00002 NN-HD-NI-NN-HD-NG-H… SWEET11p…  6.53 +        135   149 TGCAT…     3
#> 6 ROI_00002 NN-HD-NI-NN-HD-NG-H… SWEET11p…  6.53 +        130   144 TGCAT…     4
# }
```
