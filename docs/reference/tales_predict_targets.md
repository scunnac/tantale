# Predict TALE target boxes

Runs a TALE target (EBE) prediction tool over a set of TALE RVD
sequences and a set of DNA sequences, and returns the predictions as a
tibble.

## Usage

``` r
tales_predict_targets(x, subj_file, method = c("talvez", "preditale"), ...)
```

## Arguments

- x:

  TALE RVD sequences: a
  [`tales`](https://scunnac.github.io/tantale/reference/tales.md)
  object, the path to a fasta file, or a `BStringSet`. A `tales` is
  rendered with
  [`tales_rvd_strings`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md),
  which drops the termini.

- subj_file:

  Path to a fasta file of subject DNA sequence(s).

- method:

  Which backend to use, `"talvez"` (default) or `"preditale"`.

- ...:

  Passed to the chosen backend. See
  [`talvez`](https://scunnac.github.io/tantale/reference/talvez.md)
  (`opt_param`, `talvez_dir`, `conda_bin`) and
  [`preditale`](https://scunnac.github.io/tantale/reference/preditale.md)
  (`opt_param`, `predictor_path`); both accept `output_dir`. Note the
  two take different `opt_param` defaults, since the options are the
  tools' own.

## Value

A tibble of EBE predictions, with column names homogenised across
backends, plus a `method` column recording which tool produced them.

## Details

`talvez` and `preditale` are independent programs, but the package
normalises their outputs to a shared set of column names, so they can
serve as interchangeable backends of one operation. This is that
operation;
[`talvez`](https://scunnac.github.io/tantale/reference/talvez.md) and
[`preditale`](https://scunnac.github.io/tantale/reference/preditale.md)
remain available and documented individually, and are where each tool's
own options and citation live.

## See also

[`talvez`](https://scunnac.github.io/tantale/reference/talvez.md),
[`preditale`](https://scunnac.github.io/tantale/reference/preditale.md)

Other target prediction:
[`plot_target_preds()`](https://scunnac.github.io/tantale/reference/plot_target_preds.md),
[`preditale()`](https://scunnac.github.io/tantale/reference/preditale.md),
[`talvez()`](https://scunnac.github.io/tantale/reference/talvez.md)

## Examples

``` r
# \donttest{
# Needs the tantale conda environment (for talvez) and a Java runtime
# (for preditale), both built on first use.
x <- tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
                                     package = "tantale"))
subj <- system.file("extdata", "cladeIII_sweet_promoters.fasta",
                    package = "tantale")
head(tales_predict_targets(x, subj_file = subj))
#> Invoking Talvez using the following command:
#> '/home/cunnac/mamba/envs/tantale/bin/perl' TALVEZ_3.2.pl -t 0 -l 19 -e mat1 -z
#> mat2 'rvdSeqsTalvez_67778479ceeac.tsv' 'cladeIII_sweet_promoters.fasta'
#> # A tibble: 6 × 10
#>   taleId    rvds          subjSeqId score strand start   end ebeSeq  rank method
#>   <chr>     <chr>         <chr>     <dbl> <chr>  <dbl> <dbl> <chr>  <dbl> <chr> 
#> 1 ROI_00001 NN-NG-NN-HD-… SWEET11p… 10.0  -        980  1006 TGTAC…     1 talvez
#> 2 ROI_00001 NN-NG-NN-HD-… SWEET14p…  6.35 -        160   186 TGTTT…     2 talvez
#> 3 ROI_00002 NN-HD-NI-NN-… SWEET11p…  6.96 +       1395  1409 TGTAC…     1 talvez
#> 4 ROI_00002 NN-HD-NI-NN-… SWEET11p…  6.66 +       1489  1503 TCCAG…     2 talvez
#> 5 ROI_00002 NN-HD-NI-NN-… SWEET11p…  6.53 +        135   149 TGCAT…     3 talvez
#> 6 ROI_00002 NN-HD-NI-NN-… SWEET11p…  6.53 +        130   144 TGCAT…     4 talvez
head(tales_predict_targets(x, subj_file = subj, method = "preditale"))
#> # A tibble: 1 × 10
#>   subjSeqId          start   end strand score ebeSeq    pval rvds  taleId method
#>   <chr>              <dbl> <dbl> <chr>  <dbl> <chr>    <dbl> <chr> <chr>  <chr> 
#> 1 SWEET11p_93-11_Se…   210   236 +      0.242 TATAA… 5.86e-5 NN-N… ROI_0… predi…
# }
```
