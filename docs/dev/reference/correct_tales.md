# Correct TALE ORFs in error-prone sequences

As a way faster alternative to run
[`tell_tales`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)
in correction mode on error prone sequences such as ONT assembled
genomes, we provide a wrapper around the java binaries from this gitHub
[page](https://github.com/Jstacs/Jstacs/tree/master/projects/talecorrect).
It takes an input fasta file and output a file with corrected indels in
TALE coding sequences using an approach described in the Erkes et al.
[paper](https://doi.org/10.1186/s12864-023-09228-1).
[`tell_tales`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)
can subsequently be run on the corrected sequences in no correction
mode.

## Usage

``` r
correct_tales(
  uncorrected_path,
  corrected_path = file.path(getwd(), "correctedTALEs.fa"),
  hmm_path = system.file("tools", "talecorrect", "HMMs", "Xoo", package = "tantale",
    mustWork = T),
  return_corrections = FALSE,
  conda_bin = "auto"
)
```

## Arguments

- uncorrected_path:

  Path to the input sequence file

- corrected_path:

  Path of the ouput file

- hmm_path:

  Path the folder containning the profile HMM files. The default value
  points to the ones build from Xanthomonas oryzae pv. oryzae templates.
  Xox ones are also available in the parent directory. Please see the
  gitHub
  [page](https://github.com/Jstacs/Jstacs/tree/master/projects/talecorrect)
  for instructions on building custom profiles.

- return_corrections:

  Specify `TRUE` if you want the list of executed operations on the
  sequence as a tibble.

- conda_bin:

  Path to your Conda binary file if you need to specify a path different
  from the one that is automatically searched by the reticulate package
  functions.

## Value

A tibble if `return_corrections` is `TRUE` or the path to the corrected
sequences file.

## See also

Other TALE discovery:
[`tales_from_telltale()`](https://scunnac.github.io/tantale/dev/reference/tales_from_telltale.md),
[`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)

## Examples

``` r
# \donttest{
# Needs nhmmer and a Java runtime.
subj <- system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                    package = "tantale")
out_fa <- tempfile(fileext = ".fa")
correct_tales(uncorrected_path = subj, corrected_path = out_fa)
#> Running nHMMER
#> Performing TALEs cds correction on provided sequences.
#> [1] "/tmp/RtmpZLcluY/file1edd63299b180d.fa"
# }
```
