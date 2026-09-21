# Correct TALE ORFs in error-prone sequences

As a way faster alternative to run
[`tellTale`](https://scunnac.github.io/tantale/reference/tellTale.md) in
correction mode on error prone sequences such as ONT assembled genomes,
we provide a wrapper around the java binaries from this gitHub
[page](https://github.com/Jstacs/Jstacs/tree/master/projects/talecorrect).
It takes an input fasta file and output a file with corrected indels in
TALE coding sequences using an approach described in the Erkes et al.
[paper](https://doi.org/10.1186/s12864-023-09228-1).
[`tellTale`](https://scunnac.github.io/tantale/reference/tellTale.md)
can subsequently be run on the corrected sequences in no correction
mode.

## Usage

``` r
correcTales(
  uncorrectedAssemblyPath,
  correctedAssemblyPath = file.path(getwd(), "correctedTALEs.fa"),
  pathToHMMs = system.file("tools", "talecorrect", "HMMs", "Xoo", package = "tantale",
    mustWork = T),
  returnCorrectionsTble = FALSE,
  condaBinPath = "auto"
)
```

## Arguments

- uncorrectedAssemblyPath:

  Path to the input sequence file

- correctedAssemblyPath:

  Path of the ouput file

- pathToHMMs:

  Path the folder containning the profile HMM files. The default value
  points to the ones build from Xanthomonas oryzae pv. oryzae templates.
  Xox ones are also available in the parent directory. Please see the
  gitHub
  [page](https://github.com/Jstacs/Jstacs/tree/master/projects/talecorrect)
  for instructions on building custom profiles.

- returnCorrectionsTble:

  Specify `TRUE` if you want the list of executed operations on the
  sequence as a tibble.

- condaBinPath:

  Path to your Conda binary file if you need to specify a path different
  from the one that is automatically searched by the reticulate package
  functions.

- tellTaleOutDir:

  Path to a
  [`tellTale`](https://scunnac.github.io/tantale/reference/tellTale.md)
  run output directory

## Value

A tibble if `returnCorrectionsTble` is `TRUE` or the path to the
corrected sequences file.
