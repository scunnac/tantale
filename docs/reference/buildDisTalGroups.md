# Run Alvaro's Perl scripts for TALEs grouping

take distal output and classify groups

## Usage

``` r
buildDisTalGroups(path, num.groups, overwrite = F)
```

## Arguments

- path:

  directory containing DisTal output files, or the same object as the
  output of `outdir` for
  [`runDistal`](https://scunnac.github.io/tantale/reference/runDistal.md)

- num.groups:

  an integer indicating the number of TALEs groups you want to classify.

- overwrite:

  logical indicating whether to rerun the Perl scripts or only load the
  existing results.

## Value

A list containing:

- SeqOfRvdAlignments: a list of matrices of TALES alignment with RVD
  sequences

- SeqOfDistancesAlignments: a list of matrices of TALEs alignment with
  repeat distance

- repeatUnitsDistanceMatrix: a matrix of pairwise similarity scores
  between repeats

- SeqOfRepsAlignments: a list of matrices of TALEs alignment with repeat
  codes

- TALgroups: a data frame of TALEs names and their groups
