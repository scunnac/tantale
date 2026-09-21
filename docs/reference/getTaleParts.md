# Fetch Annotale parts from a tellTale output directory.

This function get sequences from tellTale "rvdSequences.fas" and
AnnoTALE "TALE_Protein_parts.fasta" and "TALE_DNA_parts.fasta" files
from a SINGLE
[`tellTale`](https://scunnac.github.io/tantale/reference/tellTale.md)
run output directory and returns a tibble. Each row describes a domain
from a tale array and includes the 'arrayID', the id of the sequence
where this array was found, the 'domainType' (type of domain, repeat,
N-term or C-term), the position of the domain inside the array, the DNA
of the corresponding domain and the RVD and amino acid sequences if
relevant.

\*\*IMPORTANT\*\*: telltale MUST have been run with the
appendExtremityCodes = TRUE

## Usage

``` r
getTaleParts(tellTaleOutDir)
```

## Arguments

- tellTaleOutDir:

  Path to a
  [`tellTale`](https://scunnac.github.io/tantale/reference/tellTale.md)
  run output directory

## Value

A tibble.
