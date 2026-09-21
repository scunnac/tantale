# Split strings of TALE sequences (\`sep\`-separated rvd or distal repeat IDs)

Load the content of fasta file containing TALE sequences (either RVD or
Distal repeat code) and return a list of vectors each one composed of
the individual elements of the sequence.

## Usage

``` r
toListOfSplitedStr(atomicStrings, sep = "-")
```

## Arguments

- atomicStrings:

  Either, the path to a fasta file, an AAStringSet or "BStringSet" or a
  list. In all cases, each element of these objects is a string of a
  tale sequence (\`sep\`-separated rvd or distal repeat IDs)

- sep:

  Separator of the elements of the sequence

## Value

A list of named vectors representing the 'splited' sequence.
