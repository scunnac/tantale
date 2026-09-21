# Plot TALE RVD sequences along a potential DNA target region

This function enable visual inspection of TALE target predictions
results that are located **within** a specified subject DNA sequence
region

## Usage

``` r
plotTaleTargetPred(predResults, subjDnaSeqFile, filterRange)
```

## Arguments

- predResults:

  A tibble of prediction results obtained with
  [`preditale`](https://scunnac.github.io/tantale/reference/preditale.md)
  or [`talvez`](https://scunnac.github.io/tantale/reference/talvez.md)
  or a custom table in this format.

- subjDnaSeqFile:

  The fasta file of subject DNA sequences that was used to predict DNA
  binding elements.

- filterRange:

  A length one genomic ranges in the form of a properly formatted
  character string (eg. "chr2:56-125") or an atomic GenomicRanges
  object. This argument specify the DNA region that will be plotted
  together with predicted binding TALEs RVD sequences whose predicted
  EBE lies **entirely whithin**.

## Value

Returns a ggplot object that can be further altered using ggplot2
package functions.

## Details

RVDs sequences predicted to target an EBE on the sense strand of the DNA
sequence are plotted on top of the double stranded DNA sequence in
parallel to its cognate EBE which is highlighted on the corresponding
strand of the DNA sequence. Those preicted to target an EBE on the
opposite strand are displayed below.

Individual RVDs are printed inside colored boxes. The color of the boxes
indicate to which degree the RVD is predicted to have affinity with the
corresponding nucleotide on the DNA sequence at that position relatively
to other nucleotides. The RVDs labelled "OO" correspond to the first
non-canonical repeat also called repeat zero in TALE protein squences.

Numeric values inside the boxes located immediately to the right of the
TALE labels reflect prediction scores.
