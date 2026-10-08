# Plot TALE RVD sequences along a potential DNA target region

This function enables visual inspection of TALE target predictions
results that are located **within** a specified subject DNA sequence
region

## Usage

``` r
plot_target_preds(preds, subj_file, filter_range)
```

## Arguments

- preds:

  A tibble of prediction results obtained with
  [`preditale`](https://scunnac.github.io/tantale/reference/preditale.md)
  or [`talvez`](https://scunnac.github.io/tantale/reference/talvez.md)
  or a custom table in this format.

- subj_file:

  The fasta file of subject DNA sequences that was used to predict DNA
  binding elements.

- filter_range:

  A length one genomic ranges in the form of a properly formatted
  character string (eg. "chr2:56-125") or an atomic GenomicRanges
  object. This argument specifies the DNA region that will be plotted
  together with predicted binding TALEs RVD sequences whose predicted
  EBE lies **entirely within**.

## Value

Returns a ggplot object that can be further altered using ggplot2
package functions.

## Details

RVDs sequences predicted to target an EBE on the sense strand of the DNA
sequence are plotted on top of the double stranded DNA sequence in
parallel to its cognate EBE which is highlighted on the corresponding
strand of the DNA sequence. Those predicted to target an EBE on the
opposite strand are displayed below.

Individual RVDs are printed inside colored boxes. The color of the boxes
indicates to which degree the RVD is predicted to have affinity with the
corresponding nucleotide on the DNA sequence at that position relatively
to other nucleotides. The RVDs labelled "OO" correspond to the first
non-canonical repeat also called repeat zero in TALE protein sequences.

Numeric values inside the boxes located immediately to the right of the
TALE labels reflect prediction scores.

## See also

Other target prediction:
[`preditale()`](https://scunnac.github.io/tantale/reference/preditale.md),
[`tales_predict_targets()`](https://scunnac.github.io/tantale/reference/tales_predict_targets.md),
[`talvez()`](https://scunnac.github.io/tantale/reference/talvez.md)

Other TALE plots:
[`plot.tales()`](https://scunnac.github.io/tantale/reference/plot.tales.md),
[`plot.tales_msa()`](https://scunnac.github.io/tantale/reference/plot.tales_msa.md),
[`talomes_heatmap()`](https://scunnac.github.io/tantale/reference/talomes_heatmap.md)

## Examples

``` r
# \donttest{
# Needs a Java runtime, for preditale().
x <- tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
                                      package = "tantale"))
rvds <- tales_rvd_strings(x)
subj <- system.file("extdata", "cladeIII_sweet_promoters.fasta",
                    package = "tantale")
preds <- preditale(rvd_seqs = rvds, subj_file = subj)
best <- preds[order(-preds$score), ][1, ]
# The site with 10 bp on either side
plot_target_preds(preds = best, subj_file = subj,
                  filter_range = paste0(best$subjSeqId, ":",
                                        best$start - 10, "-", best$end + 10))

# }
```
