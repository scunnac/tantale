# Search and report on the features of TALE protein domains potentially encoded in subject DNA sequences

`tell_tales` has been primarily written to report on 'corrected' TALE
RVD sequences in indels prone, noisy DNA sequences (suboptimally
polished genomes assembly, raw reads of long read sequencing
technologies such as PacBio or ONT) that would otherwise be missed by
conventional tools (eg AnnoTALE).

## Usage

``` r
tell_tales(
  subject_file,
  output_dir = getwd(),
  hmm_dir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T),
  nterm_min_score = 300,
  repeat_min_score = 20,
  cterm_min_score = 200,
  min_domain_hits = 4,
  min_array_length = 0,
  merge_hits = TRUE,
  min_gap = 35,
  extremity_codes = TRUE,
  rvd_sep = "-",
  hmmer_path = NULL,
  extend_len = 300,
  correct_array = FALSE,
  correction_ref = system.file("extdata", "tale_correction_ref.fa.gz", package =
    "tantale", mustWork = T),
  max_comparisons = NULL,
  frameshift = -11,
  ...
)
```

## Arguments

- subject_file:

  Fasta file with DNA sequence(s) to be searched for the presence of
  TALE coding sequences (CDS).

- output_dir:

  Path of the output directory. If not specified, results will be
  written to current working folder.

- hmm_dir:

  Specify the path to a folder holding the hmmfiles if you do not want
  to use the ones provided with tantale.

- nterm_min_score:

  Minimal nhmmer score cut_off value to consider the hit as genuine

- repeat_min_score:

  Minimal nhmmer score cut_off value to consider the hit as genuine

- cterm_min_score:

  Minimal nhmmer score cut_off value to consider the hit as genuine

- min_domain_hits:

  Minimum number of nhmmer hits for a subject sequence (a contig, a
  chromosome) to be considered further. A cheap way to discard whole
  sequences that carry nothing but stray matches, before any expensive
  work is done on them. It says nothing about the length of the TALE
  arrays found within a sequence that passes – see `min_array_length`
  for that.

- min_array_length:

  Minimum number of **repeat units** for a TALE array to be kept.
  Defaults to `0`, which keeps everything.

  Counting repeats rather than all hits means an array is not penalised
  for having had its termini missed, and matches what "array length"
  usually means for a TALE: the number of repeats is what determines how
  long a target box it recognises.

  Whether a short array is noise or a genuinely truncated TALE is a
  judgement about the biology, which is why nothing is discarded unless
  you ask. A pseudogene with three surviving repeats is real, and may be
  what you are looking for.

- merge_hits:

  Perform overlapping hits merging per domain type. Should not be
  modified.

- min_gap:

  Minimum gap in base pairs between two tale domain hits for them to be
  considered distinct. If the length of the gap is below this value,
  domains are considered "contiguous" and grouped in the same array.

- extremity_codes:

  Set this to `FALSE` if you do not want the N- and C-TREM anchor codes
  in the output sequences of RVD

- rvd_sep:

  Symbol acting as a separator in RVD sequences

- hmmer_path:

  Specify the path to a directory holding the HMMER executable if you do
  not want to use the ones provided with tantale.

- extend_len:

  number of nucleotides to extend in 3'-end at the tal ORF prediction
  stage.

- correct_array:

  True or False

- correction_ref:

  Fasta of reference TALE proteins to correct against. Two are shipped,
  both built by `data-raw/make_correction_references.R` from the same
  source:

  - `tale_correction_ref.fa.gz` (default, 494 sequences) – every
    distinct TALE protein of at least 300 aa found across 70
    *Xanthomonas oryzae* genomes.

  - `tale_correction_ref_representative.fa.gz` (136) – a
    diversity-sampled subset, for a smaller footprint.

  Both keep the pseudogenes. Their frameshifts came from high-quality
  genomes and so are real biology, and correction is meant to recover a
  sequence as it exists in nature rather than reshape every array into
  an intact TALE. A reference of only intact TALEs risks "repairing" a
  genuine pseudogene into an ORF no strain carries.

- max_comparisons:

  How many reference proteins each array may be aligned against during
  frameshift correction, and **the main control on how long correction
  takes**. `NULL`, the default, allows all of them.

  [`DECIPHER::CorrectFrameshifts()`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html)
  scores every reference with a cheap distance first, sorts them, and
  only then aligns against the closest `max_comparisons` of them –
  stopping sooner if one is close enough. So the references never
  reached cost almost nothing, and lowering this does not change *which*
  references are preferred, only how deep the search goes before
  settling for the best seen.

  Measured on four arrays against the 1057-sequence source set, all
  giving byte-identical corrected sequences:

  |                     |             |
  |---------------------|-------------|
  | **max_comparisons** | **seconds** |
  | all (1057)          | 252         |
  | 400                 | 179         |
  | 100                 | 47          |
  | 50                  | 23          |
  | 20                  | 10          |

  **The trade-off.** A cap risks a divergent array whose only good
  reference lies outside the closest `max_comparisons` by the cheap
  pre-screen. That pre-screen is an approximation, so a low cap trusts
  it to rank the truly best reference near the top.

  When it fails it does not fail by leaving an array uncorrected – it
  fails by correcting it against a poor reference, which is worse,
  because the result still looks like a corrected ORF. Against a
  deliberately small 20-sequence reference, the same four arrays give:

  |                     |                             |
  |---------------------|-----------------------------|
  | **max_comparisons** | **indels called per array** |
  | all (20), 20, 10    | 2, 2, 0, 1                  |
  | 5                   | 2, 2, 0, 2                  |
  | 2                   | 9, 11, 0, 15                |

  At 2 the aligner cannot reach a decent reference and invents indels
  wholesale. What matters is therefore not the ratio to the reference
  set but whether the closest `max_comparisons` are genuinely close: 20
  of 1057 is ample, 5 of 20 is not. With a large reference set a cap in
  the tens is safe and very much faster; with a small or a poorly
  matched one, prefer the default and pay for the full search.

- frameshift:

  This is an internal parameter of the
  [`CorrectFrameshifts`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html)
  function. The default is 11 and fiddle with this at your own risk...

- ...:

  Additional parameters for the
  [`CorrectFrameshifts`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html)
  function.

## Value

This functions has only side effects (writing files, mostly). However,
if everything ran smoothly, it will invisibly return the path of the
directory where output files were written.

List of output files:

- all_ranges.gff: gff file of all Tal arrays detected by HMMer

- array_report.tsv: report of all Tal arrays. In array_report.tsv,
  column *predicted_dels_count*/*predicted_ins_count* shows the number
  of putative deletions/insertions in the raw sequences that have been
  corrected in the corrected sequences with the function
  [`CorrectFrameshifts`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html).

- hits_report.tsv: report of all hits detected by HMMer

- hits_report.gff: gff file of all hits detected by HMMer

- domains_report.tsv: report of all Tal amino acid domains detected by
  AnnoTALE analyze

- putative_tal_orf.fasta: Tal putative ORFs

- pseudo_tal_cds.fasta: pseudo Tal CDS, putative Tal array ORFs detected
  by HMMer for whch AnnoTALE analyze failed to find RVD(s).

- rvd_sequences.fas: Sequence of RVDs (separated by rvd_sep) predicted
  to be encoded in the Tal array ORFs by AnnoTALE. Note that if
  extremity_codes is `TRUE` (by default), the N- and C-TREM anchor codes
  will be appended at the beginning and end of the sequences if the
  corresponding domain coding sequence was wound by HMMer at the DNA
  level. If no such HMMer hits were found, the "XXXXX" string will be
  appended to denote that AA sequences outside of the RVD array are
  likely to be atypical.

- c_terminus_aa_alignment.html: protein alignment of all C-termini

- c_terminus_dna_alignment.html: DNA alignment of all C-termini

- n_terminus_aa_alignment.html: protein alignment of all N-termini

- n_terminus_dna_alignment.html: DNA alignment of all N-termini

- tale_cds_all_diagnostic_regions_hmmfile.out: HMMER profile used for
  tale cds search.

- hmmer_search_out.txt: ignore

- nhmmer_human_readable_output_of_last_run.txt: primary HMMER output
  file.

- tell_tales.log: a log file

- annotale folder: folder containing result of AnnoTALE analyze for all
  Tal arrays

- correction_alignment_aa folder: folder containing protein alignment of
  Tal array detected by HMMer and corrected Tal array if `correct_array`
  = TRUE

- correction_alignment_dna folder: folder containing DNA alignment of
  Tal array detected by HMMer and corrected Tal array if
  `correct_array = TRUE`

## Details

The approach is first to use [HMMER](http://hmmer.org/) to find and
categorize regions in the input DNA sequence that are related to the
coding sequence of canonical TALE protein domains (N-Term, repeats,
C-term). Hits that are (nearly – see the `min_gap` parameter) adjacent
are grouped in "taleArrays" which are considered as potential tal genes.

If the `correct_array` parameter is turned off, the longest predicted
open reading frame (+extend_len) for each talArray is fed to
[AnnoTALE](http://www.jstacs.de/index.php/AnnoTALE) to detect TALE
domains in the predicted translation product. The Results should hence
be very similar to what would be obtained with AnnoTALE, plus many
additional informative output files such as tabular reports.

If `correct_array` is turned on, these talearrays are passed to the
[`CorrectFrameshifts`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html)
function that attemps to 'correct' potential frameshifts in the
taleArray sequences. This conveniently removes many artefactual indels
but bear in mind that this may also **erroneously** 'correct' genuine
frame shifts which can be highly relevant especially for truncTALEs or
iTALES. The resulting 'corrected' taleArray open reading frames are then
passed to AnnoTALE.

Note that occasionally, when a putative open reading frame does not
encode a canonical TALE protein (early frame shift, incomplete ORF,
etc...), the "analyze" module of AnnoTALE outputs DNA parts but no
protein parts and/or RVD sequence. This should be detected and reported
in the tell_tales log.

## See also

Other TALE discovery:
[`correct_tales()`](https://scunnac.github.io/tantale/dev/reference/correct_tales.md),
[`tales_from_telltale()`](https://scunnac.github.io/tantale/dev/reference/tales_from_telltale.md)

## Examples

``` r
# \donttest{
# Needs nhmmer and AnnoTALE, resolved from the tantale conda environment
# (and a Java runtime) on first use.
subj <- system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                    package = "tantale")
out <- tempfile("tell_tales_example")
tell_tales(subject_file = subj, output_dir = out)
#> HMMER is very picky about forbiden characters in sequence name. Renaming
#> sequences in
#> /tmp/RtmpeEXpWy/temp_libpath1edbda181ac0b4/tantale/extdata/bai3_sample_tal_genomic_regions.fasta.
#> Original seq names : talRegion5 ; talRegion6
#> Dummy seq names : seq1 ; seq2
#> HMMER 3.3.2 (Nov 2020); http://hmmer.org/
#> Copyright (C) 2020 Howard Hughes Medical Institute.
#> Warning: Some HMMER hits overlap, so the inferred RVD sequences may carry artefactual
#> insertions.
#> ℹ Check these regions: "ROI_00001", "ROI_00002", "ROI_00003", and "ROI_00004"
#> Now running AnnoTALE analyze for ROI_00001
#> Now running AnnoTALE analyze for ROI_00002
#> Now running AnnoTALE analyze for ROI_00003
#> Now running AnnoTALE analyze for ROI_00004
#> #**************************************** #** tell_tales analysis done **
#> Current date: Fri Sep 18 20:42:47 2026 #_________Provided I/O parameters
#> __________ File of subject DNA sequences:
#> /tmp/RtmpeEXpWy/temp_libpath1edbda181ac0b4/tantale/extdata/bai3_sample_tal_genomic_regions.fasta
#> TALE N-term CDS region detection HMM file:
#> /tmp/RtmpeEXpWy/temp_libpath1edbda181ac0b4/tantale/extdata/hmmProfile/Xo_TALE_Nterm_CDS_profile.hmm
#> TALE repeat unit CDS detection HMM file:
#> /tmp/RtmpeEXpWy/temp_libpath1edbda181ac0b4/tantale/extdata/hmmProfile/Xo_TALE_repeat_CDS_profile.hmm
#> TALE C-term CDS region detection HMM file:
#> /tmp/RtmpeEXpWy/temp_libpath1edbda181ac0b4/tantale/extdata/hmmProfile/Xo_TALE_Cterm_CDS_profile.hmm
#> Output directory: /tmp/RtmpZLcluY/tell_tales_example1edd6366dce955
#> #____________Other parameters________________ nterm_min_score: 300
#> repeat_min_score: 20 cterm_min_score: 200 min_domain_hits: 4 min_array_length:
#> 0 merge_hits: TRUE min_gap: 35 extend_len: 300 correct_array: FALSE
#> correction_ref:
#> /tmp/RtmpeEXpWy/temp_libpath1edbda181ac0b4/tantale/extdata/tale_correction_ref.fa.gz
#> max_comparisons: all frameshift: -11 #__________Summary measures of TALE search
#> outcome__________ Number of analysed subject sequences : 2 Total number of TALE
#> repeat DNA coding sequence motif hits found with the nhmmer approach: 88 Total
#> number of subject seqs with TALE motif hits after low hit number filtering: 2
#> Total number of distinct regions (repeat arrays) with adjacent TALE motifs : 4
#> Total number of 'complete' arrays (with both N- and C-term flanking motifs): 4
#> Minimum array length (number of TALE domain hits): 16 Maximum array length: 28
#> Median array length: 26 Number of gaps of size below 500nt between TALE motifs
#> arrays: 1.5 First quartile of size of gaps (below 500nt) between TALE motifs
#> arrays: 108 Median size of gaps (below 500nt) between TALE motifs arrays: 108
#> Upper quartile of size of gaps (below 500nt) between TALE motifs arrays: 108
#> #__________Noteworthy AnnoTale issues__________ # #*************************
tales_from_telltale(out)
#> <tales> 4 arrays, 96 parts
#>   layers: rvd   |   6 other columns
#>              rvd
#>   ROI_00001  NTERM NN NG NN HD HD NI N* NG HD NI NG NN HD NI NG NI NG NN N ...
#>   ROI_00002  NTERM NN HD NI NN HD NG HD HD NG NG NI NG NI NG CTERM
#>   ROI_00003  NTERM NN ND NN NI NK NN HD NN NG NG N* HD N* HD NI NN HD NG H ...
#>   ROI_00004  NTERM NI HD NN NS NN NG HD NG HD NG NN NG HD NS HD NI NG HD H ...
# }
```
