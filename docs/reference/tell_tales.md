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
  output_dir = tempfile("tell_tales_"),
  hmm_dir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T),
  nterm_min_score = 300,
  repeat_min_score = 20,
  cterm_min_score = 200,
  terminus_max_evalue = 1e-05,
  min_dna_hits = 4,
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
  max_comparisons = 50,
  frameshift = -11,
  ...
)
```

## Arguments

- subject_file:

  Fasta file with DNA sequence(s) to be searched for the presence of
  TALE coding sequences (CDS).

- output_dir:

  Path of the output directory, created if it does not exist. The
  default is a new directory under
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html), which R deletes
  when the session ends: give a path to keep the results. The directory
  is returned, so `tales_from_telltales(tell_tales(genome))` reads the
  run straight back.

- hmm_dir:

  Folder holding the profile HMMs, if you do not want the ones provided
  with tantale. It must hold files with the same names: the three DNA
  profiles of the nhmmer search (`Xo_TALE_Nterm_CDS_profile.hmm`,
  `Xo_TALE_repeat_CDS_profile.hmm`, `Xo_TALE_Cterm_CDS_profile.hmm`) and
  the two protein profiles of the terminus check
  (`Xo_TALE_Nterm_AA_profile.hmm`, `Xo_TALE_Cterm_AA_profile.hmm`).

- nterm_min_score:

  Minimal nhmmer score cut_off value to consider the hit as genuine

- repeat_min_score:

  Minimal nhmmer score cut_off value to consider the hit as genuine

- cterm_min_score:

  Minimal nhmmer score cut_off value to consider the hit as genuine

- terminus_max_evalue:

  Maximum `hmmsearch` E-value for the segment AnnoTALE reports on either
  side of the repeats to count as a TALE N- or C-terminus. The segments
  are searched with the TALE terminal-domain protein profiles of
  `hmm_dir`; this decides the `NTERM`, `CTERM` and `XXXXX` codes (see
  [`tales_anchor_codes`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md)).
  Genuine termini truncated to about 40 residues still match with
  E-values below 1e-18. The match must also reach, within 10 positions,
  the end of the profile that adjoins the repeats: a terminus whose
  repeat-side part is in another reading frame after a frameshift
  matches only up to the frameshift, and is coded `XXXXX`. A terminus
  shorter than the canonical one that matches is coded `NTERM`/`CTERM`,
  though it has probably lost part of its function; see
  [`tales_anchor_codes`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md).

- min_dna_hits:

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

  Merge overlapping nhmmer hits of the same domain type, since nhmmer
  can report one repeat as two overlapping hits. With `FALSE`, such a
  repeat is counted twice in `n_dna_hits` and by `min_array_length`, and
  a warning names the arrays concerned.

- min_gap:

  Minimum gap in base pairs between two tale domain hits for them to be
  considered distinct. If the length of the gap is below this value,
  domains are considered "contiguous" and grouped in the same array.

- extremity_codes:

  Set this to `FALSE` if you do not want the terminus codes in the RVD
  strings of `rvd_sequences.fas` and `array_report.tsv`.

- rvd_sep:

  Symbol acting as a separator in RVD sequences

- hmmer_path:

  Specify the path to a directory holding the HMMER executable if you do
  not want to use the ones provided with tantale.

- extend_len:

  number of nucleotides to extend in 3'-end at the tal ORF prediction
  stage.

- correct_array:

  Whether to pass each array through
  [`CorrectFrameshifts`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html)
  before AnnoTALE sees it. `FALSE` by default; see Details for what
  turning it on buys (removing artefactual indels) and risks
  (erroneously "correcting" a genuine frameshift).

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
  takes**. `NULL` allows all of them.

  [`DECIPHER::CorrectFrameshifts()`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html)
  scores every reference with a cheap distance first, sorts them, and
  only then aligns against the closest `max_comparisons` of them –
  stopping sooner if one is close enough. So the references never
  reached cost almost nothing, and lowering this does not change *which*
  references are preferred, only how deep the search goes before
  settling for the best seen.

  The default, 50, gives the same result as the full search on the four
  genomes shipped with the package, corrected against the default
  reference (494 proteins): identical corrected sequences and RVD
  strings for every array of MAI1 (10 candidate arrays), BAI3 (10),
  PXO86 (19) and BAI3-1-1 (9), at about a quarter of the time or less.
  On BAI3-1-1, an error-prone assembly, the full search took 455 s and
  50 took 61 s.

  A smaller cap risks a divergent array whose only good reference lies
  outside the closest `max_comparisons` by the cheap pre-screen. The
  array is then corrected against a poor reference, which is worse than
  leaving it uncorrected, because the result still looks like a
  corrected ORF. On BAI3-1-1, one array of the nine needs more than 20
  references:

  |                     |             |                                            |
  |---------------------|-------------|--------------------------------------------|
  | **max_comparisons** | **seconds** | **that array**                             |
  | 2 to 5              | 19-22       | N-terminus unmatched, 19 of its 26 repeats |
  | 10, 20              | 26, 34      | not parsed by AnnoTALE, absent             |
  | 50                  | 61          | N-terminus, 26 repeats, C-terminus         |

  [`tales_anomalies`](https://scunnac.github.io/tantale/reference/tales_anomalies.md)
  reports the first outcome, and
  [`tales_from_telltales`](https://scunnac.github.io/tantale/reference/tales_from_telltales.md)
  warns about the second. With a reference of your own, especially a
  small or a distant one, compare a run at the default with one at
  `NULL` before relying on the cap.

- frameshift:

  Frameshift penalty passed to
  [`CorrectFrameshifts`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html)'s
  `frameShift`. tantale's default is `-11`, overriding DECIPHER's own
  `-15`; change it only with a reason.

- ...:

  Additional parameters for the
  [`CorrectFrameshifts`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html)
  function.

## Value

Called for its side effects (writing files). If everything runs
smoothly, it invisibly returns the path of the directory where output
files were written.

List of output files:

- all_ranges.gff: gff file of all Tal arrays detected by HMMer

- array_report.tsv: one row per candidate TALE array. Columns named
  `*_dna_*` describe the nhmmer search of the subject DNA; columns named
  `*_aa_*` describe the protein segments AnnoTALE extracted from the
  array's longest ORF.

  - *array_id*, *seqnames*, *start*, *end*, *strand*: the array's
    identifier and the span of its nhmmer hits on the subject sequence.

  - *n_dna_hits*: number of nhmmer hits (N-terminus, repeats and
    C-terminus profiles together) grouped in the array. A terminus hit
    usually overlaps the adjacent repeat hit by a few nucleotides, at
    the boundary between the two domains.

  - *array_seq*: DNA sequence of that span.

  - *nterm_dna_hit*, *cterm_dna_hit*: whether an nhmmer hit of the N-
    (C-) terminus DNA profile is part of the array, anywhere in it.

  - *rvd_string*: the RVDs AnnoTALE read, separated by `rvd_sep`, with
    the terminus codes described under *rvd_sequences.fas*. Empty when
    AnnoTALE found no RVD.

  - *has_aberrant_repeat*: whether AnnoTALE flagged a repeat of
    non-canonical length (a lowercase letter in its RVD).

  - *nterm_aa_evalue*, *cterm_aa_evalue*: E-value of the `hmmsearch`
    match between the segment AnnoTALE reported upstream (downstream) of
    the repeats and the TALE N- (C-) terminal protein profile. `NA` when
    there is no segment, or no match with an E-value up to 10.

  - *nterm_aa_profile_gap*, *cterm_aa_profile_gap*: number of profile
    positions between the end of that match and the end of the profile
    that adjoins the repeats (the last position of the N-terminal
    profile, the first of the C-terminal one). `0` for a match that
    reaches the repeats, `NA` when there is no match.

  - *nterm_aa_hit*, *cterm_aa_hit*: `TRUE` when that E-value is at most
    `terminus_max_evalue` and the profile gap at most 10, `FALSE` for a
    segment that does not match, `NA` when AnnoTALE reported no segment
    on that side.

  - *nterm_aa_length*, *cterm_aa_length*: length of those segments in
    amino acid residues, excluding a stop codon, as in the `tales`
    object's `aa_seq`.

  - *longest_orf_length*, *longest_orf_seq*: the longest ORF found in
    the array region extended by `extend_len` nucleotides at its 3' end.

  - *orf_coverage*: that ORF's length as a percentage of the extended
    region's length.

  - *predicted_dels_count*, *predicted_ins_count* (with
    `correct_array = TRUE`): number of putative deletions/insertions in
    the raw sequence that
    [`CorrectFrameshifts`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html)
    corrected.

- hits_report.tsv: report of all hits detected by HMMer

- hits_report.gff: gff file of all hits detected by HMMer

- domains_report.tsv: report of all Tal amino acid domains detected by
  AnnoTALE analyze

- putative_tal_orf.fasta: for each array in which AnnoTALE found RVDs,
  the DNA of the longest ORF of the array region extended by
  `extend_len` nucleotides at its 3' end (after frameshift correction
  with `correct_array = TRUE`). This is the putative TALE coding
  sequence AnnoTALE analysed.

- pseudo_tal_cds.fasta: for each array in which AnnoTALE found no RVD,
  the DNA of the array region extended by `extend_len` nucleotides at
  its 3' end, as found in `subject_file`: candidate pseudogenes,
  assembly errors or false detections, kept for inspection.

- rvd_sequences.fas: the RVDs (separated by `rvd_sep`) of each array for
  which AnnoTALE found at least one. With `extremity_codes = TRUE` (the
  default), each string is bracketed by terminus codes (see
  [`tales_anchor_codes`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md)):
  `NTERM` (`CTERM`) when the segment AnnoTALE reported upstream
  (downstream) of the repeats matches the TALE N- (C-) terminal protein
  profile, `XXXXX` when it does not, and no code when AnnoTALE reported
  no segment on that side.

- c_terminus_aa_alignment.html: protein alignment of all C-termini (only
  written when at least 2 were found; skipped with a warning otherwise)

- c_terminus_dna_alignment.html: DNA alignment of all C-termini (same
  condition)

- n_terminus_aa_alignment.html: protein alignment of all N-termini (same
  condition)

- n_terminus_dna_alignment.html: DNA alignment of all N-termini (same
  condition)

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
function that attempts to 'correct' potential frameshifts in the
taleArray sequences. This conveniently removes many artefactual indels
but bear in mind that this may also **erroneously** 'correct' genuine
frame shifts which can be highly relevant especially for truncTALEs or
iTALEs. The resulting 'corrected' taleArray open reading frames are then
passed to AnnoTALE.

Note that occasionally, when a putative open reading frame does not
encode a canonical TALE protein (early frame shift, incomplete ORF,
etc...), the "analyze" module of AnnoTALE outputs DNA parts but no
protein parts and/or RVD sequence. This should be detected and reported
in the tell_tales log.

Each sequence of `subject_file` is treated as linear. A *tal* gene that
spans the junction of a circular molecule (the two ends of an assembled
chromosome or plasmid) is cut in two, and is reported, if at all, as two
partial arrays at the ends of the sequence. Rotating the sequence so
that it starts elsewhere avoids this.

## See also

Other TALE discovery:
[`correct_tales()`](https://scunnac.github.io/tantale/reference/correct_tales.md),
[`tale_annotations`](https://scunnac.github.io/tantale/reference/tale_annotations.md),
[`tales_from_annotale()`](https://scunnac.github.io/tantale/reference/tales_from_annotale.md),
[`tales_from_telltales()`](https://scunnac.github.io/tantale/reference/tales_from_telltales.md)

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
#> /home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/bai3_sample_tal_genomic_regions.fasta.
#> Original seq names : talRegion5 ; talRegion6
#> Dummy seq names : seq1 ; seq2
#> HMMER 3.3.2 (Nov 2020); http://hmmer.org/
#> Copyright (C) 2020 Howard Hughes Medical Institute.
#> Now running AnnoTALE analyze for ROI_00001
#> Now running AnnoTALE analyze for ROI_00002
#> Now running AnnoTALE analyze for ROI_00003
#> Now running AnnoTALE analyze for ROI_00004
#> #****************************************
#> #**   tell_tales analysis done     **
#> Current date:    Mon Oct  5 01:07:32 2026
#> #_________Provided I/O parameters __________
#> File of subject DNA sequences:   /home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/bai3_sample_tal_genomic_regions.fasta
#> TALE N-term CDS region detection HMM file:   /home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/hmmProfile/Xo_TALE_Nterm_CDS_profile.hmm
#> TALE repeat unit CDS detection HMM file: /home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/hmmProfile/Xo_TALE_repeat_CDS_profile.hmm
#> TALE C-term CDS region detection HMM file:   /home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/hmmProfile/Xo_TALE_Cterm_CDS_profile.hmm
#> Output directory:    /tmp/Rtmp8rsnUq/tell_tales_example1e3a23258bfb0c
#> #____________Other parameters________________
#> nterm_min_score: 300
#> repeat_min_score:    20
#> cterm_min_score: 200
#> terminus_max_evalue: 1e-05
#> min_dna_hits:    4
#> min_array_length:    0
#> merge_hits:  TRUE
#> min_gap: 35
#> extend_len:  300
#> correct_array:   FALSE
#> correction_ref:  /home/cunnac/Lab-Related/MyScripts/tantale/inst/extdata/tale_correction_ref.fa.gz
#> max_comparisons: 50
#> frameshift:  -11
#> #__________Summary measures of TALE search outcome__________
#> Number of analysed subject sequences :   2
#> Total number of TALE repeat DNA coding sequence motif hits found with the nhmmer approach:   88
#> Total number of subject seqs with TALE motif hits after low hit number filtering:    2
#> Total number of distinct regions (repeat arrays) with adjacent TALE motifs : 4
#> Number of arrays with nhmmer DNA hits for both termini:  4
#> Number of arrays whose AnnoTALE N-terminus matches the TALE N-terminal protein profile:  4
#> Number of arrays whose AnnoTALE C-terminus matches the TALE C-terminal protein profile:  4
#> Minimum array length (number of nhmmer DNA hits):    16
#> Maximum array length:    28
#> Median array length: 26
#> Number of gaps of size below 500nt between TALE motifs arrays:   1.5
#> First quartile of size of gaps (below 500nt) between TALE motifs arrays: 108
#> Median size of gaps (below 500nt) between TALE motifs arrays:    108
#> Upper quartile of size of gaps (below 500nt) between TALE motifs arrays: 108
#> #__________Noteworthy AnnoTale issues__________
#> # 
#> #*************************
tales_from_telltales(out)
#> <tales> 4 arrays, 96 parts
#>   layers: rvd   |   6 other columns
#>              rvd
#>   ROI_00001  NTERM NN NG NN HD HD NI N* NG HD NI NG NN HD NI NG NI NG NN N ...
#>   ROI_00002  NTERM NN HD NI NN HD NG HD HD NG NG NI NG NI NG CTERM
#>   ROI_00003  NTERM NN ND NN NI NK NN HD NN NG NG N* HD N* HD NI NN HD NG H ...
#>   ROI_00004  NTERM NI HD NN NS NN NG HD NG HD NG NN NG HD NS HD NI NG HD H ...
# }
```
