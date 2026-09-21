# Search and report on the features of TALE protein domains potentially encoded in subject DNA sequences

`tellTale` has been primarily written to report on 'corrected' TALE RVD
sequences in indels prone, noisy DNA sequences (suboptimally polished
genomes assembly, raw reads of long read sequencing technologies \[eg
PacBio, ONT\]) that would otherwise be missed by conventional tools (eg
AnnoTALE).

## Usage

``` r
tellTale(
  subjectFile,
  outputDir = getwd(),
  hmmFilesDir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T),
  TALE_NtermDNAHitMinScore = 300,
  repeatDNAHitMinScore = 20,
  TALE_CtermDNAHitMinScore = 200,
  minDomainHitsPerSubjSeq = 4,
  mergeHits = TRUE,
  minGapWidth = 35,
  appendExtremityCodes = TRUE,
  rvdSep = "-",
  hmmerpath = system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = T),
  extendedLength = 300,
  talArrayCorrection = FALSE,
  refForTalArrayCorrection = system.file("extdata", "decipher_ref_tales_aa.fa.gz",
    package = "tantale", mustWork = T),
  frameShiftCorrection = -11,
  ...
)
```

## Arguments

- subjectFile:

  Fasta file with DNA sequence(s) to be searched for the presence of
  TALE coding sequences (CDS).

- outputDir:

  Path of the output directory. If not specified, results will be
  written to current working folder.

- hmmFilesDir:

  Specify the path to a folder holding the hmmfiles if you do not want
  to use the ones provided with tantale.

- TALE_NtermDNAHitMinScore:

  Minimal nhmmer score cut_off value to consider the hit as genuine

- repeatDNAHitMinScore:

  Minimal nhmmer score cut_off value to consider the hit as genuine

- TALE_CtermDNAHitMinScore:

  Minimal nhmmer score cut_off value to consider the hit as genuine

- minDomainHitsPerSubjSeq:

  Minimum number of nhmmer hits for a subject sequence to be reported as
  having TALE diagnostic regions. This is a way to simplify output a
  little by getting ride of uninformative sequences

- mergeHits:

  Perform overlapping hits merging per domain type. Should not be
  modified.

- minGapWidth:

  Minimum gap in base pairs between two tale domain hits for them to be
  considered distinct. If the length of the gap is below this value,
  domains are considered "contiguous" and grouped in the same array.

- appendExtremityCodes:

  Set this to `FALSE` if you do not want the N- and C-TREM anchor codes
  in the output sequences of RVD

- rvdSep:

  Symbol acting as a separator in RVD sequences

- hmmerpath:

  Specify the path to a directory holding the HMMER executable if you do
  not want to use the ones provided with tantale.

- extendedLength:

  number of nucleotides to extend in 3'-end at the tal ORF prediction
  stage.

- talArrayCorrection:

  True or False

- refForTalArrayCorrection:

  Reference AA sequences for tal array predicted ORF correction if you
  do not want to use the ones provided with tantale.

- frameShiftCorrection:

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

- allRanges.gff: gff file of all Tal arrays detected by HMMer

- arrayReport.tsv: report of all Tal arrays. In the arrayReport.tsv,
  column *predicted_dels_count*/*predicted_ins_count* shows the number
  of putative deletions/insertions in the raw sequences that have been
  corrected in the corrected sequences with the function
  [`CorrectFrameshifts`](https://rdrr.io/pkg/DECIPHER/man/CorrectFrameshifts.html).

- hitsReport.tsv: report of all hits detected by HMMer

- hitsReport.gff: gff file of all hits detected by HMMer

- domainsReport.tsv: report of all Tal amino acid domains detected by
  AnnoTALE analyze

- putativeTalOrf.fasta: Tal putative ORFs

- pseudoTalCds.fasta: pseudo Tal CDS, putative Tal array ORFs detected
  by HMMer for whch AnnoTALE analyze failed to find RVD(s).

- rvdSequences.fas: Sequence of RVDs (separated by rvdSep) predicted to
  be encoded in the Tal array ORFs by AnnoTALE. Note that if
  appendExtremityCodes is `TRUE` (by default), the N- and C-TREM anchor
  codes will be appended at the beginning and end of the sequences if
  the corresponding domain coding sequence was wound by HMMer at the DNA
  level. If no such HMMer hits were found, the "XXXXX" string will be
  appended to denote that AA sequences outside of the RVD array are
  likely to be atypical.

- C-terminusAAAlignment.html: protein alignment of all C-termini

- C-terminusDNAAlignment.html: DNA alignment of all C-termini

- N-terminusAAAlignment.html: protein alignment of all N-termini

- N-terminusDNAAlignment.html: DNA alignment of all N-termini

- TALE_CDS_all_diagnostic_regions_hmmfile.out: HMMER profile used for
  tale cds search.

- hmmerSearchOut.txt: ignore

- nhmmerHumanReadableOutputOfLastRun.txt: primary HMMER output file.

- tellTale.log: a log file

- annotale folder: folder containing result of AnnoTALE analyze for all
  Tal arrays

- CorrectionAlignmentAA folder: folder containing protein alignment of
  Tal array detected by HMMer and corrected Tal array if
  `talArrayCorrection` = TRUE

- CorrectionAlignmentDNA folder: folder containing DNA alignment of Tal
  array detected by HMMer and corrected Tal array if
  `talArrayCorrection = TRUE`

## Details

The approach is first to use [HMMER](http://hmmer.org/) to find and
categorize regions in the input DNA sequence that are related to the
coding sequence of canonical TALE protein domains (N-Term, repeats,
C-term). Hits that are (nearly \[see the minGapWidth parameter\])
adjacent are grouped in "taleArrays" which are considered as potential
tal genes.

If the `talArrayCorrection` parameter is turned off, the longest
predicted open reading frame (+extendedLength) for each talArray is fed
to [AnnoTALE](http://www.jstacs.de/index.php/AnnoTALE) to detect TALE
domains in the predicted translation product. The Results should hence
be very similar to what would be obtained with AnnoTALE, plus many
additional informative output files such as tabular reports.

If `talArrayCorrection` is turned on, these talearrays are passed to the
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
in the tellTale log.
