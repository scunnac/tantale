# Search and report on the features of TALE protein domains potentially encoded in subject DNA sequences

This function is maintained in the package temporarily and should not be
used for other purpose than curiosity. You are better off using
[tellTale](https://scunnac.github.io/tantale/reference/tellTale.md).

## Usage

``` r
tellTaleLegacy(
  subjectFile,
  outputDir = getwd(),
  hmmFilesDir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T),
  minRatioOfGapForColMasking = 0.8,
  TALE_NtermDNAHitMinScore = 300,
  repeatDNAHitMinScore = 20,
  TALE_CtermDNAHitMinScore = 200,
  minDomainHitsPerSubjSeq = 4,
  mergeHits = TRUE,
  repMsaMethod = "decipher",
  minGapWidth = 35,
  minDomainHitsPerArrayForAssembl = 5,
  taleArrayStartAnchorCode = "N-TERM",
  taleArrayEndAnchorCode = "C-TERM",
  appendExtremityCodes = TRUE,
  rvdSep = " ",
  hmmerpath = system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = T)
)
```

## Arguments

- subjectFile:

  Fasta file with DNA sequence(s) to be searched for the presence of
  TALE coding sequences (CDS).

- outputDir:

  Path of the output directory. If not specified, resluts will be
  written to working folder.

- hmmFilesDir:

  Specify the path to a folder holding the hmmfiles if you do not want
  to use the ones provided with tantale.

- minRatioOfGapForColMasking:

  Columns of the tale repeat CDS alignment that contain a gap in a
  fraction of sequence higher than this value (betwen 0 and 1) will be
  masked from the alignment when translating the DNA sequences to
  protein.

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

  Perform overlapping hits merging per domain type.

- repMsaMethod:

  Value is a character string being either "decipher" or "hmmalign".
  This parameter defines which method is used to compute tale repeat CDS
  multiple alignments as a prerequisite for repeat CDS indel correction.
  hmmaling tends to remove a couple of nucleotide at the extremities of
  the repeats when they do not match the hmm profile. In contrast, the
  DECIPHER package function tends to leave spurious sequences at the end
  of repeat regions (notably for the last half repeat). In both cases,
  this can bias subsequent tale analysis based on DisTALE that relies on
  repeat sequences.

- minGapWidth:

  Minimum gap in base pairs between two tale domain hits for them to be
  considered distinct. If the length of the gap is below this value,
  domains are considered "contiguous" and grouped in the same array.

- minDomainHitsPerArrayForAssembl:

  DEPRECATED argument. Used to speficy the Minimum number of repeat in
  an array for its seq of RVD to be considered for assembly. This is a
  way to get ride of sequences that are too short reasonably be of any
  help for assembly

- taleArrayStartAnchorCode:

  This scalar character vector will symbolize a TALE N-TERM CDS hit in
  the RVD sequence

- taleArrayEndAnchorCode:

  This scalar character vector will symbolize a TALE C-TERM CDS hit in
  the RVD sequence

- appendExtremityCodes:

  Set this to `FALSE` if you do not want the N- and C-TREM anchor codes
  in the output sequences of RVD

- rvdSep:

  Symbol acting as a separator in RVD sequences

- hmmerpath:

  Specify the path to a directory holding the HMMER executable if you do
  not want to use the ones provided with tantale.

## Value

This functions has only side effects (writing files, mostly). However,
if everything ran smoothly, it will invisibly return the path of the
directory where output files were written.

## Details

`tellTaleLegacy`, the predecessor of the current
[tellTale](https://scunnac.github.io/tantale/reference/tellTale.md)
function and has been primarily written to report on 'corrected' TALE
RVD sequences in noisy DNA sequences (suboptimally polished genomes
assembly, raw reads of long read sequencing technologies \[eg PacBio,
ONT\]) that would otherwise be missed by conventional tools (eg
AnnoTALE).

It works but is far from optimal for 'corrected' repeat CDS and **it
currently tends to remove 'unconventional' portions of repeat CDS**.
This may be problematic for downstream analysis with DisTAL, especially
for 'aberrant' repeats.

N- and C-terminal domains CDS are not currently 'corrected'.

It will work also if the subject DNA sequences are of high quality **but
may 'correct' genuine frame shifts in TALE repeat array CDS**.

Will output a number of files including tabular reports on the TALE
arrays (ie full length or partial tal gene CDS) and on the identified
coding sequences for TALE domains, including repeat CDS.
