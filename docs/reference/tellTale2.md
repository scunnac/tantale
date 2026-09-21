# This function name is deprecated and will ultimately be removed. It corresponds to the [tellTale](https://scunnac.github.io/tantale/reference/tellTale.md) which should be used instead.

This function name is deprecated and will ultimately be removed. It
corresponds to the
[tellTale](https://scunnac.github.io/tantale/reference/tellTale.md)
which should be used instead.

## Usage

``` r
tellTale2(
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
