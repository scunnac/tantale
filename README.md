<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
![Coverage: 90.03%](https://img.shields.io/badge/coverage-90.03%25-brightgreen.svg)
<!-- badges: end -->

Test coverage measured locally with `covr::package_coverage()` on
2026-09-24; this repo has no CI yet, so the badge is a manual snapshot.
<p align="right">
  <img src="./man/figures/tantale_logo_small.gif">


#### ⚠️ This release is a major revision of the package with breaking changes to the previous interface. Extensive refactoring, bug fixes, and performance improvements bring the package to its first stable release. The previous release (prototypical version) remains available as a release asset:[tantale-full-history-2026-09-22.bundle](https://github.com/scunnac/tantale/releases/download/v0.9.9004).





## An integrated collection of functions for [TALE](https://en.wikipedia.org/wiki/Transcription_activator-like_effector) mining and analysis with the R language


Analyzing TALEs in  *Xanthomonas* genomes (mostly) typically means coordinating several concurrent and complementary tools running on different platforms (Java, Perl), which is cumbersome to script and automate. Making sense of the output is harder still, since there is no easy way to graphically represent the various objects of the analysis.

With `tantale`, we compiled and extended our previous code wrapping TALE analysis tools into an integrated R interface that further provides an extensive list of utilities for easy plotting. This enables a moderately proficient R programmer to perform entire analysis pipelines directly in R and access result objects for custom manipulations.

Here is a snapshot of the topics that are covered:


- A TALE-oriented OOP framework:
    - `tales`/`tales_msa` S3 classes, with subsetting, coercion, and plotting methods -- see [the `tales` class article](https://scunnac.github.io/tantale/articles/tales_class.html) and [the `tales_msa` class article](https://scunnac.github.io/tantale/articles/tales_msa_class.html)


- TALE mining in bacterial sequences:
    - Wrappers around [AnnoTALE](https://doi.org/10.1038/srep21077) and [correcTALE](https://doi.org/10.1186/s12864-023-09228-1); AnnoTALE's predictions load as a `tales` object
    - A TALE finder built for error-prone assemblies: it locates TALE arrays from DNA evidence and can correct frameshifts before AnnoTALE reads them
    - Analysis tools for RVD inventory, repeat length.
    - Compact 'talome' plots


- TALEs classification, phylogeny:
    - R reimplementations of [DisTAL](https://doi.org/10.3389/fpls.2015.00545) and [functal](https://doi.org/10.3389/fpls.2015.00545) comparisons, plus a wrapper around AnnoTALE
    - TALE groups inference with several methods
    - Easily build multiple alignments and generate nice plots


- TALE targets predictions:
    - Wrappers around target predictors ([Talvez](https://doi.org/10.1371/journal.pone.0068464) and [PrediTALE](https://doi.org/10.1371/journal.pcbi.1007206))
    - General parser for results aggregation
    - Connector with [daTALbase](https://doi.org/10.1094/MPMI-06-17-0153-FI) (to be done)



## Installation

For further details, take a look at the package
[website](https://scunnac.github.io/tantale).

### 1. Install the R package

```r
# install.packages("pak")  # if you don't already have it
pak::pkg_install("scunnac/tantale", dependencies = TRUE, upgrade = FALSE)
```

### 2. Make sure conda is available

MAFFT, HMMER, mmseqs2 and the Perl dependencies of the target predictors come
from a conda environment the package builds for itself, so **conda (or mamba, or micromamba) is a
prerequisite of the main workflow**.

You do not have to install it by hand or know anything about it. If you have
no conda at all, `reticulate` will install one from inside R:

```r
install.packages("reticulate")
reticulate::install_miniconda()
```

Note this installs **miniconda**, not mamba. If you already have conda, mamba
or micromamba, tantale finds it through `reticulate::conda_binary()` and uses
that instead, with nothing else to do.

### 3. Check everything is in place

```r
tantale::tantale_setup()
```

This reports what is present and what is missing, and changes nothing.
`tantale_setup(install = TRUE)` then builds or repairs the environment and
downloads the Java programs and example genomes described below.

Run it once with `install = TRUE` after installing tantale:

- **The Java wrappers and the example genomes need it.** They are
  downloaded by `tantale_setup()` only.
- **The first real call would otherwise be the slow one.** The conda
  environment is built on first use if it is absent, so that first call
  needs the network and takes a few minutes.
- **It checks versions as well as presence.** A worry free setup and you are good to go.

> **If you have both conda and micromamba**, `tantale_setup()` prints the binary, the default root and the
> environment actually in use, precisely so this is visible.

### Also needed

- **Java and Perl on the PATH.** Several wrappers (AnnoTALE, PrediTALE, TALE
  correction, the target predictors) are written in other languages.
  `tantale_setup()` checks for these too.
- **Linux.** tantale has been written with only Linux in mind. Its functionality with other operating systems has not been tested.

- **The Java programs and example genomes are downloaded once.**
  AnnoTALE, PrediTALE and TALE correction (about 60 MB, with no conda
  package) and the four example genomes of the articles (about 20 MB) are
  attached to [releases](https://github.com/scunnac/tantale/releases) of
  this repository. `tantale_setup(install = TRUE)` downloads them, checks
  them against checksums recorded in the package, and keeps them in your R
  data directory; set the environment variable `TANTALE_DATA_DIR` to use
  another one. On a machine without internet access, download the two
  archives yourself and pass their folder as
  `tantale_setup(install = TRUE, archive_dir = ...)`.

---

**NOTE** :

- The interface of tantale is stable. From version 1.0.0 on, any further
  change to it follows the conventions of the
  [lifecycle](https://lifecycle.r-lib.org/articles/stages.html) package: a
  function or argument to be removed or renamed is first deprecated, with
  a warning that names its replacement, and only removed in a later
  release.
- If you feel like contributing, that is great, please send me an email: sebastien.cunnac@ird.fr

---

## Licence

tantale's own code is under the MIT licence. The programs and data it
bundles or downloads from other projects keep their own terms, listed file by file,
with upstream sources, in
[`inst/COPYRIGHTS`](https://github.com/scunnac/tantale/blob/main/inst/COPYRIGHTS):

- AnnoTALE, PrediTALE and TALEcorrection, from the
  [Jstacs](https://www.jstacs.de) project, are under the GNU GPL,
  version 3 or later. Their source code is at
  <https://github.com/Jstacs/Jstacs>; the licence text travels with them
  in the tools archive (`LICENSES/COPYING.GPL-3`).
- TALVEZ 3.2 and the QueTAL FuncTAL table behind `rvd_dna_specificity`
  carry no licence and are redistributed by permission of their author,
  Alvaro L. Pérez-Quintero.

## Use of large language models

The authors used large language models (Claude, Anthropic -- including
Claude Sonnet 5) to assist with code development, debugging, and
documentation writing throughout this package. The authors affirm that they
are fully responsible for the content of the codebase and its documentation.



