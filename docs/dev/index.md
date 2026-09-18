[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

![](./reference/figures/tantale_logo_small.gif)

## An integrated collection of functions for [TALE](https://en.wikipedia.org/wiki/Transcription_activator-like_effector) mining and analysis with the R language

Analyzing TALEs in (mostly *Xanthomonas*) genomes typically means
coordinating several concurrent and complementary tools running on
different platforms (Java, Perl), which is cumbersome to script and
automate. Making sense of the output is harder still, since there is no
easy way to graphically represent the various objects of the analysis.

With `tantale`, we compiled and extended our previous code wrapping TALE
analysis tools into an integrated R interface that further provides an
extensive list of utilities for easy plotting. This enables a moderately
proficient R programmer to perform entire analysis pipelines directly in
R and access result objects for custom manipulations.

Here is a snapshot of the topics that are or will (hopefully) be covered
in the future:

- A TALE-oriented OOP framework:
  - `tales`/`tales_msa` S3 classes, with subsetting, coercion, and
    plotting methods – see [the `tales` class
    article](https://scunnac.github.io/tantale/articles/tales_class.html)
    and [the `tales_msa` class
    article](https://scunnac.github.io/tantale/articles/tales_msa_class.html)
- TALE mining in bacterial sequences:
  - Wrapper around [AnnoTALE](https://doi.org/10.1038/srep21077) and
    [correcTALE](https://doi.org/10.1186/s12864-023-09228-1)
  - tell_tales, an R function similar to AnnoTALE
  - Analysis tools for RVD inventory, repeat lenght
- TALEs classification, phylogeny:
  - Wrappers around [DisTAL](https://doi.org/10.3389/fpls.2015.00545),
    [functal](https://doi.org/10.3389/fpls.2015.00545), AnnoTALE
  - TALE groups inference
  - Easily build Multiple alignments and generate nice plots
- TALE targets predictions:
  - Wrappers around target predictors
    ([Talvez](https://doi.org/10.1371/journal.pone.0068464) and
    [PrediTALE](https://doi.org/10.1371/journal.pcbi.1007206))
  - General parser for results aggregation
  - Connector with
    [daTALbase](https://doi.org/10.1094/MPMI-06-17-0153-FI) (to be done)

## Installation

For further details, take a look at the package
[website](https://scunnac.github.io/tantale).

### 1. Install the R package

``` r
remotes::install_github("scunnac/tantale",
                        type = "source",
                        dependencies = TRUE,
                        upgrade = "never")
```

### 2. Make sure conda is available

tantale does not bundle the programs it drives. MAFFT, HMMER, mmseqs2
and the Perl dependencies of the target predictors come from a conda
environment the package builds for itself, so **conda (or mamba, or
micromamba) is a prerequisite of the main workflow**, not just of
optional extras.

You do not have to install it by hand or know anything about it. If you
have no conda at all, `reticulate` will install one from inside R:

``` r
install.packages("reticulate")
reticulate::install_miniconda()
```

Note this installs **miniconda**, not mamba. If you already have conda,
mamba or micromamba, tantale finds it through
[`reticulate::conda_binary()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
and uses that instead – nothing else to do.

### 3. Check everything is in place

``` r
tantale::tantale_setup()
```

This reports what is present and what is missing, and changes nothing.
`tantale_setup(install = TRUE)` then builds or repairs the environment.

Running it is optional – the environment is built on first use if it is
absent – but it is worth doing once, for two reasons:

- **The first real call would otherwise be the slow one.** The
  environment is built on first use rather than at install time, so it
  needs the network and takes a few minutes.
- **It checks versions, not just presence.** MAFFT changed its `--text`
  mode gap handling after 7.4x, and later versions align TALE repeat
  strings differently. An environment left over from an older version of
  tantale produces different alignments from the same input, and nothing
  else would tell you.

> **If you have both conda and micromamba**, note they keep **separate
> roots**. An environment named `tantale` in one is not the one in the
> other, and a rebuild can report success while the package goes on
> using the other copy.
> [`tantale_setup()`](https://scunnac.github.io/tantale/dev/reference/tantale_setup.md)
> prints the binary, the default root and the environment actually in
> use, precisely so this is visible.

### Also needed

- **Java and Perl on the PATH.** Several wrappers (AnnoTALE, PrediTALE,
  TALE correction, the target predictors) are written in other
  languages.
  [`tantale_setup()`](https://scunnac.github.io/tantale/dev/reference/tantale_setup.md)
  checks for these too.
- **Linux.** tantale has been written with only Linux in mind and will
  very likely not work on other operating systems.

tantale ships about 60 MB of Java programs (AnnoTALE, PrediTALE and TALE
correction) that have no conda package, which is most of its footprint.

------------------------------------------------------------------------

**NOTE** :

- tantale is under active development ahead of publication; interfaces
  may still change.
- Documentation could be improved and extended.
- If you feel like contributing, that is great, please send me an email:
  <sebastien.cunnac@ird.fr>

------------------------------------------------------------------------

## Use of large language models

The authors used large language models (Claude, Anthropic – including
Claude Sonnet 5) to assist with code development, debugging, and
documentation writing throughout this package. Where LLM assistance
extends to a manuscript describing this work, it is limited to the
copy-editing stage; the manuscript itself is written entirely by the
authors. Any figures are prepared by the authors, with LLMs used only to
help write the scripts that generate them. The authors affirm that they
are fully responsible for the content of the codebase, its
documentation, and any accompanying manuscript.
