# tantale: Transcription Activator-Like Effectors (TALEs) tools

An integrated collection of functions for TALE mining and analysis in R.

Please take a look at the package
[website](https://scunnac.github.io/tantale) for further details.

## Details

![tantale_logo](figures/tantale_logo_small.gif)

## Note

CAUTIONARY NOTES:

- tantale has been written with only Linux systems in mind and will very
  likely **not work on other OS** (eg Windows)

- Some of tantale wrappers use code written in other languages : **Java
  and Perl must be on the PATH** in your system.
  [`tantale_setup()`](https://scunnac.github.io/tantale/dev/reference/tantale_setup.md)
  checks for both.

## A TALE-oriented OOP framework

- `tales`/`tales_msa` S3 classes, with subsetting, coercion, and
  plotting methods

## TALE mining in bacterial sequences

- Wrapper around annotale_jar and correcTALE

- tell_tales, an R function similar to annotale_jar

- Analysis tools for RVD inventory, repeat lenght

## TALEs classification, phylogeny

- Wrappers around distal, functal, annotale_jar

- TALE groups inference

- Easily build Multiple alignments and generate nice plots

## TALE targets mining

- Wrappers around target predictors

- General parser for results aggregation

- Connector with daTALbase (to be done)

## Setting up

tantale does not bundle the programs it drives. MAFFT, HMMER, mmseqs2
and the Perl dependencies of the target predictors come from a conda
environment the package builds for itself, so **conda (or mamba, or
micromamba) is a prerequisite of the main workflow**, not only of
optional extras. The environment is built on first use rather than at
install time, so the first call needs a network connection and takes a
few minutes.

Start with
[`tantale_setup()`](https://scunnac.github.io/tantale/dev/reference/tantale_setup.md).
Called bare it reports what is present and changes nothing;
`tantale_setup(install = TRUE)` builds or repairs the environment. It
checks **versions**, not merely presence, which matters because MAFFT
changed its `--text` mode gap handling after 7.4x and later versions
align TALE repeat strings differently – an environment left over from an
older tantale gives different alignments from the same input, and
nothing else would report it.

If you have no conda at all, `reticulate` will install one from inside R
with
[`reticulate::install_miniconda()`](https://rstudio.github.io/reticulate/reference/install_miniconda.html)
(miniconda, note, not mamba). An existing conda, mamba or micromamba is
found automatically through
[`reticulate::conda_binary()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
and used instead.

## See also

Useful links:

- <https://scunnac.github.io/tantale>

- <https://github.com/scunnac/tantale>

- Report bugs at <https://github.com/scunnac/tantale/issues>

## Author

**Maintainer**: Sebastien Cunnac <sebastien.cunnac@ird.fr>
([ORCID](https://orcid.org/0000-0002-3695-491X))

Authors:

- Sebastien Cunnac <sebastien.cunnac@ird.fr>
  ([ORCID](https://orcid.org/0000-0002-3695-491X))

- Bao Tram Vi <vbt576@gmail.com>
