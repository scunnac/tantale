# Contributing to tantale

Contributions are welcome: bug reports, corrections to the documentation,
and new features. This page explains how to report a problem and how to
prepare a change.

## Reporting a problem

Open an issue at <https://github.com/scunnac/tantale/issues>. Please
include:

- the output of `sessionInfo()` and of `tantale::tantale_setup()`, which
  reports the conda environment in use and the versions of the external
  programs;
- a minimal example that reproduces the problem, preferably on the
  sequences shipped in `inst/extdata` or on a small FASTA file. The
  [reprex](https://reprex.tidyverse.org) package helps.

For a result that looks biologically wrong (a missed TALE, a wrong RVD, an
unexpected terminus code), name the genome and the locus, and attach the
`tell_tales()` output directory if you can.

## Proposing a change

A typo or a small documentation fix can go straight to a pull request.
For anything larger, open an issue first so that the change can be
discussed before you spend time on it.

1. Fork the repository and create a branch from `main`.
2. Build the development environment with
   `tantale::tantale_setup(install = TRUE)`. The external programs are
   pinned on purpose: later versions of MAFFT (7.453) and HMMER (3.3.2)
   change the results. Please do not update these pins in a pull request.
3. Document functions with roxygen2 (markdown is enabled), then run
   `devtools::document()`.
4. Add tests with testthat. A test that needs an external program fails
   when the program is missing; it does not skip. The full suite
   (`devtools::test()`) takes several minutes.
5. `tests/testthat/test_golden.R` records what the pipeline currently
   produces. If your change alters that output on purpose, explain each
   changed snapshot row in the pull request.
6. Raise errors, warnings and messages with `cli::cli_abort()`,
   `cli::cli_warn()` and `cli::cli_inform()`, each with a condition class
   such as `c("tantale_error_<what>", "tantale_error")`.
7. Name functions and arguments after the
   [rOpenSci conventions](https://devguide.ropensci.org/pkg_building.html#function-and-argument-naming):
   `object_verb()` function names, snake_case, data as the first argument.
8. Add an entry to `NEWS.md` for anything a user would notice.
9. Run `devtools::check()` before opening the pull request.
