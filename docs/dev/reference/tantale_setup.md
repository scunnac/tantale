# Check, and optionally build, tantale's external dependencies

Reports whether the external programs tantale drives are present and at
the versions it expects, and can build or repair the conda environment
that provides most of them.

Called bare it changes nothing – it is a diagnostic. Pass
`install = TRUE` to act on what it finds.

## Usage

``` r
tantale_setup(install = FALSE, conda = FALSE, conda_bin = "auto")
```

## Arguments

- install:

  Build the `tantale` environment if it is missing, and repair it if a
  pinned version is wrong. `FALSE` by default, so the function reports
  and changes nothing.

- conda:

  Install a conda distribution if none is found. `FALSE` by default; see
  the section above.

- conda_bin:

  Passed to `reticulate`. `"auto"` lets it choose.

## Value

Invisibly, a list with `conda` and `system` data frames of the checks,
and `prefix`, so the result can be tested as well as read.

## Details

**Why the versions are checked and not just the presence.** MAFFT
changed its `--text` mode gap handling after 7.4x, and later versions
align TALE repeat-code strings differently, leaving the N- and C-termini
unanchored. An environment built by an older version of this package can
therefore produce different alignments from the same input, with nothing
to indicate it. The pins in `tantale_conda_env.yaml` exist for that
reason, and this function is what makes them real rather than
aspirational.

**Why three paths are reported.** The conda binary, its default root,
and the environment actually in use are three different things, and on a
machine with any history they diverge – `reticulate` scans several known
locations, so two roots can each hold an environment named `tantale`.
When that happens, a rebuild can honestly report success while the
package goes on using the other one. Everything here therefore operates
on the environment's prefix rather than its name.

**Java and Perl** are checked too. They are hard requirements of the
AnnoTALE, PrediTALE and TALE-correction wrappers, they are not conda's
business, and otherwise they fail deep inside a
[`system()`](https://rdrr.io/r/base/system.html) call.

## Installing conda itself

`conda = TRUE` installs a conda distribution if none is found. It is
deliberately opt-in and separate from `install`: putting a package
manager on someone's machine is a larger side effect than building an
environment in one that already exists. Note that this installs
**miniconda**, via
[`reticulate::install_miniconda()`](https://rstudio.github.io/reticulate/reference/install_miniconda.html),
not mamba.

## See also

[`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)
and
[`tales_align()`](https://scunnac.github.io/tantale/dev/reference/tales_align.md),
the two entry points that need these tools.
