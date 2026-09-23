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
and `prefix`, so the result can be tested as well as read. `conda` is
`NULL` and `prefix` is `NA` when no conda/mamba installation, or no
`tantale` environment, was found at all; `system` is always a data
frame, since Java and Perl are checked regardless.

## Details

**Why versions are checked as well as presence.** MAFFT changed its
`--text` mode gap handling after 7.4x, and later versions align TALE
repeat-code strings differently, leaving the N- and C-termini
unanchored. An environment built by an older version of this package can
therefore produce different alignments from the same input, with nothing
to indicate it. The pins in `tantale_conda_env.yaml` exist for that
reason, and this function checks them.

**Why three paths are reported.** The conda binary, its default root,
and the environment actually in use are three different things, and on a
machine with any history they diverge – `reticulate` scans several known
locations, so two roots can each hold an environment named `tantale`.
When that happens, a rebuild can honestly report success while the
package goes on using the other one. Everything here therefore works
with the environment's prefix (its path), never its name.

**Java and Perl** are checked too. They are hard requirements of the
AnnoTALE, PrediTALE and TALE-correction wrappers, they come from outside
conda, and without this check they fail deep inside a
[`system()`](https://rdrr.io/r/base/system.html) call.

## Installing conda itself

`conda = TRUE` installs a conda distribution if none is found. It is
deliberately opt-in and separate from `install`: putting a package
manager on someone's machine is a larger side effect than building an
environment in one that already exists. It installs **miniconda**, via
[`reticulate::install_miniconda()`](https://rstudio.github.io/reticulate/reference/install_miniconda.html).

## See also

[`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md),
[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md),
[`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md),
[`talvez()`](https://scunnac.github.io/tantale/reference/talvez.md) and
[`correct_tales()`](https://scunnac.github.io/tantale/reference/correct_tales.md),
the entry points that need these tools.

## Examples

``` r
# \donttest{
# A read-only diagnostic report; nothing is installed or changed.
tantale_setup()
#> ✔ conda binary    /home/cunnac/bin/micromamba
#> ℹ default root    /home/cunnac/micromamba
#> ✔ tantale env     /home/cunnac/mamba/envs/tantale
#> ✔ mmseqs2             14.7e284
#> ✔ perl-statistics-r   0.34
#> ✔ hmmer               3.3.2
#> ✔ clustalo            1.2.4
#> ✔ igvtools            2.16.2
#> ✔ mafft               7.453
#> ✔ perl-list-moreutils 0.430
#> ✔ java                /usr/bin/java
#> ✔ perl                /usr/bin/perl
#> ✔ Everything tantale needs is present.
# }
```
