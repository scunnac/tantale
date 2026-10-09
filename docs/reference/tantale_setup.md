# Check, and optionally build, tantale's external dependencies

Reports whether the external programs tantale drives are present and at
the versions it expects, and can build or repair the conda environment
that provides most of them. It also downloads the Java programs tantale
wraps (AnnoTALE, PrediTALE, TALEcorrection) and the four example genomes
of the articles
([`tantale_genome()`](https://scunnac.github.io/tantale/reference/tantale_genome.md)),
which are too large to be part of the package.

Called bare it changes nothing – it is a diagnostic. Pass
`install = TRUE` to act on what it finds.

## Usage

``` r
tantale_setup(
  install = FALSE,
  conda = FALSE,
  conda_bin = "auto",
  archive_dir = NULL
)
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

- archive_dir:

  A directory holding the archives `tantale-tools-1.tar.gz` and
  `tantale-genomes-1.tar.gz`, downloaded beforehand, to install them
  from instead of downloading them: for a machine without internet
  access, or to download once for several machines. `NULL` (default)
  downloads them.

## Value

Invisibly, a list with `conda`, `system` and `archives` data frames of
the checks, and `prefix`, so the result can be tested as well as read.
`conda` is `NULL` and `prefix` is `NA` when no conda/mamba installation,
or no `tantale` environment, was found at all; `system` is always a data
frame, since Java is checked regardless.

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

**Where the downloads go.** The function checks two archives, attached
to releases of tantale's GitHub repository, against a sha256 recorded in
the package and unpacks them into `tools::R_user_dir("tantale", "data")`
(`~/.local/share/R/tantale` on Linux), in one folder per archive version
(`tools-1/`, `genomes-1/`). Set the environment variable
`TANTALE_DATA_DIR` to use another directory, for instance one shared by
a team. The tools take about 60 MB, the genomes 20 MB.

**Java** is checked too. It is a hard requirement of the AnnoTALE,
PrediTALE and TALE-correction wrappers, it comes from outside conda, and
without this check they fail deep inside a
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
#> ✔ tools-1  /home/cunnac/snap/codium/495/.local/share/R/tantale/tools-1
#> ✔ genomes-1  /home/cunnac/snap/codium/495/.local/share/R/tantale/genomes-1
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
#> ✔ Everything tantale needs is present.
# }
```
