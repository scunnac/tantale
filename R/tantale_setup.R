#### Checking and building the package's external dependencies ####
#
# tantale drives seven programs it does not ship: MAFFT, HMMER, mmseqs2,
# clustalo and igvtools from a conda environment, plus Java and Perl from the
# system. Before 7.4 most of these were bundled; now they are not, so "is the
# environment right?" became a question a user can lose a day to.
#
# The reason this file exists is correctness rather than convenience.
# .create_tantale_env() tests only whether an environment *named* tantale
# exists, and reports success without looking inside it. An environment built
# by an older version of this package holds MAFFT 7.520, which aligns TALE
# repeat-code strings differently -- the termini come out unanchored and the
# column count changes -- and nothing anywhere would say so. That happened
# three times during 7.4 and was caught only because a golden baseline
# existed to compare against. A user has no such thing.
#
# So the pins in tantale_conda_env.yaml were aspirational. Checking them is
# the point.


#' The version pins, read from the environment file
#'
#' The yaml is the single source of truth. A second hardcoded list here would
#' drift from it, and the drift would be silent.
#' @return A named character vector, version by package name. Packages listed
#'   without a pin get `NA`.
#' @noRd
.tantale_pins <- function(yaml = NULL) {
  if (is.null(yaml)) {
    yaml <- system.file("tools", "tantale_conda_env.yaml", package = "tantale",
                        mustWork = TRUE)
  }
  lines <- readLines(yaml, warn = FALSE)
  # dependency entries only: "  - name=version", with comments and the
  # channels block excluded. Parsed by hand rather than with a yaml package
  # to avoid a dependency for six lines of text.
  start <- grep("^dependencies:", lines)
  if (!length(start)) return(character())
  lines <- lines[seq.int(start + 1L, length(lines))]
  lines <- sub("#.*$", "", lines)
  entries <- trimws(sub("^\\s*-\\s*", "", grep("^\\s*-\\s", lines, value = TRUE)))
  entries <- entries[nzchar(entries)]
  name <- sub("[=<>].*$", "", entries)
  version <- ifelse(grepl("=", entries, fixed = TRUE),
                    sub("^[^=]*=", "", entries), NA_character_)
  stats::setNames(version, name)
}


#' What is actually installed in a conda environment
#'
#' Read from `conda-meta/`, whose file names are `name-version-build.json`.
#' That is a documented part of a conda prefix's layout and costs nothing;
#' shelling out to `conda list` for the same answer would be slower and would
#' need conda to be working in order to tell you that conda is not working.
#' @return A named character vector, version by package name.
#' @noRd
.tantale_installed <- function(prefix) {
  meta <- file.path(prefix, "conda-meta")
  if (!dir.exists(meta)) return(character())
  files <- sub("\\.json$", "", list.files(meta, pattern = "\\.json$"))
  # split off the trailing -build and -version, leaving the name: package
  # names may themselves contain hyphens (perl-statistics-r), so this counts
  # from the right rather than splitting.
  version <- sub("^.*-([^-]+)-[^-]+$", "\\1", files)
  name <- sub("^(.*)-[^-]+-[^-]+$", "\\1", files)
  ok <- name != files
  stats::setNames(version[ok], name[ok])
}


#' Compare what is pinned against what is there
#' @return A data frame with one row per pinned package.
#' @noRd
.tantale_check_conda <- function(pins, installed) {
  found <- unname(installed[names(pins)])
  data.frame(
    tool = names(pins),
    required = unname(pins),
    found = found,
    ok = !is.na(found) & (is.na(pins) | found == pins),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}


#' The system tools conda does not provide
#'
#' Java and Perl are hard requirements of the AnnoTALE, PrediTALE and
#' TALEcorrection wrappers. They are not conda's business, and today they
#' fail deep inside a `system()` call with nothing useful said. This is the
#' only place they are ever checked.
#' @noRd
.tantale_check_system <- function() {
  needs <- c(
    java = "AnnoTALE, PrediTALE and TALE correction",
    perl = "the target-prediction wrappers"
  )
  paths <- Sys.which(names(needs))
  data.frame(
    tool = names(needs),
    required = "on PATH",
    found = ifelse(nzchar(paths), unname(paths), NA_character_),
    ok = nzchar(paths),
    needed_by = unname(needs),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}


#' Check, and optionally build, tantale's external dependencies
#'
#' @description
#' Reports whether the external programs tantale drives are present and at the
#' versions it expects, and can build or repair the conda environment that
#' provides most of them.
#'
#' Called bare it changes nothing -- it is a diagnostic. Pass
#' `install = TRUE` to act on what it finds.
#'
#' @details
#' **Why the versions are checked and not just the presence.** MAFFT changed
#' its `--text` mode gap handling after 7.4x, and later versions align TALE
#' repeat-code strings differently, leaving the N- and C-termini unanchored.
#' An environment built by an older version of this package can therefore
#' produce different alignments from the same input, with nothing to indicate
#' it. The pins in `tantale_conda_env.yaml` exist for that reason, and this
#' function is what makes them real rather than aspirational.
#'
#' **Why three paths are reported.** The conda binary, its default root, and
#' the environment actually in use are three different things, and on a
#' machine with any history they diverge -- `reticulate` scans several known
#' locations, so two roots can each hold an environment named `tantale`. When
#' that happens, a rebuild can honestly report success while the package goes
#' on using the other one. Everything here therefore operates on the
#' environment's prefix rather than its name.
#'
#' **Java and Perl** are checked too. They are hard requirements of the
#' AnnoTALE, PrediTALE and TALE-correction wrappers, they are not conda's
#' business, and otherwise they fail deep inside a `system()` call.
#'
#' @section Installing conda itself:
#' `conda = TRUE` installs a conda distribution if none is found. It is
#' deliberately opt-in and separate from `install`: putting a package manager
#' on someone's machine is a larger side effect than building an environment
#' in one that already exists. Note that this installs **miniconda**, via
#' [reticulate::install_miniconda()], not mamba.
#'
#' @param install Build the `tantale` environment if it is missing, and
#'   repair it if a pinned version is wrong. `FALSE` by default, so the
#'   function reports and changes nothing.
#' @param conda Install a conda distribution if none is found. `FALSE` by
#'   default; see the section above.
#' @param conda_bin Passed to `reticulate`. `"auto"` lets it choose.
#' @return Invisibly, a list with `conda` and `system` data frames of the
#'   checks, and `prefix`, so the result can be tested as well as read.
#' @seealso [tell_tales()] and [tales_align()], the two entry points that
#'   need these tools.
#' @export
tantale_setup <- function(install = FALSE, conda = FALSE, conda_bin = "auto") {

  ## conda itself ------------------------------------------------------------
  bin <- tryCatch(reticulate::conda_binary(conda_bin), error = function(e) NULL)
  if (is.null(bin)) {
    if (!isTRUE(conda)) {
      cli::cli_alert_danger("No conda or mamba installation found.")
      cli::cli_alert_info("Run {.run tantale_setup(install = TRUE, conda = TRUE)} to install one, or install mamba yourself.")
      return(invisible(list(conda = NULL, system = .tantale_check_system(),
                            prefix = NA_character_)))
    }
    cli::cli_alert_info("Installing miniconda (this is a one-off, and takes a few minutes)...")
    reticulate::install_miniconda()
    bin <- reticulate::conda_binary(conda_bin)
  }
  # computed into a variable: cli reads a leading dot inside braces as inline
  # markup, so {.tantale_conda_root(bin)} is parsed as a class, not a call
  root <- .tantale_conda_root(bin)
  cli::cli_alert_success("conda binary    {.file {bin}}")
  cli::cli_alert_info("default root    {.file {root}}")

  ## the environment ---------------------------------------------------------
  prefix <- tryCatch(.tantale_env_prefix(conda_bin = conda_bin),
                     error = function(e) NULL)

  if (is.null(prefix)) {
    if (!isTRUE(install)) {
      cli::cli_alert_danger("No {.val tantale} environment.")
      cli::cli_alert_info("Run {.run tantale_setup(install = TRUE)} to build it.")
      return(invisible(list(conda = NULL, system = .tantale_check_system(),
                            prefix = NA_character_)))
    }
    prefix <- .tantale_env_prefix(conda_bin = conda_bin)
  }
  cli::cli_alert_success("tantale env     {.file {prefix}}")

  ## versions ----------------------------------------------------------------
  pins <- .tantale_pins()
  checks <- .tantale_check_conda(pins, .tantale_installed(prefix))

  if (isTRUE(install) && any(!checks$ok)) {
    checks <- .tantale_repair(prefix, checks, conda_bin = conda_bin)
  }

  sys <- .tantale_check_system()
  width <- max(nchar(c(checks$tool, sys$tool)))
  for (i in seq_len(nrow(checks))) {
    r <- checks[i, ]
    label <- format(r$tool, width = width)
    if (isTRUE(r$ok)) {
      cli::cli_alert_success("{label} {r$found}")
    } else if (is.na(r$found)) {
      cli::cli_alert_danger("{label} missing (need {r$required})")
    } else {
      cli::cli_alert_danger("{label} {r$found} -- need {r$required}")
    }
  }

  ## java and perl -----------------------------------------------------------
  for (i in seq_len(nrow(sys))) {
    r <- sys[i, ]
    label <- format(r$tool, width = width)
    if (isTRUE(r$ok)) cli::cli_alert_success("{label} {r$found}")
    else cli::cli_alert_danger("{label} not on PATH -- needed by {r$needed_by}")
  }

  ## what to do next ---------------------------------------------------------
  if (any(!checks$ok) && !isTRUE(install)) {
    cli::cli_alert_info("Run {.run tantale_setup(install = TRUE)} to repair the environment.")
  } else if (all(checks$ok) && all(sys$ok)) {
    cli::cli_alert_success("Everything tantale needs is present.")
  }

  invisible(list(conda = checks, system = sys, prefix = prefix))
}


#' Bring an existing environment up to the pinned versions
#'
#' Separate from creation because the failure it addresses is specific:
#' `conda create` and `micromamba create` against an environment that already
#' exists will not *downgrade* a package, so an environment holding MAFFT
#' 7.520 stays on 7.520 and the create reports success. An explicit install of
#' the pinned version does downgrade it.
#'
#' Always addressed by prefix. `-n <name>` resolves against the binary's
#' default root, which is not necessarily the root the environment lives in.
#' @noRd
.tantale_repair <- function(prefix, checks, conda_bin = "auto") {
  wrong <- checks[!checks$ok, ]
  spec <- ifelse(is.na(wrong$required), wrong$tool,
                 paste0(wrong$tool, "=", wrong$required))
  cli::cli_alert_info("Repairing {nrow(wrong)} package{?s}: {.val {spec}}")
  reticulate::conda_install(
    envname = prefix,
    packages = spec,
    channel = c("conda-forge", "bioconda"),
    conda = conda_bin
  )
  .tantale_check_conda(.tantale_pins(), .tantale_installed(prefix))
}


#' Warn once per session if the environment does not match the pins
#'
#' The lazy path has to keep working -- users cannot be made to call
#' `tantale_setup()` -- so the check rides along with the first use of the
#' environment instead. Once per session, because it is advice rather than an
#' error, and repeating it on every call to an internal would be noise.
#' @noRd
.tantale_warn_if_unpinned <- function(prefix) {
  if (isTRUE(.tantale_state$version_checked)) return(invisible(NULL))
  .tantale_state$version_checked <- TRUE
  checks <- .tantale_check_conda(.tantale_pins(), .tantale_installed(prefix))
  bad <- checks[!checks$ok, ]
  if (!nrow(bad)) return(invisible(NULL))
  cli::cli_warn(c(
    "The {.val tantale} conda environment does not match the versions this package pins.",
    "!" = "{.val {bad$tool}}: found {.val {bad$found}}, need {.val {bad$required}}",
    "i" = "MAFFT in particular aligns TALE repeat strings differently after 7.4x, so results may differ.",
    "i" = "Run {.run tantale_setup(install = TRUE)} to repair it."
  ), class = "tantale_warning_conda_versions")
  invisible(NULL)
}

.tantale_state <- new.env(parent = emptyenv())


#' Where the conda binary would create an environment named with `-n`
#'
#' Not `dirname(dirname(bin))`: micromamba's binary usually sits in `~/bin`
#' while its root is wherever `MAMBA_ROOT_PREFIX` points, and the difference
#' between those two is the whole reason 7.4's rebuilds appeared to succeed
#' while the package kept using an environment in the other root.
#' @noRd
.tantale_conda_root <- function(bin) {
  if (grepl("mamba", basename(bin), fixed = TRUE)) {
    root <- Sys.getenv("MAMBA_ROOT_PREFIX", unset = "")
    return(if (nzchar(root)) root else "(MAMBA_ROOT_PREFIX unset)")
  }
  dirname(dirname(bin))
}
