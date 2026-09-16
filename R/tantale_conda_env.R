
.run_in_conda <- function(env_name, command,
                             conda_bin = "auto",
                             cwd = getwd(),
                             ...) {
  conda_bin <- reticulate::conda_binary(conda_bin)
  # activateEnvCmd <- glue::glue("eval \"$({conda_bin} shell hook -s posix)\"",
  #                              "; micromamba activate {env_name}")
  # fullCommand <- glue::glue_collapse(c(activateEnvCmd, command), sep = "; ")
  fullCommand <- glue::glue("eval \"$({conda_bin} shell hook -s posix)\"",
                            "{conda_bin} run --cwd {cwd} -n {env_name} {command}",
                            .sep = "; ")
  system(command = fullCommand, ...)
}

# .run_in_conda <- function(env_name, command,
#                              conda_bin = "auto",
#                              intern = FALSE) {
#   logger::log_debug("Starting the following command in the '{env_name}' conda env :
#                    {command}")
#   reticulate::conda_run2( conda = conda_bin,
#                           envname = env_name,
#                           cmd_line = command,
#                           intern = intern,
#                           echo = FALSE)
# }





.create_tantale_env <- function(conda_bin = "auto") {
  env_name <- "tantale"
  if (!env_name %in% (reticulate::conda_list(conda = conda_bin)["name"] %>% unlist())) {
    cli::cli_inform("A custom conda env will be installed on your system to run external dependencies...")
    condayml <- system.file("tools", "tantale_conda_env.yaml", package = "tantale", mustWork = T)
    res <- reticulate::conda_create(envname = env_name,
                                    environment = condayml)
    if (!is.character(res)) {
      cli::cli_warn("Installation of the conda environment failed.")
      return(invisible(res))
    }
    return(invisible(0L))
  } else {
    # Deliberately silent. This used to announce that the environment "can be
    # used for analysis" -- on every run, and without having looked inside it.
    # Both halves were wrong: it is noise when true, and a false assurance
    # when the environment holds the wrong MAFFT (7.4a).
    return(invisible(0L))
  }
}

# reticulate::condaenv_exists(envname = env_name, conda = conda_bin)
# reticulate::conda_remove(envname = env_name, conda = conda_bin)
# reticulate::conda_list(conda = conda_bin)

# reticulate::conda_binary()
# reticulate::conda_list(conda = "/home/cunnac/bin/miniconda3/condabin/conda")["name"] %>% unlist()
# .create_tantale_env(conda_bin = "/home/cunnac/bin/miniconda3/condabin/conda")
# #perl-data-dumper

# use warnings;
# use strict;
# use Getopt::Std;
# use Statistics::R;
# use List::MoreUtils qw(uniq);
# use List::Util qw( min max );
# use Algorithm::NeedlemanWunsch;
# use Bio::Perl;
# use Statistics::Basic qw(:all);
# use List::Util qw( min max );
# use POSIX qw(ceil);


#' Locate the tantale conda environment, creating it if necessary
#'
#' The environment is where the package's external tools live -- MAFFT,
#' HMMER, mmseqs2 and the Perl dependencies of the target predictors. They are
#' no longer shipped inside the package, so this is how they are found.
#'
#' Versions are pinned in \code{inst/tools/tantale_conda_env.yaml} and the
#' pins matter: MAFFT in particular aligns TALE repeat strings differently
#' after 7.4x (restructuring-notes.md 7.4). An environment created by an older
#' version of this package may hold the wrong ones, so the pins are checked
#' here rather than trusted.
#'
#' @param conda_bin Passed to \code{reticulate}.
#' @return The environment's prefix directory.
#' @noRd
.tantale_env_prefix <- function(conda_bin = "auto") {
  if (as.logical(.create_tantale_env(conda_bin = conda_bin))) {
    cli::cli_abort(
      c("Could not create the {.val tantale} conda environment.",
        "i" = "It provides MAFFT, HMMER and mmseqs2, which this package needs.",
        "i" = "Check that conda or mamba is installed and that you are online.",
        "i" = "{.run tantale_setup()} reports what is missing."),
      class = c("tantale_error_conda_env", "tantale_error"))
  }
  envs <- reticulate::conda_list(conda = conda_bin)
  python <- envs$python[envs$name == "tantale"]
  if (length(python) == 0L) {
    cli::cli_abort("No conda environment named {.val tantale} was found.",
                   class = c("tantale_error_conda_env", "tantale_error"))
  }
  if (length(python) > 1L) {
    # conda and micromamba keep separate roots, so the same environment name
    # can exist in both. Prefer the one belonging to the binary in use.
    root <- dirname(dirname(reticulate::conda_binary(conda_bin)))
    owned <- python[startsWith(python, root)]
    if (length(owned) >= 1L) {
      python <- owned[1]
    } else {
      cli::cli_warn(c("{length(python)} conda environments are named {.val tantale}.",
                      "i" = "Using {.file {dirname(dirname(python[1]))}}."))
      python <- python[1]
    }
  }
  prefix <- dirname(dirname(python))
  # The pins are only real if something checks them. tantale_setup() is the
  # place to do it deliberately; this is the safety net for users who never
  # call it, which is most of them (7.4a).
  .tantale_warn_if_unpinned(prefix)
  prefix
}


#' Paths to the MAFFT executables
#'
#' MAFFT's text mode needs three programs, not one: the aligner itself and the
#' two converters that move sequences in and out of its hexadecimal encoding.
#' They sit in different places depending on how MAFFT was installed, so the
#' layout is resolved here rather than assumed.
#'
#' @param mafft_path \code{NULL} to use the conda environment, or the root of
#'   a standalone MAFFT directory (one holding \code{mafft.bat} with the
#'   helpers under \code{mafftdir/libexec}).
#' @param conda_bin Passed to \code{reticulate}.
#' @return A list of the three executable paths.
#' @noRd
.mafft_binaries <- function(mafft_path = NULL, conda_bin = "auto") {
  if (is.null(mafft_path)) {
    prefix <- .tantale_env_prefix(conda_bin = conda_bin)
    bins <- list(mafft = file.path(prefix, "bin", "mafft"),
                 hex2text = file.path(prefix, "libexec", "mafft", "hex2maffttext"),
                 text2hex = file.path(prefix, "libexec", "mafft", "maffttext2hex"))
  } else {
    # a standalone MAFFT directory, laid out as the distributed archives are
    bins <- list(mafft = file.path(mafft_path, "mafft.bat"),
                 hex2text = file.path(mafft_path, "mafftdir", "libexec", "hex2maffttext"),
                 text2hex = file.path(mafft_path, "mafftdir", "libexec", "maffttext2hex"))
  }
  missing <- names(bins)[!vapply(bins, file.exists, logical(1))]
  if (length(missing) > 0L) {
    cli::cli_abort(
      c("Cannot find {length(missing)} of MAFFT's executable{?s}.",
        "x" = "Missing: {.file {unlist(bins[missing])}}",
        "i" = "Text-mode alignment needs the two hex converters as well as {.file mafft} itself."),
      class = c("tantale_error_mafft_missing", "tantale_error"))
  }
  bins
}
