
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
    python <- .tantale_pick_env(dirname(dirname(python)))
  }
  prefix <- if (length(python) == 1L) dirname(dirname(python)) else python
  # The pins are only real if something checks them. tantale_setup() is the
  # place to do it deliberately; this is the safety net for users who never
  # call it, which is most of them (7.4a).
  .tantale_warn_if_unpinned(prefix)
  prefix
}


#### Running the environment's programs ####
#
# One mechanism: resolve the executable to an absolute path inside the
# environment, then run it. No activation, no PATH, no `conda run`.
#
# That is not a stylistic preference, it is the only approach that is
# actually safe here, for three measured reasons.
#
#   1. This machine carries /usr/bin/mafft 7.505 and /usr/bin/nhmmer 3.4,
#      against pins of 7.453 and 3.3.2. Anything PATH-based picks the wrong
#      one whenever the environment is not on the front of PATH.
#   2. `conda run` only puts the *first* command of a compound string inside
#      the environment: the `;` and `&&` are eaten by the outer shell before
#      conda ever sees them. talecorrection() joins three nhmmer calls with
#      "; ", so two of them were running against system HMMER 3.4 with
#      nothing said. Absolute paths make shell splitting irrelevant.
#   3. `conda run` costs ~1.3 s per call against ~36 ms for a direct call,
#      and mmseqs2 alone makes four calls in a row.
#
# Verified that nothing needs the environment's variables: `env -i` runs of
# mafft, mmseqs and perl all work, and the env's perl resolves its own @INC
# (conda bakes the prefix into its binaries at build time).


# Executables that do not live in <prefix>/bin. MAFFT's text mode needs the
# two hex converters, which conda installs under libexec.
.TANTALE_BIN_SUBDIR <- c(hex2maffttext = "libexec/mafft",
                         maffttext2hex = "libexec/mafft")


#' Absolute paths to executables inside the tantale environment
#'
#' @param tools Character vector of program names.
#' @param prefix The environment prefix; resolved if not supplied.
#' @param conda_bin Passed to \code{reticulate}.
#' @return A named character vector of absolute paths, one per tool.
#' @noRd
.tantale_bin <- function(tools, prefix = NULL, conda_bin = "auto") {
  if (is.null(prefix)) prefix <- .tantale_env_prefix(conda_bin = conda_bin)
  subdir <- ifelse(tools %in% names(.TANTALE_BIN_SUBDIR),
                   .TANTALE_BIN_SUBDIR[tools], "bin")
  paths <- file.path(prefix, subdir, tools)
  names(paths) <- tools
  missing <- tools[!file.exists(paths)]
  if (length(missing) > 0L) {
    cli::cli_abort(
      c("Cannot find {length(missing)} of the programs tantale needs.",
        "x" = "Missing: {.file {unname(paths[missing])}}",
        "i" = "{.run tantale_setup()} reports what the environment holds."),
      class = c("tantale_error_tool_missing", "tantale_error"))
  }
  paths
}


#' Run a command built from .tantale_bin() paths
#'
#' @param command A shell command. Every executable in it must already be an
#'   absolute path from \code{.tantale_bin()} -- this function does nothing
#'   to the environment, which is the point.
#' @param cwd Directory to run in, or \code{NULL}.
#' @param stderr_file Capture stderr here rather than letting it through.
#'   Replayed in the error message if the command fails.
#' @param check Abort on a non-zero exit status.
#' @param what What to call the command in messages.
#' @return The exit status, invisibly.
#' @noRd
.tantale_exec <- function(command, cwd = NULL, stderr_file = NULL,
                          check = TRUE, what = "command") {
  if (!is.null(cwd)) {
    command <- paste0("cd ", shQuote(cwd), " && ", command)
  }
  if (!is.null(stderr_file)) {
    command <- paste0("{ ", command, " ; } 2> ", shQuote(stderr_file))
  }
  status <- system(command = command, intern = FALSE)
  if (check && !identical(as.integer(status), 0L)) {
    saidWhy <- if (!is.null(stderr_file) && file.exists(stderr_file)) {
      utils::tail(readLines(stderr_file, warn = FALSE), 20)
    } else character()
    cli::cli_abort(
      c("{what} failed.",
        "x" = "Exit status {status}.",
        if (length(saidWhy)) c("i" = "It said:"),
        stats::setNames(saidWhy, rep(" ", length(saidWhy)))),
      class = c("tantale_error_exec_failed", "tantale_error"))
  }
  invisible(status)
}


#' The three MAFFT executables, from the environment or a standalone install
#'
#' MAFFT's text mode needs three programs, not one: the aligner and the two
#' converters that move sequences in and out of its hexadecimal encoding.
#'
#' The conda case is just \code{.tantale_bin()}. This exists for the other
#' case -- \code{tales_align(mafft_path = )}, where the user points at a
#' MAFFT outside the environment, which is laid out differently from the
#' conda one and cannot be resolved from a prefix.
#'
#' @param mafft_path \code{NULL} for the tantale environment, or the root of
#'   a standalone MAFFT directory (one holding \code{mafft.bat} with the
#'   helpers under \code{mafftdir/libexec}).
#' @return A named character vector of three absolute paths.
#' @noRd
.mafft_paths <- function(mafft_path = NULL, conda_bin = "auto") {
  tools <- c("mafft", "hex2maffttext", "maffttext2hex")
  if (is.null(mafft_path)) return(.tantale_bin(tools, conda_bin = conda_bin))

  paths <- c(mafft = file.path(mafft_path, "mafft.bat"),
             hex2maffttext = file.path(mafft_path, "mafftdir", "libexec", "hex2maffttext"),
             maffttext2hex = file.path(mafft_path, "mafftdir", "libexec", "maffttext2hex"))
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing) > 0L) {
    cli::cli_abort(
      c("Cannot find {length(missing)} of MAFFT's executable{?s} under {.file {mafft_path}}.",
        "x" = "Missing: {.file {unname(paths[missing])}}",
        "i" = "Text-mode alignment needs the two hex converters as well as {.file mafft} itself."),
      class = c("tantale_error_mafft_missing", "tantale_error"))
  }
  paths
}


#' Choose between several environments all named "tantale"
#'
#' conda and micromamba keep separate roots and \code{reticulate} scans
#' several of them, so the same environment name can exist more than once.
#'
#' The previous rule here was "prefer the one under the conda binary's
#' root", computed as \code{dirname(dirname(binary))}. That is not a root:
#' micromamba's binary usually lives in \code{~/bin}, so the expression
#' yields the home directory and matches every candidate, leaving the choice
#' to whatever order \code{conda_list()} returned -- silently, because the
#' warning only fired when nothing matched.
#'
#' The criterion used instead is the one that actually matters: which
#' environment holds the tools at the versions this package pins. That is
#' not a heuristic, it is the requirement. Picking a \code{tantale}
#' environment with MAFFT 7.520 in it would produce different alignments
#' from the same input with nothing to indicate it (restructuring-notes.md
#' 7.4a).
#'
#' @param prefixes Candidate environment prefixes.
#' @return One prefix.
#' @noRd
.tantale_pick_env <- function(prefixes) {
  override <- getOption("tantale.env_prefix")
  if (!is.null(override)) return(override)

  pins <- .tantale_pins()
  satisfies <- vapply(prefixes, function(p) {
    all(.tantale_check_conda(pins, .tantale_installed(p))$ok)
  }, logical(1))

  if (sum(satisfies) == 1L) return(prefixes[satisfies])

  if (sum(satisfies) > 1L) {
    # Interchangeable as far as this package is concerned; say which, so a
    # surprising choice is at least visible.
    cli::cli_warn(
      c("{sum(satisfies)} environments named {.val tantale} match the pinned versions.",
        "i" = "Using {.file {prefixes[satisfies][1]}}.",
        "i" = "Set {.code options(tantale.env_prefix = )} to choose."),
      class = "tantale_warning_multiple_envs")
    return(prefixes[satisfies][1])
  }

  cli::cli_abort(
    c("{length(prefixes)} environments are named {.val tantale}, and none matches the pinned versions.",
      "x" = "Candidates: {.file {prefixes}}",
      "i" = "Using one of them anyway risks results from the wrong tool versions.",
      "i" = "Run {.run tantale_setup(install = TRUE)}, or set {.code options(tantale.env_prefix = )}."),
    class = c("tantale_error_conda_env", "tantale_error"))
}


#' The environment could not be created
#'
#' Four call sites raised this with four near-identical hand-written
#' sentences, none of which said what to do about it. They now share one,
#' which points at the function written for exactly this situation.
#'
#' @param what The tool that was about to be run.
#' @noRd
.abort_no_env <- function(what) {
  cli::cli_abort(
    c("Could not create the {.val tantale} conda environment, needed to run {what}.",
      "i" = "{.run tantale_setup()} reports what is present and what is missing.",
      "i" = "{.run tantale_setup(install = TRUE)} builds or repairs it."),
    class = c("tantale_error_conda_env", "tantale_error"))
}
