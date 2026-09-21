# Retired: functal()
#
# Superseded by tales_compare_functal() (R/functal.R), an R reimplementation
# built on universalmotif::compare_motifs() rather than a wrapper around the
# vendored QueTAL Perl script this function called. The Perl script cannot
# work in this environment: it needs Bio::Perl, dropped by BioPerl in the
# 1.7 reorganisation, and no available perl/bioperl combination restores it
# (dev/restructuring-notes.md ledger section 12b has the full diagnosis).
# tales_compare_functal()'s numbers are not a port of this function's --
# compare_motifs()'s Pearson correlation over matched columns is a
# different statistic from FuncTAL's single correlation over the whole
# flattened, padded region (ledger 12b again, the 0.360 vs 0.489 finding).
#
# The vendored Perl tool itself moved alongside this file, from
# inst/tools/QueTAL_v1.1/ to inst/legacy/QueTAL_v1.1/.
#
# Kept rather than deleted per the project's standing "never delete code
# that looks dead" rule -- this is the only record of how the original
# QueTAL FuncTAL wrapper invoked the script (working directory and output
# directory quirks included), should anyone need it as a reference.
#
# Moved out of R/ on 2026-09-20. See dev/restructuring-notes.md section 12b.

#' Run functal from QueTAL to build a phylogenetic tree of TALE RVD sequences.
#'
#' A R wrapper around the \href{https://doi.org/10.3389/fpls.2015.00545}{QueTAL} 'functal' perl script.
#'
#' @param tal_file Path to a QueTAL-formatted file of TALE RVD sequences.
#' @param tree_format Tree layout passed to functal's `-n` option (default `"fan"`).
#' @param output_prefix Prefix used for functal's output file names.
#' @param output_dir Directory where output will be copied (default: current working directory).
#' @param functal_path Path to the functal perl script if you want to use another
#'   version than the one provided with tantale.
#' @param conda_bin Path to your Conda binary file if you need to specify a
#'   non-standard location, otherwise leave to "auto".
#' @return Returns invisibly the exit code of the shell call to functal (ie '0' if successful).
functal <- function(tal_file,
                    tree_format = "fan",
                    output_prefix = "FuncTALE",
                    output_dir = getwd(),
                    functal_path = system.file("legacy", "QueTAL_v1.1", "FuncTAL", "FuncTAL_v.1.1.pl", package = "tantale", mustWork = T),
                    conda_bin = "auto") {
  # Running functal from QueTAL_v1.1
  # A few observations:
  # Refuse to use another output directory than the "Ouputs" one in the program folder
  # Cannot invoke the program from another working directory than the one where the pl script is located
  # Crashes on some real-world CDS inputs; not yet isolated to a specific cause
  # So here is a caller function to get around these issues:

  functal_dir <- dirname(functal_path)

  # Assembling the command to be run
  # -I functal_dir is required for perl to find Statistics.pm, which ships
  # alongside the script rather than as an installed module. List::MoreUtils
  # and Bio::Perl come from the 'tantale' conda environment.
  functal_cmd <- paste(shQuote(.tantale_bin("perl", conda_bin = conda_bin)),
                       "-I", shQuote(functal_dir), shQuote(functal_path),
                      "-n", tree_format,
                      shQuote(tal_file),
                      shQuote(output_prefix))

  # Run the command inside the tantale conda environment
  cli::cli_inform(c("Running functal", " " = "{functal_cmd}"))
  envReady <- !as.logical(.create_tantale_env(conda_bin = conda_bin))
  if (envReady) {
    exitCom <- .tantale_exec(functal_cmd, cwd = functal_dir,
                             check = FALSE, what = "FuncTAL")
  } else {
    .abort_no_env("FuncTAL")
  }
  if (exitCom != 0) {
    cli::cli_abort(
      c("FuncTAL failed with perl exit status {exitCom}.",
        "i" = "See the console output above for what it reported."),
      class = c("tantale_error_functal_failed", "tantale_error"))
  }

  # Transferring the ouput to the output dir and deleting it in the functal "Outputs" directory
  functal_output_files <-
    list.files(file.path(functal_dir, "Outputs"), full.names = TRUE)
  file.copy(
    from = functal_output_files,
    to = output_dir,
    overwrite = FALSE,
    recursive = FALSE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  unlink(functal_output_files, recursive = TRUE, force = FALSE)
  return(invisible(exitCom))
}
