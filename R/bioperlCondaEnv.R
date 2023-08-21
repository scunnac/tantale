#' Run a bash command in a specified Conda environment
#'
#'
#' @param envName A character string with the name of an installed Conda environment
#' @param command A Bash command to execute in this env with system() as a character string
#' @param condaBinPath Path the conda binary to be used
#' @param ... Additional parameters for \code{system()}
#' @return The \code{system()} return value
#' @export
systemInCondaEnv <- function(envName, command,
                             condaBinPath = "auto",
                             ...) {
  # Found this 'conda shell.bash hook ' trick for the 'conda activate' command to run smoothly
  # in a bash call at the end of this post:
  # https://stackoverflow.com/questions/55854189/how-to-activate-anaconda-environment-within-r/55854475
  activateEnvCmd <- glue::glue("eval \"$({reticulate::conda_binary(conda = condaBinPath)} shell.bash hook)\"; conda activate {envName}")
  fullCommand <- glue::glue_collapse(c(activateEnvCmd, command), sep = "; ")
  logger::log_debug("Starting the following command in the '{envName}' conda env :
                   {fullCommand}")
  system(command = fullCommand, ...)
}

#' Install Conda environments on a system
#'
#'
#' @param customCondaEnvYaml A named vector of Conda env recipe yaml file paths. Names correspond to
#' env names in the recipes.
#' @param condaBinPath Path to the conda binary to be used
#' @param ... Additional parameters for \code{system()}
#' @return The \code{system()} return value
#' @export
condaEnvCheckInstall <- function(customCondaEnvYaml, condaBinPath = "auto")
{
  condaEnvs <- reticulate::conda_list(conda = condaBinPath)
  if (!all(names(customCondaEnvYaml) %in% condaEnvs$name)) {
    condaToInstall <- customCondaEnvYaml[!names(customCondaEnvYaml) %in% condaEnvs$name]
    logger::log_warn("This custom conda env recipes are not installed on your system : {condaToInstall}")
    res <- lapply(condaToInstall, function(condayml) {
      systemInCondaEnv(envName = "base",
                       condaBinPath = reticulate::conda_binary(conda = condaBinPath),
                       command = glue::glue("mamba env create -f {condayml}"))
    })
    if (!all(sapply(res, function(x) x == 0L))) {
      logger::log_warn("Installation of a conda environment failed.")
    }
  } else {
    logger::log_info("A Conda environment with the name '{names(customCondaEnvYaml)}' has been found on your system and can be used for analysis.")
  }
}
