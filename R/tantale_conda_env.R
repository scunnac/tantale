
systemInCondaEnv <- function(envName, command,
                             condaBinPath = "auto",
                             cwd = getwd(),
                             ...) {
  condaBinPath <- reticulate::conda_binary(condaBinPath)
  # activateEnvCmd <- glue::glue("eval \"$({condaBinPath} shell hook -s posix)\"",
  #                              "; micromamba activate {envName}")
  # fullCommand <- glue::glue_collapse(c(activateEnvCmd, command), sep = "; ")
  fullCommand <- glue::glue("eval \"$({condaBinPath} shell hook -s posix)\"",
                            "{condaBinPath} run --cwd {cwd} -n {envName} {command}",
                            .sep = "; ")
  logger::log_debug("Starting the following command in the '{envName}' conda env :
                   {fullCommand}")
  system(command = fullCommand, ...)
}

# systemInCondaEnv <- function(envName, command,
#                              condaBinPath = "auto",
#                              intern = FALSE) {
#   logger::log_debug("Starting the following command in the '{envName}' conda env :
#                    {command}")
#   reticulate::conda_run2( conda = condaBinPath,
#                           envname = envName,
#                           cmd_line = command,
#                           intern = intern,
#                           echo = FALSE)
# }





createTantaleEnv <- function(condaBinPath = "auto") {
  envName <- "tantale"
  if (!envName %in% (reticulate::conda_list(conda = condaBinPath)["name"] %>% unlist())) {
    logger::log_warn("A custom conda env will be installed on your system to run external dependencies...")
    condayml <- system.file("tools", "tantale_conda_env.yaml", package = "tantale", mustWork = T)
    res <- reticulate::conda_create(envname = envName,
                                    environment = condayml)
    if (!is.character(res)) {
      logger::log_warn("Installation of the conda environment failed.")
      return(invisible(res))
    }
    return(invisible(0L))
  } else {
    logger::log_info("A Conda environment with the name '{envName}' has been found on your system and can be used for analysis.")
    return(invisible(0L))
  }
}

# reticulate::condaenv_exists(envname = envName, conda = condaBinPath)
# reticulate::conda_remove(envname = envName, conda = condaBinPath)
# reticulate::conda_list(conda = condaBinPath)

# reticulate::conda_binary()
# reticulate::conda_list(conda = "/home/cunnac/bin/miniconda3/condabin/conda")["name"] %>% unlist()
# createTantaleEnv(condaBinPath = "/home/cunnac/bin/miniconda3/condabin/conda")
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





