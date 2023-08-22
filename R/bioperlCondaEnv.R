
systemInCondaEnv <- function(envName, command,
                             condaBinPath = "auto",
                             ...) {
  activateEnvCmd <- glue::glue("eval \"$({reticulate::conda_binary(conda = condaBinPath)} shell.bash hook)\"; conda activate {envName}")
  fullCommand <- glue::glue_collapse(c(activateEnvCmd, command), sep = "; ")
  logger::log_debug("Starting the following command in the '{envName}' conda env :
                   {fullCommand}")
  system(command = fullCommand, ...)
}


createBioPerlEnv <- function(condaBinPath = "auto") {
  envName <- "perlforal"
  if (!envName %in% (reticulate::conda_list(conda = condaBinPath)["name"] %>% unlist())) {
    logger::log_warn("A custom conda env will be installed on your system to run external dependencies...")
    condayml <- system.file("tools", "bioperl_conda_env.yaml", package = "tantale", mustWork = T)
    res <- systemInCondaEnv(envName = "base",
                            condaBinPath = condaBinPath,
                            command = glue::glue("mamba env create -f {condayml}"))
    if (!res == 0L) {
      logger::log_warn("Installation of the conda environment failed.")
    }
    return(invisible(res))
  } else {
    logger::log_info("A Conda environment with the name '{envName}' has been found on your system and can be used for analysis.")
    return(invisible(0L))
  }
}

# reticulate::conda_binary()
# reticulate::conda_list(conda = "/home/cunnac/bin/miniconda3/condabin/conda")["name"] %>% unlist()
# createBioPerlEnv(condaBinPath = "/home/cunnac/bin/miniconda3/condabin/conda")
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





