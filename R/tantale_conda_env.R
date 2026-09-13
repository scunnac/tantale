
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
    cli::cli_inform("A Conda environment with the name {.val {env_name}} has been found on your system and can be used for analysis.")
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





