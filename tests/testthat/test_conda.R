
require_tantale_conda_env <- function() {
  testthat::skip_on_cran()
  condaYaml <- system.file("tools", "tantale_conda_env.yaml", package = "tantale", mustWork = FALSE)
  installHint <- paste(
    "The 'tantale' conda environment is required to run this test but was not found.",
    "To set it up, either run:",
    "    tantale:::.create_tantale_env()",
    "(needs conda or mamba on your PATH), or create it manually with:",
    sprintf("    conda env create -n tantale -f %s", condaYaml),
    sep = "\n"
  )
  conda <- tryCatch(reticulate::conda_binary("auto"), error = function(e) NA_character_)
  if (is.na(conda)) {
    testthat::fail(paste("No conda installation found on this machine.", installHint, sep = "\n"))
  }
  envs <- tryCatch(unlist(reticulate::conda_list(conda = conda)["name"]), error = function(e) character())
  if (!"tantale" %in% envs) {
    testthat::fail(installHint)
  }
}

test_that(".run_in_conda does not error with a test command",
          {
            require_tantale_conda_env()
            expect_true(0 == .run_in_conda(env_name = "tantale", command = "mmseqs createdb --help", intern = FALSE))
          }
)

test_that(".run_in_conda returns a character string corresponding to stdoutput with a test command",
          {
            require_tantale_conda_env()
            expect_true(is.character(.run_in_conda(env_name = "tantale", command = "mmseqs createdb --help", intern = TRUE)))
          }
)

