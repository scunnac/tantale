
test_that("systemInCondaEnv does not error with a test command",
          {expect_true(0 == systemInCondaEnv(envName = "tantale", command = "mmseqs createdb --help", intern = FALSE))}
)

test_that("systemInCondaEnv returns a character string corresponding to stdoutput with a test command",
          {expect_true(is.character(systemInCondaEnv(envName = "tantale", command = "mmseqs createdb --help", intern = TRUE)))}
)



