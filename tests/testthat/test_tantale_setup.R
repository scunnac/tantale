# tantale_setup() and the version checking behind it (ledger 7.4a).
#
# The parsing and comparison are pure, so they are tested against fixtures
# written here rather than against whatever happens to be installed. The one
# test that touches the real environment is at the bottom, and it fails
# rather than skips if conda is missing -- a check that silently does not run
# is worse than none.


#### reading the pins out of the yaml ####

test_that(".tantale_pins() parses the shipped environment file", {
  pins <- .tantale_pins()
  expect_true(all(c("mafft", "hmmer", "mmseqs2") %in% names(pins)))
  # the pin that matters: MAFFT after 7.4x aligns TALE strings differently
  expect_identical(unname(pins[["mafft"]]), "7.453")
})

test_that(".tantale_pins() ignores comments and the channels block", {
  f <- withr::local_tempfile(fileext = ".yaml")
  writeLines(c(
    "name: tantale",
    "channels:",
    "  - conda-forge",
    "  - bioconda",
    "dependencies:",
    "  - hmmer=3.3.2",
    "  # a comment explaining the pin below",
    "  - mafft=7.453",
    "  - somethingunpinned"
  ), f)
  pins <- .tantale_pins(f)
  # conda-forge and bioconda are channels, not dependencies
  expect_setequal(names(pins), c("hmmer", "mafft", "somethingunpinned"))
  expect_identical(unname(pins[["mafft"]]), "7.453")
  expect_true(is.na(pins[["somethingunpinned"]]))
})

test_that(".tantale_pins() strips a trailing comment from a dependency line", {
  f <- withr::local_tempfile(fileext = ".yaml")
  writeLines(c("dependencies:", "  - mafft=7.453  # not something newer"), f)
  expect_identical(unname(.tantale_pins(f)[["mafft"]]), "7.453")
})


#### reading what is installed out of conda-meta ####

test_that(".tantale_installed() reads name and version from conda-meta", {
  d <- withr::local_tempdir()
  dir.create(file.path(d, "conda-meta"))
  file.create(file.path(d, "conda-meta", c(
    "mafft-7.453-h516909a_1.json",
    "hmmer-3.3.2-hdbdd923_4.json",
    # a hyphenated package name, which a naive split would mangle
    "perl-statistics-r-0.34-pl5321_1.json",
    # a version that is not all digits
    "mmseqs2-14.7e284-pl5321h6a68c12_2.json"
  )))
  got <- .tantale_installed(d)
  expect_identical(unname(got[["mafft"]]), "7.453")
  expect_identical(unname(got[["perl-statistics-r"]]), "0.34")
  expect_identical(unname(got[["mmseqs2"]]), "14.7e284")
})

test_that(".tantale_installed() is empty rather than an error with no prefix", {
  expect_length(.tantale_installed(tempfile()), 0L)
})


#### comparing the two ####

test_that(".tantale_check_conda() flags a wrong version, not just a missing one", {
  pins <- c(mafft = "7.453", hmmer = "3.3.2", clustalo = "1.2.4")
  # the exact situation 7.4a exists for: an environment built by an older
  # version of the package, holding a MAFFT that aligns differently
  installed <- c(mafft = "7.520", hmmer = "3.3.2")
  out <- .tantale_check_conda(pins, installed)

  expect_false(out$ok[out$tool == "mafft"])      # present but wrong
  expect_true(out$ok[out$tool == "hmmer"])
  expect_false(out$ok[out$tool == "clustalo"])   # absent
  expect_true(is.na(out$found[out$tool == "clustalo"]))
  expect_identical(out$found[out$tool == "mafft"], "7.520")
})

test_that(".tantale_check_conda() accepts any version for an unpinned package", {
  out <- .tantale_check_conda(c(thing = NA_character_), c(thing = "1.0"))
  expect_true(out$ok)
})


#### the conda root ####

test_that(".tantale_conda_root() does not assume the binary sits under the root", {
  # micromamba usually lives in ~/bin while its root is elsewhere; treating
  # dirname(dirname(bin)) as the root is what made 7.4's rebuilds land in the
  # wrong place.
  withr::local_envvar(MAMBA_ROOT_PREFIX = "/somewhere/else")
  expect_identical(.tantale_conda_root("/home/u/bin/micromamba"), "/somewhere/else")
  # plain conda does keep its binary under the root
  expect_identical(.tantale_conda_root("/opt/miniconda3/condabin/conda"),
                   "/opt/miniconda3")
})

test_that(".tantale_conda_root() says so when the root is unset", {
  withr::local_envvar(MAMBA_ROOT_PREFIX = NA)
  expect_match(.tantale_conda_root("/home/u/bin/micromamba"), "unset")
})


#### the system tools ####

test_that(".tantale_check_system() looks for java and perl and says what needs them", {
  out <- .tantale_check_system()
  expect_setequal(out$tool, c("java", "perl"))
  expect_true(all(nzchar(out$needed_by)))
})


#### against the real environment ####

test_that("tantale_setup() reports the real environment and changes nothing", {
  # Fails rather than skips: this is the function whose whole purpose is to
  # notice a wrong environment, so a run where it cannot look is a failure.
  out <- tantale_setup()

  expect_type(out, "list")
  expect_setequal(names(out), c("conda", "system", "prefix"))
  expect_true(dir.exists(out$prefix))
  # every pinned package present and at the pinned version
  expect_true(all(out$conda$ok),
              info = paste("unmet pins:",
                           paste(out$conda$tool[!out$conda$ok], collapse = ", ")))
  expect_setequal(out$conda$tool, names(.tantale_pins()))
})


#### Resolving and running the environment's programs (ledger 12) ####

test_that(".tantale_bin() resolves into the environment, not the PATH", {
  # This machine carries /usr/bin/mafft 7.505 and /usr/bin/nhmmer 3.4
  # against pins of 7.453 and 3.3.2, so "which binary" is not academic.
  bins <- .tantale_bin(c("mafft", "nhmmer", "mmseqs", "perl"))
  prefix <- .tantale_env_prefix()
  expect_true(all(startsWith(unname(bins), prefix)))
  expect_true(all(file.exists(bins)))
  expect_named(bins, c("mafft", "nhmmer", "mmseqs", "perl"))
})

test_that(".tantale_bin() knows the tools that are not in bin/", {
  bins <- .tantale_bin(c("hex2maffttext", "maffttext2hex"))
  expect_true(all(grepl("libexec/mafft/", bins, fixed = TRUE)))
  expect_true(all(file.exists(bins)))
})

test_that(".tantale_bin() names everything that is missing", {
  expect_error(.tantale_bin(c("mafft", "no_such_tool", "other_missing_tool")),
               class = "tantale_error_tool_missing")
})

test_that("the resolved binaries are the pinned versions, not the system ones", {
  v <- system2(.tantale_bin("mafft"), "--version", stdout = TRUE, stderr = TRUE)
  expect_match(paste(v, collapse = " "), .tantale_pins()[["mafft"]], fixed = TRUE)
})

test_that(".tantale_exec() reports a failure instead of returning quietly", {
  expect_error(.tantale_exec("exit 3", what = "a deliberately failing command"),
               class = "tantale_error_exec_failed")
  expect_invisible(.tantale_exec("true"))
})

test_that(".tantale_exec() runs every command of a compound string", {
  # The failure mode `conda run` has: with it, only the first command lands
  # inside the environment and the rest fall back to the system PATH.
  out <- withr::local_tempfile()
  .tantale_exec(paste("echo one >", shQuote(out), "&& echo two >>", shQuote(out)))
  expect_identical(readLines(out), c("one", "two"))
})

test_that(".tantale_exec() honours cwd and captures stderr", {
  d <- withr::local_tempdir()
  .tantale_exec("pwd > here.txt", cwd = d)
  expect_identical(normalizePath(readLines(file.path(d, "here.txt"))),
                   normalizePath(d))

  err <- withr::local_tempfile()
  expect_error(.tantale_exec("echo 'the reason' >&2; exit 1", stderr_file = err),
               "the reason")
})


#### Choosing between environments of the same name ####

test_that(".tantale_pick_env() chooses by the pins, not by order", {
  good <- .tantale_env_prefix()
  bad <- file.path(tempdir(), "empty_env")   # no conda-meta -> satisfies nothing
  dir.create(bad, showWarnings = FALSE)
  # the failing candidate first, so a "take the first" rule would pick wrong
  expect_identical(.tantale_pick_env(c(bad, good)), good)
})

test_that(".tantale_pick_env() refuses to guess when none matches the pins", {
  a <- file.path(tempdir(), "empty_a"); dir.create(a, showWarnings = FALSE)
  b <- file.path(tempdir(), "empty_b"); dir.create(b, showWarnings = FALSE)
  expect_error(.tantale_pick_env(c(a, b)), class = "tantale_error_conda_env")
})

test_that(".tantale_pick_env() lets the user override", {
  withr::local_options(tantale.env_prefix = "/chosen/by/hand")
  expect_identical(.tantale_pick_env(c("/a", "/b")), "/chosen/by/hand")
})
