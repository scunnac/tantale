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
