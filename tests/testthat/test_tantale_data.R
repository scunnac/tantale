# The downloaded archives (ledger §50). These tests build a small archive of
# their own and install it into a temporary TANTALE_DATA_DIR, so they touch
# neither the network nor the real installation.

fake_archive <- function(dir) {
  stage <- file.path(withr::local_tempdir(.local_envir = parent.frame()), "fake-1")
  dir.create(stage)
  writeLines("hello", file.path(stage, "a.txt"))
  writeLines(paste(digest::digest(file.path(stage, "a.txt"), file = TRUE, algo = "sha256"),
                   "a.txt", sep = "  "),
             file.path(stage, "MANIFEST"))
  archive <- file.path(dir, "tantale-fake-1.tar.gz")
  withr::with_dir(dirname(stage),
                  utils::tar(archive, files = "fake-1", compression = "gzip", tar = "internal"))
  list(fake = list(version = "fake-1", file = basename(archive),
                   sha256 = digest::digest(archive, file = TRUE, algo = "sha256"),
                   what = "a test archive"))
}

test_that("an archive installs from a local directory and checks out", {
  withr::local_envvar(TANTALE_DATA_DIR = withr::local_tempdir())
  src <- withr::local_tempdir()
  spec <- fake_archive(src)
  expect_false(.tantale_archive_ok("fake", spec))
  dest <- .tantale_install_archive("fake", archive_dir = src, archives = spec)
  expect_identical(normalizePath(dest), normalizePath(file.path(Sys.getenv("TANTALE_DATA_DIR"), "fake-1")))
  expect_true(.tantale_archive_ok("fake", spec))
  # an altered file is caught
  writeLines("changed", file.path(dest, "a.txt"))
  expect_false(.tantale_archive_ok("fake", spec))
})

test_that("an archive with the wrong sha256 is refused", {
  withr::local_envvar(TANTALE_DATA_DIR = withr::local_tempdir())
  src <- withr::local_tempdir()
  spec <- fake_archive(src)
  spec$fake$sha256 <- strrep("0", 64)
  expect_error(.tantale_install_archive("fake", archive_dir = src, archives = spec),
               class = "tantale_error_archive_checksum")
})

test_that("a missing archive file or download is a classed error", {
  withr::local_envvar(TANTALE_DATA_DIR = withr::local_tempdir())
  spec <- fake_archive(withr::local_tempdir())
  expect_error(.tantale_install_archive("fake", archive_dir = withr::local_tempdir(), archives = spec),
               class = "tantale_error_missing_file")
  expect_error(suppressMessages(.tantale_install_archive("fake", archives = spec,
                                                         base_url = "https://127.0.0.1:9/none")),
               class = "tantale_error_download")
})

test_that("tools and genomes not installed are classed errors", {
  withr::local_envvar(TANTALE_DATA_DIR = withr::local_tempdir())
  expect_error(.tantale_tool("annotale"), class = "tantale_error_tool_missing")
  expect_error(tantale_genome("MAI1"), class = "tantale_error_genome_missing")
  expect_error(tantale_genome("nope"))
})

test_that("the installed tools and genomes are found", {
  # the real installation, which the rest of the suite needs anyway
  for (tool in c("annotale", "preditale", "talecorrection", "talecorrection_hmm")) {
    expect_true(file.exists(.tantale_tool(tool)))
  }
  for (g in c("MAI1", "BAI3", "BAI3-1-1", "PXO86")) {
    expect_true(file.exists(tantale_genome(g)))
  }
  expect_true(.tantale_archive_ok("tools"))
  expect_true(.tantale_archive_ok("genomes"))
})
