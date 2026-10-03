#### The tools and genomes tantale downloads ####
#
# The three Java programs tantale wraps (58 MB) and the four example genomes
# (19 MB) made the source package 57 MB, against the 5 MB a package is
# expected to stay under. They now live in two archives attached to GitHub
# releases of the package's repository (tools-1, genomes-1), built by
# dev/make_archives.R. tantale_setup(install = TRUE) downloads them, checks
# their sha256 against the pins below and unpacks them into the user's data
# directory (ledger §50). A published archive never changes: a new version
# of any file means a new archive (tools-2, ...) and a new pin here.

.tantale_archives <- list(
  tools = list(
    version = "tools-1",
    file = "tantale-tools-1.tar.gz",
    sha256 = "6ebdfbfe66e269fcdb1285b64f06dffaf30f466ef3c234f6a90e81d089627501",
    what = "AnnoTALE, PrediTALE and TALEcorrection"
  ),
  genomes = list(
    version = "genomes-1",
    file = "tantale-genomes-1.tar.gz",
    sha256 = "f7eb357025fbabec836b6b58323bcfa692319b94de989d303f79f0b963105ff4",
    what = "the four example genomes"
  )
)


#' Where the downloaded tools and genomes live
#'
#' `TANTALE_DATA_DIR` when set (a shared or larger disk, say), the user's R
#' data directory otherwise.
#' @noRd
.tantale_data_dir <- function() {
  dir <- Sys.getenv("TANTALE_DATA_DIR")
  if (nzchar(dir)) dir else tools::R_user_dir("tantale", "data")
}

.tantale_archive_path <- function(which, archives = .tantale_archives) {
  file.path(.tantale_data_dir(), archives[[which]]$version)
}


#' Is an unpacked archive complete and unaltered?
#'
#' Checks every file listed in its `MANIFEST` against the recorded sha256.
#' @return `TRUE` or `FALSE`.
#' @noRd
.tantale_archive_ok <- function(which, archives = .tantale_archives) {
  dir <- .tantale_archive_path(which, archives)
  manifest <- file.path(dir, "MANIFEST")
  if (!file.exists(manifest)) return(FALSE)
  lines <- readLines(manifest, warn = FALSE)
  hashes <- sub("  .*$", "", lines)
  files <- file.path(dir, sub("^[^ ]*  ", "", lines))
  all(file.exists(files)) &&
    identical(unname(vapply(files, digest::digest, character(1),
                            file = TRUE, algo = "sha256")),
              hashes)
}


#' Download, check and unpack one archive
#'
#' @param which `"tools"` or `"genomes"`.
#' @param archive_dir A directory holding the archive file already, used
#'   instead of downloading it, or `NULL`.
#' @param archives The archive specifications (tests pass their own).
#' @param base_url Where the release files are served.
#' @return The unpacked directory, invisibly.
#' @noRd
.tantale_install_archive <- function(which, archive_dir = NULL,
                                     archives = .tantale_archives,
                                     base_url = "https://github.com/scunnac/tantale/releases/download") {
  spec <- archives[[which]]
  if (is.null(archive_dir)) {
    src <- tempfile(fileext = ".tar.gz")
    on.exit(unlink(src), add = TRUE)
    url <- paste(base_url, spec$version, spec$file, sep = "/")
    cli::cli_alert_info("Downloading {.url {url}}...")
    oldTimeout <- options(timeout = max(1200, getOption("timeout")))
    on.exit(options(oldTimeout), add = TRUE)
    ok <- tryCatch(utils::download.file(url, src, mode = "wb", quiet = TRUE) == 0L,
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) {
      cli::cli_abort(
        c("Could not download {.file {spec$file}} ({spec$what}).",
          "i" = "Check the network, or download it from {.url {url}} yourself and pass its directory as {.arg archive_dir} to {.fn tantale_setup}."),
        class = c("tantale_error_download", "tantale_error"))
    }
  } else {
    src <- file.path(archive_dir, spec$file)
    if (!file.exists(src)) {
      cli::cli_abort("{.file {spec$file}} is not in {.file {archive_dir}}.",
                     class = c("tantale_error_missing_file", "tantale_error"))
    }
  }
  if (!identical(digest::digest(src, file = TRUE, algo = "sha256"), spec$sha256)) {
    cli::cli_abort(
      c("{.file {spec$file}} does not have the expected sha256.",
        "i" = "The file is damaged or is not the one this version of tantale expects."),
      class = c("tantale_error_archive_checksum", "tantale_error"))
  }
  dest <- .tantale_archive_path(which, archives)
  unlink(dest, recursive = TRUE)
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  utils::untar(src, exdir = dirname(dest), tar = "internal")
  if (!.tantale_archive_ok(which, archives)) {
    cli::cli_abort("{.file {spec$file}} unpacked into {.file {dest}} with missing or altered files.",
                   class = c("tantale_error_archive_checksum", "tantale_error"))
  }
  invisible(dest)
}


#' The path of one of the downloaded Java tools
#'
#' The default of every wrapper's tool argument, so a user can still point a
#' wrapper to a copy of their own.
#' @param tool One of the names below.
#' @return A path that exists.
#' @noRd
.tantale_tool <- function(tool = c("annotale", "preditale", "talecorrection",
                                   "talecorrection_hmm")) {
  tool <- match.arg(tool)
  rel <- switch(tool,
                annotale = "AnnoTALEcli-1.5.jar",
                preditale = "PrediTALE.jar",
                talecorrection = file.path("talecorrect", "TALEcorrection.jar"),
                talecorrection_hmm = file.path("talecorrect", "HMMs", "Xoo"))
  path <- file.path(.tantale_archive_path("tools"), rel)
  if (!file.exists(path)) {
    cli::cli_abort(
      c("{.file {basename(rel)}} is not installed.",
        "i" = "Run {.run tantale::tantale_setup(install = TRUE)} to download the Java tools tantale wraps."),
      class = c("tantale_error_tool_missing", "tantale_error"))
  }
  path
}


#' Path to one of the example genomes
#'
#' Four *Xanthomonas oryzae* pv. *oryzae* genome assemblies are used
#' throughout the articles and examples. They are not part of the package
#' itself: [tantale_setup()] downloads them, with the Java tools, when called
#' with `install = TRUE`. This function returns where one of them is.
#'
#' MAI1 (GenBank CP025609.1), BAI3 (CP025610.1) and PXO86 (RefSeq
#' NZ_CP007166.1) are complete genomes of African and Asian strains;
#' the sequences are those of the records, with shortened FASTA headers.
#' BAI3-1-1 is an unpublished assembly of a BAI3 derivative that carries
#' sequencing and assembly errors in its TALE loci, distributed with tantale
#' to illustrate frameshift correction.
#'
#' @param strain Which genome.
#' @return The path of a FASTA file.
#' @seealso [tantale_setup()], which downloads them.
#' @export
#' @examples
#' \donttest{
#' # Needs the genomes downloaded by tantale_setup(install = TRUE).
#' tantale_genome("MAI1")
#' }
tantale_genome <- function(strain = c("MAI1", "BAI3", "BAI3-1-1", "PXO86")) {
  strain <- match.arg(strain)
  path <- file.path(.tantale_archive_path("genomes"), paste0(strain, ".fa"))
  if (!file.exists(path)) {
    cli::cli_abort(
      c("The example genome {.val {strain}} is not installed.",
        "i" = "Run {.run tantale::tantale_setup(install = TRUE)} to download the example genomes."),
      class = c("tantale_error_genome_missing", "tantale_error"))
  }
  path
}
