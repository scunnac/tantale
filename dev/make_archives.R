# Builds the two archives tantale_setup() downloads (ledger §50):
#
#   tantale-tools-1.tar.gz    the Java tools tantale wraps, with their
#                             licence, provenance and checksums
#   tantale-genomes-1.tar.gz  the four example genomes of the articles
#
# Each unpacks to one folder (tools-1/, genomes-1/) holding README.md and
# MANIFEST (sha256 of every other file). The archives are attached to the
# GitHub releases `tools-1` and `genomes-1` of scunnac/tantale, and
# R/tantale_data.R pins their sha256. A new version of any file means a new
# archive version (tools-2, ...) and a new pin; a published archive never
# changes.
#
# Run from the package root:  Rscript dev/make_archives.R <source> <out>
# <source> holds the files listed below (the package's inst/tools/ and
# inst/extdata/ before they left the package, commit 20755db); <out> receives
# the archives. Prints each archive's sha256.

args <- commandArgs(trailingOnly = TRUE)
src <- args[1]
out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

sha256 <- function(f) digest::digest(file = f, algo = "sha256")

write_manifest <- function(dir) {
  files <- sort(list.files(dir, recursive = TRUE))
  files <- setdiff(files, "MANIFEST")
  writeLines(paste(vapply(file.path(dir, files), sha256, character(1)), files,
                   sep = "  "),
             file.path(dir, "MANIFEST"))
}

pack <- function(name, version, files, readme) {
  stage <- file.path(tempfile(), version)
  for (f in names(files)) {
    dest <- file.path(stage, f)
    dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
    stopifnot(file.copy(file.path(src, files[[f]]), dest))
  }
  writeLines(readme, file.path(stage, "README.md"))
  write_manifest(stage)
  archive <- normalizePath(file.path(out, paste0("tantale-", name, "-1.tar.gz")),
                           mustWork = FALSE)
  # relative paths, from inside the staging folder's parent, so the archive
  # unpacks to <version>/
  old <- setwd(dirname(stage))
  on.exit(setwd(old))
  utils::tar(archive, files = version, compression = "gzip", tar = "internal")
  cat(basename(archive), sha256(archive), "\n")
}

talecorrect <- file.path("talecorrect", list.files(file.path(src, "talecorrect"),
                                                  recursive = TRUE))
tool_files <- c(`AnnoTALEcli-1.5.jar` = "AnnoTALEcli-1.5.jar",
                `PrediTALE.jar` = "PrediTALE.jar",
                `LICENSES/COPYING.GPL-3` = "COPYING.GPL-3",
                stats::setNames(talecorrect, talecorrect))

pack("tools", "tools-1", as.list(tool_files), c(
  "# tantale tools, version 1",
  "",
  "The Java programs that the R package tantale",
  "(https://github.com/scunnac/tantale) wraps. `tantale::tantale_setup()`",
  "downloads this archive, checks its sha256 and unpacks it; `MANIFEST`",
  "gives the sha256 of every file. The jars and the HMMs are byte-identical",
  "to their upstream copies (checked 2026-09-23).",
  "",
  "| file | program | upstream |",
  "|---|---|---|",
  "| `AnnoTALEcli-1.5.jar` | AnnoTALE 1.5 | https://www.jstacs.de/downloads/AnnoTALEcli-1.5.jar |",
  "| `PrediTALE.jar` | PrediTALE | https://www.jstacs.de/downloads/PrediTALE.jar |",
  "| `talecorrect/` | TALEcorrection (jar, HMMs, sources) | `TALECorrection_scripts.zip` at https://www.jstacs.de/downloads/ |",
  "",
  "All three are part of Jstacs (J. Grau, J. Keilwagen and co-authors) and",
  "are distributed under the GNU General Public License, version 3 or",
  "(at your option) any later version",
  "(`LICENSES/COPYING.GPL-3`). Their source code is at",
  "https://github.com/Jstacs/Jstacs (`projects/xanthogenomes/`,",
  "`projects/tals/prediction/`, `projects/talecorrect/`).",
  "",
  "References:",
  "",
  "- AnnoTALE: Grau J. et al. (2016). Scientific Reports 6, 21077. https://doi.org/10.1038/srep21077",
  "- PrediTALE: Erkes A. et al. (2019). PLoS Computational Biology 15, e1007206. https://doi.org/10.1371/journal.pcbi.1007206",
  "- TALEcorrection: Erkes A. et al. (2023). BMC Genomics 24. https://doi.org/10.1186/s12864-023-09228-1"
))

pack("genomes", "genomes-1",
     list(`MAI1.fa` = "MAI1.fa", `BAI3.fa` = "BAI3.fa",
          `BAI3-1-1.fa` = "BAI3-1-1.fa", `PXO86.fa` = "PXO86.fa"), c(
  "# tantale example genomes, version 1",
  "",
  "Four *Xanthomonas oryzae* pv. *oryzae* genome assemblies used by the",
  "articles and examples of the R package tantale",
  "(https://github.com/scunnac/tantale). `tantale::tantale_setup()`",
  "downloads this archive; `tantale::tantale_genome()` returns the path of",
  "one genome. `MANIFEST` gives the sha256 of every file.",
  "",
  "| file | strain | source |",
  "|---|---|---|",
  "| `MAI1.fa` | MAI1 | GenBank CP025609.1 (RefSeq NZ_CP025609.1) |",
  "| `BAI3.fa` | BAI3 | GenBank CP025610.1 (RefSeq NZ_CP025610.1) |",
  "| `PXO86.fa` | PXO86 | RefSeq NZ_CP007166.1 |",
  "| `BAI3-1-1.fa` | BAI3-1-1 | unpublished assembly, distributed with tantale |",
  "",
  "The sequences of MAI1, BAI3 and PXO86 are identical to the cited records",
  "(checked 2026-10-04); only the FASTA headers are shortened. BAI3-1-1 is",
  "an assembly of a BAI3 derivative that carries sequencing and assembly",
  "errors in its TALE loci, used to illustrate frameshift correction."
))
