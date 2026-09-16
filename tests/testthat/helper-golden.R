# Fingerprinting for the regression baseline in test_golden.R.
#
# Snapshotting whole tables would produce files too large to read and too noisy
# to diff. A fingerprint keeps one row per column of each artefact, carrying
# what a refactor must not change -- type, length, distinct count, missingness,
# and a digest of the values. When something moves, the snapshot diff names the
# artefact and the column, which is the part that saves time; the value itself
# is then reproduced in a couple of lines at the console.
#
# Small artefacts whose contents are worth reading are snapshotted whole
# instead, by test_golden.R directly.

.fingerprint_column <- function(x) {
  # Digest the values only. Names live in the `names` column of the fingerprint
  # so a rename shows up as its own change rather than as a value change.
  v <- unname(x)
  if (is.factor(v)) v <- as.character(v)
  # Doubles are rounded before digesting: an alignment or a distance recomputed
  # on another machine can differ in the last bits without differing in any way
  # that matters here.
  if (is.double(v)) v <- round(v, 8)
  list(
    type       = typeof(v),
    n          = length(v),
    n_distinct = length(unique(v)),
    n_missing  = sum(is.na(v)),
    digest     = digest::digest(v, algo = "md5")
  )
}

fingerprint <- function(x) {
  if (inherits(x, "XStringSet")) {
    x <- data.frame(name = names(x), seq = as.character(unname(x)),
                    stringsAsFactors = FALSE)
  }
  if (is.matrix(x)) {
    x <- data.frame(value = as.vector(x), stringsAsFactors = FALSE)
  }
  if (!is.data.frame(x)) {
    x <- data.frame(value = x, stringsAsFactors = FALSE)
  }
  x <- as.data.frame(x)
  rows <- lapply(names(x), function(nm) {
    c(list(column = nm), .fingerprint_column(x[[nm]]))
  })
  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}

# expect_snapshot_value() skips on CRAN by default. This baseline exists to be
# run, and a check that silently does not run is worse than no check, so it is
# forced on. Everything it needs (MAFFT, arlem) ships in inst/tools.
expect_golden <- function(x) {
  testthat::expect_snapshot_value(x, style = "json2", cran = TRUE)
}


#### tell_tales() ####

# tell_tales() returns invisible(output_dir): its real output is the ~36 files
# it writes, so that directory is what has to be pinned.
#
# Some of them record the run rather than the result -- AnnoTALE's
# protocol_analyze.txt and HMMER's echoed command line embed temp paths, the
# log and both GFFs stamp the date. Rather than keep a list of which files are
# affected, every file is digested with those lines removed. Keeping a list
# was in fact wrong: the GFFs were not on it, and the omission only showed up
# when a session ran past midnight.

# Lines that say where or when the run happened, rather than what it found.
.RUN_SPECIFIC <- paste(
  "^##date",            # rtracklayer's GFF header
  "^# Date:",           # HMMER's own header
  "^# Version:",        # HMMER's build, not its findings: a version that
                        # changed results would show up in the hits instead
  "Current date",       # tell_tales.log
  "^# CPU time:",       # HMMER, genuinely varies run to run
  "^# Mc/sec:",         # HMMER throughput, likewise
  "/tmp/",              # any absolute temp path
  "Rtmp",               # R's per-session temp directory
  "file[0-9a-f]{10,}",  # tempfile() basenames
  sep = "|"
)

.telltale_file_digest <- function(path) {
  txt <- readLines(path, warn = FALSE)
  keep <- txt[!grepl(.RUN_SPECIFIC, txt)]
  list(n_lines = length(txt),
       n_dropped = length(txt) - length(keep),
       digest = digest::digest(keep, algo = "md5"))
}

telltale_fingerprint <- function(dir) {
  files <- sort(list.files(dir, recursive = TRUE))
  rows <- lapply(files, function(f) {
    d <- .telltale_file_digest(file.path(dir, f))
    data.frame(file = f, n_lines = d$n_lines, n_dropped = d$n_dropped,
               digest = d$digest, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# tell_tales() takes ~9 seconds, and three separate expectations want to look
# at the same run. Run it once per test file and hand the same output
# directory to all of them.
telltale_run <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      out <- file.path(tempdir(), "golden_telltale")
      unlink(out, recursive = TRUE)
      ret <- suppressWarnings(suppressMessages(tantale::tell_tales(
        subject_file = system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                                   package = "tantale", mustWork = TRUE),
        output_dir = out)))
      cache <<- list(dir = out, returned = ret)
    }
    cache
  }
})
