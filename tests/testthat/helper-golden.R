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
# Seven of them carry the run rather than the result -- AnnoTALE's
# protocol_analyze.txt, HMMER's echoed command line, and the log all embed
# temp paths and a timestamp. Verified by running twice: the other 29 are
# byte-identical, and those seven differ in nothing else. They are compared by
# line count after the volatile lines are dropped, which still catches a stage
# that stops writing or starts writing more.

.TELLTALE_VOLATILE <- c("hmmerSearchOut.txt",
                        "nhmmerHumanReadableOutputOfLastRun.txt",
                        "tell_tales.log",
                        "protocol_analyze.txt")

.is_volatile <- function(path) {
  any(vapply(.TELLTALE_VOLATILE, function(v) endsWith(path, v), logical(1)))
}

# Content digest for a file, ignoring anything that encodes where or when the
# run happened.
.telltale_file_digest <- function(path) {
  if (.is_volatile(path)) {
    txt <- readLines(path, warn = FALSE)
    drop <- grepl("/tmp/|Rtmp|[0-9]{4}$|Current date|file[0-9a-f]{8,}", txt)
    return(list(kind = "volatile", n_lines_kept = sum(!drop)))
  }
  list(kind = "stable", digest = digest::digest(file = path, algo = "md5"))
}

telltale_fingerprint <- function(dir) {
  files <- sort(list.files(dir, recursive = TRUE))
  rows <- lapply(files, function(f) {
    d <- .telltale_file_digest(file.path(dir, f))
    data.frame(file = f, kind = d$kind,
               value = if (d$kind == "stable") d$digest else as.character(d$n_lines_kept),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
