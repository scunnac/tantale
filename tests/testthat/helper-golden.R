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
# forced on. Everything it needs is bundled or in the conda environment.
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
#
# The temp-directory rule cannot be a generic "looks temp-y" pattern -- a
# bare "/tmp/" and a bare "Rtmp" were both tried, in that order, and both
# over-matched the same way: neither can tell "this run's own scratch output
# directory" from "some unrelated path that merely happens to live under a
# tempdir too", which is exactly what happens to the four reference-file
# lines below whenever the *package itself* is installed under a tempdir --
# not this test session's tempdir, some other process's, at a different
# time -- as covr::package_coverage()'s own install step does by default
# (confirmed: it reintroduced the identical failure the "/tmp/" fix was
# for, byte-identical digest and all, until this was root-caused properly).
#
# Matching this session's *actual, current* `tempdir()` value instead is
# precise by construction: it can only ever match this specific run's own
# output directory, whatever it happens to be named or wherever `tempdir()`
# decided to put it -- and cannot match some other process's tempdir just
# because both happen to start with "Rtmp", the same way two people sharing
# a surname are not the same person.
.run_specific_pattern <- function() {
  paste(
    "^##date",            # rtracklayer's GFF header
    "^# Date:",           # HMMER's own header
    "^# Version:",        # HMMER's build, not its findings: a version that
                          # changed results would show up in the hits instead
    "Current date",       # tell_tales.log
    "^# Current dir:",    # HMMER echoes the working directory, not a finding
    "^# CPU time:",       # HMMER, genuinely varies run to run
    "^# Mc/sec:",         # HMMER throughput, likewise
    "file[0-9a-f]{10,}",  # tempfile() basenames -- a property of the name
                          # itself (a real reference file is never bare hex),
                          # not of where it sits, so this one is safe as-is
    gsub("([.\\+*?\\[\\]^$(){}|\\\\])", "\\\\\\1", tempdir(), perl = TRUE),
    sep = "|"
  )
}

# Absolute directories say where this machine keeps a file, not what the run
# found. tell_tales.log echoes four of them -- the three HMM profiles and
# correction_ref -- and they differ between a load_all() run from the source
# tree and a run from an installed package, so digesting them verbatim made
# this baseline reproduce only on the machine that recorded it (ledger 8.1d).
#
# Rewriting rather than dropping the line is what keeps the useful half. The
# basename says *which* reference was used, and that is exactly the signal
# that caught the 8.1b rename; only the directory is noise.
# The pattern has to be narrow. A first attempt at "anything between two
# slashes" also ate `</head>`, HMMER's `//` record separators and the `//` in
# `http://hmmer.org/` -- all of which are content, and silently digesting
# them away would have made the baseline weaker, not more portable.
#
# So: a slash that starts a token (preceded by nothing or by a delimiter),
# followed by one or more non-empty `segment/` groups. `//` cannot match
# because a segment must be non-empty, and `</head>` cannot because `<` is
# not a delimiter. Only the directory part is consumed; the basename stays.
.PATH_PREFIX <- "(^|[[:space:]=,;'\"(\\[])/(?:[A-Za-z0-9._+~-]+/)+"

.normalise_paths <- function(txt) {
  gsub(.PATH_PREFIX, "\\1<path>/", txt, perl = TRUE)
}

# The log names the package's own files (the subject file of the golden run,
# the HMM profiles, correction_ref) by absolute path. Where tantale is
# installed is not a finding, and pkgcheck's R CMD check installs it under the
# test session's own tempdir(), where the tempdir rule above dropped those five
# lines and changed the digest. The installation directory is therefore
# rewritten first, to a path that .normalise_paths() reduces to "<path>/"
# exactly as it reduces the real one.
.telltale_file_digest <- function(path, pkg_dir = system.file(package = "tantale")) {
  txt <- readLines(path, warn = FALSE)
  if (nzchar(pkg_dir)) txt <- gsub(pkg_dir, "/pkg", txt, fixed = TRUE)
  keep <- txt[!grepl(.run_specific_pattern(), txt, perl = TRUE)]
  list(n_lines = length(txt),
       n_dropped = length(txt) - length(keep),
       digest = digest::digest(.normalise_paths(keep), algo = "md5"))
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
