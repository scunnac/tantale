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
