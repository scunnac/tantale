# The consumer requirements table must match reality, or it is worse than
# nothing: a stale table is a promise the code does not keep. This checks it by
# ablation -- drop each column and see whether the function still works.

fixture <- function() {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  tales(readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))$tale_parts)
}

works <- function(f, x) {
  isTRUE(tryCatch({suppressWarnings(suppressMessages(f(x))); TRUE}, error = function(e) FALSE))
}

callers <- list(
  tales_coded_strings    = function(x) tales_coded_strings(x),
  tales_domain_codes     = function(x) tales_domain_codes(x),
  tales_rvd_strings      = function(x) tales_rvd_strings(x),
  plot.tales = function(x) plot(x),
  tales_compare          = function(x) tales_compare(x)
)

test_that("every all_of column really is required", {
  x <- fixture()
  for (fn in names(callers)) {
    req <- tantale:::TALES_REQUIREMENTS[[fn]]
    for (col in intersect(req$all_of, names(x))) {
      expect_false(works(callers[[fn]], x[, setdiff(names(x), col)]),
                   label = paste0(fn, " is documented as needing ", col,
                                  " but works without it"))
    }
  }
})

test_that("every optional column really is optional", {
  x <- fixture()
  for (fn in names(callers)) {
    req <- tantale:::TALES_REQUIREMENTS[[fn]]
    for (col in intersect(req$optional, names(x))) {
      expect_true(works(callers[[fn]], x[, setdiff(names(x), col)]),
                  label = paste0(fn, " is documented as treating ", col,
                                 " as optional but fails without it"))
    }
  }
})

test_that("an any_of group needs at least one member, and any one suffices", {
  x <- fixture()
  for (fn in names(callers)) {
    req <- tantale:::TALES_REQUIREMENTS[[fn]]
    grp <- intersect(req$any_of, names(x))
    if (length(grp) < 2L) next
    # dropping all of them must fail
    expect_false(works(callers[[fn]], x[, setdiff(names(x), grp)]),
                 label = paste0(fn, " works with none of ", paste(grp, collapse = "/")))
    # keeping any single one must suffice
    for (keep in grp) {
      expect_true(works(callers[[fn]], x[, setdiff(names(x), setdiff(grp, keep))]),
                  label = paste0(fn, " fails with only ", keep))
    }
  }
})

test_that("tales_requirements() renders the table", {
  tbl <- tales_requirements()
  expect_s3_class(tbl, "data.frame")
  expect_named(tbl, c("fn", "requirement", "columns"))
  expect_setequal(unique(tbl$fn), names(tantale:::TALES_REQUIREMENTS))
})
