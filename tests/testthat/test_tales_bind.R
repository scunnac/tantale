# Tests for tales_bind(). See ledger §5.2, "Implementation plan, as of
# 2026-09-19", for the design these pin.

# A minimal N-term/repeat/C-term array, aa_seq given explicitly so namespace
# behaviour is controllable.
bind_df <- function(array_id, aa_seq) {
  n <- length(aa_seq)
  tibble::tibble(
    array_id = array_id,
    position_in_array = seq_len(n),
    domain_type = c("N-terminus", rep("repeat", n - 2L), "C-terminus"),
    rvd = c("NTERM", paste0("R", seq_len(n - 2L)), "CTERM"),
    aa_seq = aa_seq
  )
}


#### type checks ####

test_that("every argument must be a tales object", {
  a <- tales(bind_df("A1", c("M", "AAA", "P")))
  expect_error(tales_bind(a, "nope"), class = "tantale_error_bind_type")
  expect_error(tales_bind(), class = "tantale_error_bind_type")
})

test_that("a tales_msa input is rejected, pointing at as_tales()", {
  a <- tales(bind_df("A1", c("M", "AAA", "P")))
  aligned <- data.frame(
    array_id = "A1", position_in_array = 1:3,
    alignment_position = 1:3, rvd = c("NTERM", "R1", "CTERM")
  )
  msa <- tales_msa(aligned)
  expect_error(tales_bind(a, msa), class = "tantale_error_bind_tales_msa")
  expect_error(tales_bind(a, msa), regexp = "as_tales")
})


#### array_id disjointness ####

test_that("colliding array_id across inputs is an error naming the ids", {
  a <- tales(bind_df("A1", c("M", "AAA", "P")))
  b <- tales(bind_df("A1", c("M", "BBB", "P")))
  err <- expect_error(tales_bind(a, b), class = "tantale_error_bind_array_id_collision")
  expect_match(conditionMessage(err), "A1")
})

test_that("disjoint array_id across inputs binds cleanly", {
  a <- tales(bind_df("A1", c("M", "AAA", "P")))
  b <- tales(bind_df("A2", c("M", "BBB", "P")))
  out <- tales_bind(a, b)
  expect_s3_class(out, "tales")
  expect_setequal(out$array_id, c("A1", "A2"))
  expect_equal(nrow(out), nrow(a) + nrow(b))
})


#### group column ####

test_that("group is dropped, with a message, when present on either input", {
  a <- tales(bind_df("A1", c("M", "AAA", "P")))
  a$group <- 1L
  b <- tales(bind_df("A2", c("M", "BBB", "P")))
  # disjoint labels -- still dropped, not just the colliding ones
  expect_message(out <- tales_bind(a, b), class = "tantale_message_bind_group_dropped")
  expect_false("group" %in% names(out))
})

test_that("group is dropped the same way when both inputs carry colliding labels", {
  a <- tales(bind_df("A1", c("M", "AAA", "P")))
  a$group <- 1L
  b <- tales(bind_df("A2", c("M", "BBB", "P")))
  b$group <- 1L
  expect_message(out <- tales_bind(a, b), class = "tantale_message_bind_group_dropped")
  expect_false("group" %in% names(out))
})

test_that("no message when neither input carries group", {
  a <- tales(bind_df("A1", c("M", "AAA", "P")))
  b <- tales(bind_df("A2", c("M", "BBB", "P")))
  expect_no_message(tales_bind(a, b))
})


#### dom_code namespace ####

test_that("a shared namespace binds without recoding", {
  # Coding both from one combined object, then splitting, gives them the same
  # namespace -- carried through subsetting, not recomputed (tales_class.R's
  # .tales_regrade()).
  both <- tales(dplyr::bind_rows(
    bind_df("A1", c("M", "AAA", "P")),
    bind_df("A2", c("M", "BBB", "P"))
  ))
  coded <- tales_assign_domain_codes(both)
  a <- coded[coded$array_id == "A1", ]
  b <- coded[coded$array_id == "A2", ]
  expect_identical(tales_namespace(a), tales_namespace(b))

  expect_no_message(out <- tales_bind(a, b))
  expect_identical(tales_namespace(out), tales_namespace(a))
  # codes unchanged from the shared coding -- a and b are disjoint row
  # subsets of coded, in the same relative order, so binding them back
  # reproduces it exactly.
  expect_equal(as.data.frame(out), as.data.frame(coded))
})

test_that("mismatched namespaces recode dom_code by default, with a message", {
  a <- tales_assign_domain_codes(tales(bind_df("A1", c("M", "AAA", "P"))))
  b <- tales_assign_domain_codes(tales(bind_df("A2", c("M", "BBB", "P"))))
  expect_false(identical(tales_namespace(a), tales_namespace(b)))

  expect_message(out <- tales_bind(a, b), class = "tantale_message_bind_namespace_recoded")
  expect_false(identical(tales_namespace(out), tales_namespace(a)))
  expect_false(identical(tales_namespace(out), tales_namespace(b)))
  # a fresh, one-to-one code <-> sequence correspondence over the union
  expect_equal(dplyr::n_distinct(out$dom_code), dplyr::n_distinct(out$aa_seq))
})

test_that("on_namespace_mismatch = 'error' aborts instead of recoding", {
  a <- tales_assign_domain_codes(tales(bind_df("A1", c("M", "AAA", "P"))))
  b <- tales_assign_domain_codes(tales(bind_df("A2", c("M", "BBB", "P"))))
  expect_error(
    tales_bind(a, b, on_namespace_mismatch = "error"),
    class = "tantale_error_bind_namespace_mismatch"
  )
})

test_that("two all-NULL namespaces agree and never trigger a recode", {
  a <- tales(bind_df("A1", c("M", "AAA", "P")))
  b <- tales(bind_df("A2", c("M", "BBB", "P")))
  expect_null(tales_namespace(a))
  expect_null(tales_namespace(b))
  expect_no_message(out <- tales_bind(a, b))
  expect_null(tales_namespace(out))
})


#### distances untouched ####

test_that("tales_bind() never touches tale_distances/domain_distances", {
  # No distance-table argument to feed or reconcile -- a caller who needs one
  # re-runs tales_compare_distal() on the bound result (documented, not enforced).
  expect_false("tal_sim" %in% names(formals(tales_bind)))
})


#### round trip against the manual workaround ####

test_that("tales_bind() matches the manual paste0(strain, '_', array_id) workaround", {
  genomes <- list(
    g1 = tales(bind_df("A1", c("M", "AAA", "P"))),
    g2 = tales(bind_df("A1", c("M", "CCC", "P"))), # same array_id, different genome
    g3 = tales(bind_df("A2", c("M", "DDD", "P")))
  )
  manual <- dplyr::bind_rows(
    lapply(names(genomes), function(nm) {
      x <- tibble::as_tibble(genomes[[nm]])
      x$array_id <- paste0(nm, "_", x$array_id)
      x
    })
  )
  prefixed <- Map(function(nm, x) {
    x$array_id <- paste0(nm, "_", x$array_id)
    tales(tibble::as_tibble(x))
  }, names(genomes), genomes)

  out <- do.call(tales_bind, unname(prefixed))

  common <- setdiff(names(manual), "dom_code")
  expect_equal(
    as.data.frame(out[order(out$array_id), common]),
    as.data.frame(manual[order(manual$array_id), common]),
    ignore_attr = TRUE
  )
})
