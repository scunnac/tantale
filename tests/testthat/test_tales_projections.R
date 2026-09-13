# Tests for the 'tales' projections. See dev/restructuring-notes.md §1.

coded_tales <- function() {
  tales(tibble::tibble(
    array_id = c(rep("a1", 3), rep("a2", 2)),
    position_in_array = c(1:3, 1:2),
    dom_code = c("7", "2", "7", "2", "13"),
    aa_seq = c("MDP", "LTPA", "MDP", "LTPA", "SIVA"),
    rvd = c("NTERM", "NI", "NTERM", "NI", "CTERM")
  ))
}

test_that("tales_coded_strings() renders one string per array, in part order", {
  out <- tales_coded_strings(coded_tales())
  expect_s4_class(out, "BStringSet")
  expect_setequal(names(out), c("a1", "a2"))
  expect_identical(as.character(out[["a1"]]), "7 2 7")
  expect_identical(as.character(out[["a2"]]), "2 13")
})

test_that("tales_coded_strings() orders by position_in_array, not row order", {
  x <- coded_tales()
  shuffled <- x[rev(seq_len(nrow(x))), ]
  expect_identical(as.character(tales_coded_strings(shuffled)),
                   as.character(tales_coded_strings(x)))
})

test_that("tales_domain_codes() is one row per distinct code", {
  out <- tales_domain_codes(coded_tales())
  expect_setequal(names(out), c("dom_code", "aa_seq", "rvd"))
  expect_identical(nrow(out), 3L)               # codes 7, 2, 13
  expect_identical(out$aa_seq[out$dom_code == "7"], "MDP")
})

test_that("both projections require dom_code", {
  x <- tales(tibble::tibble(
    array_id = "a", position_in_array = 1L, rvd = "NI"
  ))
  expect_error(tales_coded_strings(x), class = "tantale_error_projection_column")
  expect_error(tales_domain_codes(x), class = "tantale_error_projection_column")
})

test_that("an NA dom_code is refused rather than pasted as the text 'NA'", {
  # the failure mode restructuring-notes.md §1 warns about
  x <- coded_tales()
  x$dom_code[1] <- NA_character_
  expect_error(tales_coded_strings(x), class = "tantale_error_projection_na")
})

test_that("tales_domain_codes() needs aa_seq and rvd", {
  x <- coded_tales()
  expect_error(tales_domain_codes(x[setdiff(names(x), "aa_seq")]),
               class = "tantale_error_projection_column")
})


#### Fidelity to what distalr() used to store ####

test_that("the projections reproduce the stored slots of a real distalr run", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(out$tale_parts)

  # coded.repeats.str: same names, same sequences
  rebuilt <- tales_coded_strings(x)
  stored <- out$coded.repeats.str
  expect_setequal(names(rebuilt), names(stored))
  expect_identical(as.character(rebuilt[names(stored)]), as.character(stored))

  # repeats.code: same code -> sequence/rvd mapping
  codes <- tales_domain_codes(x)
  expect_identical(nrow(codes), nrow(out$repeats.code))
  joined <- merge(
    transform(as.data.frame(codes), code = as.integer(dom_code)),
    as.data.frame(out$repeats.code), by = "code"
  )
  expect_identical(nrow(joined), nrow(codes))
  expect_true(all(joined$aa_seq == joined$`AA Seq`))
  expect_true(all(joined$rvd.x == joined$rvd.y))
})
