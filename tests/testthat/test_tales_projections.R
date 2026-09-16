# Tests for the 'tales' projections. See dev/restructuring-notes.md §1.

# a1 is a complete array (both termini), a2 an incomplete one. Code 2 recurs
# within a1, which is what the rendering tests below rely on.
#
# This fixture used to carry rvd = NTERM, NI, NTERM for a1 -- two N-termini,
# the second of them mid-array. That is not a TALE, and it went unnoticed
# because without a domain_type column the anomaly checks cannot run. Adding
# the column made tales() flag it immediately (terminus_duplicated,
# terminus_misplaced), which is the class working as intended.
coded_tales <- function() {
  tales(tibble::tibble(
    array_id = c(rep("a1", 4), rep("a2", 2)),
    position_in_array = c(1:4, 1:2),
    dom_code = c("7", "2", "2", "13", "2", "13"),
    aa_seq = c("MDP", "LTPA", "LTPA", "SIVA", "LTPA", "SIVA"),
    rvd = c("NTERM", "NI", "NI", "CTERM", "NI", "CTERM"),
    domain_type = c("N-terminus", "repeat", "repeat", "C-terminus",
                    "repeat", "C-terminus")
  ))
}

test_that("tales_coded_strings() renders one string per array, in part order", {
  out <- tales_coded_strings(coded_tales())
  expect_s4_class(out, "BStringSet")
  expect_setequal(names(out), c("a1", "a2"))
  expect_identical(as.character(out[["a1"]]), "7 2 2 13")
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


#### RVD strings ####

test_that("tales_rvd_strings() drops termini by default", {
  x <- tales(tibble::tibble(
    array_id = rep("a1", 4),
    position_in_array = 1:4,
    rvd = c("NTERM", "NI", "HD", "CTERM")
  ))
  expect_identical(as.character(tales_rvd_strings(x)[["a1"]]), "NI-HD")
  expect_identical(as.character(tales_rvd_strings(x, rvd_only = FALSE)[["a1"]]),
                   "NTERM-NI-HD-CTERM")
})

test_that("tales_rvd_strings() honours sep and part order", {
  x <- tales(tibble::tibble(
    array_id = rep("a1", 3),
    position_in_array = 3:1,
    rvd = c("NG", "HD", "NI")
  ))
  expect_identical(as.character(tales_rvd_strings(x, sep = " ")[["a1"]]), "NI HD NG")
})

test_that("tales_rvd_strings() matches the format of the shipped sample fasta", {
  shipped <- Biostrings::readBStringSet(
    system.file("extdata", "Sample_TALEs_RVDSeqs_AnnoTALE.fasta",
                package = "tantale", mustWork = TRUE)
  )
  x <- suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
  built <- tales_rvd_strings(x)
  # same shape: dash-separated RVDs, no anchor codes
  expect_false(any(grepl(paste(tales_anchor_codes(), collapse = "|"),
                         as.character(built))))
  expect_true(all(grepl("^[A-Z*]+(-[A-Z*]+)+$", as.character(built))))
  expect_true(all(grepl("^[A-Z*]+(-[A-Z*]+)+$", as.character(shipped))))
})


test_that("tales_domain_codes() includes rvd when present but does not require it", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  withRvd <- tales_domain_codes(x)
  expect_named(withRvd, c("dom_code", "aa_seq", "rvd"))
  # a dom_code + aa_seq object is valid and must not be blocked: the
  # correspondence between them is the substance, and rvd is decoration
  noRvd <- tales_domain_codes(x[, setdiff(names(x), "rvd")])
  expect_named(noRvd, c("dom_code", "aa_seq"))
  expect_identical(nrow(noRvd), nrow(withRvd))
})

test_that("tales_domain_codes() still requires aa_seq", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  expect_error(tales_domain_codes(x[, setdiff(names(x), "aa_seq")]),
               class = "tantale_error_projection_column")
})

#### sep and repeats_only (ledger 8.2b) ####

test_that("tales_coded_strings() honours sep", {
  x <- coded_tales()
  expect_identical(
    as.character(tales_coded_strings(x, sep = "-")),
    gsub(" ", "-", as.character(tales_coded_strings(x)), fixed = TRUE))
})

test_that("tales_coded_strings() defaults to a space, as ARLEM and MAFFT expect", {
  # Pinned deliberately: this default differs from tales_rvd_strings()'s "-",
  # and the difference is load-bearing rather than an oversight.
  expect_match(as.character(tales_coded_strings(coded_tales()))[1], " ")
})

test_that("tales_coded_strings(repeats_only = TRUE) drops the termini", {
  x <- coded_tales()
  full <- tales_coded_strings(x)
  bare <- tales_coded_strings(x, repeats_only = TRUE)
  nRepeat <- table(x$array_id[x$domain_type == "repeat"])
  expect_identical(unname(lengths(strsplit(as.character(bare), " ", fixed = TRUE))),
                   as.integer(nRepeat[names(bare)]))
  # and it really is a subset: fewer codes than the unfiltered rendering
  expect_true(all(nchar(bare) < nchar(full)))
})

test_that("tales_coded_strings(repeats_only = TRUE) needs domain_type", {
  x <- coded_tales()
  x$domain_type <- NULL
  expect_error(tales_coded_strings(x, repeats_only = TRUE),
               class = "tantale_error_projection_column")
})
