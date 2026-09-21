# summary() reports what print() is too cheap to compute. These pin the
# numbers and the vocabulary, not the layout.

fx <- function() {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  tales(readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))$tale_parts)
}
shown <- function(x) paste(utils::capture.output(print(x)), collapse = "\n")

test_that("summary() returns an object whose numbers can be used, not just read", {
  s <- summary(fx())
  expect_s3_class(s, "summary.tales")
  expect_type(s, "list")
  expect_identical(s$n_arrays, 44L)
  expect_identical(s$n_parts, 955L)
})

test_that("distinct domains counts domains, not repeats", {
  # A dom_code names any distinct part sequence, and the termini are parts.
  # Calling the total a repeat count overstates it by the number of termini.
  x <- fx()
  s <- summary(x)
  expect_identical(s$n_distinct_domains, length(unique(x$dom_code)))
  expect_gt(s$n_distinct_domains, s$n_distinct_by_type[["repeat"]])
  expect_identical(sum(s$n_distinct_by_type), s$n_distinct_domains)
  # and the printed line says "domains"
  expect_match(shown(s), "distinct domains")
  expect_no_match(shown(s), "distinct repeats")
})

test_that("distinct RVDs counts repeats only, and ignores missing values", {
  x <- fx()
  before <- summary(x)$n_distinct_rvds
  gapped <- x
  gapped$rvd[which(gapped$domain_type == "repeat")[1:3]] <- NA
  # NA is a missing RVD, not another kind of RVD
  expect_identical(summary(suppressWarnings(tales(gapped)))$n_distinct_rvds, before)
})

test_that("repeats per array excludes the termini", {
  x <- fx()
  s <- summary(x)
  perArray <- tapply(x$domain_type == "repeat", x$array_id, sum)
  expect_equal(unname(s$repeats_per_array[["max"]]), unname(max(perArray)))
  # strictly fewer than the part count, since every array here has two termini
  expect_lt(s$repeats_per_array[["max"]], max(table(x$array_id)))
})

test_that("summary() surfaces anomalies that print() does not", {
  x <- fx()
  expect_identical(nrow(summary(x)$anomalies), 0L)
  expect_match(shown(summary(x)), "anomalies\\s+none")

  bad <- x
  bad$rvd[c(3, 40)] <- NA
  s <- summary(suppressWarnings(tales(bad)))
  expect_gt(nrow(s$anomalies), 0L)
  expect_match(shown(s), "missing_rvd")
  expect_match(shown(s), "tales_anomalies")
})

test_that("summary.tales_msa() reports gaps and per-layer consensus", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:4], ]
  al <- suppressWarnings(suppressMessages(tales_align(sub, residue_col = "dom_code")))
  s <- summary(al)

  expect_s3_class(s, "summary.tales_msa")
  expect_identical(s$width, tales_width(al))
  expect_identical(s$n_gaps, sum(is.na(as.matrix(al, value = "dom_code"))))
  expect_named(s$layers, c("dom_code", "rvd"))

  # repeats that differ can share an RVD, so the coarser layer agrees at least
  # as often as the finer one
  expect_gte(s$layers$dom_code$n_no_consensus, s$layers$rvd$n_no_consensus)
})

test_that("both summaries print and return invisibly", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleTalesMsa.rds")))
  for (s in list(summary(fx()),
                 summary(readRDS(test_path("data_for_tests", "sampleTalesMsa.rds"))))) {
    expect_invisible(print(s))
    expect_gt(length(utils::capture.output(print(s))), 3L)
  }
})
