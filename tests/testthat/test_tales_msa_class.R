# Tests for the 'tales_msa' class. See dev/class-design.md §4.

# Two arrays of 3 parts each, aligned over 4 columns. Array a2 has a gap at
# alignment position 2 -- represented by the *absence* of a row.
minimal_msa_df <- function() {
  tibble::tibble(
    array_id = c(rep("a1", 3), rep("a2", 3)),
    position_in_array = c(1:3, 1:3),
    alignment_position = c(1L, 2L, 3L, 1L, 3L, 4L),
    rvd = c("NTERM", "NI", "HD", "NTERM", "HD", "CTERM")
  )
}

test_that("tales_msa() builds a subclass of tales", {
  x <- tales_msa(minimal_msa_df())
  expect_s3_class(x, "tales_msa")
  expect_s3_class(x, "tales")
  expect_true(is_tales(x))
  expect_true(is_tales_msa(x))
})

test_that("every tales invariant still applies", {
  df <- minimal_msa_df()
  df$position_in_array[2] <- 1L            # duplicate key
  expect_error(tales_msa(df), class = "tantale_error_tales_duplicate_key")
})

test_that("alignment_width defaults to the largest alignment_position", {
  x <- tales_msa(minimal_msa_df())
  expect_identical(tales_width(x), 4L)
})

test_that("an explicit alignment_width is kept", {
  x <- tales_msa(minimal_msa_df(), alignment_width = 10L)
  expect_identical(tales_width(x), 10L)
})


#### Invariants specific to an alignment ####

test_that("alignment_position is required", {
  df <- minimal_msa_df()
  df$alignment_position <- NULL
  expect_error(tales_msa(df), class = "tantale_error_tales_missing_column")
})

test_that("alignment_position must be positive and non-NA", {
  df <- minimal_msa_df()
  df$alignment_position[1] <- NA_integer_
  expect_error(tales_msa(df), class = "tantale_error_msa_position")
})

test_that("alignment_position must be unique within an array", {
  df <- minimal_msa_df()
  df$alignment_position <- c(1L, 1L, 3L, 1L, 3L, 4L)
  expect_error(tales_msa(df), class = "tantale_error_msa_duplicate")
})

test_that("an alignment may not reorder parts", {
  df <- minimal_msa_df()
  df$alignment_position[1:3] <- c(3L, 2L, 1L)   # reversed against position_in_array
  expect_error(tales_msa(df), class = "tantale_error_msa_order")
})

test_that("alignment_position may not exceed the declared width", {
  df <- minimal_msa_df()
  expect_error(tales_msa(df, alignment_width = 3L),
               class = "tantale_error_msa_width")
})

test_that("gaps are absent rows, so arrays may have different row counts", {
  df <- minimal_msa_df()
  df <- df[-2, ]                                # drop a part of a1
  expect_s3_class(tales_msa(df), "tales_msa")
})


#### as.matrix() ####

test_that("as.matrix() materialises the grid with NA gaps", {
  x <- tales_msa(minimal_msa_df())
  m <- as.matrix(x)
  expect_equal(dim(m), c(2L, 4L))
  expect_setequal(rownames(m), c("a1", "a2"))
  expect_identical(m["a1", ], c("1" = "NTERM", "2" = "NI", "3" = "HD", "4" = NA))
  expect_identical(m["a2", "2"], NA_character_)   # the implicit gap
})

test_that("as.matrix() honours the gap argument and the chosen layer", {
  df <- minimal_msa_df()
  df$dom_code <- c("1", "2", "3", "1", "3", "4")
  x <- tales_msa(df)
  expect_identical(as.matrix(x, gap = "-")["a2", "2"], "-")
  expect_identical(as.matrix(x, value = "dom_code")["a1", "2"], "2")
})

test_that("as.matrix() errors on a layer the object does not have", {
  x <- tales_msa(minimal_msa_df())
  expect_error(as.matrix(x, value = "dom_code"),
               class = "tantale_error_msa_layer")
})


#### Graded degradation ####

test_that("dropping alignment_position demotes to tales, not to a tibble", {
  x <- tales_msa(minimal_msa_df())
  out <- dplyr::select(x, -alignment_position)
  expect_false(is_tales_msa(out))
  expect_true(is_tales(out))
  expect_null(tales_width(out))
})

test_that("dropping a tales key column degrades all the way to a tibble", {
  x <- tales_msa(minimal_msa_df())
  out <- dplyr::select(x, -array_id)
  expect_false(is_tales(out))
  expect_s3_class(out, "tbl_df")
})

test_that("row subsetting keeps the class and the declared width", {
  x <- tales_msa(minimal_msa_df())
  out <- dplyr::filter(x, array_id == "a1")
  expect_s3_class(out, "tales_msa")
  expect_identical(tales_width(out), 4L)   # width is carried, not recomputed
})


#### tales_align() ####

test_that("tales_align() returns a tales_msa that preserves every layer", {
  x <- suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
  msa <- suppressWarnings(tales_align(x, residue_col = "rvd"))

  expect_s3_class(msa, "tales_msa")
  # all input columns survive, plus the new coordinate
  expect_true(all(names(x) %in% names(msa)))
  expect_true("alignment_position" %in% names(msa))
  # no parts gained or lost
  expect_identical(nrow(msa), nrow(x))
})

test_that("tales_align() output is internally consistent", {
  x <- suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
  msa <- suppressWarnings(tales_align(x, residue_col = "rvd"))

  expect_silent(validate_tales_msa(msa))
  # the matrix view is as wide as the declared alignment
  m <- as.matrix(msa)
  expect_identical(ncol(m), tales_width(msa))
  expect_setequal(rownames(m), unique(x$array_id))
  # every non-gap cell of the matrix corresponds to exactly one part
  expect_identical(sum(!is.na(m)), nrow(msa))
})

test_that("tales_align() round-trips the residues it aligned on", {
  x <- suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
  msa <- suppressWarnings(tales_align(x, residue_col = "rvd"))
  m <- as.matrix(msa, value = "rvd")

  # reading a row left to right, ignoring gaps, gives back that array's RVDs
  # in position_in_array order -- the property the positional back-map relies on
  for (a in rownames(m)) {
    original <- x$rvd[x$array_id == a][order(x$position_in_array[x$array_id == a])]
    expect_identical(unname(m[a, !is.na(m[a, ])]), original)
  }
})

test_that("tales_align() refuses an incomplete tales", {
  x <- suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
  repeats_only <- dplyr::filter(x, domain_type == "repeat")
  expect_error(tales_align(repeats_only),
               class = "tantale_error_tales_incomplete")
})

test_that("tales_align() refuses a layer the object does not have", {
  x <- suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
  expect_error(tales_align(x, residue_col = "dom_code"),
               class = "tantale_error_msa_layer")
})

test_that("tales_align() carries the dom_code namespace", {
  x <- suppressWarnings(
    tales_from_telltale(test_path("data_for_tests", "tellTaleExampleOutput"))
  )
  attr(x, "dom_code_namespace") <- "run-xyz"
  msa <- suppressWarnings(tales_align(x, residue_col = "rvd"))
  expect_identical(tales_namespace(msa), "run-xyz")
})


#### plot() method ####

test_that("plot() on a tales_msa produces the same object as the direct call", {
  # plot_tales_msa() itself has no assertions anywhere in the suite, so this
  # also serves as the first executable check that the plotting path runs.
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(out$tale_parts)
  msa <- suppressWarnings(suppressMessages(
    tales_align(dplyr::filter(x, array_id %in% unique(x$array_id)[1:4]),
                residue_col = "rvd")
  ))

  viaMethod <- suppressWarnings(suppressMessages(
    plot(msa, fill = "rvd", label = NULL)
  ))
  direct <- suppressWarnings(suppressMessages(
    plot_tales_msa(repeat_align = as.matrix(msa, value = "rvd"))
  ))
  expect_s3_class(viaMethod, class(direct)[1])
})

test_that("plot() errors on a layer the alignment lacks", {
  x <- tales_msa(minimal_msa_df())
  expect_error(plot(x, fill = "dom_code"), class = "tantale_error_msa_layer")
})

test_that("plot() accepts similarity tables in either vocabulary", {
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(out$tale_parts)
  keep <- unique(x$array_id)[1:4]
  msa <- suppressWarnings(suppressMessages(
    tales_align(dplyr::filter(x, array_id %in% keep), residue_col = "rvd")
  ))
  legacy <- out$tal.similarity[out$tal.similarity$TAL1 %in% keep &
                                 out$tal.similarity$TAL2 %in% keep, ]

  fromLegacy <- suppressWarnings(suppressMessages(
    plot(msa, fill = "rvd", label = NULL, tal_sim = legacy)
  ))
  fromTyped <- suppressWarnings(suppressMessages(
    plot(msa, fill = "rvd", label = NULL, tal_sim = tale_sim(legacy))
  ))
  expect_s3_class(fromLegacy, class(fromTyped)[1])
})
