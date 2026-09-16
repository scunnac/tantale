# print() has to say what the class knows and a tibble does not. These pin the
# facts, not the layout -- the arrangement is free to change.

fixture <- function() {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  tales(readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))$tale_parts)
}

shown <- function(x, ...) paste(utils::capture.output(print(x, ...)), collapse = "\n")

test_that("print.tales() names the class, the arrays and the parts", {
  x <- fixture()
  out <- shown(x)
  expect_match(out, "tales")
  # arrays and parts are different numbers and both belong there: a tibble
  # only ever showed the second
  expect_match(out, paste0(length(unique(x$array_id)), " arrays"))
  expect_match(out, paste0(nrow(x), " parts"))
})

test_that("print.tales() reports which residue layers are present", {
  x <- fixture()
  expect_match(shown(x), "layers:.*rvd")
  expect_match(shown(x), "layers:.*dom_code")
  # and only the ones that are there
  bare <- as_tales(as.list(c(t1 = "NI-HD-NG")), sep = "-", residue_col = "rvd")
  expect_match(shown(bare), "layers: rvd")
  expect_no_match(shown(bare), "dom_code")
})

test_that("print.tales() shows the namespace only when the object is stamped", {
  x <- fixture()
  expect_no_match(shown(x), "namespace:")
  stamped <- tales(readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))$tale_parts,
                   dom_code_namespace = "abc123")
  expect_match(shown(stamped), "namespace: abc123")
})

test_that("print.tales() elides the middle, and does not when everything fits", {
  x <- fixture()
  ids <- unique(x$array_id)
  expect_match(shown(x, n = 2), "\\.\\.\\.")
  # first and last arrays are shown, the ones between are not
  expect_match(shown(x, n = 2), ids[1], fixed = TRUE)
  expect_match(shown(x, n = 2), ids[length(ids)], fixed = TRUE)
  expect_no_match(shown(x, n = 2), ids[10], fixed = TRUE)
  # four arrays with n = 2 is exactly the boundary: nothing to elide
  small <- x[x$array_id %in% ids[1:4], ]
  expect_no_match(shown(small, n = 2), "^\\s+\\.\\.\\.\\s*$")
})

test_that("print.tales() previews dom_code when both layers exist", {
  x <- fixture()
  first <- unique(x$array_id)[1]
  codes <- x$dom_code[x$array_id == first][order(x$position_in_array[x$array_id == first])]
  expect_match(shown(x), paste(utils::head(codes, 4), collapse = " "), fixed = TRUE)
})

test_that("print.tales() survives an object with no rows", {
  x <- fixture()
  expect_no_error(shown(x[0, ]))
  expect_match(shown(x[0, ]), "0 arrays")
})

test_that("print.tales_msa() reports the alignment width and draws gaps", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleTalesMsa.rds")))
  msa <- readRDS(test_path("data_for_tests", "sampleTalesMsa.rds"))
  out <- shown(msa)
  expect_match(out, "tales_msa")
  expect_match(out, paste0(tales_width(msa), " alignment positions"))
})

test_that("print.tales_msa() aligns its columns so gaps line up", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  d <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales(d$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:4], ]
  al <- suppressWarnings(suppressMessages(tales_align(sub, residue_col = "dom_code")))
  out <- shown(al, gap = "-")
  expect_match(out, "-")
  # every previewed row is the same printed length, which is what makes an
  # alignment readable down the page
  rows <- grep("ROI_", strsplit(out, "\n")[[1]], value = TRUE)
  expect_gt(length(rows), 1L)
  expect_length(unique(nchar(rows)), 1L)
})

test_that("print() returns its object invisibly", {
  x <- fixture()
  expect_invisible(print(x))
  expect_identical(utils::capture.output(y <- print(x)), utils::capture.output(print(x)))
  expect_identical(y, x)
})
