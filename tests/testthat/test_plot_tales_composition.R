# plot() had no test and was broken: it calls mutate() and
# ggplot() unqualified, and neither dplyr nor ggplot2 was imported into the
# package namespace, so it only worked if the *user* had attached them.
# R CMD check's "no visible global function definition" NOTE was pointing at
# real breakage, not the usual NSE false positive.

test_that("plot() on a tales returns a ggplot", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  p <- plot(tales_quietly(out$tale_parts))
  expect_s3_class(p, "ggplot")
})

test_that("plot() works with only the package attached", {
  # The regression this guards: unqualified calls resolving through the user's
  # search path rather than the package namespace.
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  expect_no_error(plot(tales_quietly(out$tale_parts)))
})

test_that("seqnames is optional: the facet is added only when present", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales_quietly(out$tale_parts)
  withFacet <- plot(x)
  withoutFacet <- plot(x[, setdiff(names(x), "seqnames")])
  expect_s3_class(withFacet, "ggplot")
  expect_s3_class(withoutFacet, "ggplot")
  expect_s3_class(withFacet$facet, "FacetGrid")
  expect_s3_class(withoutFacet$facet, "FacetNull")
})

test_that("it reports the columns it needs rather than failing obscurely", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales_quietly(out$tale_parts)
  expect_error(plot(x[, setdiff(names(x), "aa_seq")]),
               class = "tantale_error_projection_column")
})

test_that("a legacy tale_parts data frame is still accepted", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  expect_s3_class(suppressWarnings(plot(tales_quietly(out$tale_parts))), "ggplot")
})

test_that("plot() dispatches to the composition plot for a tales", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales_quietly(out$tale_parts)
  p <- plot(x)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "Overview of TALE composition by genome")
})

test_that("a tales_msa still dispatches to plot.tales_msa, not plot.tales", {
  # tales_msa inherits from tales, so the more specific method must win. If
  # plot.tales ever shadowed it, an alignment would silently render as a
  # composition plot.
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales_quietly(out$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:3], ]
  msa <- suppressWarnings(suppressMessages(tales_align(sub, residue_col = "rvd")))
  p <- suppressWarnings(suppressMessages(plot(msa)))
  expect_false(identical(p$labels$title, "Overview of TALE composition by genome"))
})

test_that("position = 'alignment' lays parts out on the alignment coordinate", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales_quietly(out$tale_parts)
  sub <- x[x$array_id %in% unique(x$array_id)[1:6], ]
  msa <- suppressWarnings(suppressMessages(tales_align(sub, residue_col = "rvd")))
  back <- suppressWarnings(suppressMessages(as_tales(msa)))

  arrayLayout <- plot(back, position = "array")
  alignLayout <- plot(back, position = "alignment")
  # the aligned layout spans the alignment width; the array one only the longest array
  expect_equal(max(alignLayout$data$.x), tales_msa_width(msa))
  expect_equal(max(arrayLayout$data$.x), max(back$position_in_array))
  expect_gt(max(alignLayout$data$.x), max(arrayLayout$data$.x))
})

test_that("position = 'alignment' needs an aligned object", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales_quietly(out$tale_parts)
  expect_error(plot(x, position = "alignment"),
               class = "tantale_error_projection_column")
})

test_that("the default layout is unchanged", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales_quietly(out$tale_parts)
  expect_identical(plot(x)$data$.x, x$position_in_array)
})

test_that("facet_by chooses the panel columns", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales_quietly(out$tale_parts)
  x$strain <- sub("_ROI_.*", "", x$array_id)
  p <- plot(x, facet_by = "strain")
  expect_identical(names(p$facet$params$rows), "strain")
  p <- plot(x, facet_by = c("strain", "seqnames"))
  expect_identical(names(p$facet$params$rows), c("strain", "seqnames"))
  expect_s3_class(plot(x, facet_by = NULL)$facet, "FacetNull")
})

test_that("facet_by refuses an absent column and one that varies within an array", {
  skip_if_not(file.exists(test_path("data_for_tests", "sampleDistalrOutput.rds")))
  out <- readRDS(test_path("data_for_tests", "sampleDistalrOutput.rds"))
  x <- tales_quietly(out$tale_parts)
  expect_error(plot(x, facet_by = "strain"), class = "tantale_error_plot_facet")
  expect_error(plot(x, facet_by = "domain_type"), class = "tantale_error_plot_facet")
  expect_error(plot(x, facet_by = 1), class = "tantale_error_plot_facet")
  # an explicit request for an absent seqnames is an error, the default is not
  y <- x[, setdiff(names(x), "seqnames")]
  expect_error(plot(y, facet_by = "seqnames"), class = "tantale_error_plot_facet")
  expect_s3_class(plot(y)$facet, "FacetNull")
})


#### Colours (R/palette.R, ledger §46) ####

colour_fixture <- function() {
  aa <- function(n) strrep("A", n)
  tales(tibble::tibble(
    array_id = rep(c("a1", "a2"), each = 5),
    position_in_array = rep(1:5, 2),
    domain_type = rep(c("N-terminus", "repeat", "repeat", "repeat", "C-terminus"), 2),
    rvd = rep(c("NTERM", "HD", "NI", "NG", "CTERM"), 2),
    aa_seq = c(aa(288), aa(34), aa(33), aa(20), aa(278),
               aa(264), aa(34), aa(34), aa(20), aa(278))
  ))
}

test_that("black or white text is chosen by the fill's luminance", {
  expect_identical(.text_colour_on(c("#FFFFFF", "#000000", "#DDCC77", "#332288")),
                   c("black", "white", "black", "white"))
})

test_that("plot.tales() colours parts by role and length", {
  p <- plot(colour_fixture(), facet_by = NULL)
  fills <- stats::setNames(ggplot2::layer_data(p, 1)$fill,
                           paste(p$data$domain_type, nchar(p$data$aa_seq)))
  expect_identical(unname(fills["repeat 34"]), .tol_muted[["sand"]])
  expect_identical(unname(fills["repeat 20"]), .tol_muted[["cyan"]])
  expect_identical(unname(fills["repeat 33"]), .tol_muted[["rose"]])
  expect_identical(unname(fills["C-terminus 278"]), .tol_muted[["teal"]])
  # the longer N-terminus takes the full wine, the shorter a lighter shade
  expect_identical(unname(fills["N-terminus 288"]), .tol_muted[["wine"]])
  expect_false(fills[["N-terminus 264"]] %in% c(.tol_muted[["wine"]], fills[["N-terminus 288"]]))
  expect_identical(levels(p$data$part)[1:2], c("N-terminus, 264 aa", "N-terminus, 288 aa"))
})

test_that("plot.tales() writes dark text on light fills", {
  p <- plot(colour_fixture(), facet_by = NULL)
  d <- p$data
  expect_true(all(d$label_colour[d$colour == .tol_muted[["sand"]]] == "black"))
})
