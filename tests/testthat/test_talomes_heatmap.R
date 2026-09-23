# talomes_heatmap() draws with base graphics, so these tests look at what
# reaches the device: the strings in its display list, the file written to
# save_path, and the devices left open afterwards.

ann <- data.frame(
  group = c("G1", "G1", "G1", "G2", "G2", "G3"),
  strain = c("S1", "S2", "S3", "S1", "S2", "S3"),
  rvdseq = c("NI-HD-NG", "NI-HD-NG", "NN-HD-NG",
             "HD-NI-NG-NG", "HD-NI-NG-NG", "NG-NG-HD"),
  trunc = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
  origin = c("Mali", "Mali", "Burkina", "Mali", "Mali", "Burkina")
)

# The character arguments of every drawing operation on the current device.
drawn_strings <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grDevices::dev.control("enable")
  force(expr)
  ops <- grDevices::recordPlot()[[1]]
  unlist(lapply(ops, function(op) Filter(is.character, as.list(op[[2]]))))
}

test_that("both plot types draw and return NULL invisibly", {
  for (type in c("all", "single")) {
    grDevices::pdf(NULL)
    expect_invisible(res <- talomes_heatmap(ann, "group", "strain", "rvdseq",
                                            plot_type = type))
    grDevices::dev.off()
    expect_null(res)
  }
})

test_that("title and axis labels are the ones supplied, in both plot types", {
  for (type in c("all", "single")) {
    txt <- drawn_strings(talomes_heatmap(ann, "group", "strain", "rvdseq",
                                         plot_type = type, title = "My title",
                                         x_lab = "My groups", y_lab = "My strains"))
    expect_true(all(c("My title", "My groups", "My strains") %in% txt), info = type)
    expect_false("RVD sequences variants" %in% txt, info = type)
  }
})

test_that("column labels count the distinct RVD sequences per group", {
  txt <- drawn_strings(talomes_heatmap(ann, "group", "strain", "rvdseq"))
  # G1 has two distinct sequences, G2 one, G3 one
  expect_true(all(c("G1 #2", "G2 #1", "G3 #1") %in% txt))
})

test_that("a numeric group is shown with a G prefix", {
  num <- transform(ann, group = as.integer(sub("G", "", group)))
  txt <- drawn_strings(talomes_heatmap(num, "group", "strain", "rvdseq"))
  expect_true("G1 #2" %in% txt)
})

test_that("truncated TALEs are marked and the extra column gets a legend", {
  txt <- drawn_strings(talomes_heatmap(ann, "group", "strain", "rvdseq",
                                       trunc_tales_col = "trunc",
                                       extra_col = "origin"))
  expect_true("T" %in% txt)
  expect_true(all(c("Burkina", "Mali") %in% txt))
})

test_that("save_path writes the file and closes its device", {
  for (type in c("all", "single")) {
    f <- tempfile(fileext = ".png")
    devs <- grDevices::dev.list()
    talomes_heatmap(ann, "group", "strain", "rvdseq", plot_type = type,
                    save_path = f)
    expect_true(file.exists(f) && file.size(f) > 0, info = type)
    expect_identical(grDevices::dev.list(), devs, info = type)
  }
})

test_that("an explicit save_path = NULL draws on the current device", {
  txt <- drawn_strings(talomes_heatmap(ann, "group", "strain", "rvdseq",
                                       save_path = NULL))
  expect_true("G1 #2" %in% txt)
})

test_that("an unknown plot_type is an error", {
  expect_error(talomes_heatmap(ann, "group", "strain", "rvdseq",
                               plot_type = "heatmap"),
               "should be one of")
})
