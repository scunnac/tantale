# Records which functions take or return a matrix, a plain list of atomic
# vectors, or a plain list of matrices, for the "Matrices and lists of
# vectors" section of dev/function-graph.qmd (ledger §29).
#
# Every function in the namespace, internals included, is traced while the
# test suite runs. Each call logs, on entry, the shape of every supplied
# argument (each element of `...` too) and, on exit, the shape of the
# return value, or of each element of a returned named, unclassed list
# that is not itself one of the shapes. Functions the tests never call are
# written as rows with direction "never called". Run from the package root:
#
#   Rscript dev/function-graph-shapes.R
#
# Takes as long as the full test suite (several minutes) and needs the same
# external tools. An optional argument is passed to testthat's `filter`, to
# try the script on a few test files: Rscript dev/function-graph-shapes.R msa

filter <- commandArgs(TRUE)[1]
if (is.na(filter)) filter <- NULL

pkgload::load_all(".", quiet = TRUE)
ns <- asNamespace("tantale")

fns <- Filter(function(f) is.function(get(f, ns)), ls(ns, all.names = TRUE))
fns <- fns[!startsWith(fns, ".__") & !fns %in% c(".onLoad", ".onAttach")]

# Only unclassed lists count as lists here: a classed list (a ggplot, a
# summary object) is an object with its own interface.
shape_of <- function(x) {
  if (is.data.frame(x) || is.function(x) || is.environment(x)) return(NA_character_)
  if (is.matrix(x)) return(sprintf("matrix<%s>", typeof(x)))
  if (identical(class(x), "list") && length(x) > 0) {
    types <- paste(sort(unique(vapply(x, typeof, ""))), collapse = "/")
    if (all(vapply(x, function(e) is.atomic(e) && !is.null(e) && is.null(dim(e)), logical(1)))) {
      # A named record of single values (file paths, settings) is told
      # apart from a list holding one vector per array or sequence.
      if (all(lengths(x) <= 1)) return(sprintf("list of scalars<%s>", types))
      return(sprintf("list of vectors<%s>", types))
    }
    if (all(vapply(x, is.matrix, logical(1))))
      return(sprintf("list of matrices<%s>", types))
  }
  NA_character_
}
safe_shape <- function(x) tryCatch(shape_of(x), error = function(e) NA_character_)

log_env <- new.env()
log_env$rows <- list()
called <- new.env()
add <- function(fn, direction, arg, shape) {
  if (is.na(shape)) return(invisible())
  key <- paste(fn, direction, arg, shape, sep = "\t")
  log_env$rows[[key]] <- (if (is.null(log_env$rows[[key]])) 0L else log_env$rows[[key]]) + 1L
}

# On entry, before the body can reassign an argument. Only supplied
# arguments are evaluated, so defaults are left to the body. No function
# in the package captures its arguments unevaluated (substitute(),
# enquo(), match.call()), so forcing them early changes nothing.
# R turns tracing off while a tracer runs. Forcing an argument can run
# package code (validate_tales(new_tales(x)) runs new_tales() here), so
# tracing is turned back on for the duration, or such calls go unrecorded.
on_entry <- function(fn, env) {
  called[[fn]] <- TRUE
  old <- tracingState(TRUE)
  on.exit(tracingState(old))
  for (a in names(formals(get(fn, ns)))) {
    if (a == "...") {
      dots <- tryCatch(eval(quote(list(...)), env), error = function(e) list())
      for (d in dots) add(fn, "in", "...", safe_shape(d))
      next
    }
    supplied <- tryCatch(!eval(call("missing", as.name(a)), env), error = function(e) FALSE)
    if (supplied) add(fn, "in", a, tryCatch(safe_shape(get(a, envir = env)),
                                            error = function(e) NA_character_))
  }
}
on_exit <- function(fn, value) {
  shape <- safe_shape(value)
  add(fn, "out", "<value>", shape)
  if (is.na(shape) && identical(class(value), "list") && !is.null(names(value))) {
    for (n in names(value)) add(fn, "out", paste0("$", n), safe_shape(value[[n]]))
  }
}
assign(".shape_entry", on_entry, envir = .GlobalEnv)
assign(".shape_exit", on_exit, envir = .GlobalEnv)

# One trace() per function: a second call would replace the first.
for (fn in fns) {
  suppressMessages(trace(
    fn, where = ns, print = FALSE,
    tracer = bquote(.GlobalEnv$.shape_entry(.(fn), environment())),
    exit = bquote({
      .v <- tryCatch(returnValue(), error = function(e) NULL)
      if (!is.null(.v)) .GlobalEnv$.shape_exit(.(fn), .v)
    })
  ))
}

# A generic in another package (dplyr_reconstruct(), print()) finds a
# method through the S3 registry, which still holds the untraced function.
s3 <- getNamespaceInfo(ns, "S3methods")
# The generic is looked up from tantale's namespace, or from whichever
# loaded namespace defines it (dplyr_col_modify() is not imported).
generic_home <- function(g) {
  if (exists(g, envir = ns, mode = "function")) return(ns)
  for (p in loadedNamespaces())
    if (exists(g, envir = asNamespace(p), mode = "function", inherits = FALSE))
      return(asNamespace(p))
  NULL
}
for (i in seq_len(nrow(s3))) {
  home <- generic_home(s3[i, 1])
  if (!s3[i, 3] %in% fns || is.null(home)) next
  registerS3method(s3[i, 1], s3[i, 2], get(s3[i, 3], ns), envir = home)
}

testthat::test_dir(file.path("tests", "testthat"), package = "tantale",
                   load_package = "none", stop_on_failure = FALSE,
                   reporter = "summary", filter = filter)

rows <- names(log_env$rows)
out <- as.data.frame(do.call(rbind, strsplit(rows, "\t")), stringsAsFactors = FALSE)
names(out) <- c("fn", "direction", "arg", "shape")
out$calls <- unlist(log_env$rows, use.names = FALSE)
never <- setdiff(fns, ls(called, all.names = TRUE))
out <- rbind(out, data.frame(fn = never, direction = rep("never called", length(never)),
                             arg = "", shape = "", calls = 0L))
out <- out[order(out$fn, out$direction, out$arg), ]

header <- sprintf("# Recorded %s at commit %s by dev/function-graph-shapes.R%s",
                  format(Sys.Date()), system("git rev-parse --short HEAD", intern = TRUE),
                  if (is.null(filter)) "" else sprintf(" (test filter \"%s\")", filter))
path <- file.path("dev", "function-graph-shapes.tsv")
writeLines(header, path)
suppressWarnings(write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE,
                             append = TRUE))
cat("Traced", length(fns), "functions;", length(ls(called, all.names = TRUE)), "called;",
    length(unique(out$fn[out$direction != "never called"])),
    "with a matrix or list of vectors; written to", path, "\n")
