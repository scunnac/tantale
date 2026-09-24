# Records which object classes the exported functions take and return, for
# the data-flow view in dev/function-graph.qmd (ledger §29).
#
# Every exported function and S3 method is traced while the test suite
# runs. Each call logs, on entry, the class of its first argument (the
# first element of `...` when that comes first) and of any other supplied
# argument of a class listed in `tracked`; on exit, the class of its return
# value and of the elements of a returned list. The distinct combinations are written
# to dev/function-graph-dataflow.tsv. Run from the package root:
#
#   Rscript dev/function-graph-dataflow.R
#
# Takes as long as the full test suite (several minutes) and needs the same
# external tools. Rerun it when exported functions are added or change what
# they take or return; function-graph.qmd lists exports missing from the TSV.

pkgload::load_all(".", quiet = TRUE)
ns <- asNamespace("tantale")

exports <- getNamespaceExports(ns)
exports <- exports[vapply(exports, function(f) is.function(get(f, ns)), logical(1))]
s3 <- as.data.frame(getNamespaceInfo(ns, "S3methods")[, 1:3, drop = FALSE])
names(s3) <- c("generic", "class", "method")
traced <- sort(unique(c(exports, s3$method)))

describe <- function(x) {
  cl <- class(x)
  if (is.function(x)) return("function")
  if (inherits(x, "ggplot")) return("ggplot")
  cl[1]
}

.flow_log <- new.env()
.flow_log$rows <- list()

tracked <- c("tales", "tales_msa", "pairwise_distances", "tale_distances",
             "domain_distances", "BStringSet", "DNAStringSet", "AAStringSet")

# On entry, before the body can reassign an argument. Only supplied
# arguments are evaluated, so defaults are left to the body.
# R turns tracing off while a tracer runs. Evaluating an argument can run
# another traced function (tales_align(tales_compare_distal(x)$tales)), so
# tracing is turned back on for the duration, or that call goes unrecorded.
entry_args <- function(fn, env) {
  old <- tracingState(TRUE)
  on.exit(tracingState(old))
  args <- names(formals(get(fn, ns)))
  if (!length(args)) return(c(first = NA_character_, other = ""))
  what <- if (args[1] == "...") quote(..1) else as.name(args[1])
  first <- tryCatch(describe(eval(what, env)), error = function(e) NA_character_)
  other <- character()
  for (a in setdiff(args[-1], "...")) {
    supplied <- tryCatch(!eval(call("missing", as.name(a)), env), error = function(e) FALSE)
    if (!supplied) next
    cl <- tryCatch(describe(get(a, envir = env)), error = function(e) NA_character_)
    if (cl %in% tracked) other <- c(other, sprintf("%s=%s", a, cl))
  }
  c(first = first, other = paste(other, collapse = ";"))
}

record <- function(fn, entry, value) {
  parts <- ""
  if (is.list(value) && !is.data.frame(value) && !is.null(names(value))) {
    parts <- paste(sprintf("%s=%s", names(value), vapply(value, describe, "")),
                   collapse = ";")
  }
  key <- paste(fn, entry[["first"]], entry[["other"]], describe(value), parts,
               sep = "\t")
  .flow_log$rows[[key]] <- (if (is.null(.flow_log$rows[[key]])) 0L else .flow_log$rows[[key]]) + 1L
}

for (fn in traced) {
  suppressMessages(trace(
    fn, where = ns, print = FALSE,
    tracer = bquote(.tantale_flow_in <- .GlobalEnv$.tantale_flow_entry(.(fn), environment())),
    exit = bquote({
      # A function that returns NULL (talomes_heatmap() draws and returns
      # invisible(NULL)) is recorded too; only an exit by error is skipped.
      .v <- returnValue(default = quote(.tantale_no_value))
      if (!identical(.v, quote(.tantale_no_value)))
        .GlobalEnv$.tantale_flow_record(.(fn), .tantale_flow_in, .v)
    })
  ))
}
# A generic in another package (print(), dplyr_reconstruct()) can find a
# method through the S3 registry, which still holds the untraced function.
generic_home <- function(g) {
  if (exists(g, envir = ns, mode = "function")) return(ns)
  for (p in loadedNamespaces())
    if (exists(g, envir = asNamespace(p), mode = "function", inherits = FALSE))
      return(asNamespace(p))
  NULL
}
for (i in seq_len(nrow(s3))) {
  home <- generic_home(s3$generic[i])
  if (is.null(home)) next
  registerS3method(s3$generic[i], s3$class[i], get(s3$method[i], ns), envir = home)
}
assign(".tantale_flow_entry", entry_args, envir = .GlobalEnv)
assign(".tantale_flow_record", record, envir = .GlobalEnv)

testthat::test_dir(file.path("tests", "testthat"), package = "tantale",
                   load_package = "none", stop_on_failure = FALSE,
                   reporter = "summary")

rows <- names(.flow_log$rows)
out <- as.data.frame(do.call(rbind, strsplit(paste0(rows, "\t "), "\t")),
                     stringsAsFactors = FALSE)[, 1:5]
names(out) <- c("fn", "arg_class", "other_args", "return_class", "return_parts")
out$return_parts <- trimws(out$return_parts)
out$calls <- unlist(.flow_log$rows, use.names = FALSE)
out <- out[order(out$fn, -out$calls), ]

header <- sprintf("# Recorded %s at commit %s by dev/function-graph-dataflow.R",
                  format(Sys.Date()), system("git rev-parse --short HEAD", intern = TRUE))
path <- file.path("dev", "function-graph-dataflow.tsv")
writeLines(header, path)
suppressWarnings(write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE,
                             append = TRUE))
cat("Traced", length(traced), "functions;", length(unique(out$fn)), "seen;",
    nrow(out), "rows written to", path, "\n")
cat("Never called by the tests:", paste(setdiff(traced, out$fn), collapse = ", "), "\n")
