---
name: golden-rebaseline
description: Run the tantale golden regression baseline and, when it reports changes, decide whether to accept them. Use whenever test_golden.R fails, after any change that could alter what tell_tales() writes or what the tales/tales_msa classes produce, or when asked to re-baseline or accept snapshots.
---

# Re-baselining the golden tests

`tests/testthat/test_golden.R` records what the pipeline currently produces so
that a refactor meant to change nothing can be shown to have changed nothing.
It asserts nothing about what the values *should* be — which is exactly why
accepting a snapshot carelessly is how it stops being worth anything.

**The rule this skill exists to enforce: never accept before you have
explained every changed row.** An accepted snapshot is indistinguishable from
a correct one. Once wrong values are in `_snaps/golden.md`, every later run
agrees with them.

## Procedure

### 1. Run, and do not accept yet

```bash
Rscript -e 'suppressMessages(devtools::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test_golden.R")'
```

If it passes, stop. There is nothing to accept.

### 2. Identify what changed

The failure prints row indices, not file names. For a `telltale_fingerprint`
diff, rows are files in `sort(list.files(dir, recursive = TRUE))` order, so
resolve an index to a name:

```r
f <- sort(list.files(out, recursive = TRUE)); f[36]
```

For a `fingerprint()` diff, rows are columns of the artefact, and the
`column` field names them directly.

### 3. Explain every changed row before accepting

For each one, state what changed and why, and satisfy yourself it is
intended. Useful things to know:

- **`tell_tales.log` alone changed** — usually a parameter was added,
  renamed, or its default changed. Check `n_lines` in the diff: a line count
  that moved by exactly one matches one new parameter echoed into the log.
- **Many files changed at once** — suspect the fingerprint machinery rather
  than the pipeline. `.RUN_SPECIFIC` (dropped lines) and `.normalise_paths()`
  (rewritten paths) in `helper-golden.R` apply to every file, so a change
  there moves everything.
- **Nothing should have changed** — then something did, and that is the
  finding. Do not accept it.

When you touch `helper-golden.R`, audit what the change actually rewrites
before trusting it:

```r
txt  <- readLines(file.path(out, f), warn = FALSE)
keep <- txt[!grepl(.RUN_SPECIFIC, txt)]
cbind(keep, .normalise_paths(keep))[keep != .normalise_paths(keep), ]
```

This matters. A normaliser that quietly rewrites content produces a perfectly
stable digest while making the baseline *weaker* — a first version of
`.normalise_paths()` also ate `</title></head>`, HMMER's `//` record
separators and the `//` in `http://hmmer.org/`, and every test would still
have gone green.

### 4. Accept, then re-run to confirm

`snapshot_accept()` only accepts what is already on disk from a previous
failing run, so accept and re-run as two separate steps:

```bash
Rscript -e 'suppressMessages(devtools::load_all(".", quiet=TRUE)); testthat::snapshot_accept("golden", path="tests/testthat")'
Rscript -e 'suppressMessages(devtools::load_all(".", quiet=TRUE)); testthat::test_file("tests/testthat/test_golden.R")'
```

Running accept and the test in one command accepts the *previous* state and
then fails again on the new one.

### 5. Say what moved

In the commit message and in `dev/restructuring-notes.md`, record which
artefacts changed and why. "Re-baselined golden" is not a description.

## Notes

- The suite takes a few minutes; run it in the background and do other work.
- It needs MAFFT, HMMER and arlem. If they are missing these tests **fail
  rather than skip** — deliberately, since a baseline that quietly does not
  run is worse than none. Run `tantale_setup()` to check the environment.
- `tests/testthat/_snaps/golden.html` is a `snapshot_review()` artefact and
  is gitignored. Do not commit it.
