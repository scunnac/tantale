# Working on tantale

An R package for analysing TALEs (transcription activator-like effectors) of
*Xanthomonas*. Current work is pre-publication cleanup on the `dev` branch.

**`dev/restructuring-notes.md` is the ledger** and the main source of
context: what has been done, what is deferred, and why. It is ~3200 lines
and is a record, not a reading list — **start at its `START HERE` block**,
which lists the open items, then read only the sections bearing on the task
in hand.

Markers: `[V]` verified/done, `[A]` agreed but not executed, `[P]` parked
pending a judgement call, `[superseded]` a record of a replaced plan and
therefore *not* work.

Write outcomes back into the relevant section — including what was tried
and rejected.

## Standing rules

**Never delete code that looks dead.** Park it in
`R/unused_pending_review.R` (or `inst/legacy/` for whole retired classes)
with a comment saying what superseded it. This is the maintainer's explicit
instruction; "obsolete now, plausibly useful later" is a real category here.

**pkgdown articles are written as Quarto (`.qmd`), not R Markdown.**
Agreed direction for the future website. `vignettes/articles/*.qmd` --
pkgdown 2.2.0 supports quarto vignettes natively (see its NEWS). R's own
build machinery never descends into `vignettes/articles/`
(`tools::pkgVignettes()`), so these are never built by `R CMD build`/`check`
and need no `VignetteBuilder` entry. The numbered `vignettes/*.Rmd` files
are a separate, older thing -- §7.5, deliberately last -- not yet migrated.

**Vignettes come last, and never constrain the code.** They will be rebuilt
from the finished API, not the other way round. Do not let an existing
vignette dictate a signature, and do not spend effort keeping them building
mid-refactor.

**Conda: always `-p <prefix>`, never `-n <name>`.** conda and micromamba
keep separate roots, and `reticulate::conda_list()` can return two
environments both named `tantale`. Using the name let three rebuilds report
success while the package went on using the other copy. `tantale_setup()`
prints the binary, the root and the environment in use for exactly this
reason.

**Roxygen markdown is ON** (`Roxygen: list(markdown = TRUE)`). So
`[square brackets]` in prose become link targets — write `(parentheses)`
instead. Existing `\code{}`/`\strong{}` macros still work and need not be
converted.

**Run `pkgdown::check_pkgdown()` after adding or retiring an exported
topic.** A dangling entry in `_pkgdown.yml` is a hard error in
`build_site()`, and `R CMD check` does not look at that file.

**Use cli conditions, never `stop()`, `warning()` or `message()`.**
`cli::cli_abort()` / `cli_warn()` / `cli_inform()`, each with a condition
class (`c("tantale_error_<what>", "tantale_error")`) so tests can assert on
the class rather than the wording. Re-audit with R's parser, not grep, so
comments and strings cannot give false positives:

```r
pd <- getParseData(parse(f, keep.source = TRUE))
pd[pd$token == "SYMBOL_FUNCTION_CALL" &
   pd$text %in% c("stop", "warning", "message"), ]
```

Legitimate exceptions: `cat()` inside `print`/`format` methods,
`packageStartupMessage()` in `startup.R`, and `stopifnot()` for internal
invariants that are not addressed to the user. `classification.R` still has
nine unconverted sites, deliberately — see ledger §11.

**Run the environment's programs by absolute path, never via `PATH` or
`conda run`.** `.tantale_bin(tools)` resolves them inside the conda prefix
and `.tantale_exec()` runs them with the exit status checked. This is not
style: this machine carries `/usr/bin/mafft` 7.505 and `/usr/bin/nhmmer`
3.4 against pins of 7.453 and 3.3.2, and `conda run` puts only the *first*
command of a compound string inside the environment. See ledger §12.

**MAFFT is pinned to 7.453 on purpose.** Later versions changed `--text`
mode gap handling and align TALE repeat-code strings differently, leaving
the termini unanchored. Do not "update" it. Same care for HMMER (3.3.2).

## Tests

- **Run targeted test files**, not the whole suite, unless the change
  ripples broadly. The full suite takes several minutes.
- **Tests fail rather than skip when an external tool is missing.** A check
  that silently does not run is worse than none.
- **The golden baseline is the safety net for refactors** meant to change
  nothing. Use the `golden-rebaseline` skill when it reports a change — the
  discipline is to explain every changed row *before* accepting.
- **Beware vacuous assertions.** `x$missing_column` is `NULL`, and
  `expect_identical(NULL, NULL)` passes. This has bitten twice. When a test
  reads a column that may not exist, assert the column exists first.
- **A fixture that omits the columns the validators key on is not
  exercising the validators.** Adding `domain_type` to one minimal fixture
  immediately exposed an array with two N-termini that had sat there
  unnoticed.

## API conventions

rOpenSci guidelines (ledger §9.0): `object_verb()` naming, **data first**,
snake_case for arguments and columns, no clashes with base or tidyverse.

Two S3 class families: `tales` → `tales_msa`, and `pairwise_distances` →
`tale_distances` / `domain_distances`.

**`dom_code` names a distinct *domain* sequence, not a repeat.** The N- and
C-termini are parts like the repeats are and get codes too — that is why the
word is "domain". Codes come from `dplyr::cur_group_id()`, so they are
meaningful only within the call that minted them; the `dom_code_namespace`
stamp exists to catch cross-run mixing, and is enforced, not advisory.

## Documentation

- **Readers are biologists too.** Embed the TALE biology in the docs, not
  just the R mechanics.
- **No implementation archaeology in user docs.** What changed and why
  belongs in code comments or the ledger, not in `@details`.
- **Put user-facing prose in the exported function's block**, not in an
  internal's `@noRd` — easy to get wrong when the real work lives in a
  helper.
- `@param ...` must name the arguments it accepts if it forwards to
  something the reader cannot open (ledger §8.1c).

## Where things stand

`dev/restructuring-notes.md` is ~3200 lines. **Read its `START HERE` block,
not the whole file.** It lists the six open items and what each is blocked
on. Three of them are waiting on the maintainer's decision, not on work.

## Commits

Commit when a piece of work is coherent and its tests pass. Say what moved
and why; reference the ledger section. Do not push without being asked.
