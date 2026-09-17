# Working on tantale

An R package for analysing TALEs (transcription activator-like effectors) of
*Xanthomonas*. Current work is pre-publication cleanup on the `dev` branch.

**`dev/restructuring-notes.md` is the ledger** and the main source of
context: what has been done, what is deferred, and why. Sections are marked
`**[A]**` (actionable), `**[V]**` (verified/done) or `**[ ]**` (open). Read
the relevant section before starting on something, and write the outcome
back into it — including what was tried and rejected.

## Standing rules

**Never delete code that looks dead.** Park it in
`R/unused_pending_review.R` (or `inst/legacy/` for whole retired classes)
with a comment saying what superseded it. This is the maintainer's explicit
instruction; "obsolete now, plausibly useful later" is a real category here.

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

## Commits

Commit when a piece of work is coherent and its tests pass. Say what moved
and why; reference the ledger section. Do not push without being asked.
