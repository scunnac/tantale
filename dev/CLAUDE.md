# Working on tantale

An R package for analysing TALEs (transcription activator-like effectors) of
*Xanthomonas*. Current work is pre-publication cleanup on the `dev` branch.

New to TALE biology? [The Wikipedia
page](https://en.wikipedia.org/wiki/Transcription_activator-like_effector)
is worth a look before touching code that deals with repeats, RVDs or
target prediction.

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

**Never delete code that looks dead.** Park it in `R/unused_pending_review.R`
(create it fresh if it does not currently exist -- its last batch was
reviewed and moved out to `inst/legacy/` on 2026-09-21, see ledger §17) or
`inst/legacy/` directly for whole retired classes, with a comment saying
what superseded it. This is the maintainer's explicit instruction;
"obsolete now, plausibly useful later" is a real category here. Before
parking or retiring anything, check every similarly-named sibling for its
own, separate test dependency -- two near-identical functions in the same
file can each be kept alive by a different test file (exactly what
happened when this convention's own file was last reviewed). Grep is not
enough on its own: a retired function can have a *dedicated* test file
that exercises only it and nothing else (`test_conda.R` for
`.run_in_conda()`, missed on the same 2026-09-21 pass, caught only by a
full `devtools::test()` run afterwards) -- run the full suite before
calling a retirement done, don't rely on spot-checking the files you
already suspect.

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

**`devtools::load_all()` is not enough before a `pkgdown::build_site()`
or `docs/` rebuild after changing `R/`.** Quarto's article-rendering step
spawns a *separate* R subprocess that runs the article's own `library(tantale)`
call, resolving against whatever is *installed*, not the live source
tree `load_all()` patched into the calling session. A stale install
fails with a plain "could not find function" deep inside a `quarto
render` error that pkgdown reports only as a generic "Error running
quarto CLI" -- `quiet = FALSE` is needed to see the real cause. Run
`devtools::install(quick = TRUE, upgrade = FALSE)` first, and confirm in
a fresh `Rscript` session (`exists("fn", where = asNamespace("tantale"))`),
not the same session that ran the install. Cost two failed rebuilds to
find, 2026-09-20 -- see ledger's `START HERE` header.

**Four articles share one cached analysis; prime it before a full
build, every time, not just when something breaks.**
`tale_classification.qmd`, `tale_msa.qmd`, `tales_msa_class.qmd` and
`tale_target_prediction.qmd` all run the same three-genome
discovery/comparison/grouping pipeline. Only `tale_classification.qmd`
actually computes it (visibly, as the teaching example) and caches the
result to `vignettes/articles/_cache/` (gitignored); the other three
silently `readRDS()` it and abort with a clear message if it is missing.
`pkgdown::build_site()`/`build_articles()` renders every `.qmd` under
`vignettes/` in one `quarto render` project pass, and **that pass's
internal order is not alphabetical, not `_pkgdown.yml`'s navbar order,
and not controllable via a committed `vignettes/_quarto.yaml`
`project: render:` list -- all three tried and falsified empirically,
2026-09-20.** A clean-cache full build reliably renders
`tale_target_prediction.qmd` (position 2 of 7, observed) before
`tale_classification.qmd`, so it reliably aborts unless primed first.
**Don't rebuild articles via `pkgdown::build_articles()`/`build_site()`
at all -- call `pkgdown::build_article()` once per file instead.** A
priming call before the batch (render `tale_classification.qmd` alone,
then let the whole-project pass render everything, including it, again)
avoids the abort but still re-renders the canonical article's own
uncached content (the small-fixture backend comparison, the functal/
motif section) a second time for nothing. `build_article(name, pkg =
".")` renders exactly one file directly -- confirmed by reading its
`build_quarto_articles(pkg, article = ...)` branch -- so looping it over
all seven names (`articles/tale_classification` first, the rest and
`getting_started` in any order) plus one call to
`build_articles_index(pkg)` for the listing page builds every article
exactly once, in a chosen order, with no double rendering and no
dependence on `quarto`'s own (uncontrollable) project order. Verified
end-to-end: identical `docs/articles/*.html` output to a normal build,
10m20s total for a cold cache, zero errors, zero re-renders. See ledger
§15 for the full record, including the `_quarto.yaml` experiment so it
is not retried.

**`docs/` publishes in release mode directly -- `docs/dev/` no longer
exists (dropped 2026-09-20).** `_pkgdown.yml`'s `development: mode:` is
hardcoded to `release` (changed from `auto` on 2026-09-21, after testing
directly that it builds cleanly on this pkgdown/quarto version -- see
ledger §15) for exactly this reason: with `mode: auto`, the version's
third component being `>= 9000` would route any build into `docs/dev/`
unless every single `build_*` call passed an override. No override is
needed for `build_*()` calls any more; a bare call now writes straight
to `docs/` as intended.

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

**Branch is `main`, not `dev`.** A full repository history reset was
executed 2026-09-22 (ledger §26): `master`/`dev` and all pre-reset history
are retired, replaced by a single orphan commit on `main`. Old history
(327 commits, all branches/tags) survives as a `git bundle` attached to
the `v0.1.9553` GitHub release, and the pre-reset working directory was
kept on disk as `tantale-old-before-reset`, not deleted. Nothing about
any tracked file's *content* changed in the reset itself.

`dev/restructuring-notes.md` is ~7750 lines. **Read its `START HERE`
block for the pre-2026-09-21 history; for everything since, read §17-§27
directly** (numbered, in order, at the end of the file) -- this note is
the short pointer, not a re-summary of either. §27 (2026-09-22) closed
clean: all three test-suite findings from §26's post-reset `R CMD check`
(`reshape2` leftover dependency, a hardcoded `ncores`, and a golden-
baseline mismatch in the frameshift-correction test) are fixed and
verified; `test_golden.R` passes clean under both `load_all()` and a real
install. Two unrelated findings surfaced by that same check run and are
flagged, not yet investigated: `tales_group_kmedoids()`'s own `@examples`
fails under a real check (does not reproduce under `load_all()`), and
`inst/tools/arlem/arlem` triggers an "undeclared executable file"
`R CMD check` warning (plausibly an accepted cost of bundling a
third-party binary, not confirmed).

**§24 is a maintainer triage of the whole open-items list, 2026-09-21,
same day as §17-23 but a later session -- read it before assuming any
"still open" item below needs a decision.** It also prompted a ledger
cleanup pass (fixed ~10 stale status markers found by a full read-through,
compressed the "-original"/"-superseded" historical subsections) -- see
§24 for the priorities that came out of it and what got dropped from
active tracking. Not yet pushed as of this note; re-check `git log
--oneline origin/dev..dev` directly rather than trusting this note if any
doubt.

**Session (2026-09-21, earlier) ends at §21 mid-item, on purpose --
an intermediate pause, not a close, and is written up as one.**

**§17, the in-depth documentation review of every exported function and
S3 method -- DONE.** All 51 exports and 19 S3 methods (re-verified count,
not assumed -- see §17's own closing note for how the original "49" was
itself wrong) read against their actual current behaviour, not just
structurally checked. One real, live functional bug found and fixed with
the maintainer's sign-off: `correct_tales()` was feeding two nHMMER
result files to `TALEcorrection.jar`'s flags swapped, confirmed against
the tool's own printed usage -- corrected sequences may now differ from
before. Version bumped to `0.9.9003`.

**§18 correction:** the earlier retirement of `unused_pending_review.R`'s
last batch had missed an orphaned test file (`test_conda.R`, testing
`.run_in_conda()` directly, nothing else) -- only caught by running the
*full* suite once, at the very end of §17, not the targeted files the
retirement itself checked. Removed; lesson folded into this file's own
"never delete code that looks dead" rule above.

**§19, the `repeat_sims`/`tal_sim`/`domain_sim`/`fill_type` naming
question -- acted on, see §23.** §20 is still parked, explicitly not to
be started without going back to the maintainer first: `tale_parts_to_rvd()`
as a rename/refactor candidate.

**§21, `plot.tales_msa()`'s matrix round trips -- in progress, not
done.** The function round-trips its already-long `tales_msa` input
through an array-by-position matrix and back at least nine times before
handing long data to `ggplot()`, which wants it long anyway; only two of
those (the consensus/match computation) have been converted so far, via
two new private, `tales_msa`-native functions
(`.tales_consensus_long()`/`.tales_consensus_match_long()` in
`tales_consensus.R`) that are candidates to become public drop-in
replacements for `tales_consensus()`/`tales_consensus_match()` later --
documented to that standard already, `@noRd` for now. Internal "repeat_*"
identifiers throughout `tales_plot.R` renamed to "domain_*" while in
there; the public API (`fill_type`, `domain_sim` -> `domain_distances`,
`tal_sim` -> `tale_distances`) was untouched at the time but has since
been renamed too, see §23. **Read §21's own closing note before touching
this function again** -- it lists the three remaining matrix-shaped helpers in
the order they'd naturally get done, what falls out for free once they
are (`domain_align`/`rvd_align` disappearing from the function entirely),
and two separate, only-just-noticed findings not to conflate with the
main thread: a mislabelled `position_in_array` column (actually alignment
position), and `dev/class-design.md` §4.6 -- an older, stale design note
that already called for exactly this `tales_consensus_match()`-as-method
direction back on 2026-09-13 and was never revisited until this session
rediscovered the same idea independently.

**§21's own remaining items (the three matrix-shaped helpers,
`.pick_ref_name()` included) are the maintainer's to do, not an
assistant's -- do not start them without being asked, same footing as
§19/§20 above.** Two tests (`test_plot_tales_msa.R`,
`test_error_conditions.R`) call `.pick_ref_name()`/`.rvd_to_match_align()`
directly and assert their current matrix-in/matrix-out shape; converting
these helpers to the long-tibble form the rest of §21 already moved to
means rewriting those test assertions too, and the maintainer wants that
part done personally. See the ledger's §21 closing note, 2026-09-21, for
the full reasoning.

**§22, the `position_in_array` mislabel §21 flagged but left alone --
fixed.** `domainAlignLong`/`rvdAlignLong`'s per-position column (and
everything joined on it, and the base plot's `aes()`/`scale_x_discrete()`
`limits`) renamed to `alignment_position`, matching what
`as.matrix.tales_msa()`'s columns actually are. Label only, no values
changed -- golden baseline confirmed a single column-name diff and
nothing else. The visible x-axis title, still "Position in array", was
deliberately left as-is -- user-facing plot output, a separate decision
from the internal rename. See §22 for the full record.

**§23, §19's naming question -- acted on.** `repeat_sims` ->
`domain_distances` (`tales_align()`); `tal_sim` -> `tale_distances` and
`domain_sim` -> `domain_distances` (`plot.tales_msa()`); `tal_sim` ->
`tale_distances` (`tales_group_hclust()`/`tales_group_kmedoids()`);
`fill_type`'s `"repeat_clust"`/`"repeat_sim"` -> `"domain_clust"`/
`"domain_sim"`. A real breaking public-API change across four exported
functions -- man pages regenerated, `pkgdown::check_pkgdown()` clean, full
suite `FAIL 0 | WARN 35 | SKIP 0 | PASS 756`. Version bumped to
`0.9.9004`. **`vignettes/articles/tale_msa.qmd` and
`tales_msa_class.qmd`, initially left calling the old names, were fixed
the same session on request** -- package reinstalled first (this file's
own quarto-subprocess rule), `docs/` rebuilt in full via the established
decomposed sequence (never `build_site()`/`build_articles()` directly),
zero errors, `pkgdown::check_pkgdown()` clean. §21's own reserved matrix
helpers (`.domain_to_sim_align()` etc.) still take a parameter literally
named `domain_sim`, untouched, still the maintainer's to do -- see §23's
own record for exactly how the two renames meet at that call site without
colliding.

`tales_bind()` (§5.2), `tales_group()`'s split into
`tales_group_hclust()`/`tales_group_kmedoids()` (§11), `tales_compare()`'s
rename to `tales_compare_distal()`, and a new
`tales_compare_functal()`/`tales_to_universalmotif()` pair replacing the
unfixable Perl `functal()` for that one path (§12b) remain built and
unaffected by anything in this session.

Sections are kept in ascending numeric order within each chapter (fixed
2026-09-17, after §7 and §8 had drifted into add-order). Content only, never
renumbered -- ~35 code comments cite specific section numbers.

**Treat `[V]`/DONE markers as "verified once," not "still true."**
Repeatedly confirmed again this session (§17's own "49 exports" turned
out wrong on re-count; `dev/class-design.md` §4.6 sat stale for over a
week). Spot-check against the actual code before trusting a closed
section, especially before building on top of it.

## Commits

Commit when a piece of work is coherent and its tests pass. Say what moved
and why; reference the ledger section. Do not push without being asked.
