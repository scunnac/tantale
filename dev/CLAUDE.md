# Working on tantale

An R package for analysing TALEs (transcription activator-like effectors) of
*Xanthomonas*. Current work is pre-publication cleanup on the `main` branch.

New to TALE biology? [The Wikipedia
page](https://en.wikipedia.org/wiki/Transcription_activator-like_effector)
is worth a look before touching code that deals with repeats, RVDs or
target prediction.

**`dev/restructuring-notes.md` is the ledger** and the main source of
context: what has been done, what is deferred, and why. It was compacted
on 2026-09-23 to ~1600 lines (full earlier text: `git show
7ef1fe9:dev/restructuring-notes.md`). **Start with "Where things stand"
below, then the ledger's START HERE block**, which holds the ranked list
of pending issues; read only the sections bearing on the task in hand.

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
and need no `VignetteBuilder` entry. The old numbered `vignettes/*.Rmd`
set is gone, merged into the articles (§7.5b/c). The only traditional
vignette is `vignettes/getting_started.qmd`, which `R CMD build` does
render (`VignetteBuilder: quarto`).

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
invariants that are not addressed to the user. No other call remains in
`R/` (re-checked 2026-09-23).

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

**Delete `docs/`'s entire contents before a full site build, every time
-- maintainer's explicit instruction, 2026-09-22.** Building article by
article (above) only ever writes or overwrites the files a given render
produces; it never removes a file that source no longer accounts for --
a renamed or retired article's old `.html`/`.md` (and duplicated/stale
paths like an old `docs/articles/articles/...` nesting) simply stays on
disk otherwise, silently shipped alongside the real site. `rm -rf docs/*`
(or equivalent) first, then rebuild from an empty directory, every time a
*full* site build is done -- not needed for a single `build_article()`
check of one file while drafting.

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

**The exported interface is declared stable (maintainer, 2026-09-24), but
no deprecation cycle before 1.0.0.** Until the 1.0.0 release, renames and
removals stay hard (no alias, no `lifecycle` warnings), recorded in
`NEWS.md` as before. From 1.0.0 on, they follow the lifecycle package's
conventions: keep the old name working with a warning that names the
replacement (`lifecycle::deprecate_warn()`), remove it in a later release.

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
- **Writing tone (maintainer's repeated correction).** "Try to have a
  neutral, scientific tone. Humour is not discouraged if relevant." Two
  separate habits to cut: (1) the "A, not B" contrast in any position;
  (2) parallel/paired phrasing in general, with or without a "not"
  (balanced doublets, rhetorical pairs and triplets). Both show in "two
  readings of one geometry, not two independent alignments". State the
  fact once; if a contrast really matters, give it its own plain
  sentence. Factual enumerations are fine. Check all prose, ledger
  entries and commit messages included, before calling it done.
- **Every statement must match the rendered output.** Claims about a
  table, figure or return value in an article are checked against the
  render (`docs/articles/<name>.md` carries prose plus chunk output;
  figures read as images). §30 found several that did not match.
- **No implementation archaeology in user docs.** What changed and why
  belongs in code comments or the ledger, not in `@details`.
- **Put user-facing prose in the exported function's block**, not in an
  internal's `@noRd` — easy to get wrong when the real work lives in a
  helper.
- `@param ...` must name the arguments it accepts if it forwards to
  something the reader cannot open (ledger §8.1c).

## Where things stand

*Updated 2026-09-24, end of session. Re-check with `git log --oneline -5`
and `git status` before trusting any of it.*

### State

- Branch `main`, version **0.9.9010**, everything pushed at the end of
  2026-09-24 (check `git status -sb`). The installed
  `tantale` may be older: reinstall before rendering articles.
- **The GitHub repository is private** since 2026-09-24 (maintainer's
  choice, to change things without users). The README's install command
  works only with access; GitHub Pages for a private repository needs a
  paid plan, so the published site may be offline until it is public.
- Full `devtools::check()` on 0.9.9010, 2026-09-24: tests 0 failures,
  examples (with `--run-donttest`) and vignette OK. Its one remaining
  NOTE (no news entries in `NEWS.md`) is fixed: see "NEWS.md" below.
  The two check fixes made afterwards were re-checked with a quick,
  no-tests `R CMD check`.
- The site in `docs/` is current with 0.9.9010 as of 2026-09-24 (partial
  rebuilds: reference, home, news, llm docs, search, articles re-rendered
  one by one). It has not been wiped and rebuilt in full recently.
- The Claude Code pointer file is `.claude/CLAUDE.md` (it imports this
  file). It was moved out of the repo root because pkgdown publishes every
  root `.md` as a page, whatever `.Rbuildignore` says.

### Open items

The full, ranked list is in the ledger's START HERE block (reviewed
against the code 2026-09-23). Headlines:

**Pending issues deserving urgent action** (all done 2026-09-24):
1. Done: `inst/COPYRIGHTS` lists every bundled third-party file with its
   licence and source; GPL-3 text in `inst/tools/COPYING.GPL-3`; TALVEZ
   and QueTAL redistributed with their author's permission.
3. Done 2026-09-24: `NEWS.md`'s top heading now carries the version
   number, which R's news parser needs (rule under "NEWS.md" below).
5. Done 2026-09-24: `inst/legacy/docs_temp/` deleted.
7. Done 2026-09-24: §32.2 option 2 (identical rare RVDs score 1 in the
   `rvd_sim` fill); `rvd_dna_specificity`'s `NA` row fixed.
8. Done: README declares the interface stable; lifecycle conventions
   from 1.0.0 on.

**Decisions to make before 1.0.0:** `tales_rvd_strings(rvd_only =)` ->
`repeats_only`; §2; §20; the "Position in array" axis title (§22); §21
option (c); `tell_tales()`'s argument list; exposing ARLEM's
duplication/insertion costs; the distribution channel (§34).

**Reserved for the maintainer; do not start unasked:** §21 items 1-4 (the
matrix helpers of `plot.tales_msa()` and their tests), §20, §2.

**Deferred by the maintainer:** §30's parallel-phrasing sweep; an Rcpp
ARLEM (§33); rOpenSci and the one-archive plan (§34); §5.2.

### Done recently (details in the ledger)

- **§29.1/§29.2** `dev/function-graph.qmd`: a data-flow view (exported
  functions and the classes they take and return, recorded from the test
  suite by `dev/function-graph-dataflow.R` into a committed TSV) and a
  call graph of all 187 functions parsed from `R/`, folded by file.
  Render with `quarto render dev/function-graph.qmd`; the HTML is
  gitignored. Rerun the recorder after changing what an export takes or
  returns. Possible later reuse on the site (§29.2).
- **§29.3** `dev/function-graph.qmd` also lists every function,
  internals included, that takes or returns a matrix or a list of
  vectors, from `dev/function-graph-shapes.R` (traces all functions
  during the test suite, ~9 min, writes a dated TSV). Rerun it after
  changing what a function takes or returns.
- **§30** site prose review: 8 articles, README, `pkgdown/index.md` and
  all 63 published reference pages checked against rendered output or the
  code; many content errors fixed along with tone.
- **§32.1** `as_tales()` on a `tales_msa` now drops `alignment_width`.
- **§32.3** domain distances: DECIPHER backend counts gaps
  (`penalizeGapLetterMatches = TRUE`), Biostrings backend uses free end
  gaps; all three backends follow DisTAL's definition. Golden
  re-baselined, grouping unchanged.
- **§34** distribution findings recorded (sizes, licences, upstream
  checksums); `test_correct_tales.R` no longer reads an absolute
  `/home/...` path or writes to `~`; sweep found no other such case.
- **§32.4** `array_report.tsv` terminus lengths no longer count the stop
  codon.
- **§33** ARLEM computed in R (`R/arlem.R`), identical scores; executable
  removed (`inst/legacy/arlem_binary.R` keeps the old driver).
- **§25/§25b/§31** truncTALE article published; `tale_mining.qmd`
  correction sections rewritten.
- Earlier, still true: §17 documentation review and the `correct_tales()`
  flag-swap fix; §23 public renames (`domain_distances`,
  `tale_distances`, `"domain_clust"`/`"domain_sim"`); §26 history reset
  (old history in a git bundle on the `v0.1.9553` release).

**Treat `[V]`/DONE markers as "verified once"** and spot-check against the
code before building on a closed section.
Ledger sections are in ascending order and never renumbered (~35 code
comments cite them).

### Practices learned this session

- **Two Claude sessions may share this checkout.** Run `ListAgents` at the
  start. Stage explicit paths only, and before staging a shared file
  (`dev/restructuring-notes.md`, `dev/CLAUDE.md`, `NEWS.md`,
  `DESCRIPTION`) check its diff for hunks you did not write; stage only
  yours by building your version from HEAD and writing it to the index
  (`git hash-object -w` + `git update-index --cacheinfo`), since
  `git add -p` is interactive and does not work here. Message the other
  session before installing, wiping `docs/`, or rendering. A separate git
  worktree per session avoids most of this.
- **An untracked file in `R/` is picked up** by `load_all()` and
  `devtools::install()`. If another session has work in progress there,
  test and install from a scratchpad copy of the tracked files instead.
- **When a change alters `domain_distances` or `tale_distances`**, delete
  `vignettes/articles/_cache/compare.rds` and `group.rds` (keep
  `discovery.rds` unless `tell_tales()` output changed), reinstall, and
  re-render `articles/tale_classification` first, then `tale_msa`,
  `tales_msa_class` and `tale_target_prediction`. Check that the group
  numbers the articles hard-code (`group == 6`) still point at the same
  locus.
- **Prose must match the render.** After re-rendering, read
  `docs/articles/<name>.html` (or `.md` after `build_llm_docs()`) and the
  figure PNGs, and fix any quoted number that moved.
- **`build_home()` needs the network** (it queries CRAN for a link). If
  DNS fails, it aborts the build script; rerun the remaining steps once
  the network is back.
- **Golden changes:** use the `golden-rebaseline` skill and explain every
  changed row before accepting.

### Site builds

`pkgdown::build_site()`/`build_articles()` must not be called directly
(§15's quarto ordering bug). For a **full** rebuild, delete `docs/`'s
contents first, then run this sequence, read out of
`pkgdown:::build_site_local()` so nothing is skipped:

```r
pkgdown::init_site(".")
pkgdown::build_home(".")
pkgdown::build_reference(".")           # no `quiet` argument
pkgdown::build_articles_index(".")
# then pkgdown::build_article(name, pkg = ".") once per article,
# "articles/tale_classification" first to prime the shared cache
pkgdown::build_tutorials(".")
pkgdown::build_news(".")                # no `quiet` argument either
pkgdown:::build_sitemap(pkgdown::as_pkgdown("."))
pkgdown::build_llm_docs(".")            # bs_version 5; skip if bs3
pkgdown::build_redirects(".")
pkgdown::build_search(".")              # bs_version 5; build_docsearch_json() if bs3
pkgdown:::check_built_site(pkgdown::as_pkgdown("."))
```

For a **partial** update (no page added or removed), skip the wipe and
run only what changed: `build_article()` for edited articles,
`build_reference()` after roxygen changes (reinstall first), then
`build_home()` if README or `pkgdown/index.md` changed, `build_news()`,
`build_llm_docs()`, `build_search()` and the checks. Deleting `docs/`
before a full build matters: it is the only thing that removes pages whose
source is gone. `docs/articles/articles/<name>.html` files are
deliberate redirect stubs from `build_redirects()`.

## NEWS.md

**The top heading carries the current version number: `# tantale
0.9.9010`, never `# tantale (development version)`.** R's news parser
(`tools:::.build_news_db_from_package_NEWS_md()`) ignores every heading
without a version number, so a development-version heading makes
`R CMD check` report "No news entries found in NEWS.md". Whenever
`Version:` in `DESCRIPTION` changes, change the top heading to match, in
the same commit. Entries since the last release stay under that one
heading; at a release, the heading keeps the released number and the
next version starts a new one above it. Entries are `##` sections with
prose; each counts as one news entry, so bullets are not needed.

## Commits

Commit when a piece of work is coherent and its tests pass. Say what moved
and why; reference the ledger section. Do not push without being asked.
