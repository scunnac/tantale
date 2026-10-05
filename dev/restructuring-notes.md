# tantale — restructuring notes and action ledger

Working record of the pre-publication overhaul: findings, decisions,
deferred questions. Branch `main` (the only branch since the history reset
of §26).

**Compacted 2026-09-23.** Closed sections were cut down to their decisions
and to the facts still worth knowing; statements found wrong against the
code were corrected or removed. The full text before compaction (9661
lines) is in git: `git show 7ef1fe9:dev/restructuring-notes.md`. Commit
hashes quoted in that older text predate the history reset and no longer
resolve; the pre-reset history is in the git bundle attached to the
`v0.1.9553` GitHub release (§26).

Status markers:

- **[V]** verified against the code/data at the time of writing
- **[A]** agreed direction, not yet executed
- **[P]** parked: needs a judgement call or a dedicated pass
- **[superseded]** a record of a replaced plan; nothing to do

Section numbers are never reused or renumbered: code comments, tests and
articles cite them (`ledger §N`). Treat every `[V]` as
"verified once" and spot-check it against the code before building on it.

---

## START HERE

`dev/CLAUDE.md` ("Where things stand") is the short briefing; this block
is the detailed list it points to. Updated 2026-09-23, after a full review
of this file against the code.

### Pending issues deserving urgent action

Ranked by risk; re-checked against the code on 2026-09-23. All eight
were done on 2026-09-24 (outcomes kept below).

1. **DONE 2026-09-24: licence notices for the bundled tools (§34).**
   - `inst/COPYRIGHTS`: every bundled third-party file with its
     copyright holder, licence, upstream URL, sha256 and reference
     (citations checked against doi.org). Covers the three Jstacs jars
     and `talecorrect/`, TALVEZ 3.2, and QueTAL FuncTAL in
     `inst/legacy/` plus the `rvd_dna_specificity` dataset taken from it.
   - `inst/tools/COPYING.GPL-3`: the GPL-3 text, byte-identical to the
     `COPYING.txt` inside each jar (it differs from R's own copy only in
     `http`/`https` URLs).
   - DESCRIPTION `Copyright:` field pointing to `inst/COPYRIGHTS`;
     README "Licence" section. `LICENSE`/`LICENSE.md` untouched (CRAN's
     MIT template; GitHub's licence detection).
   - **Maintainer, 2026-09-24:** permission to redistribute TALVEZ and
     QueTAL obtained from A. L. Pérez-Quintero (the maintainer co-authored
     both papers); the notice says so. No `cph` entries in `Authors@R`
     (maintainer's choice).
   - Checked: all three tools' sources are in github.com/Jstacs/Jstacs
     (`projects/xanthogenomes/`, `projects/tals/prediction/`,
     `projects/talecorrect/`); the jars hold no `.java` files.
   - **Two residual points, both closed by the maintainer 2026-09-24**
     (no further action: (a) "no need to go further"; (b) due diligence
     considered done):
     (a) the source pointer is Jstacs' `master`, not the exact revision
     each jar was built from (Jstacs publishes no versioned source for
     them). Usual practice for unmodified redistribution; a source
     snapshot could go into the one-archive release if wanted.
     (b) TALVEZ's `simplescancode/*.class` are, per the script's header,
     "java code from Matzieu and Hatzigeorgiu 2010" (DIANA PlantTFBS).
     The permission obtained is the TALVEZ author's; whether it covers
     that third-party code is not established.
   Record of the problem:
   **Licence notices for the bundled tools (§34).** The three jars
   (AnnoTALE, PrediTALE, TALEcorrection) are GPL-3 and may be
   redistributed only with the licence and a pointer to their source
   (github.com/Jstacs/Jstacs). Nothing in the package gives either, and
   `LICENSE`/`DESCRIPTION` present the whole package as MIT. TALVEZ has
   no licence at all; the maintainer is asking its author. GitHub
   distribution is already distribution, so this applies today, whatever
   channel §34 ends up choosing. Cheap interim fix: a third-party notice
   (licence text, upstream URL, source link per tool) shipped with the
   package and referenced from `LICENSE`/README.
   **2026-09-24: the GitHub repository is private for now** (maintainer,
   to change things freely without users). While it stays private this is
   not distribution, so the notices must be in place before it goes
   public again, or before anyone outside gets a copy.
2. **DONE 2026-09-24: every external program's exit status is checked.**
   nHMMER (`.run_nhmmer_search()`, `correct_tales()`), PrediTALE and
   TALEcorrection go through `.tantale_exec()`; AnnoTALE's three stages
   through a new `.annotale_exec()` (`R/annotale.R`), which keeps the
   `tantale_error_annotale_failed` class and, in quiet mode, replays
   AnnoTALE's stderr. Inside `tell_tales()` an AnnoTALE failure still
   skips the array with a warning, which now carries AnnoTALE's error as
   its parent. `correct_tales()`'s three nHMMER searches were joined with
   `"; "`, so only the last one's status was ever seen; now `&&`.
   `run_annotale_predict()`/`run_annotale_build()` also had unquoted paths
   (contradicting §8.0b's "every command quotes its paths"); fixed.
   Checked: golden baseline unchanged; PXO86 still gives 18 arrays;
   `test_external_exit_status.R` (8 tests) fails on the old code for
   every path the old code left unchecked. `.check_hmmer()` already
   checked its status and is unchanged.
3. **DONE 2026-09-24: full `devtools::check()`** (0.9.9010, `--as-cran`,
   17 min): tests pass in the isolated library (0 failures), all examples
   including `\donttest{}` pass, vignette builds. It found:
   - a real bug: `.functal_pwm()` read `rvd_dna_specificity` by its bare
     name, which resolves only when tantale is attached, so
     `tantale::tales_to_universalmotif()` (and `tales_compare_functal()`)
     failed without `library(tantale)`. Now `tantale::rvd_dna_specificity`;
     a test in `test_tales_compare_functal.R` fails if any package
     function reads a dataset by bare name;
   - a broken Rd link (`universalmotif-class`), four undeclared globals
     (`:=` now imported from rlang; `alignment_position`, `sq_len` in
     `R/globals.R`). Fixed; a quick re-check shows neither any more;
   - **still open, a decision:** NOTE "No news entries found in NEWS.md".
     The first diagnosis here (the parser wants bullets) was wrong.
     **Re-checked 2026-09-24** by reading
     `tools:::.build_news_db_from_package_NEWS_md()`: it skips every
     heading until one contains a version number, and our only top
     heading is `# tantale (development version)`, so it finds nothing.
     Under a version heading, each `##` section becomes one entry
     (category = heading, text = its paragraphs); bullets are not needed.
     With the first line changed to `# tantale 0.9.9010`, the parser
     returns 21 entries, one per `##` section. The fix is the heading;
     bullets are a separate, purely stylistic choice.
     **DONE 2026-09-24 (maintainer's decision):** top heading changed
     to `# tantale 0.9.9010`; prose sections kept. Standing rule in
     `dev/CLAUDE.md` ("NEWS.md"): the heading follows `DESCRIPTION`'s
     version at every bump.
   - INFO only: 38 non-default Imports; installed size 88 MB (§34).
4. **DONE 2026-09-24: ggplot2 deprecations.** `label.size` and the line
   `size` replaced by `linewidth`; `ggplot2 (>= 3.5.0)` declared (the
   first version where `geom_label()` takes `linewidth`). Five
   representative figures pixel-identical before and after; no
   deprecation warning left in the plot tests.
5. **DONE 2026-09-24: `inst/legacy/docs_temp/` deleted** (maintainer's
   decision). An identical copy (plus an `.Rhistory`) is in
   `tantale-old-before-reset/docs_temp/`. The `.gitignore`/`.Rbuildignore`
   entries are kept. Record of the problem:
   **`inst/legacy/docs_temp/`** (3.8 MB of old notebooks, untracked and
   gitignored) sat inside `inst/`. **Confirmed 2026-09-24:** a tarball
   built from this checkout contains its 18 files. The maintainer's
   files: delete them, or move them out of the package tree. Note that
   `extra/` does not exist in this checkout (untracked since §26); it
   survives only in `tantale-old-before-reset/extra` (274 MB).
6. **DONE 2026-09-24: tests for the untested exports.**
   `test_talomes_heatmap.R`, `test_annotale.R`, `test_tales_consensus.R`,
   and `preditale()`/`plot_target_preds()` tests in
   `test_target_predictions.R`, written against known answers where one
   exists (a site covers position 0 plus one base per RVD; the toy
   regions' frameshifted copy is AnnoTALE's pseudogene; the
   `tales_msa`-native consensus equals the matrix one). Writing them found
   three `talomes_heatmap()` bugs, fixed with the default figure
   pixel-identical: with `plot_type = "all"` the `title`/`x_lab`/`y_lab`
   arguments were ignored; `par()` read before the output device was
   chosen left an extra device open after `save_path`; an explicit
   `save_path = NULL` failed. An unknown `plot_type` is now an error.
   Coverage re-measured 2026-09-24: 90.03% (was 73.96%); README badge
   updated. Lowest files: `tantale_conda_env.R` 59% and
   `tantale_setup.R` 67% (install branches, which would modify the
   machine), `requirements.R` 73%, `tales_ingest.R` 76%, `conversion.R`
   78%.
7. **DONE 2026-09-24: §32.2, option 2**, plus a fix for the nameless
   `NA` row of `rvd_dna_specificity` found on the way. Record of the
   question: **§32.2 decision** (`rvdSimDf`, 17 RVDs). Narrower than first thought:
   the alignment side already scores an unknown RVD pair as neutral (0)
   and any RVD against itself as 1 (`.rvd_score_table()`). Only
   `plot.tales_msa(fill_type = "rvd_sim")` still greys out identical RVDs.
8. **DONE 2026-09-24, maintainer's decision: the interface is stable.**
   README now says so, and that from 1.0.0 on changes follow the
   lifecycle package's conventions (deprecate with a warning first,
   remove later). **No deprecation cycle before 1.0.0** (maintainer): until
   then renames stay hard. `NEWS.md` announces it; rule in
   `dev/CLAUDE.md`'s API conventions.

### Decisions to make before 1.0.0 (breaking or user-visible)

Not urgent in themselves, but cheaper before the release: until 1.0.0 a
rename is a hard rename (no deprecation cycle, maintainer 2026-09-24);
from 1.0.0 on it needs a lifecycle deprecation (START HERE item 8).

- **DONE (§37):** `tales_rvd_strings(rvd_only =)` renamed
  `repeats_only`, like its sibling (§8.2b).
- **DONE (§37, D1):** `repeat_to_rvd_map()` (§2) and
  `tale_parts_to_rvd()` (§20) retired to
  `inst/legacy/conversion_retired.R`, with `repeat_to_rvd_map_distalr()`.
- **DONE 2026-10-02 (§22):** `plot.tales_msa()`'s x-axis title now
  reads "Position in alignment".
- **Decided 2026-10-02 (§43):** `tell_tales()` keeps its flat argument
  list; its two `TODO` blocks are settled.
- **Moved after 1.0.0 (§43):** §21 option (c) and ARLEM's costs. Both
  can be added without breaking a call.
- **Distribution channel** and the one-archive plan (§34).
- **DONE 2026-10-03 (§48): should functions that currently consume
  tales strings (e.g. output of `tales_rvd_strings()`) also accept a
  `tales` object directly, calling `tales_rvd_strings()` internally?**
  Answer: only `talomes_heatmap()` needed it; done.
  Raised by the maintainer. Touches argument signatures and therefore
  cheaper to decide before 1.0.0. Relevant functions to identify: any
  exported function whose first argument is a character vector of RVD
  strings rather than a `tales` object. Decide before coding.
  **Survey of the 51 exports, 2026-10-03:** three take RVD strings.
  `talvez()` and `preditale()` (`rvd_seqs`: fasta path or `XStringSet`)
  already have a `tales`-aware front end, `tales_predict_targets(x)`,
  which takes a `tales`, a fasta path or a `BStringSet` and renders a
  `tales` with `tales_rvd_strings()`. The one real gap is
  `talomes_heatmap(tale_annotation, group_col, strain_col, rvd_col)`:
  `tale_classification.qmd` spends five lines building that data frame
  from a grouped `tales` (RVD strings, then a join on `group` and
  `strain`). Possible answers for the maintainer: (a) `talomes_heatmap()`
  also takes a `tales` carrying the group and strain columns and computes
  the RVD strings itself; (b) leave it, the data frame keeps it usable on
  tables from other tools; (c) as (a), and ask whether `talvez()` and
  `preditale()` should stay exported next to `tales_predict_targets()`.
- **DONE (§42, §49): the terminus check's coverage.** The match must
  reach the profile end next to the repeats (§42). No rule on the far
  end (maintainer, 2026-10-04, Q63): a short terminus that matches stays
  `NTERM`/`CTERM`, and the docs explain that its function is probably
  impaired (`?tales_anchor_codes`, `trunctale_correction.qmd`).
- **A beta of 1.0.0** (decided 2026-10-02, §43: yes, as a release
  candidate when the repository goes public or the rOpenSci submission is
  made): number it **0.99.0** (then 0.99.1, ...).
  R versions are numeric only (`1.0.0-beta` is invalid, `1.0.0-1` sorts
  after 1.0.0); 0.9.9010 < 0.99.0 < 1.0.0; Bioconductor also requires
  0.99.z for new submissions. Mark it with a `# tantale 0.99.0` NEWS
  heading, a `v0.99.0` tag and a GitHub *pre-release*. Hard renames are
  still allowed at 0.99.x (lifecycle rules start at 1.0.0).

### Housekeeping pending (2026-09-24) -- both DONE the same day

- **DONE:** partial rebuild (reinstall; reference, home, news, llm docs,
  search; `trunctale_correction` re-rendered, only its printed timings
  changed). `pkgdown/index.md` got README's Licence section, since the
  home page is built from it and showed no licence information.
  **Correction, same evening:** that partial rebuild left six articles
  and the articles index rendered with 0.9.9005-0.9.9008 (the check only
  looked for 0.9.9009). A **full rebuild** followed (`docs/` and
  `_cache/` emptied, every step of "Site builds", 16 min): every page is
  now 0.9.9010. Compared with the previous build: same file list, all 25
  figures byte-identical, article outputs differ only in printed run
  times and in source line numbers inside printed error messages
  (`tales_class`). **Finding:** on a cold cache,
  `tale_classification.qmd`'s `discover_all` chunk prints DECIPHER's
  progress bars and timings (21 lines; printed output, which
  `message = FALSE` does not hide). A warm-cache render does not, so the
  article was re-rendered once more and is identical to the previous
  one. A lasting fix (e.g. `results = "hide"` on that chunk, or
  `capture.output()`) is left for later. To check stamps after a partial
  rebuild: `grep -rL "0.9.9010</small>" docs --include=*.html` (the
  redirect stubs in `docs/reference/` carry no stamp).
  Record of the item: **`docs/` is behind 0.9.9010**: reference pages (the `run_annotale_*()`
  return values, `talomes_heatmap()`), home page (README's stability
  note) and news. Partial rebuild per `dev/CLAUDE.md` (reinstall first;
  no page added or removed, so no wipe). The site may be offline anyway
  while the repository is private.
- **DONE** (see §29.3, "Done 2026-09-24"). Record of the item:
  **Re-run `dev/function-graph-dataflow.R`** (now ~10 min, the suite is
  longer) so the §29.2 data-flow view records the functions tested since:
  `talomes_heatmap()`, `preditale()`, `plot_target_preds()`,
  `run_annotale_*()`.

### Worth investigating, no deadline

- `correct_array = TRUE` fabricated a ~489 nt region on a two-array PXO86
  excerpt that does not exist in the genome; the fixture
  `dev/fixtures/pxo86_roi18_19_excerpt.fa` reproduces it in under a
  minute (§25). Did not reproduce on the full genome.
- **Settled (§49):** the article cache builds BAI3-1-1 at `max_comparisons = 50`, now the default (§25b).
- Validate the 136-sequence correction reference against the 494 one
  (maintainer's own, §8.1b).
- Each toy region yields a spurious single-hit array at its 3' end,
  because `min_domain_hits` filters subject sequences (§8.1).
- **Closed:** `.rvds_from_annotale_file()` is called by
  `.tale_parts_assemble()`, which `tales_from_telltales()` and
  `tales_from_annotale()` both use (§35, §37; §29.1).
- Consumers of `tales_rvd_strings()` other than
  `tales_to_universalmotif()` were never checked for silently dropping
  an array with no repeats (§12b).
- **DONE 2026-10-02 (Q50):** the last errors with only the generic
  `tantale_error` class were five, not seven (`conversion.R` 2,
  `talecorrection_java.R` 1, `tales_ingest.R` 1, `tales_plot.R` 1;
  found with the parser). They now carry `tantale_error_bad_argument`,
  `_missing_file` or the new `_annotale_file`; their messages were
  rewritten, and `test_error_conditions.R` asserts the classes (the
  `tales_plot.R` one sits in an unreachable `else`, untested). Still
  open: 13 `cli_warn()` calls carry no class at all, so tests can
  only match their wording (§9.5).
- Rendered error messages from an installed package show a source path
  ("at tantale/R/tales_class.R:818:3"); cosmetic (§30).
- **DONE 2026-10-02 (Q50):** commented-out developer snippets in `R/`
  no longer carry `/home/cunnac/...` paths (`system.file()`,
  `reticulate::conda_binary()`, `tempdir()` instead). Left alone: the
  recorded run in `inst/extdata/tellTaleExampleOutput/` (its
  `tell_tales.log` and `hmmer_search_out.txt` print the paths of the
  machine that produced them) and the upstream ant file
  `inst/tools/talecorrect/TALEcorrection.xml` (§34).
- **DONE 2026-10-02 (Q50):** the stale `pipeline.svg`/`.png` moved from
  `man/figures/` to `dev/figures/` (§7.6).

### Parked, reserved or deferred

- **Reserved for the maintainer, do not start unasked:** §21 items 1-4
  (the matrix helpers of `plot.tales_msa()` and the tests pinning their
  shape). §2 and §20 were settled by the §37 retirements.
- **Deferred by the maintainer:** §30's parallel-phrasing sweep; an Rcpp
  ARLEM (§33); rOpenSci (§34); §5.2 (talome-wide MSA plot, leaning
  shape A).
- **Dropped from tracking, unless raised again:** §9.3, §8.6b, §12b's
  `universalmotif` follow-ups.

### Cross-cutting lessons

Rules already in `dev/CLAUDE.md` are not repeated here.

- **An `R CMD check` finds what `load_all()` hides**: undeclared
  dependencies still installed in the ambient library, the core limit
  (`_R_CHECK_LIMIT_CORES_`, at most 2), relative `test_path()` results
  (§27).
- **A golden digest pins a code path and cannot tell a right answer from
  a wrong one.** Where a known answer exists, assert it (§8.1's toy
  fixture). When diagnosing golden mismatches, give the control run the
  same output-directory basename, and keep diagnostic output inside
  `tempdir()` (§27).
- **Read the source when the docs are silent** (§8.1b: DECIPHER's
  `maxComparisons` ranks before truncating).
- **Genomes are gold-quality unless flagged** (BAI3-1-1 is flagged): a
  frameshift found in one is biology (§25).
- **`Filter()` over a predicate that can return `logical(0)`** silently
  shifts every later keep/drop decision (§7.7).
- **An `Edit` whose `old_string` ends on a heading can delete the
  heading**; re-read the section boundaries afterwards (§7.5c).
- **cli pluralisation needs a quantity in the same bullet**
  (`{cli::qty(n)}`); an inline style span such as `{.fn x}` is not one
  (§8).
- Headless Chromium for screenshots: the snap build only reads and writes
  inside `$HOME`, and does not open a debugging port under `chromote`;
  the Playwright build in `~/.cache/ms-playwright/` works (§7.5c, §29.1).

---

## 0. Framing

`tale_parts` (now the `tales` class) arrived late, while `distalr()` was
being written, to replace a family of conversion functions that reshaped
`runDistal()` output for MSA and plotting. Two cautions that apply
throughout:

- **Absence from the old workflow figure says nothing about a function's
  value.** The package also ships standalone utilities.
- **Zero internal call sites says nothing either.** Exported functions
  are meant for users, and an unused internal may be dormant.

---

## 1. Shrink `distalr()`'s returned list -- OBSOLETE, `distalr()` no longer exists

The analysis showed that three of `distalr()`'s six returned elements were
projections of `tale_parts`, and that only `tale_parts` and the two
similarity tables were irreducible. It motivated the class design (§5);
`distalr()` itself is gone. The projections exist as
`tales_coded_strings()`/`tales_rvd_strings()`/`tales_domain_codes()`.

---

## 2. Conversion functions -- retirements done in §37 **[V]**

- **`repeat_to_rvd_map()` is redundant with `tales`** (same 251 rows as
  the `dom_code`/`rvd` columns of the fixture). Its one-RVD-per-repeat
  assertion is pinned by `test_untested_exports.R`, so retiring it cannot
  drop that check silently; the `aa_seq` -> `rvd` anomaly check of
  `tales()` (§6) covers the same ground. **Retirement not executed**; the
  function is still exported. Reserved for the maintainer ("we will see
  that later").
- `repeat_to_rvd_map_distalr()`'s name is misleading: it needs a
  `dom_code` column, and `distalr()` no longer exists.
- **Dormant internals repaired and wired up [V]**: `.rvd_to_match_align()`
  is `fill_type = "rvd_sim"` (RVD-level specificity similarity to the
  reference, signed on [-1, 1], diverging scale; termini have no RVD and
  render grey); `rvdSimDf` also scores RVD alignments
  (`tales_align(domain_distances = "rvd")`). Lesson: they were only
  unwired, and the capability existed nowhere else. See §32.2
  for `rvdSimDf`'s coverage.

---

## 3. Legacy cemetery (`inst/legacy/`) — DONE **[V]**

The AnnoTALE <-> QueTAL shims are in `inst/legacy/annotale_quetal_shims.R`
beside `tellTaleLegacy.R`: uncalled, kept because the file conventions they
encode are recorded nowhere else. `inst/legacy/` ships but is never
sourced. It now also holds `annout_class.R`, `arlem_binary.R` (§33),
`functal.R` + `QueTAL_v1.1/` (§12b), `msa_heatmap.R` (§4) and
`unused_pending_review.R` (§18).

---

## 4. Plotting: `msa_heatmap()` — RETIRED **[V]**

Retired to `inst/legacy/msa_heatmap.R` once `plot.tales_msa()` could draw
a consensus. Every `plot_type` of the old function maps onto separate
arguments of the new one (`fill_type`, supplying an RVD layer, a
reference pattern, `consensus = TRUE`).

Implementation note worth keeping: the consensus is a **separate `aplot`
panel**. `aplot::insert_left()` reorders the main plot's y axis onto the
tree's leaves and silently drops a y level with no matching leaf, so a
consensus "row" inside the alignment cannot survive a tree. The panel's
height is a ratio of the main plot, set to `1 / n_arrays` (clamped) so it
stays about one row high.

---

## 5. OOP restructuring

Class definitions (identity, invariants, constructors, method policy) live
in **`dev/class-design.md`**.

Settled and built **[V]**: additive S3 classes over tibbles (dplyr keeps
working); per-stage classes and no session/project container; the long
table as the canonical alignment form, with matrices as `as.matrix()`
views; snake_case columns; two families, `tales` -> `tales_msa` and
`pairwise_distances` -> `tale_distances`/`domain_distances`. The S4
`annout` class was retired to `inst/legacy/` (§5.3).

The structural problem this solved: every core shape existed as both a
wide matrix and a long table, converted ad hoc at each call site (17
`melt`/`acast`/`dcast` calls, 16 positional `colnames<-`).

**`dom_code` run-dependence** is the invariant with the most bite:
`cur_group_id()` over `aa_seq` mints `1..N` on every run, so a join across
runs succeeds and silently maps domains to the wrong sequences. Enforced
by the `dom_code_namespace` stamp (a content hash) on `tales` and
`domain_distances`, checked by every method that consumes two of them.
`tale_distances` carries no stamp on purpose: `array_id` does not collide
across runs.

### 5.1 Does `diagnose_tale_parts()` survive the `tales` class? — RETIRED **[V]**

Gone. `tales_anomalies()` reports the same conditions and
`tales(sanitize = TRUE)` drops the offending arrays.

### 5.2 Talome-wide MSA summary plot **[P]** — parked for much later, leaning A

Idea: every group's alignment in one faceted figure,
`plot.tales(position = "alignment")`. Each group aligns independently, with
its own width and coordinate system, so a concatenated `tales_msa` would
be valid in structure and meaningless in content. Two shapes:

- **A**, a list of `tales_msa`, one per group (maintainer's leaning,
  2026-09-19);
- **B**, a plain `tales` carrying `group` and `alignment_position` as
  ordinary columns, faceted with `scales = "free_x"`.

A group-aware `tales_align()` would follow from A with no bind method.

Settled along the way **[V]**:
- `tales_group_*()` take the `tales` and return it with `group` filled;
  taking `x` is the one point where "these distances came from this
  object" can be checked (`tantale_error_group_mismatch`).
- **`tales_bind()`** (built, `R/tales_class.R`):
  `tales_bind(..., on_namespace_mismatch = c("recode", "error"),
  sanitize = FALSE)`. Refuses a `tales_msa` (demote with `as_tales()`
  first); hard error on colliding `array_id`s; **drops `group`** (a
  clustering result of one distance matrix, meaningless across runs);
  on a namespace mismatch, recodes `dom_code` over the union and informs
  that any companion distance table is stale; never touches distance
  tables (re-run `tales_compare_distal()` on the result). Named
  `tales_bind()` because base `c()` dispatch drops attributes silently on
  mixed inputs and has no room for the mismatch argument. Tests:
  `test_tales_bind.R`.
- **Which operations invalidate a companion distance table** (worth a
  user-facing home): subsetting never does (filter the table to surviving
  ids); binding does. `tale_distances` then needs only the new cross
  pairs, while `domain_distances` needs rekeying to the new `dom_code`s
  first. Attaching the distance tables to the `tales` object as slots was
  rejected: every tibble verb would have to keep them honest.

### 5.3 `tell_tales()` refactoring — MECHANICAL PASS DONE **[V]**

745 -> 160 lines, 17 named internals, the body reads as a pipeline. Each
extraction was checked against a baseline in both directions. Two defects
no test would have caught: the GFF export tested
`exists("reducedOlapGr")`, which silently became `FALSE` once the merge
moved into a function (199 -> 103 records); removing the `annout` class
removed an `@import Biostrings` four calls depended on.

**Not done, on purpose:** the 17 arguments are unchanged (grouping them
changes the user interface), and two `TODO` blocks remain in the body
(circular molecules; what two output files should contain), both
questions for the maintainer.

### 5.4 Reassemble whole-TALE sequences from ordered domain parts — DONE **[V]**

`tales_get_protein_seq()`/`tales_get_dna_seq()` (`R/tales_projections.R`):
one `AAStringSet`/`DNAStringSet` per array, parts pasted in
`position_in_array` order, names = `array_id`. A `tales_msa` gives the
same result, since gaps are never rows. No terminus fallback (anomaly
checks own that); an NA guard stops a missing sequence being pasted as
"NA".

---

## 6. Correctness review backlog -- DONE **[V]**

All items closed. Kept for reference:

- **The `as.dist()` inversion**: only `.repeat_to_cluster_align()` (now
  `.domain_to_cluster_align()`) clustered a similarity as if it were a
  distance; fixed (45 clusters instead of 59 on the reference output),
  `h_cut` default moved 90 -> 10.
- **ggplot2 API drift**: `plot_tales_msa()` aborted on every call under
  ggplot2 4.0.3 (`palette =` passed to `scale_fill_manual()`); replaced by
  `discrete_scale()`. The same kind of drift (deprecated `label.size`/
  `size`) was fixed on 2026-09-24 (START HERE item 4).
- `.build_repeat_msa()` takes the residue type from `tales_align()`
  instead of guessing from six frequent RVDs.
- `aa_seq` -> `rvd` consistency is a soft anomaly
  (`aa_seq_rvd_inconsistent`); `aa_seq` <-> `dom_code` is a hard bijection
  in `validate_tales()`.
- `tell_tales.log`: the subject file is logged as given (it used to show
  the renamed temp copy), and every parameter row is two tab-separated
  fields.
- The `tidytree::MRCA()` "Invalid edge matrix" noise is a cli message, not
  a failure; silenced at its one call site.
- The hclust dendrogram cut line is drawn at `cutOff` (was `cutOff / 2`,
  no reason found).
- `plot.tales_msa()`'s collected legend is placed through
  `aplot::as.patchwork() & theme(...)`, the only level patchwork reads;
  the returned value is still the `aplot`.
- `plyr` removed (three sites in `telltale.R`).
- `tales_consensus()` ties: see §8.7.

**ARLEM score semantics.** `arlem_score` is a cost; `tale_distances`
stores the normalised cost (`arlem_score / max_length`) as `dissim`
directly (the old `Sim = 100 - cost` round trip is gone, §9.6). Two facts
still hold:
- The range is compressed: aligned TALEs mostly match, so `dissim` spans
  about 0-8 of its nominal 0-100 (0-7.04 on the article data, §32.3).
- **Duplication and insertion costs are fixed at 10**
  (`.arlem_dup_cost`, `.arlem_indel_cost`, `R/distalr.R:336`) against
  substitution costs up to 99. Undocumented to users and not a parameter;
  a biological modelling choice (see START HERE, decisions before 1.0.0).

---

## 7. Documentation and artifacts

- Self-deprecatory language in README, `?tantale` and comments rewritten
  (2026-09-18); `?tantale`'s four `@section` tags fixed on the way
  (roxygen2 8.0 misparsed them).
- `talomes_heatmap()`'s rows/columns docs corrected.
- `DESCRIPTION` carries `Config/roxygen2/version: 8.0.0` (roxygen2 8.0
  replaced `RoxygenNote`).
- `pipeline.svg`/`.png` (now in `dev/figures/`): see §7.6.

### 7.1 Vignettes 1-4 cannot be built **[V]** — SUPERSEDED, see §7.5b/§7.5c

The old numbered vignettes were chained through `save.image()` to a
hardcoded home path and could not build anywhere else. Rewritten, then
merged into the articles.

### 7.2 `R CMD check` results **[V]**

First full check (2026-09): undeclared `::` imports added, six unused
Imports removed (`GenomicFeatures`, `RColorBrewer`, `dichromat`,
`optparse`, `scales`, `msa`), broken Rd links and stale `@param`s fixed.
The NSE false positives are declared in `R/globals.R`; ~40 "no visible
global function" notes were all real and fixed. Over-long fixture paths
(65 over the 100-byte tar limit) fixed by renaming four test fixture
directories (`err_missing_dna`, `err_array_count`, `err_missing_nterm`,
`example_output`; the shipped `inst/extdata/tellTaleExampleOutput` kept
its name). The executable-file warning went with ARLEM's binary (§33).
Last full check: §27. See START HERE item 3.

### 7.3 A systemic habit: unqualified calls to non-imported packages **[V]**

Three defects of one shape (`mutate()`/`ggplot()`, `Quote()`,
`countMatches()` called bare). `Imports` makes a package installable;
only `@import`/`@importFrom` or `pkg::` makes it visible. **In this
package, a bare call to anything outside base is suspect**; `R CMD
check`'s "no visible global function definition" finds them.

### 7.4 Shrink the payload — DONE **[V]**

`inst/tools` 122 -> 63 MB: MAFFT and HMMER come from the conda
environment. HMMER 3.3 -> 3.3.2 changed nothing but the banner. **MAFFT
7.520 changes results** (different column counts, termini no longer
anchored at the first and last columns; no gap-penalty setting fixes it);
bioconda's 7.453 is byte-identical to the old bundled 7.450, hence the
pin. The three jars stayed (no conda package; §34).

### 7.4a `tantale_setup()` -- DONE **[V]**

Checks, and on request builds or repairs, everything outside R. Versions
are compared with the pins parsed from
`inst/tools/tantale_conda_env.yaml` (read from `conda-meta/`, no
subprocess); repair installs the pinned specs explicitly (`create` does
not downgrade); Java and Perl are checked; three paths are reported
(conda binary, default root, environment in use), because on a machine
with history they diverge; installing conda is opt-in. The lazy path
warns once per session when the environment is off its pins.

### 7.4b Tell users how to get conda, and that they now need it -- DONE **[V]**

README and `?tantale` explain the three-step setup (install, make conda
available, `tantale_setup()`). `_pkgdown.yml` lost a dangling
`annout-class` entry, which had made `build_site()` impossible.

### 7.5 Worked examples: vignettes and `@examples` **[A]**

Superseded by §7.5b/§7.5c and §17: every exported function has a run and
verified example.

### 7.5a An article on the `tales` class, and what a `dom_code` is -- DONE **[V]**

`vignettes/articles/tales_class.qmd`. Key points it carries: a `dom_code`
names a distinct **domain** sequence (N-terminus, repeat or C-terminus);
coding repeats as symbols is what makes arrays alignable; it is finer
than an RVD; the codes mean nothing outside the call that minted them.
Cross-references in plain `.qmd` are inline code (`` `fn()` ``), which
downlit autolinks; roxygen link syntax is dropped silently.

### 7.5z Articles are Quarto, not R Markdown -- agreed **[V]**

In `dev/CLAUDE.md`'s standing rules.

### 7.5b The four numbered vignettes rebuilt; a `tales_msa` article; `@examples` on the core API -- DONE **[V]**

Rebuilt self-contained, then merged into the articles (§7.5c). Found on
the way: `plot.tales_msa()` returned visibly (every figure printed twice),
now `invisible()`.

### 7.5c The numbered-vignette/article duality abrogated; genome choice and figure sizing fixed; one real vignette reinstated -- DONE **[V]**

- All walkthroughs are `vignettes/articles/*.qmd` (pkgdown only; R's
  build never descends there). One real vignette,
  `vignettes/getting_started.qmd` (`VignetteBuilder: quarto`), built from
  a shipped fixture so it needs no external tool at install time.
- **Why the other articles cannot become real vignettes as written**
  (relevant to §34): the quarto vignette engine renders a minimal,
  unthemed page; building vignettes needs the quarto CLI at install time;
  `tale_mining.qmd` runs `tell_tales()` for minutes and needs the conda
  and Java tools, which `tantale_setup()` only builds on first use.
- PXO86 is left out of the classification/alignment/prediction articles:
  with the African strains it turns one clean group per locus into
  singletons and paralog clusters.
- `_pkgdown.yml`: an `articles:` section needs `navbar: ~` on its group to
  keep the Articles dropdown (`pkgdown:::navbar_articles()`); quarto
  cannot resolve `@sec-` cross-references across separate articles. Theme
  `sandstone`.

### 7.6 README/index/docs follow-ups **[V]** five of five done, plus a sixth

Done: README section on the use of large language models; a dedicated
`pkgdown/index.md`; DisTAL and functal named on the index page; `pak`
install instructions; README and `?tantale` descriptions reconciled on the
two concrete errors (tool-name styling left to drift, maintainer's
choice); lifecycle badge `stable` (maintainer's call); a static coverage
badge with a dated caveat in README (73.96% on 2026-09-21, 90.03% on
2026-09-24).

Open:
- The "stable" badge against README's own "interfaces may still change"
  (`README.md:116`). START HERE item 8.
- **`pipeline.svg`/`.png`** show pre-restructuring function names and
  were referenced from nowhere outside `dev/`. Moved from `man/figures/`
  to `dev/figures/` on 2026-10-02 (maintainer's go-ahead, Q50), so they
  no longer ship (~740 KB). Redraw or replace with the §29.2 data-flow
  view if a workflow figure is wanted. Re-exporting the PNG needs
  Inkscape (librsvg changes the typography).

A one-off: `extra/tantale_logo.png` vanished during one `covr` run and
never again; restored from git, cause unknown.

### 7.7 Version scheme decided: 0.9.x pre-publication, 1.0.0 at release; docs/ now tracks the latest build — DONE **[V]**

`0.9.x` until publication, `1.0.0` at release; current `0.9.9009`.
`docs/` tracks the latest build, published in release mode
(`_pkgdown.yml` `development: mode: release`, §15). Background:
`pkgdown:::dev_mode_auto()` treats a version as "devel" only when its
third component is >= 9000, which is why `0.9.1` once wrote a release
build straight into `docs/`; with `mode: release` hardcoded this no
longer matters.

pkgdown publishes every root `.md` except README/LICENSE/NEWS as a site
page (`pkgdown:::package_mds()`, no exclude option); the Claude Code
pointer file therefore lives in `.claude/CLAUDE.md`. A post-build script
that edited `search.json` was abandoned after it silently dropped 59
entries (the `Filter()` trap in START HERE).

### 7.8 Audit docs and website for "repeat" used where "domain" is meant -- DONE **[V]**

Seven instances fixed across `distalr.R`, `pairwise_distances_class.R`,
`tales_msa_class.R`, `tales_plot.R`, `functal.R` and `tale_msa.qmd`:
domain-level objects (`domain_distances`, anything keyed by `dom_code`)
described as repeat-only. FuncTAL's side is genuinely repeat-specific
(termini carry no RVD). The same check found that `fill_type`'s
`"repeat_clust"`/`"repeat_sim"` scored termini too; renamed in §23.

---

## 8. Tests — error conditions now covered **[V]**

`test_error_conditions.R` exists because 11 condition classes had no
test. It found two handlers raising cli's "Cannot pluralize without a
quantity" instead of their own class (fixed with `{cli::qty()}`).

### 8.0 `tell_tales()` has an unguarded filter — FIXED **[V]**

- `min_domain_hits` counts hits **per subject sequence** (a cheap
  pre-filter); `>=`, as documented. When it removes everything the run
  now stops with a classed error naming the argument.
- `min_array_length` (default 0) drops arrays with fewer repeat units
  after grouping; 0 because a short array may be a real truncated TALE.
- All three give-up points of `tell_tales()` are tested
  (`test_tell_tales_guards.R`).

### 8.0b Quoting paths in shell commands — DONE **[V]**

Every shell command built by the package quotes its interpolated paths
(`shQuote()`); paths are built with `file.path()` before quoting. (The
AnnoTALE wrappers were missed and fixed on 2026-09-24.)

### 8.1 A purpose-built fixture for `tell_tales()` -- DONE **[V]**

- The DECIPHER correction dominates runtime and scales with the number of
  arrays times the reference size; `correction_ref_20.fa.gz` (20
  sequences at even intervals) makes the path testable in ~19 s.
- `data-raw/make_toy_tale_regions.R` -> `toy_tal_regions.fasta`
  (`toy_intact`, `toy_frameshift` = the same TALE with one inserted base,
  `toy_no_tale`) with `toy_tal_regions_truth.tsv`.
  `test_tell_tales_correction.R` asserts the known answer: the insertion
  truncates the ORF to 57% coverage, correction restores the same ORF
  length and RVD string with one more insertion charged.
- Incidental, still true: each toy region yields a spurious single-hit
  array at its 3' end (START HERE, worth investigating).

### 8.1b Curating the shipped correction reference -- DONE **[V]**

The original 1057 references were uncurated `tell_tales()` output over 70
*X. oryzae* genomes (555 duplicates, 8 fragments). Shipped:
`tale_correction_ref.fa.gz` (494, the default) and
`tale_correction_ref_representative.fa.gz` (136, pseudogenes kept on
purpose), built by `data-raw/make_correction_references.R`.

**Reference size is not the lever; `max_comparisons` is.** DECIPHER ranks
all references by a cheap distance, truncates to `maxComparisons`, then
aligns: capping at 20 gave byte-identical corrections in 10 s instead of
252 s. The cap fails in the worse direction when the closest references
are not close (a poor correction still looks like an ORF), documented in
`@param max_comparisons` and pinned by a test. `processors` gives ~10%
and is not exposed.

**Deferred to the maintainer:** validate the 136 set against the 494 on
real frameshifted arrays.

### 8.1c `...` must not hide arguments behind an internal **[V]** — audited

Rule in `dev/CLAUDE.md`. The one violation (`tales_align()` forwarding to
`.build_repeat_msa()`) was fixed by promoting `mafft_opts`/`mafft_path`.

### 8.1d Golden baseline records machine-specific paths -- DONE **[V]**

`helper-golden.R` rewrites absolute directories to `<path>/` (keeping the
basename, so a change of reference file still shows) and drops lines
containing **this session's own `tempdir()`**. Generic patterns such as
`/tmp/` or `Rtmp` broke under `covr`, which installs the package under a
tempdir; the pattern was also narrowed once already after "anything
between two slashes" rewrote HMMER's `//` separators and URLs.

### 8.2 Silence MAFFT by default — DONE **[V]**

`mafft_verbose = FALSE` on `tales_align()`: stderr goes to a temp file,
replayed only on failure (`tantale_error_mafft_failed`).

### 8.2b `tales_coded_strings()` needs a `sep` argument -- DONE **[V]**

| | `tales_rvd_strings()` | `tales_coded_strings()` |
|---|---|---|
| `sep` | `"-"` (AnnoTALE) | `" "` (ARLEM, MAFFT) |
| filter | `rvd_only = TRUE` | `repeats_only = FALSE` |

Different defaults on purpose: different consumers. **Naming debt:**
`rvd_only` means "repeats only" and should become `repeats_only` before
1.0.0 (breaking). The fixture fix found here generalises: a fixture that
omits the columns the validators key on does not exercise the validators.

### 8.3 Regression baseline — DONE **[V]**

`test_golden.R` + `helper-golden.R`: snapshots of the column contract, the
anomaly report, the requirements table, the projections,
`tales_compare_distal()`, both alignment layers, grouping, `plot()` data
and consensus, plus two `tell_tales()` runs digested file by file. Tables
are reduced to per-column fingerprints (type, length, distinct count,
missingness, md5 of rounded values). `expect_golden()` forces
`cran = TRUE` so the baseline never silently skips. Use the
`golden-rebaseline` skill to accept a change.

### 8.4 `print()` methods for `tales` and `tales_msa` — DONE **[V]**

`format()`, `print()` and `summary()` for both classes (`R/tales_print.R`,
`R/tales_summary.R`). Print writes to stdout through `cat()` of
`cli::format_inline()` output; `format(x)` equals
`capture.output(print(x))` and no element contains a newline (tested).
Left out on purpose: an RVD frequency table, per-array breakdowns.

## 8.5 Internals audit — function census **[V]**

Superseded by `dev/function-graph.qmd` (§29), which lists every function
with its callers and callees. Rule kept: a single-caller helper earns its
name when a reader needs the name; it then lives next to its caller.

### 8.5b Break `tales_compare()` into three composable steps -- DONE **[V]**

`tales_compare_distal()` = `tales_assign_domain_codes()` ->
`tales_domain_distances()` -> `tales_tale_distances()`. The chain states
the model: two TALEs are compared by aligning their arrays, and the cost
of substituting one domain for another is how different those domains are
as proteins. Step 1 stamps a namespace; step 3 refuses mismatched
namespaces (`tantale_error_namespace_mismatch`). The ARLEM cost matrix
relies on `dom_code`s being exactly `1..n` in row order.

### 8.6 Legacy preconditions leaking through class methods — DONE **[V]**

`plot_tales_msa()` folded into `plot.tales_msa()` and unexported: three of
its four argument checks could not fire through the method.
`plot_tales_composition()` is internal behind `plot.tales()`.

### 8.6b Co-locating single-caller internals — PARTLY DONE **[V]**

Files now follow responsibilities (plotting in `tales_plot.R`, the
`tales_msa` class and MAFFT in `tales_msa_class.R`, ingestion in
`tales_ingest.R`). Remaining row dropped from tracking (maintainer,
2026-09-21).

### 8.6c Plot methods collected into `tales_plot.R` **[V]**

Method-family layout (print, summary and plot each in one file across
both classes), maintainer's choice; `msa.R`'s remainder became
`tales_consensus.R`.

### 8.7 `tales_consensus()` depended on row order **[V]** — FIXED

A tie was won by whichever array came first. Now a tie (or a gap
majority) returns `NA`: no strict majority, no consensus.

## 9. Long-term systematic passes -- all DONE, including 9.2b **[V]**

### 9.0 Governing convention — rOpenSci package API guidelines **[reference]**

<https://devguide.ropensci.org/pkg_building.html#package-api>:
1. `object_verb()` naming for functions sharing a data type; methods on
   generics for verbs applied to a class;
2. data first;
3. snake_case throughout;
4. no name clashes with base or popular packages (re-checked 2026-09-23:
   none of the 51 exports clashes with base, stats, utils, graphics,
   methods, ggplot2, dplyr, tidyr, purrr, magrittr, Biostrings or
   tibble);
5. consistent argument names and order across functions with similar
   inputs.

### 9.1 Argument names — DONE **[V]**

No argument anywhere contains a capital letter or a dot (checked with
`formals()`).

### 9.2 Column names — DONE for the classes, PARTLY for the report files **[V]**

Class and distance-table columns are snake_case (`array_id`,
`domain_type`, `position_in_crd`, `dna_seq`, `source_directory`,
`position_in_array`, `aa_seq`, `rvd`, `seqnames`, `dom_code`; distances
`id1`, `id2`, `dissim`). `.tales_rename_legacy()` still accepts camelCase
input from old output directories; `.pairwise_distances_rename_legacy()`
does the same for distance tables. **Prefer `[[ ]]` to `$` in tests**: a
missing column gives `NULL` and `expect_identical(NULL, NULL)` passes.
One intermediate plot-data column, `rvdSimVsRef` (`tales_plot.R`), is
still camelCase; not part of any returned table.

### 9.2b The sweep stopped at `array_id` in the report files -- DONE **[V]**

`array_report.tsv`: `array_id, seqnames, start, end, strand,
n_domain_hits, array_seq, has_all_domains, predicted_ins_count,
predicted_dels_count, rvd_string, has_aberrant_repeat, nterm_aa_length,
cterm_aa_length, longest_orf_length, orf_coverage, longest_orf_seq`.
The names originate as `mcols()` of one `GRangesList`, so the TSVs and
both GFFs follow from one place. No internal reader needs the old names.

### 9.2c `tell_tales()`'s output file names snake_cased -- DONE **[V]**

`hits_report.tsv`, `domains_report.tsv`, `array_report.tsv`,
`all_ranges.gff`, `rvd_sequences.fas`, `correction_alignment_*`, the
termini alignment pages, etc. All output paths are built in
`.telltale_paths()` (§17 moved the last inline ones there). AnnoTALE's own
output names inside `annotale/` are untouched. The golden history of
`tell_tales.log` has one unexplained digest change at this rename.

### 9.2d Inventory `inst/extdata/`; park what nothing uses in `extra/` -- DONE **[V]**

10 unused files (~20 MB) moved to `extra/` (untracked since §26). A
basename grep gives false positives here: several `inst/extdata/` files
have namesakes under `tests/testthat/data_for_tests/`, some with
different content; check each hit's path. `PXO86.fa` is kept because
articles use it (§25).

### 9.3 Decide `@internal` vs `@noRd` per function — PARTLY DONE **[V]**

Finding: `@keywords internal` without a title generates no Rd at all, so
it behaves exactly like `@noRd`. Only a titled `@keywords internal` block
produces a checked, hidden `man/dot-*.Rd`. Rest dropped from tracking
(maintainer, 2026-09-21).

### 9.4 Use `@family` wherever justified -- DONE **[V]**

Eight families feed `_pkgdown.yml`'s `has_concept()` reference groups.
Blocks with `@rdname` carry no tag (duplicate `\concept{}`).

### 9.5 Unify the user-messaging system -- DONE **[V]**

`logger` removed; cli is the only messaging system; the empty-message
`stop()`s (text sent to the logger, condition empty) are gone.
`logger` had also read and written the user's global logger settings.
Seven sites still carry only the generic `tantale_error` class (START
HERE, worth investigating).

### 9.6 "similarity" and "repeat" are both wrong names -- DONE **[V]**

`pairwise_sim`/`tale_sim`/`repeat_sim` -> `pairwise_distances`/
`tale_distances`/`domain_distances` (28% of domain ids are termini).
`dissim` is the stored column; `sim`/`norm_arlem_score` inputs are
converted on ingest. Entity word singular, head noun plural
(`tale_distances`), since `tales_*` means "operates on a tales object".

### 9.7 Roxygen markdown enabled -- DONE **[V]**

`Roxygen: list(markdown = TRUE)`, verified by diffing `man/` before and
after. Square brackets in prose are links (rule in `dev/CLAUDE.md`).

## 10. Explicitly ruled out

- Deleting dormant internals such as the old HMMER wrappers (now in
  `inst/legacy/unused_pending_review.R`).
- Treating `run_annotale_predict()`/`run_annotale_build()` as dead: no
  internal callers, legitimate standalone utilities.

---

## 11. Rework `tales_group()` — together, and last -- SPLIT INTO TWO **[V]**

Joint session, 2026-09-19. `tales_group()` retired without an alias:
- `tales_group_hclust(x, tale_distances, k, plot_tree = FALSE)`: clusters
  `as.dist()` of the distances (the DisTAL semantics; the Euclidean
  distance between rows that had been switched on was an experiment),
  `ward.D`, cut with `cutree(k = )` (the old height bisection failed on
  tied merges).
- `tales_group_kmedoids(x, tale_distances, k_range, k, seed = 7,
  plot_silhouette = TRUE)`: PAM, a second clustering method in its own
  right (one full clustering per candidate `k`). The stdin prompt for
  `k = NULL` runs only when `interactive()`; `k = "auto"` picks the elbow
  (`.tales_group_kmedoids_elbow()`, unchanged heuristic); `k_range` is
  validated against the number of arrays (§28).
- Plots are drawn only when asked; `seed` is an argument.

Articles rewritten for the split (2026-09-20). Not checked: whether the
two methods agree on the demo data.

---

## 12. One mechanism for running the environment's programs -- DONE **[V]**

Measured: `conda run` costs ~1.3 s per call against ~36 ms by absolute
path, and in a compound command only the first program runs inside the
environment (this machine's `/usr/bin/nhmmer` 3.4 was answering two of
TALEcorrection's three calls). Shipped: `.tantale_bin(tools)` (absolute
paths inside the prefix, naming everything missing at once) and
`.tantale_exec()` (runs and **checks the exit status**; `conda_run2()`
does so only on its micromamba branch). The environment's programs run
correctly from a scrubbed environment (`env -i`). The Java tools and the
nHMMER search of `tell_tales()` were brought under the same check on
2026-09-24 (START HERE item 2).

### 12a The environment-choosing rule was broken **[V]**

`dirname(dirname(conda_binary()))` is the home directory for micromamba
in `~/bin`, so every environment matched and the first one listed won.
`.tantale_pick_env()` picks the candidate holding the tools at the pinned
versions; one wins silently, several warn, none aborts.
`options(tantale.env_prefix = )` overrides.

### 12b `functal()` cannot work, and cannot be fixed with a package -- REBUILT ON `universalmotif`, NOT PATCHED **[V]**

FuncTAL's Perl needs `Bio::Perl`, gone since BioPerl 1.7. Rebuilt in R
(`R/functal.R`):
- `tales_to_universalmotif(x)`: one PWM per array from the repeat RVDs
  (termini excluded), looked up in **`rvd_dna_specificity`** (exported
  dataset, 404 RVDs, from QueTAL's table); aborts on an array with no
  repeats (`tantale_error_no_repeats`).
- `tales_compare_functal(x, method = "PCC", ...)`: `compare_motifs()` on
  those PWMs, returns a `tale_distances` usable by both `tales_group_*()`.
  `min.overlap = 1` and `normalise.scores = TRUE` depart from
  `compare_motifs()`'s defaults (arrays differ widely in length).
- **The scores do not reproduce FuncTAL's**: FuncTAL correlates the whole
  padded, flattened matrix once, `compare_motifs()` combines per-column
  scores (0.360 against 0.489 on the same pair, pinned in a test).
  Maintainer's choice: the standard PCC.
- `functal()` and the vendored Perl moved to `inst/legacy/`.
- Follow-ups (`scan_sequences()`, `merge_motifs()`, significance with a
  TALE-specific null via `make_DBscores()`, keeping the raw `score`,
  stamping the `method` used) dropped from tracking; `motif_tree()`/
  `view_motifs()` are shown in `tale_classification.qmd`.

---

## 13. `stop()`/`message()`/`warning()` converted to cli **[V]**

Re-checked 2026-09-23 with the parser audit in `dev/CLAUDE.md`: **no
`stop()`, `warning()` or `message()` call remains anywhere in `R/`**,
`classification.R` included (the nine sites once left for §11 went with
the rewrite). One `.abort_no_env(what)` replaces four wordings of "could
not create the conda environment".

---

## 14. `DESCRIPTION` Imports audit -- checked, mostly clean; `XVector` dropped **[V]**

- Every package in `Imports` had a live call site, apart from `XVector`
  (used only by the parked `.extract_seqs_from_hits()`, now in
  `inst/legacy/`): dropped. Re-add it if that function is ever revived.
- `reshape2` replaced (18 sites) by three helpers verified with
  `identical()` against the calls they replace: `.pairwise_long_to_matrix()`
  (reproduces `acast()`'s sorted-level order), `.matrix_to_long()`
  (reproduces `melt()` returning numeric margins when all dimnames parse
  as numbers, which ARLEM's `id1 + 1` relies on), `.dcast_count_matrix()`
  (`fun.aggregate` called on empty cells, via `summarise(.drop = FALSE)`).
- `ggcorrplot`/`corrr` dropped from `Suggests` (unused).
- A grep for `pkg::` cannot see dependencies used only through S4
  dispatch or class inheritance; confirm removals with a real `R CMD
  check`.

---

## 15. Cache the repeated discovery/compare/group pipeline across articles **[V]**

`tale_classification.qmd` computes the three-genome discovery,
`tales_compare_distal()` and `tales_group_kmedoids()` and caches them in
`vignettes/articles/_cache/{discovery,compare,group}.rds` (gitignored);
`tale_msa.qmd`, `tales_msa_class.qmd` and `tale_target_prediction.qmd`
`readRDS()` them and abort naming the canonical article if missing.
Checking and computing in every article was rejected as upkeep.

A whole-project `quarto render` (what `build_articles()`/`build_site()`
run) uses an order that follows neither the alphabet nor `_pkgdown.yml`
nor a `_quarto.yaml` `render:` list (all three tested). Hence the
per-article `build_article()` sequence in `dev/CLAUDE.md`. Also here:
`_pkgdown.yml` `template: opengraph: image:` (logo copy in
`man/figures/`), and `dev/dev-notes.Rmd` (copy-paste snippets for dev
sessions).

---

## 16. Preserving the early prototype (`v0.1.9553`) ahead of an eventual repo-bloat cleanup -- DONE, see §26 **[V]**

A branch or tag pointing at old commits keeps them from being collected,
so a real size reduction needs a copy of the history outside the repo
first. Decided: fresh single-commit history, a full `git bundle` kept as
a release asset, one branch `main`. Executed in §26.

---

## 17. In-depth review of every exported function's documentation -- DONE **[V]**

All 51 exports and 19 S3 methods read against their behaviour, one file
at a time; every example extracted with `tools::Rd2ex()` and run
(examples needing external tools run for real). Standards: accuracy
against the code, completeness (every argument, a real `@return`, a run
example; `\dontrun{}` replaced by `try()` for error demonstrations),
biology embedded, no archaeology.

Fixes with consequences beyond wording:
- **`correct_tales()` fed TALEcorrection's nHMMER results on swapped
  flags** (`r=` repeats and `c=` C-terminus reversed). Fixed with the
  maintainer's sign-off; `test_correct_tales.R` moved from 63 to 70
  corrections. Articles updated in §25b.
- `tale_parts_to_rvd()`'s `sep` did nothing (join hardcoded); now used.
- `[.tales` and `[.pairwise_distances` were exported without any docs.
- `tell_tales()`'s termini-alignment file names moved into
  `.telltale_paths()`; `talomes_heatmap()` now returns `invisible(NULL)`.
- `AnnoTALE_QueTAL_functions_library.R` renamed `annotale.R`.

Noted here and closed on 2026-09-24 (START HERE items 2, 4, 6): no test
for `talomes_heatmap()` or `run_annotale_*()`, the analyze stage's exit
status unchecked, `plot_target_preds()`'s ggplot2 deprecations.

---

## 18. `R/unused_pending_review.R`'s last batch reviewed and retired -- DONE **[V]**

`repeat_to_rvd_align()` and `.rvd_to_repeat_align()` had test dependents
in two different test files and moved into `R/conversion.R` (still
internal, used only by tests). The other six (`.write_hmm_file()`,
`.run_hmmer_search()`, `.run_hmmalign()`, `.extract_seqs_from_hits()`,
`.tales_compare_core()`, `.run_in_conda()`) went to
`inst/legacy/unused_pending_review.R`. `.run_in_conda()` had its own test
file, found only by the full suite; removed with it. Both lessons are in
`dev/CLAUDE.md`'s parking rule.

---

## 19. `repeat_sims`/`tal_sim`/`domain_sim`/`fill_type` values -- naming may not have followed §9.6's rename **[superseded]**

The observation: argument names and `fill_type` values kept the
"repeat"/"sim" vocabulary §9.6 had retired, and `"repeat_clust"`/
`"repeat_sim"` scored termini as well. Acted on in §23.

---

## 20. `tale_parts_to_rvd()` -- candidate for a rename and a refactor/rewrite -- retired in §37 **[V]**

Maintainer's flag (2026-09-21). The name carries the legacy "tale_parts"
term, and the function builds RVD strings next to two siblings named
around "map". §17 fixed its `sep` argument and its docs. Reserved for the
maintainer.

---

## 21. `plot.tales_msa()` -- consensus computation moved off the matrix round-trip; internal "repeat_*" naming corrected to "domain_*" -- DONE, partial **[V]**

Audit: `plot.tales_msa()` round-trips its long `tales_msa` through an
array-by-position matrix about nine times; only the two `hclust()` calls
are genuinely matrix-shaped (and they read the distance tables directly).

Done: `.tales_consensus_long()`/`.tales_consensus_match_long()`
(`R/tales_consensus.R`, `@noRd`, documented to public standard) take the
`tales_msa` itself and count implicit gaps as `n_arrays - rows at the
position`; numerically identical to the matrix versions; guarded by
`.assert_tales_msa_layer()`. Internal names changed from "repeat" to
"domain" (`.domain_to_sim_align()`, `.domain_to_cluster_align()`, ...);
a dead consensus computation removed.

**Reserved for the maintainer, do not start unasked** (the tests in
`test_plot_tales_msa.R` and `test_error_conditions.R` pin the current
matrix shapes, and the maintainer wants to rewrite them personally):
1. `.domain_to_cluster_align()`/`.domain_to_sim_align()` as joins against
   the distance table;
2. `.rvd_to_match_align()`, the same;
3. `.pick_ref_name()` as a grouped summary;
4. then `.consensus_panel()` on `.tales_consensus_long()`, after which
   `as.matrix()` leaves `plot.tales_msa()` entirely. The two
   `.pairwise_long_to_matrix()` calls in `tales_plot.R` could become
   `as.matrix()` in the same pass (the two in `distalr.R` need its numeric
   ordering and stay).

**Open:** option (c), a public `tales_consensus()`/
`tales_consensus_match()` taking a `tales_msa` natively
(`dev/class-design.md` §4.6).

## 22. `plot.tales_msa()`'s `position_in_array` mislabel -- fixed **[V]**

The per-position column of the plot's long tables held the alignment
position; renamed `alignment_position`. The visible axis title read
"Position in array", where `plot.tales()` says "Position in alignment"
for the same coordinate. **Maintainer, 2026-10-02:** rename to match;
done the same day (`R/tales_plot.R`, NEWS).

## 23. `repeat_sims`/`tal_sim`/`domain_sim`/`fill_type` renamed -- §19 acted on **[V]**

| old | new |
|---|---|
| `tales_align(repeat_sims =)` | `domain_distances =` |
| `plot.tales_msa(tal_sim =, domain_sim =)` | `tale_distances =`, `domain_distances =` |
| `fill_type = "repeat_clust"` (default), `"repeat_sim"` | `"domain_clust"`, `"domain_sim"` |
| `tales_group_*(tal_sim =)` | `tale_distances =` |

Docs now say every fill mode scores the whole alignment, termini
included. `.tales_group_distmat()`'s parameter is `dists`, to avoid
shadowing the `tale_distances()` constructor. Version 0.9.9004.

---

## 24. Maintainer triage of the open-items list -- decisions recorded, 2026-09-21 **[V]**

Decisions of that triage, all since carried out or recorded in place:
- `tales_consensus_match(long = TRUE)`'s coordinate column renamed
  `alignment_position` (label only; the matrix interface stays).
- §16 prioritised (done in §26); §7.2's fixture paths fixed (see §7.2).
- §2 kept open at lower priority; §5.2 stays parked; §8.7 confirmed
  (ties give `NA`); §9.3, §8.6b and §12b's follow-ups dropped from
  tracking; `mode: release` confirmed as permanent (§15).

---

## 25. Proposed article: how correction handles genuine truncTALEs -- DONE, published via §31 **[V]**

`vignettes/articles/trunctale_correction.qmd`, on PXO86. Biology from Ji
et al. 2016 (doi 10.1038/ncomms13435, "iTALEs") and Read et al. 2016
(doi 10.3389/fpls.2016.01516, "truncTALEs"): TALEs lacking the activation
domain that suppress *Xa1*-mediated resistance; both papers name PXO86.

Two genuine truncTALEs (confirmed by the maintainer), two different
DNA-level events, visible before any correction runs:

| | `ROI_00019` | `ROI_00001` |
|---|---|---|
| event | clean in-frame early stop | genuine frameshift |
| C-terminus nHMMER hit | none | full length, `frameshift_count = 2` |
| `has_all_domains` | FALSE | TRUE |
| N-terminus | 230 aa (283-288 elsewhere) | 230 aa |
| `correct_array = TRUE` | untouched | C-terminus extended by 33 aa |
| `correct_tales()` | untouched | unchanged |

The DECIPHER correction extends `ROI_00001` because the span `tell_tales()`
hands it already contains a C-terminus-shaped template;
`max_comparisons` (20, 50, all 494) makes no difference. Conclusion stated
in the article: `correct_tales()` left both sequences alone,
`correct_array = TRUE` rewrote the frameshift-type one. (The old text of
this section called the full run "1057 references, the default"; the
default is the 494-sequence set, §8.1b.)

**Standing rule from this work:** treat a genome as gold-quality unless it
is flagged (BAI3-1-1 is): a frameshift in it is biology.

**Open:** on a two-array excerpt (`dev/fixtures/pxo86_roi18_19_excerpt.fa`)
`correct_array = TRUE` fabricated a ~489 nt region that does not exist in
the genome; the full-genome run did not. A robustness question for the
DECIPHER wrapper. Excerpts also shift nHMMER e-values (search-space size),
so excerpt-based correction results need confirming on the full genome.

### §25b, `tale_mining.qmd`'s correction chapters are stale after §17's `correct_tales()` fix -- DONE **[V]**

With the flags fixed, `correct_tales()` alone clears both broken BAI3-1-1
arrays; the article was rewritten and re-rendered accordingly. **Flagged,
not decided:** the shared article cache still builds on `bai311_best`
(`max_comparisons = 50`); switching to `correct_tales()`'s result would
ripple through the four cached articles.

---

## 26. §16's history reset, executed -- DONE, with two real findings surfaced along the way **[V]**

2026-09-22. Full-history `git bundle` (327 commits, all branches and tags,
185 MB) attached to the `v0.1.9553` GitHub release and verified by
cloning; orphan commit on `main`; `extra/` untracked and gitignored;
default branch switched before deleting `master`/`dev`/`v0.1.9553`;
`.git` 223 -> 46 MB. The pre-reset working directory is kept on disk as
`tantale-old-before-reset`. The two `devtools::check()` findings went to
§27.

---

## 27. Chasing "a clean, entire test suite" -- all three real bugs fixed and verified **[V]**

- `test_pairwise_distances_class.R` still used `reshape2` (installed in
  the ambient library, undeclared): replaced by a hand-computed matrix and
  `.pairwise_long_to_matrix()`.
- `ncores = 4` in a test: `R CMD check` limits BiocParallel to 2.
- Golden mismatch under `R CMD check` only: `test_path()` is relative
  there and absolute under `load_all()`, and the path normaliser only
  rewrites absolute paths. Fixed with `normalizePath(test_path(...))`.
  The apparent isolated-vs-whole-file difference seen during the
  investigation was produced by the diagnostic itself (output moved off
  `tempdir()` changes which lines are dropped). Methodology in START HERE.

---

## 28. The two §27 follow-up findings, triaged -- both resolved (the ARLEM one by §33) **[V]**

- `tales_group_kmedoids()`'s example tried `k = 4` on 4 arrays; fixed,
  and `k_range` is now validated against `n - 1`
  (`tantale_error_group_kmedoids_krange`).
- The bundled ARLEM binary: Linux x86-64 only, an undeclared executable
  for `R CMD check`, run without an exit check, and its own licence
  string ("Unauthorized commercial usage and distribution of this program
  is prohibited") does not clearly allow redistribution; no public source
  or conda package found. Resolved by §33.

---

## 29. Function dependency diagram, to rebuild the maintainer's mental map -- phase 1 (29.1), phase 2 (29.2) and the matrix/list table (29.3) built **[V]**

Maintainer's request (2026-09-22): a map of the package covering
internals, and the class objects if legible. Decisions: a dev-only `.qmd`
under `dev/`, interactive (`visNetwork`), generated from the code so it
cannot go stale, opening grouped by file. Distinct from
`dev/figures/pipeline.svg` (§7.6), a curated workflow figure.

### 29.1 Phase 1 built, 2026-09-23 **[V]**

`dev/function-graph.qmd`, rendered with `quarto render
dev/function-graph.qmd` (about 5 s) into a self-contained
`dev/function-graph.html` (gitignored).

**Extraction.** Static parse of `R/*.R`, nothing loaded or run. Nodes:
top-level `name <- function(...)` definitions, 187 in the 22 files that
define any (51 exported, 19 S3 methods, 117 internal; kind from
`NAMESPACE`). Edges from `getParseData()` tokens inside each definition:
calls (`SYMBOL_FUNCTION_CALL` or `SPECIAL` for `%||%`, skipping other
packages' `pkg::name`), functions passed by name (excluding local
bindings and `$`/`@` accessors), and S3 dispatch from `as_tales()`, the
only generic the package defines. String literals equal to a function
name were checked: all are class names or message text. Constants
(`TALES_KEY_COLS` etc.) are left out.

**View.** One node per file at first; double-click opens a file or folds
a function's file back; buttons open or fold all; a click highlights
direct callers and callees. The folding is an `htmlwidgets::onRender()`
script clustering on a `file` field: `visClusteringByGroup()` was tried
and dropped, because it clusters on `group` and repaints opened nodes in
vis.js's default palette. vis.js makes no one-node cluster, so
`talecorrection_java.R` shows `correct_tales` itself. The page also has a
function table and the list of internals with no caller in `R/`:
`repeat_to_rvd_align()`, `.rvd_to_repeat_align()` (test-only), `.onAttach()`
(called by R), `.rvds_from_annotale_file()` (called nowhere).

**Blind spots**, listed on the page: methods of base and dplyr generics
show no callers; calls built from strings are not followed.

### 29.2 Phase 2: data-flow view -- built, 2026-09-23 **[V]**

**Maintainer's decision:** a separate, smaller view (exported functions
plus the classes), first on the page; may later be reproduced on the
site's home or getting-started page, so its code is self-contained.

**Method: run-time recording.** Static reading was dropped: 17 functions,
internals included, run the same `is_tales()` check, and returned classes
are often attached indirectly; Rd `\value` sections mention classes that
are not returned. `dev/function-graph-dataflow.R` traces the 70 exported
functions and S3 methods (`trace()` with entry and exit expressions),
runs `testthat::test_dir(load_package = "none")` after
`pkgload::load_all()` (about 4 minutes, 0 failures), and writes
`dev/function-graph-dataflow.tsv` (committed): class of the first
argument on entry (`..1` when the first formal is `...`), supplied
arguments of a tracked class, the return class and the classes of a
returned list's elements. Entry capture matters: `plot.tales()` does
`x <- tales(x)`. Not recorded: the three `dplyr_*` methods (dplyr
dispatches through its own copy) and five exports no test calls (START
HERE item 6); the page lists them.

**View.** 44 nodes: the five tantale classes, the three Biostrings sets
that link functions (`BStringSet` from `tales_rvd_strings()` into
`talvez()`), and the functions taking or returning them. Plain R types
are in tooltips only (a shared "tibble" node would draw false paths). One
hand-declared link, `tell_tales()` -> `tales_from_telltale()` through the
results folder, guarded by `stopifnot()`. Validators, predicates,
`print`/`format`/`[` methods and the `as_tales` methods are left out by a
name pattern. The hierarchical layout drew one long line (cycles such as
`tales` -> `tales_group_hclust()` -> `tales`); force layout used.
Constructors share their class's name, so node ids carry `fn:`/`class:`
prefixes.

**For the site:** the chunk reads only the TSV and `NAMESPACE`; it would
need the TSV where the site build can read it, user-facing prose, and a
check that `visNetwork` works inside pkgdown's Bootstrap 5 pages (DT
tables break there through a jQuery conflict).

---

### 29.3 Matrices and lists of vectors -- built, 2026-09-24 **[V]**

Maintainer asked for a table of every function, internals included, that
takes or returns a matrix or a list of vectors, classified by input/output
and exported/internal, kept updatable. Built the same way as §29.2:
`dev/function-graph-shapes.R` traces all 187 namespace functions during
the test suite (~9 min) and writes `dev/function-graph-shapes.tsv`
(dated, with commit); the new section of `dev/function-graph.qmd` builds
the table from it and lists the functions the tests never call. Columns
at the maintainer's request: `input`, `output` (structure only, "none"
where a side has none of these shapes), `exported?`, `function_name`;
sorted by `output`, "none" last. Argument names, element types and file
are in each row's expandable detail.

Shapes: `matrix<type>`, `list of vectors<type>`, `list of matrices<type>`,
and `list of scalars<type>` (all elements of length 1). The last covers
both one repeat-code string per array (`.build_repeat_msa()`,
`as_tales()`) and records such as `tell_tales()`'s internal `paths` and
`params`; a length-based rule cannot tell them apart, so the page shows
it as its own shape and says so. Only unclassed lists count.

First recording (e30603f): 36 functions (4 exported, 4 S3 methods, 28
internal); 182 of 187 functions called. The 5 never called
(`.format_domain_distances_mat()`, `.rvds_from_annotale_file()`,
`.greet_startup_cli()`, `.abort_no_env()`, `.tantale_repair()`) were
read: none takes or returns one of these shapes.

**Three tracing pitfalls found on the way** (each gave a wrong first
recording):
- `ls(env)` hides names starting with `.`: every internal came out
  "never called". Use `all.names = TRUE`.
- R turns tracing off while a tracer runs. A tracer that forces argument
  promises runs any package code in them untraced, so
  `validate_tales(new_tales(x))` never recorded `new_tales()`. Fixed with
  `tracingState(TRUE)` inside the entry tracer.
- `trace(fn, where = ns)` replaces the namespace binding only. A generic
  from another package (`dplyr_reconstruct()`) dispatches through the S3
  registry, which keeps the untraced method: re-register each traced
  method with `registerS3method()`, looking the generic up in its own
  namespace.

**Done 2026-09-24:** `dev/function-graph-dataflow.R` (§29.2) got the same
two fixes (`tracingState(TRUE)` in the entry tracer, S3 methods
re-registered) and was rerun at 4c1b60b: 69 of 70 exports recorded. The
one missing, `talomes_heatmap()`, is a fourth pitfall: the exit tracer
skipped `NULL` return values, and it returns `invisible(NULL)`. Fixed
afterwards (`returnValue(default = )`, so only an exit by error is
skipped) but not rerun; it will appear at the next recording.

## 30. Full website prose review against `feedback_writing_tone` -- articles/README/index and reference pages DONE; parallel-phrasing sweep deferred **[V]**

2026-09-23. All 8 articles, README, `pkgdown/index.md` and the 63
published reference topics, checked against their rendered output
(`docs/articles/<name>.md`, figures) or the code, with the tone checklist
(memory `feedback_tantale_doc_language`). Many content errors fixed, among
them: `tale_msa.qmd` claimed the scoring matrix gave a more compact
alignment (it displaced a half-repeat, §32.3); `pkgdown/index.md` called a
nine-group dendrogram "three clean groups"; the motif tree figure had no
tip labels; the `pairwise_distances` pages described a similarity with a
`sim` score; README and `?tantale` said tantale bundles no programs.

**Deferred (maintainer):** a second sweep for parallel/paired phrasing
(habit 2 in `dev/CLAUDE.md`).

Open from this pass: rendered error messages of an installed package show
a source path ("at tantale/R/tales_class.R:818:3"), cosmetic. Its other
findings became §32.

---

## 31. §25/§25b's official-publish steps, executed -- DONE **[V]**

The truncTALE article registered in `_pkgdown.yml`, `NEWS.md` entries,
version 0.9.9005, a full site rebuild from an emptied `docs/` (which
removed 228 stale files from the numbered-vignette era that the
per-article loop had kept). `build_reference()` and `build_news()` take
no `quiet` argument. `docs/articles/articles/<name>.html` files are
`build_redirects()` stubs.

---

## 32. Four findings from §30's render check -- all four fixed **[V]**

### 32.1 `as_tales()` keeps `alignment_width` on a demoted object -- **FIXED 2026-09-23** **[V]**

`new_tales()` calls `tibble::as_tibble()`, which drops the `tales_msa`
class but keeps attributes; `new_tales()` now removes `alignment_width`.
`dom_code_namespace` still survives demotion, as intended. Test in
`test_tales_msa_class.R`.

### 32.2 `rvdSimDf` covers 17 RVDs, `rvd_dna_specificity` covers 404 -- **DONE 2026-09-24, option 2** **[V]**

Maintainer: known; unsure what is best; describe the options so a
decision can be made later.

Facts: internal `rvdSimDf` (289 rows = 17 x 17 RVDs, columns
`rvd1`/`rvd2`/`Cor`; Spearman correlation of TALVEZ's `mat1` base
profiles) feeds two things:
- `plot.tales_msa(fill_type = "rvd_sim")` via `.rvd_to_match_align()`:
  any RVD outside the 17 gets `NA` and plots grey, **even when identical
  to the reference** (`NV` in group 6 of `tale_msa.qmd`);
- the MAFFT matrix of `tales_align(domain_distances = "rvd")`, through
  `.rvd_score_table()` (`R/tales_msa_class.R:436`). **Checked
  2026-09-23: this side already falls back**: an unknown pair scores 0
  (neutral) and any RVD against itself scores 1. So the alignment handles
  rare RVDs; only the plot does not.

Exported `rvd_dna_specificity` has 404 RVDs with A/C/G/T preference
counts.

Options, not evaluated:
1. **Recompute the similarity from `rvd_dna_specificity`** (Spearman as
   now, or Pearson) for all RVDs, or on demand for those present. One
   source of truth, full coverage. Rare RVDs have noisy profiles, a
   4-point correlation is crude, and the 17 current values would change
   (existing `rvd_sim` plots and `"rvd"` alignments; golden impact to
   check).
2. **Keep the 17, fall back for the rest**, as `.rvd_score_table()`
   already does: identical RVDs score 1, other unknown pairs stay `NA` in
   the plot. Smallest change; fixes the misleading grey on identical RVDs.
3. **Hybrid**: the 17 from `rvdSimDf`, the rest computed from
   `rvd_dna_specificity` above a minimum count, `NA` below it.
4. **A different measure** (e.g. 1 - Jensen-Shannon divergence between
   normalised profiles), better suited to probability profiles; changes
   existing values as in option 1.

Whatever is chosen, the plot legend and the `rvd_sim` docs should say what
a grey cell means.

**Facts gathered 2026-09-24, for the decision:**
- The 17 are exactly the RVDs of TALVEZ's own `inst/tools/TALVEZ_3.2/mat1`
  (`data-raw/sysdata.R`: Spearman over its rows). So `rvd_sim` uses the
  same RVD model as TALVEZ target prediction.
- `rvd_dna_specificity` is a different table (FuncTAL's `2014mat18`).
  Spearman recomputed from it for the same 17 RVDs differs from
  `rvdSimDf$Cor` by up to 1.07; option 1 would therefore change the
  existing 17 values as well.
- Of its 404 profiles, 78 are flat (1/1/1/1, correlation undefined) and
  81 have a row sum of 4 or less; only 210 profiles are distinct. A
  4-point rank correlation on such counts is mostly noise.
- In the three article genomes (`_cache/discovery.rds`), 510 RVDs, 9
  distinct; the only one outside the 17 is `NV` (5 occurrences).
- Side finding: `rvd_dna_specificity` had
  no `"NA"` row. `readr::read_tsv()` read the RVD name `NA` (Asn-Ala) as
  a missing value (`data-raw/rvd_dna_specificity.R`, default `na =`), so
  the row survives with `rvd = NA_character_` and
  `tales_to_universalmotif()` gives NA repeats the flat fallback profile
  instead of 1/2/1/0.

**Outcome 2026-09-24 (maintainer chose option 2 and the NA fix):**
- `.rvd_to_match_align()` (`R/tales_plot.R`): a cell with no value in
  `rvdSimDf` scores 1 when its RVD equals the reference's; gaps and
  termini (`tales_anchor_codes()`) stay `NA`. `XX` needed no special
  case: `rvdSimDf` already holds `XX`-`XX` = 1. Legend title now reads
  "RVD specificity vs reference (grey: no score)"; `plot.tales_msa()`
  details and `tales_msa_class.qmd` updated (its `NV` column now scores
  1; figure re-rendered and checked). Test in `test_plot_tales_msa.R`,
  failing on the old code.
- NA row: `data-raw/rvd_dna_specificity.R` reads with `na = character()`
  (and now points at `inst/legacy/`, where the source table moved with
  FuncTAL); `.rda` rebuilt, identical except the restored name. Test in
  `test_tales_compare_functal.R`, failing on the old data. The article
  genomes carry no `NA` RVD, so no article output changes. RVDs read at
  run time come from FASTA by string splitting and were never affected.
- The docs of `.rvd_to_match_align()` said it was `rvdSimDf`'s only
  consumer; corrected (`.rvd_score_table()` reads it too).

### 32.3 The `domain_distances` matrix displaced an identical half-repeat -- **FIXED 2026-09-23 (option 1, `penalizeGapLetterMatches = TRUE`)** **[V]**

Symptom: with `domain_distances` as MAFFT's scoring matrix, BAI3's
terminal half-repeat (`dom_code` 29, a 20-aa prefix of the full repeat
31) aligned against MAI1's full repeat 31 instead of MAI1's identical 29.
Cause: `DECIPHER::DistanceMatrix()` defaults to
`penalizeGapLetterMatches = FALSE`, so an overhang counted as nothing and
`dissim(29, 31) = 0`.

DisTAL's definition (Pérez-Quintero et al. 2015, doi
10.3389/fpls.2015.00545): global alignment with free end gaps, distance =
share of residues that differ, based on the longer repeat. Half vs full
repeat = 14/34 = 41.2. Now followed by all three backends:
- DECIPHER: `penalizeGapLetterMatches = TRUE` (version 0.9.9007);
- Biostrings: `type = "overlap"` (free end gaps; internal gaps keep their
  cost; 0.9.9008);
- mmseq2 already agreed.

About 7300 of 12769 pairs still differ by a median of ~1 point between
DECIPHER (distances inside one multiple alignment) and Biostrings
(pairwise), which is inherent. Grouping unchanged (k = 9, same
partition). Test in `test_distalPairwiseAlign.R` pins 29 vs 31 for all
backends; golden re-baselined with every row explained. After the fix
the scored and unscored alignments of group 6 are identical.

### 32.4 C-terminus length differs by one between `tales` and `array_report.tsv` -- **FIXED 2026-09-23 (stop codon no longer counted)** **[V]**

AnnoTALE's C-terminal part includes the stop codon as `*` whenever the
CDS ends inside it (truncTALEs, and the 278-aa variant of normal BAI3
TALEs). `.aa_residue_count()` now counts residues excluding `*` for both
terminus lengths, the rule `tales_ingest.R` applies to `aa_seq`. Version
0.9.9009. Tests in `test_tell_tales.R` (PXO86 excerpt, now also
`data_for_tests/pxo86_roi18_19_excerpt.fa`).

---

## 33. ARLEM re-implemented in R, wired in; the executable removed -- DONE **[V]**

2026-09-23. Model from Abouelhoda, Giegerich, Behzadi & Steyaert (APBC
2008, "Alignment of minisatellite maps: a minimum spanning tree-based
approach"): per-interval duplication histories grown from the leftmost
or rightmost unit, then an alignment DP with simultaneous right growth
(O(n^3) via the paper's A' table), with a `$` sentinel as the binary had.
Two binary behaviours not in the paper, found by probing: only `# Indel
hist` of the cost file affects scores; the leading unit is explained as
an insertion.

- `R/arlem.R`: `.arlem_histories()`, `.arlem_align()`, `.arlem_scores_r()`,
  **`identical()` to the binary** on ~6060 random pairs and on all real
  data tried. `tales_tale_distances()` uses it; `.run_arlem()` and the
  cost-file writer are in `inst/legacy/arlem_binary.R`;
  `inst/tools/arlem/` deleted outright (keeping the binary anywhere would
  still redistribute it).
- The binary's answers are kept as
  `tests/testthat/data_for_tests/arlem_reference_scores.rds` (1126 pairs,
  `data-raw/make_arlem_reference_scores.R`), checked by `test_arlem_r.R`.
- `matrixStats` added to `Imports` (`colMins()` halves the run time).
  Speed: ~4.5x the binary (26 arrays: 0.70 s against 0.16 s).
- Kept: `arlem_score`/`norm_arlem_score` column names, `.arlem_*`
  internal names, one citation in `tales_tale_distances()`, past
  `NEWS.md` entries. Version 0.9.9006.

**Deferred:** an Rcpp version, likely faster than the binary.

## 34. Distribution strategy -- findings recorded; one-archive plan and rOpenSci both parked **[P]**

*2026-09-23.* Discussion with the maintainer; nothing decided on the
channel. This section records what was checked, so it need not be redone.
Current state: the tools-archive plan at the end of this section is
written up and parked by the maintainer; rOpenSci is parked until the
package matures; the channel itself is open.

### What limits the choice

**Size.** The source package, as `git ls-files` minus `docs/`, `dev/`,
`pkgdown/`, compresses to 57.8 MB. CRAN and Bioconductor both cap the
source tarball at 5 MB.
- `inst/tools`: 62 MB, of which 58 MB are three jars (TALEcorrection
  27 MB, AnnoTALEcli-1.5 16 MB, PrediTALE 15 MB). Jars do not compress.
- `inst/extdata`: 19 MB, almost all of it four whole genomes (`BAI3.fa`,
  `BAI3-1-1.fa`, `MAI1.fa`, `PXO86.fa`).
- `inst/legacy`: 4 MB.

Without the jars, TALVEZ, `talecorrect/`, the four genomes and
`inst/legacy`, the same estimate gives **1.5 MB**.

**Licences of the bundled tools.**
- All three jars carry `COPYING.txt`, GPL-3 (Jstacs). Redistribution is
  allowed with a pointer to the source (github.com/Jstacs/Jstacs), which
  nothing in the package gives. `LICENSE`/`DESCRIPTION` mention none of
  them, so the package reads as entirely MIT.
- TALVEZ 3.2 (A. Pérez-Quintero, IRD): no licence in the script, the
  zip or the web page. By default that grants no redistribution right,
  the same position ARLEM was in (§28). Its Java part ships as `.class`
  files only. The maintainer is asking the author for permission (see
  below).

**Upstream sources, verified 2026-09-23 by download and `sha256sum`.**
Every bundled file is byte-identical to its upstream copy.

| tool | upstream URL | sha256 (first 12) |
|---|---|---|
| AnnoTALEcli 1.5 | `https://www.jstacs.de/downloads/AnnoTALEcli-1.5.jar` (older versions also kept online) | `fe99d0840733` |
| PrediTALE | `https://www.jstacs.de/downloads/PrediTALE.jar` (unversioned URL; Last-Modified 2019-01-16) | `67a3ef81c2ba` |
| TALEcorrection | only inside `https://www.jstacs.de/downloads/TALECorrection_scripts.zip`, 283 MB (a 227 MB test BAM) | `9adf9de20a41` |
| TALEcorrection HMMs (Xoo, Xoc, custom fasta) | same zip, `HMMs/` | all identical |
| TALVEZ 3.2 | `https://bioinfo-web.mpl.ird.fr/xantho/talvez/downloads/TALVEZ_3.2.zip` (http times out, https works; Last-Modified 2016) | zip `5661ca5825e7`; the 19 bundled files identical, the zip adds a `tmp/` of example output |

The Jstacs GitHub releases (v2.3b to v2.4.1) carry no assets. None of
AnnoTALE, PrediTALE, Jstacs or TALVEZ is on bioconda or conda-forge.

### Channels as discussed

- **GitHub + r-universe** (`scunnac.r-universe.dev`): no review, no size
  cap known (to check), binaries for Linux/macOS, ordinary
  `install.packages()`. Zenodo gives each release a DOI.
- **Bioconductor**: the best audience fit (~15 Bioconductor imports,
  `biocViews` set). Needs the 5 MB tarball, `BiocCheck`, likely
  Bioconductor classes at the interfaces, the twice-yearly cycle, and a
  declared Windows exception (`OS_type: unix`).
- **CRAN**: poor fit. Every external-tool test would have to skip there,
  against the fail-don't-skip rule; jars need their sources in `java/`.

### rOpenSci -- parked for later, maintainer finds it tempting

Peer review of the code (open GitHub issue, editor + two reviewers)
against the rOpenSci dev guide, which tantale already follows for naming
(§9.0). On acceptance: optional transfer to the `ropensci` org,
`docs.ropensci.org/tantale`, a badge, promotion, a fast-tracked JOSS
review with a short paper, and distribution via `ropensci.r-universe.dev`.
Cross-listing on CRAN or Bioconductor stays possible.

Scope is the open question. The wrappers (AnnoTALE, PrediTALE,
TALEcorrection, TALVEZ, MAFFT, HMMER, with parsing and a managed
install) fit "scientific software wrappers"; the guide counts an
improved installation as added value. Data visualisation and
statistical/modelling libraries are out of scope, which touches
`plot.tales_msa()`, the DisTAL distances and the classification. The
licence of wrapped tools is judged case by case, so TALVEZ would come up.
A pre-submission inquiry (a short issue) settles scope before any review.
**Maintainer, 2026-09-23: revisit once the package has matured further.**

### First proposal: download each tool from its upstream **[superseded]**

Replaced the same day by the one-archive plan below. `tantale_setup()`
would have fetched each tool from its authors' page against a pinned
sha256. TALEcorrection has no standalone download (only the 283 MB zip),
TALVEZ's only upstream is a 2016 IRD server, and the maintainer objected
to depending on several third-party servers.

### Side finding: absolute path in `test_correct_tales.R` -- FIXED 2026-09-23 **[V]**

`tests/testthat/test_correct_tales.R` lines 4 and 17 read `BAI3-1-1.fa`
by an absolute path under `/home/cunnac/...`, so the test only ran on
this machine. The second test also wrote its output to
`tempfile(tmpdir = "~")`, into the user's home directory. Both now use
`system.file("extdata", "BAI3-1-1.fa", ...)` and a plain `tempfile()`.
The file passes under `load_all()` (real nhmmer + TALEcorrection run,
1 min 13 s), and nothing is left in `~`. The plan below moves the
genomes out of `inst/extdata`, so this test will need its fixture again
then.

Sweep for the same kind of bug, same day:
- `R/`, `tests/testthat/*.R`, the articles, `getting_started.qmd`,
  README, `pkgdown/index.md`, `_pkgdown.yml`, `DESCRIPTION`: no other
  absolute path in executed code or in `@examples`. No other test writes
  to `~`, calls `setwd()`, or writes into the working directory.
- Commented-out developer snippets still carry `/home/cunnac/...`
  paths: `R/target_predictions.R` (102-104, 249-253), `R/telltale.R`
  (22-23), `R/distalr.R` (620, 624), `R/tantale_conda_env.R` (44-45).
  Never executed; left alone, a cleanup candidate for the maintainer.
- Fixtures echo old absolute paths as provenance only:
  `tell_tales.log`, AnnoTALE's `protocol_analyze.txt`, HMMER output,
  and the `source_directory` column of `sampleDistalrOutput.rds` and
  `sampleTalesMsa.rds`. Nothing reads a path back out of them (checked
  with grep over `R/` and `tests/`); the golden test already masks them.

### Maintainer, 2026-09-23 (later): one archive on GitHub, TALVEZ settled

- **TALVEZ:** the maintainer will ask A. Pérez-Quintero for permission
  to redistribute; treat it as settled.
- **Concern with the first proposal:** it depends on several third-party
  servers (jstacs.de, the 2016 IRD server) staying up and keeping the
  same files. **Preferred instead:** bundle everything into one archive
  attached to a GitHub release of `scunnac/tantale`. The only server
  involved is then github.com.
- GitHub limits, checked in its docs: each asset under 2 GiB, up to 1000
  assets per release, "no limit on the total size of a release, nor
  bandwidth usage". An *immutable* release (repository setting) locks
  its assets and tag after publication.
- The plan that came out of this follows. Its last paragraph lists what
  is still to decide.

### The one-archive plan, as proposed 2026-09-23 -- parked by the maintainer

*Copied from the discussion at the maintainer's request, as a record.
Not to be acted on until the maintainer says so.*

**The release.** A dedicated GitHub release on `scunnac/tantale` with a
tag of its own, for example `tools-1`, kept separate from package
versions. The archive changes only when a tool does, while the package
changes much more often. The package pins the archive's URL and SHA-256.
A new PrediTALE, say, would mean a `tools-2` release and a new pin in the
package.

GitHub's own documentation sets these limits:
- each file must be under 2 GiB;
- "there is no limit on the total size of a release, nor bandwidth
  usage";
- a repository setting, *immutable releases*, locks a release's assets
  and tag once it is published, so the file behind the URL can never
  change.

**What the archive contains**
- the three jars, the Xoo and Xoc HMMs, and TALVEZ 3.2;
- `README`: for each tool, its upstream URL, upstream SHA-256, version,
  licence, citation and source link;
- `LICENSES/`: the GPL-3 text for the Jstacs tools and Alvaro
  Pérez-Quintero's permission for TALVEZ. GPL-3 allows redistributing
  the jars as long as the licence travels with them and the notice says
  where the source is;
- `MANIFEST`: a SHA-256 for every file, so an unpacked copy can be
  checked later.

**How it gets built.** A script in `dev/` assembles the archive, from
upstream or from the current `inst/` copies (shown identical to upstream
above), writes the README and prints the archive's SHA-256. The
maintainer uploads the archive through the GitHub web page (the `gh`
command-line tool is not installed on this machine).

**Package side.**
- `tantale_setup()` fetches the archive once, checks its SHA-256, unpacks
  it into `tools::R_user_dir("tantale", "data")/tools-1/` (`"data"`
  rather than `"cache"`, which users and operating systems clean), and
  checks the per-file hashes.
- `tantale_setup(tools_from = "path/to/archive.tar.gz")` installs from a
  local copy, for machines without internet or a cluster where one
  person downloads for everyone.
- One internal resolver (`.tantale_tool("annotale")`) returns each
  tool's path, or aborts with `tantale_error_tool_missing` pointing to
  `tantale_setup()`. The defaults of `run_annotale_predict()`,
  `run_annotale_build()`, the internal `.run_annotale_analyze()` (called
  by `tell_tales()`), `preditale()`, `talvez()` and `correct_tales()`
  (including `hmm_path`) move from `system.file(..., mustWork = TRUE)`
  to that resolver. The arguments stay, so a user can still point to
  their own copy.
- The files leave `inst/`. `inst/tools/talecorrect/`'s upstream `.java`
  sources and shell scripts, which tantale does not call, move to
  `inst/legacy/`. `inst/legacy` goes into `.Rbuildignore`.

**Sizes, compressed:** tools 50 MB, the four genomes 5.6 MB. Suggested:
the same release, two assets, since only the articles and a few
examples need the genomes.

**Still depends on a server:** github.com, the one download location,
and conda-forge/bioconda for the conda tools, unchanged. Removing the
jars from `inst/` does not shrink `.git` (~60 MB of history).

**Genomes.** Used by the articles, the `@examples` of `annotale.R`, and
`test_correct_tales.R`. Articles and examples would fetch them on first
use; tests would use a subset cut around the TALE loci, like
`bai3_sample_tal_genomic_regions.fasta`. PXO86 is `NZ_CP007166`.

**Open decisions when this is picked up:** one asset or two; whether to
deposit the same archive on Zenodo as a fallback URL with a DOI; the
genome accessions for BAI3 and MAI1, and whether BAI3-1-1 is public.

### Code hosting on the institutional GitLab (`forge.ird.fr`) -- noted, not pursued

*2026-09-23.* The maintainer asked whether moving the repository to IRD's
GitLab would hurt the options above. Not pursued for now; the forge
answered HTTP 502 that day.

- **GitLab only:** r-universe still works (its docs: packages "can be
  hosted on any public Git server"; only the `packages.json` registry
  must be on GitHub). rOpenSci's process is built around GitHub, so a
  GitLab-hosted package would be a question for the pre-submission
  inquiry. Bioconductor's submission guide assumes a GitHub repository.
  Zenodo's automatic release archiving is GitHub-only.
- **What a move entails:** new `URL`/`BugReports` in `DESCRIPTION`,
  `url:` in `_pkgdown.yml`, badges and install instructions; the site on
  GitLab Pages via CI (check it is reachable from outside IRD); GitHub
  Actions rewritten as GitLab CI (check shared runners can build the
  conda environment and run Java); the GitHub repository archived with a
  pointer, and the `v0.1.9553` history bundle copied or left there.
- **Main practical risk:** institutional instances often create accounts
  for staff only, so outside users could not open issues.
- **Suggested setup if it is ever wanted:** GitLab as the main
  repository with push mirroring to GitHub (free tier, unless the admins
  disabled it). That keeps r-universe, Zenodo DOIs, a route to rOpenSci
  or Bioconductor, and public issues.
- **Tools archive:** keep it on GitHub releases (plus Zenodo) whatever
  the code host. An institutional server is the kind of dependency the
  one-archive plan was chosen to avoid.

## 35. Catch-up review of 42deefd..b54f277 and the termini audit (`dev/notes_for_claude.md`) -- plan A+B done **[V]**

Reviewed 2026-09-30 (the maintainer's four commits since `claude-reviewed`).
`dev/notes_for_claude.md` is the maintainer's road map for this phase:
retirements, renames, `tales_from_annotale()`, the size proposal and the
`tell_tales()`/`tales_from_telltale()` termini audit.

**Full `devtools::test()` at b54f277: 2 errors**, both from the
work-in-progress tightening of `.tale_parts()` (`R/tales_ingest.R`):
- `test_tale_parts.R:7` (`err_array_count` fixture): `unmatched = "error"`
  in the RVD join now aborts in dplyr where the test expects a warning.
- `test_tell_tales.R:31` (PXO86 excerpt, real run): the re-enabled
  `stopifnot(nrow(taleProtString) == nrow(taleDnaString))` fails.

Other points in the same diff: the missing-terminus warning says "missing a
N-terminus and C-terminus domains" when both are missing (cli takes the
quantity from `{array_id}`); the array-length warning carries the error
class `tantale_error_parts_inconsistent`; `tales_anchor_codes()` now
returns a named vector (`N-`, `-C`, `??`), which the tests accept.

**The PXO86 excerpt shows both directions of the termini disagreement**
(`tests/testthat/data_for_tests/pxo86_roi18_19_excerpt.fa`, default
parameters):
- `ROI_00001` (a partial array cut by the excerpt's start): nhmmer hits are
  one N-terminus and two repeats, no C-terminus. AnnoTALE writes a DNA
  parts file with `N-terminus`, `repeat 1`, `C-terminus`, no protein parts
  file and no RVD entry. So AnnoTALE names a C-terminus that no nhmmer hit
  supports, and the protein and DNA files disagree in row count.
- `ROI_00003` (the genuine truncTALE, PXO86 ROI_00019): AnnoTALE reports a
  42-aa C-terminus; nhmmer has no C-terminus hit above `cterm_min_score`,
  so the RVD string ends in `XXXXX`. Here the terminus is real and the
  hit is missing, so "no supporting hit" cannot by itself mean "not a
  terminus" for short, truncated termini.

**Mechanism confirmed in the code.** `.telltale_finish_rvd_strings()`
writes `NTERM`/`CTERM` when a hit of that type exists anywhere in the
array's hits; `.tale_parts()` then attaches that code to whatever
AnnoTALE called the terminus of the longest ORF. Nothing compares the two
positions.

**Maintainer's decisions, 2026-10-01.** (Q1) A terminus that fails the
check is recoded `XXXXX`. The criterion is biological: whatever AnnoTALE
extracts on either side of the repeat region of a predicted ORF is a real
terminal segment by definition (the ORF is assumed translated); the open
question is whether that segment resembles a canonical TALE N- or
C-terminal domain. Test that with the *protein* profiles already shipped
in `inst/extdata/hmmProfile/` (`Xo_TALE_Nterm_AA_profile.hmm`, 288
positions; `Xo_TALE_Cterm_AA_profile.hmm`, 279; built 2020 from 357/359
X. oryzae sequences; unused by any code so far). The position overlap
with nhmmer DNA hits proposed first is dropped. (Q2) Finish the
`.tale_parts()` tightening, while the maintainer keeps thinking about it.
(Q3) Order: termini check and `array_report.tsv` renames, then the two
stray `tell_tales()` warnings, then `tales_from_annotale()`, renames and
retirements. **No coding until the plan is agreed.**

**Feasibility run, 2026-10-01** (scratchpad script, `hmmsearch` 3.3.2
from the env, `-T 0 --domT 0`, one row per terminus, best domain):
- MAI1 (clean, 10 arrays): every AnnoTALE terminus hits its profile,
  E < 1e-163, 535-632 bits, full profile length.
- PXO86 excerpt `ROI_00003` (genuine truncTALE): the 42-aa C-terminus hits
  at E = 6.6e-19 (55 bits, profile 1-37). The nhmmer DNA search missed it,
  so the current code writes `XXXXX`; the protein check would write `CTERM`.
- BAI3-1-1 raw (`cterm_min_score = 300`): six termini now coded
  `NTERM`/`CTERM` have **no hit at all**, even at bit threshold 0:
  C-termini of `ROI_00001` (226 aa), `ROI_00007` (137 aa), `ROI_00008`
  (26 aa); N-termini of `ROI_00006` and `ROI_00009` (24 aa each). The
  N-terminus of `ROI_00001` (247 aa) hits only profile 103-150 (68 bits,
  E = 4e-22). Genuine short fragments (42-45 aa) score E < 1e-18, so an
  E-value cutoff around 1e-5 separates the two groups on this data.

**First proposal, superseded by the maintainer's decision above:** place each AnnoTALE terminus on the
genome (the ORF's offset within the array region plus the part's offset
within the ORF, or a `matchPattern()` of the part's DNA in the array
sequence) and test overlap with the nhmmer hit of the same type. Record
the outcome per terminus as a column of the `tales` object and of
`array_report.tsv`, which also answers the note's request for
per-terminus nhmmer columns in place of `rvd_string`. Open question for
the maintainer: whether an unsupported terminus stays in the object
flagged, or is recoded.

**Plan review, 2026-10-01 (maintainer's answers; still no coding).**
- A2 corrected: `XXXXX` marks an AnnoTALE terminus with no protein-profile
  hit, and nothing else. When AnnoTALE reports no terminus there is
  probably nothing on that side of the repeats: `.tale_parts_from_file()`
  drops that row (today it adds an all-`NA` row).
  `position_in_array` must then be renumbered per array (repeats are
  numbered `position_in_crd + 1` on the assumption of an N-terminus row).
  In the three runs examined, AnnoTALE reported both termini for every
  array it analysed; the missing-terminus case exists only in the
  synthetic fixture `err_missing_nterm`.
- A3, A4, A7, B2, Q4 (E <= 1e-5), Q7 (`min_domain_hits` -> `min_dna_hits`)
  agreed. Q5: `rvd_string` stays; after plan A its terminus codes describe
  protein-profile matches, and the `array_report.tsv` columns get
  documented one by one.
- A5 (anomaly for `XXXXX`) postponed: `tales_anomalies()` and the checks
  in `.tale_parts()` are to be re-examined one by one.
- **B3: AnnoTALE's RVD count equals its repeat-part count** in all 21
  arrays examined (MAI1, PXO86 excerpt, BAI3-1-1 raw). With both read from
  AnnoTALE (A4), the separate length check is redundant with a join that
  errors on unmatched rows on either side.
- **B1: AnnoTALE `analyze` behaviour understood** (PXO86 excerpt
  `ROI_00001`, `protocol_analyze.txt`). It splits the ORF DNA into parts,
  then translates it. The split produced a 9-nt "repeat 1"; translation
  failed ("Both RVD positions are gaps"), then "Protein version ... could
  not be analyzed, splitting into regions failed". It wrote
  `TALE_DNA_parts.fasta` anyway, and neither `TALE_Protein_parts.fasta`
  nor `TALE_RVDs.fasta`. `.telltale_run_annotale()` warns and skips the
  array, but leaves the DNA parts file, which `.tale_parts()` then reads.
  Different case, BAI3-1-1 `ROI_00003`/`ROI_00005`: protein and DNA parts
  hold an N-terminus (285 aa) and a 3-aa "C-terminus", no repeat, no RVD
  file.
- Q8: `tales_get_protein_seq()`/`tales_get_dna_seq()` already abort with
  `tantale_error_projection_na` on an `NA` sequence, and
  `tales(x, sanitize = TRUE)` removes anomalous arrays.
- Agreed 2026-10-01: A1 in `tell_tales()`; A2 drop the row and recompute
  positions; A6 wording as drafted; Q6 abort on an old directory with a
  message to rerun `tell_tales()`; Q9 (a) `tell_tales()` removes the DNA
  parts file of an array it skips for lack of protein parts. B1: parts
  lacking either the DNA or the protein sequence are filtered out, at the
  level of the whole array (dropping a single part would leave a hole in
  the array and shift the positions). The consolidated plan A+B is in the
  reply of the same day; awaiting final approval.

**Plan A+B executed 2026-10-01 (0.9.9011).** Coding approved by the
maintainer, with one more instruction: a helper called from one place goes
back into its caller (`.telltale_finish_rvd_strings()` inlined into
`tell_tales()`).
- `R/telltale.R`: `.tale_termini_hmmsearch()` (kept as a function for the
  future `tales_from_annotale()`), called in `tell_tales()` on the termini
  already read by `.telltale_align_termini(type = "AA")`, now run before
  the RVD strings. New argument `terminus_max_evalue = 1e-5`;
  `min_domain_hits` -> `min_dna_hits`; `array_report.tsv`: `n_dna_hits`,
  `nterm_dna_hit`, `cterm_dna_hit`, `nterm_aa_evalue`, `cterm_aa_evalue`,
  `nterm_aa_hit`, `cterm_aa_hit`. `.telltale_run_annotale()` deletes
  `TALE_DNA_parts.fasta` along with an absent/empty protein parts file.
  Run log: counts of arrays whose termini match the protein profiles.
  `?tell_tales` documents every `array_report.tsv` column and what
  `hmm_dir` must hold.
- `R/tales_ingest.R`: `.tale_parts()` rewritten as planned (B1-B5):
  RVDs from `TALE_RVDs.fasta`, codes from `array_report.tsv`
  (`tantale_error_telltale_outdated` without the columns), whole arrays
  dropped on a protein/DNA disagreement (`tantale_warning_parts_inconsistent`),
  absent termini warned about (`tantale_warning_terminus_absent`) and
  positions counted per array, RVD/repeat mismatch an error
  (`tantale_error_parts_inconsistent`). `rvd_sequences.fas` is no longer
  read by the package. `.tale_parts_from_file()` no longer adds rows.
- Fixtures: `data-raw/make_telltale_test_fixtures.R` regenerates
  `inst/extdata/tellTaleExampleOutput` and `example_output` from
  `bai3_sample_tal_genomic_regions.fasta` (RVD strings byte-identical to
  the 2022 run) and derives the three `err_*` directories from it, keeping
  only what `.tale_parts()` reads. `err_array_count` now plants an RVD
  missing from `TALE_RVDs.fasta`. The old fixtures used the pre-rename
  column names (`SeqOfRVD`, `AllDomains`...); none remain in tests or
  examples, which bears on the maintainer's question whether
  `.tales_rename_legacy()`/`.pairwise_distances_rename_legacy()` are still
  needed (not examined further). New fixture
  `termini_profile_cases.fa` (four real segments from MAI1, PXO86,
  BAI3-1-1).
- Seen on BAI3-1-1 raw: `ROI_00001`'s 247-aa N-terminus matches only
  profile positions 103-150 (E = 4e-22) and is coded `NTERM`. The check
  has no coverage criterion; a partial match counts as a terminus.
- Left for the maintainer: `tale_mining.qmd` and
  `trunctale_correction.qmd` select `has_all_domains`, and will fail to
  render until updated; the articles' cache (`discovery.rds` included)
  is stale.
- Checks: full `devtools::test()` (only the golden tests changed), then
  quick `R CMD check` (no tests, no vignettes): 0 errors, 0 warnings, 0
  notes. **Golden re-baselined**, every row explained: `all_ranges.gff`
  (array attributes renamed/added), `array_report.tsv` (renamed and new
  columns; all old columns keep their digests), `tell_tales.log` (+3
  lines: `terminus_max_evalue` and two protein-profile counts), in both
  the plain and the corrected run. The first run of the golden test also
  caught a bug of the inlined code, `has_aberrant_repeat` all `NA`
  (`nzchar()` drops names), fixed before accepting. Weakness noticed: the
  fingerprint rounds doubles to 8 decimals, so the E-value columns (around
  1e-190) all digest as 0.

## 36. The two stray `tell_tales()` warnings (road map step 2) **[V]**

Investigated 2026-10-01; plan awaiting the maintainer's approval.

**36.1 "invalid seqlevels 'seq2' ignored"** (`GenomeInfoDb::renameSeqlevels()`
in `.telltale_hits_to_ranges()`). `.telltale_prepare_subject()` renames
every subject sequence `seq1..seqN` and returns the whole renaming
vector. The hit ranges only carry the sequences with hits (the hit table
was `droplevels()`'d after the `min_dna_hits` filter), and
`renameSeqlevels()` warns about every name it is given that the ranges do
not have. Harmless: BAI3-1-1's second sequence carries no hit. Fix:
pass only the names present in the ranges.

**36.2 "Some HMMER hits overlap, so the inferred RVD sequences may carry
artefactual insertions"** fires on nearly every array (PXO86 17/18, BAI3
9/10, MAI1 9/10, BAI3-1-1 raw 8/9). Measured on the merged hits of those
four genomes: 51 overlapping pairs, **all between a terminus hit and the
adjacent repeat hit, none between two repeats**. C-terminus/repeat: 43
pairs, 16-20 nt; N-terminus/repeat: 8 pairs, 1-4 nt. The merge step
(`.telltale_merge_overlapping_hits()`) merges overlaps within a domain
type only, by design (its own doc calls a terminus/repeat overlap "a real
feature of where one domain ends and the next begins"), so the check in
`.telltale_group_arrays()` contradicts it. The message's rationale is also
obsolete: RVDs are read by AnnoTALE from the ORF; the nhmmer hits only
delimit the array and feed `n_dna_hits`, the `*_dna_hit` columns and
`hits_report.tsv`. The merge step does act on same-type duplicates
(6 merged hits in PXO86, 1 in BAI3 and in MAI1).

**36.3 Side finding:** `.telltale_prepare_subject()` calls
`Rsamtools::indexFa(subject_file)`, which writes `<subject>.fai` next to
the user's input file, and fails where that directory is read-only. This
is where the `.fai` files in `inst/extdata/` and `data_for_tests/` come
from.

**Plan proposed 2026-10-01:** P1 pass `renameSeqlevels()` only the names
present (test: `toy_tal_regions.fasta`, which has a sequence without a
TALE). P2 check overlaps between hits of the same type only (possible
with `merge_hits = FALSE`), new message naming `n_dna_hits` and
`merge_hits`; correct the obsolete RVD rationale in the merge/group docs;
one sentence in `?tell_tales` on normal terminus/repeat overlaps; tests
with and without `merge_hits`. P3 build the seqinfo from the sequences in
memory, no `indexFa()`. **Open:** Q10, terminus/repeat overlaps never
reported (recommended) or reported above a threshold (e.g. > 30 nt);
Q11, include P3 and remove the 7 tracked `.fai` files.

**Maintainer, 2026-10-01:** P1 and P2 agreed. Q10: terminus/repeat
overlaps are never reported. P3 questioned ("indexing is the standard
way"), then agreed after these arguments (Q11 follows P3: the 7 `.fai`
files are removed):
- the `.fai` is written next to the input, which CRAN's policy forbids
  outside `tempdir()` in examples, tests and vignettes. Examples on
  `system.file()` data write into the installed library, which is how
  the 7 tracked `.fai` files arose. On a read-only directory the call
  fails;
- the index serves only to get the sequence lengths, and the sequences
  are already in memory one line earlier (`Biostrings::readDNAStringSet()`);
- **latent bug found while checking:** the `.fai` names a sequence by the
  first word of its header, `readDNAStringSet()` by the whole header. For
  a header with a space (`>ctg1 plasmid pXO1`), `seqinfo[...]` in
  `.telltale_hits_to_ranges()` finds no match and the ranges get
  `seqlengths` `NA`, silently. This is the very case the renaming step
  exists for. Building the seqinfo from `names()`/`width()` of the
  sequences in memory fixes both.
  Consequence (read from the code, not run): `tell_tales()` extends each
  array by `extend_len` and relies on `GenomicRanges::trim()` to clip it
  at the contig end; with `seqlengths` `NA` it cannot, so an array within
  `extend_len` of the end of such a contig would be extended past it
  before `getSeq()`.

**Done 2026-10-01 (P1-P3):**
- P1: `.telltale_hits_to_ranges()` passes `renameSeqlevels()` only the
  names present in the ranges.
- P2: `.telltale_group_arrays()` checks `isDisjoint()` per domain type
  within each array; new message names `n_dna_hits` and `merge_hits`.
  The obsolete RVD rationale is gone from the merge/group docs;
  `?tell_tales` says, under `n_dna_hits`, that terminus/repeat overlaps
  are normal; `merge_hits` documented. On the PXO86 excerpt, the default
  run gives no overlap warning, `merge_hits = FALSE` flags ROI_00003.
- P3: seqinfo from `names()`/`width()` of the sequences in memory;
  Rsamtools dropped from Imports (it had no other use). The 7 `.fai`
  files removed; test runs leave none behind.
- New guard: a subject with duplicated or empty sequence names is an
  error (`tantale_error_seqnames`). `Seqinfo()` rejects both; the old
  index accepted them, but a hit on an unnamed sequence already failed
  in `renameSeqlevels()`. Two tests wrote an unnamed random sequence;
  they now name it.
- Tests in `test_tell_tales.R`: subject preparation (full headers, no
  file written, name guard), hit ranges on a sequence without hits,
  overlap check on synthetic arrays (terminus/repeat quiet, repeat/repeat
  reported), and a `merge_hits = FALSE` run. `tell_tales`,
  `tell_tales_guards`, `tell_tales_correction`, `golden`, `annotale`,
  `external_exit_status`, `tale_parts` pass; golden unchanged. Quick
  `R CMD check` (no tests, no vignettes): 0/0/0.
- **New finding, not fixed (to discuss):** a third stray warning,
  "GRanges object contains N out-of-bound ranges", from
  `GenomicRanges::resize()` in `tell_tales()`'s array extension (seen on
  the PXO86 excerpt and `toy_tal_regions.fasta`, at HEAD too). The
  following `trim()` clips the range, so the result is right; the
  warning is emitted before the clip. Fix: compute the extended end
  clipped to the sequence length, or muffle only that warning around
  `resize()`.
  **Done 2026-10-01 (maintainer: suppress):** a calling handler around
  `resize()` muffles only warnings whose message contains "out-of-bound
  range"; the termini test in `test_tell_tales.R` asserts the PXO86
  excerpt run no longer emits it.

## 37. `tales_from_annotale()`, renames and retirements (road map step 3) **[V]**

Investigated 2026-10-01; plan awaiting the maintainer's approval.

**AnnoTALE's own output.** `run_annotale_predict()` on `MAI1.fa` (2 min)
writes `Predict/` (5.8 MB: GFF3, GenBank, TALE DNA and protein fasta) and
`Analyze/` (80 KB: `TALE_Protein_parts.fasta`, `TALE_DNA_parts.fasta`,
`TALE_RVDs.fasta`), one file each for all TALEs. Names look like
`MAI1-tempTALE1 [624136-627961:1]` (0-based start, then strand); the GFF3
gives the contig (`seqid`) and `Id=MAI1-tempTALE1`. 9 TALEs on MAI1,
where `tell_tales()` finds 10 arrays. `.tale_parts_from_file()` and
`.rvds_from_annotale_file()` already parse these files unchanged.

**Plan proposed:**
- T1 `tales_from_annotale(annotale_dir, terminus_max_evalue = 1e-5,
  sanitize = FALSE)`: finds the three `Analyze` files under
  `annotale_dir` (so `run_annotale_predict()`'s `output_dir` works), codes
  the termini with `.tale_termini_hmmsearch()` as `tell_tales()` does,
  `array_id` = AnnoTALE's name up to the space, `seqnames` from the
  `Predict` GFF3 when present. The assembly part of `.tale_parts()`
  (protein/DNA join, inconsistent arrays, absent termini, positions, RVD
  join) becomes a helper shared by both readers; the terminus codes come
  from `array_report.tsv` for one and `hmmsearch` for the other. Test
  fixture: the `Analyze` output of `bai3_sample_tal_genomic_regions.fasta`,
  written by a `data-raw/` script.
- R1 `tales_from_telltale()` -> `tales_from_telltales()`;
  R2 `tales_width()` -> `tales_msa_width()`;
  R3 `repeat_to_rvd_align()` -> `.repeat_to_rvd_align()`;
  R4 (optional, from the pre-1.0.0 list) `tales_rvd_strings(rvd_only =)`
  -> `repeats_only`, as in `tales_coded_strings()`.
- D1 retire `tale_parts_to_rvd()` (same as
  `tales_rvd_strings(rvd_only = FALSE)`), `repeat_to_rvd_map()`,
  `repeat_to_rvd_map_distalr()` to `inst/legacy/`; the
  `test_plot_tales_msa.R` fixture builds its `dom_code`/`rvd` map inline;
  the two golden rows move to `tales_rvd_strings()` or are dropped.
- D2 retire `.rvd_to_repeat_align()` with its tests in
  `test_error_conditions.R`.
- D3 `.tales_rename_legacy()`/`.pairwise_distances_rename_legacy()`: no
  code in `R/` produces the old spellings and no shipped data carries
  them; they serve only tables saved by old versions. Documented in
  `?tales`, `?pairwise_distances` (example uses `TAL1`/`TAL2`/`Sim`),
  `?tales_group_hclust`; `.as_mafft_score_table()` has its own
  `RepU1`/`RepU2`/`Sim` branch. Recommended: retire all of it.

**Maintainer, 2026-10-01:** Q14 yes (`array_id` up to the space), Q15
`tales_msa_width()`, Q16 yes (`repeats_only`), Q17 yes (retire D3),
Q18 update the articles and rebuild the site.

**Done 2026-10-01 (T1, R1-R4, D1-D3):**
- T1: `tales_from_annotale()` in `R/tales_ingest.R`. `.tale_parts()` split
  into `.tale_parts_assemble()` (reading, protein/DNA join, inconsistent
  arrays, absent termini, positions, RVD join) and `.tale_parts_finish()`
  (terminus codes, `seqnames`, missing-sequence warning), shared with
  `.tale_parts_annotale()`. `.tale_parts_from_file()` and
  `.rvds_from_annotale_file()` cut ids at the first space (tell_tales'
  `ROI_*` ids have none). New error `tantale_error_annotale_missing` when
  no parts file is found. Fixture `inst/extdata/annotaleExampleOutput`
  (36 KB; analyze's three files plus predict's GFF3, renamed without the
  parentheses R CMD check calls non-portable), from
  `data-raw/make_annotale_example_output.R`. On that fixture the RVD
  strings, termini included, are the same set as `tellTaleExampleOutput`'s.
  Tests: `test_tales_from_annotale.R`.
- R1, R2, R4 everywhere (R, tests, articles, earlier NEWS entries of this
  version); R3 `.repeat_to_rvd_align()`.
- D1, D2: the four functions in `inst/legacy/conversion_retired.R`; their
  tests removed; the XXXXX regression test now runs on
  `tales_rvd_strings()`. The plot fixture builds its map inline
  (`rvd_map_of()` in `test_plot_tales_msa.R`). §2's condition for
  retiring `repeat_to_rvd_map()` (its one-RVD-per-code assertion) is met
  by the `aa_seq_rvd_inconsistent` anomaly plus the `dom_code`/`aa_seq`
  bijection check, both tested.
- D3: converters and their name vectors in
  `inst/legacy/legacy_column_names.R`; `.as_mafft_score_table()` lost its
  `RepU1` branch; docs updated. The test fixture
  `sampleDistalrOutput.rds` carried the old names in `repeat.similarity`
  and `tal.similarity`; converted by `data-raw/rename_sample_distalr_output.R`
  (it has no generator). Name-clash and rename tests removed.
- Golden re-baselined, both rows explained: the requirements table lost
  the rows of `repeat_to_rvd_map_distalr()` and `tale_parts_to_rvd()`; in
  the projections test the fourth snapshot (`repeat_to_rvd_map_distalr()`)
  is replaced by `tales_rvd_strings(x, repeats_only = FALSE)`, whose
  digest equals the retired `tale_parts_to_rvd()`'s output on the same
  fixture (checked by sourcing the legacy file, C collation). The
  `tale_parts_to_rvd()` line never had a recorded snapshot (the old file
  held four projection entries).
- Side finding, not acted on: `tales_rvd_strings()` orders arrays with
  `split()`, so by the session's collation; under `fr_FR` "BAI3_ROI_*"
  sorts before "BAI3-1-1_ROI_*", under C (testthat) after. The retired
  `tale_parts_to_rvd()` used dplyr's C-locale ordering.
- Full `devtools::test()`: only the two golden failures above.

**Site rebuilt 2026-10-01** (wipe of `docs/` and the cache, 27 min, no
error; `check_built_site()` clean). Articles: calls renamed;
`has_all_domains` replaced by `nterm_dna_hit`/`cterm_dna_hit` in
`tale_mining.qmd` and by `cterm_dna_hit` in `trunctale_correction.qmd`
(the report tables narrowed to five columns so `cterm_aa_length` and
`orf_coverage` stay visible; its quoted numbers 42/183/216 aa, 83/86/93%
still match). `tale_classification.qmd`'s `discover_all` chunk gets
`results = "hide"`: on a cold cache DECIPHER's progress output was
printed (`verbose` is not passed by `tell_tales()`). `tale_msa` and
`tale_target_prediction` renders are unchanged.

**`tale_mining.qmd` no longer matches its render** (left to the
maintainer, who is revising it). Checked by running the article's two
BAI3-1-1 cases directly:
- raw (`cterm_min_score = 300`): `tales_anomalies()` returns 0 rows. The
  text expects `missing_rvd` for `ROI_00003`/`ROI_00005`. Since §35 those
  two arrays enter the `tales` object as two termini and no repeat (N
  coded `NTERM`, the 3-aa "C-terminus" `XXXXX`), and no anomaly check
  flags an array without repeats. `tales_from_telltales()` emits no
  warning either, so "the warnings obtained when importing" has nothing
  to show (and the article sets `warning = FALSE` globally).
- corrected (`max_comparisons = 20`): `tales_anomalies()` returns 0 rows;
  `ROI_00001` (ORF coverage 72%) failed in AnnoTALE, so `tell_tales()`
  deleted its parts and the array is simply absent from the object. The
  text says it is now flagged. Coverage of the other arrays: 90-93%.
- Candidate for the postponed A5 (anomalies re-examined one by one): an
  array with no repeat part as an anomaly.

**2026-10-01, after the maintainer's answers (Q19-Q21):** pushed
(Q20). Q21 done: `tale_mining.qmd` opens with a section on AnnoTALE as
the standard tool (`run_annotale_predict()`, `run_annotale_build()`,
`tales_from_annotale()` on the shipped example), then what `tell_tales()`
adds. Facts behind it, measured with `run_annotale_predict()`: BAI3, 9
TALEs, none flagged; BAI3-1-1, 8 TALEs, all flagged "putative pseudo
gene" (names carry "(Pseudo)"); MAI1, 9. Running it on BAI3-1-1 found a
bug of `tales_from_annotale()`: `tempTALE8` has no repeat and an empty
RVD record, and `.rvds_from_annotale_file()` built `1:0`. Fixed
(`seq_along()`, the empty-element warning muffled there); test with a
planted record. Not recorded in the `tales` object: AnnoTALE's
"(Pseudo)" flag.

## 38. Non-standard TALE structure as an anomaly (Q19) **[V]**

Maintainer, 2026-10-01: anything that is not a standard TALE,
`NTERM - RVD x n - CTERM`, should be reported by `tales_anomalies()`.

Measured (arrays not of that form / arrays): test fixture 1/44 (one
`XXXXX`); the articles' discovery cache (MAI1, BAI3, BAI3-1-1 corrected)
0/26; `tellTaleExampleOutput` 0/4; AnnoTALE predict on BAI3-1-1 7/8
(seven `XXXXX`, one without repeats); `tell_tales()` on BAI3-1-1 raw 7/8
(`ROI_00002` is the only standard one).

**Plan proposed:** three checks in `.tales_anomalies()`:
`terminus_absent` (no N- / no C-terminus part; needs `domain_type`),
`terminus_unmatched` (a terminus coded `XXXXX`; needs `rvd`),
`no_repeat` (no repeat part; needs `domain_type`). Objects without
those columns (e.g. `as_tales()` on repeat-only strings) are not
checked. Consequences: `tales()` warns about such arrays, and
`sanitize = TRUE` drops them; golden anomaly table +1 row; the
`tale_mining.qmd` BAI3-1-1 narrative gets its signal back.

**Done 2026-10-01 (Q22 approved as proposed):**
- `.tales_anomalies()`: `terminus_absent` and `no_repeat` (need
  `domain_type`), `terminus_unmatched` (needs `rvd` too). Documented in
  `tales_anomalies()`; NEWS entry. Tests in `test_tales_class.R`
  (§ "Standard TALE structure") and `test_tales_from_annotale.R`.
- Internal rebuilds no longer repeat the anomaly warning:
  `tales_assign_domain_codes()` and `tales_align()` muffle
  `tantale_warning_tales_anomalous` (their input was reported on when it
  was built).
- The shared fixture `sampleDistalrOutput.rds` holds one non-standard
  array, `PXO86_ROI_00019` (C-terminus `XXXXX`), the genuine truncTALE.
  Tests that only load it use `tales_quietly()`
  (`helper-fixture-tales.R`); the projection fixtures with a deliberately
  incomplete `a2` too. Golden: one row added to the anomaly table
  (`PXO86_ROI_00019 terminus_unmatched`), accepted.
- In the articles' own PXO86 run, `ROI_00019`'s C-terminus is coded
  `CTERM`, so `trunctale_correction.qmd`'s "Neither shows up in
  `tales_anomalies()`" still holds. Why the fixture codes it `XXXXX` is
  not investigated.
- `tales_class.qmd`: definition of an anomaly rewritten (it said a missing
  terminus is not flagged).
- `tale_mining.qmd`, BAI3-1-1: raw run flags 7 of 8 arrays (all but
  `ROI_00002`): `no_repeat` for `ROI_00003`/`ROI_00005`, `XXXXX` termini
  for the others; import warning now shown (`warning = TRUE` on that
  chunk). `correct_tales()` route: `ROI_00006` and `ROI_00009` keep an
  unmatched N-terminus (both DECIPHER runs repair them), so the text no
  longer calls that result clean, and the `sanitize` example now drops
  those two from `bai311_java` instead of nothing from `bai311_corr`.
- **Still wrong in `tale_mining.qmd`, left to the maintainer:** the
  `max_comparisons = 20` narrative says `ROI_00001` is "now flagged"; it
  is absent from the object (AnnoTALE failed on it and `tell_tales()`
  drops it), and an absent array cannot be an anomaly. The coverage
  section reads `ROI_00001`/`ROI_00005` as "the two originally-broken
  arrays" (those are `ROI_00003`/`ROI_00005`), and says all three routes
  bring both to the same coverage, while `max_comparisons = 20` leaves
  `ROI_00001` at 72%.
- An investigation of a suspected state leak between successive
  `tell_tales()` calls in one session was a misreading of the diff: no
  leak (raw then corrected in one session gives the same result as a
  fresh session; three repeated corrected runs agree).

## 39. Road map R5, R2, R3, R10 (Q25) **[V]**

Agreed 2026-10-01 (Q25, Q28-Q32). Done the same day:
- R5: `vignettes/getting_started.qmd` renamed `vignettes/tantale.qmd`;
  pkgdown links a vignette named after the package as "Get started" in
  the navbar and exempts it from the articles index, so it left the
  "Learn tantale" list. Links in `tale_msa`, `tales_msa_class`,
  `tale_target_prediction`, `pkgdown/index.md`, `dev/CLAUDE.md` and
  `dev/dev-notes.Rmd` updated. Full site wiped and rebuilt (17 min).
- R2: `plot.tales_msa()` builds every RVD column (consensus match and
  `rvdSimVsRef`) in one block. `refTaleId` comes from `rvd_align` there
  and is overridden by the domain block when domain distances are given,
  the same final value as before. Golden: plot data identical, only
  `rvdSimVsRef`'s column position moved.
- R3: `plot.tales(facet_by = "seqnames")`, a character vector of column
  names or `NULL` (Q30: name `facet_by`, vectors accepted). The default
  gives no panels when `seqnames` is absent; an explicit absent column,
  or one that varies within an array, is an error
  (`tantale_error_plot_facet`). Panels stay `facet_grid()` rows with free
  heights (Q31). Example added to `tale_classification.qmd` (Q32):
  `plot(all_tales, facet_by = "strain")`.
- R10: README "TALE mining" bullets reworded without function names
  (Q28).

## 40. Arrays in C-locale order (D7) **[V]**

Found in §37: `split()` groups arrays by a factor whose levels follow the
session's collation, so under en_US/fr_FR "BAI3_..." came before
"BAI3-1-1_..." and under C (tests, `R CMD check`) after. Maintainer chose
O1 (Q36, Q37), 2026-10-01. `.array_factor()` (`R/tales_projections.R`)
sorts the levels with `sort(method = "radix")`, which compares bytes;
used by `tales_coded_strings()`, `tales_rvd_strings()`,
`.tales_assemble_seq()` (`tales_get_protein_seq()`,
`tales_get_dna_seq()`) and `tales_align()`'s MAFFT input. Their
`order()` calls take `method = "radix"` too. `format.tales()` (order of
appearance) and the internal checks were left alone. Test sets
`en_US.UTF-8` collation with `withr::local_collate()` (fails if that
locale is absent). Golden unchanged (it runs under C). Article renders
(en_US session): only order changes, in `tale_classification`'s functal
table, tree and logos and `tale_target_prediction`'s figure; no alignment
changed. D6 (keep AnnoTALE's "(Pseudo)" flag): maintainer said no.

## 41. Road map R1, R4, R8 (Q40-Q44) **[V]**

Maintainer, 2026-10-01: R1 P1, R4 yes, R8 by hand for now (switch to
push/PR triggers when the package is mature) and test other OSes; R6
(vdiffr) asked for an opinion, then deferred by the maintainer
2026-10-02 (opinion given: low priority, five plots would do); R7 and
R9 wait, both tied to meeting rOpenSci's requirements for review and a
JOSS paper.
- R1 done: `.tidy_biostrings_msa()`'s comment now describes what it takes
  (a named XStringSet of equal-length sequences) and its one caller,
  `plot_target_preds()`. The file input it claimed was never needed.
- R4 done: `.telltale_log()` printed the log with `cli_inform()`, which
  rewraps lines into one paragraph and evaluates `{...}`: a braced path
  crashed `tell_tales()` at its last step. Now `cli::cli_verbatim()`.
  Test with an output directory `run{1}`.
- R8 written, not yet run: `.github/workflows/R-CMD-check.yaml`,
  `workflow_dispatch` with an `os` choice. Builds the environment through
  micromamba + `tantale_setup(install = TRUE)`, as a user would.
  Platforms, checked with `micromamba search` 2026-10-01: every pinned
  tool has osx-64 builds; MAFFT 7.453, HMMER 3.3.2 and mmseqs2 14 have
  none for osx-arm64 or win-64. So Windows is out (also `OS_type: unix`),
  Apple Silicon Macs cannot build the environment as pinned (an osx-64
  environment under Rosetta would need `tantale_setup()` support), and
  the workflow offers an Intel Mac. 2026-10-02: `macos-13` had been
  retired by GitHub (4 December 2025), so the option is now
  `macos-15-intel`, GitHub's last Intel image (supported until autumn
  2027; no Intel macOS runner after that). Private repository on the free
  plan: macOS minutes count 10x.
- R8 first run, 2026-10-02 (Q45), `ubuntu-latest`, run 37062861053:
  passed at the first attempt. `R CMD check` Status OK (no NOTE), tests
  `[ FAIL 0 | WARN 48 | SKIP 0 | PASS 1047 ]`, R 4.6.1. 30 min: R
  dependencies 17 min (no cache yet), environment 1 min (all pins
  resolved: MAFFT 7.453, HMMER 3.3.2, mmseqs2 14.7e284, clustalo 1.2.4,
  igvtools 2.16.2), check 10 min. The "X File existence..." annotations
  are stderr of tests that provoke tool failures on purpose. GitHub
  annotates `actions/checkout@v4` and `setup-micromamba@v2` as Node 20
  actions, forced onto Node 24. macOS not run yet.

## 42. Terminus check: a coverage rule (Q51) -- adopted (Q56) **[V]**

The terminus check (§35) codes a terminus `NTERM`/`CTERM` when its best
`hmmsearch` match against the TALE N- or C-terminal protein profile has
a full-sequence E <= `terminus_max_evalue`, whatever part of the profile
the match covers. Profiles: N-terminus 288 positions, C-terminus 279.
BAI3-1-1 raw `ROI_00001`: a 247-aa N-terminus matching positions
103-150 (E = 4e-22) is coded `NTERM`. A minimum covered fraction would
reject genuine truncated termini (PXO86 truncTALE: 42-aa C-terminus,
positions 1-37).

Proposed rule: the match must also reach the profile end that faces the
repeats, within 10 positions (N-terminus: `hmm_to >= 278`; C-terminus:
`hmm_from <= 11`), read from `hmmsearch --domtblout`. **Maintainer,
2026-10-02:** agreed to a trial first: score every AnnoTALE terminus of
the article genomes and the test fixture under both rules and list
every terminus whose code changes, for the maintainer to judge before
the rule is adopted.

Trial, 2026-10-02 (scripts kept out of the repo; rerunnable from this
description). `tell_tales()` at HEAD on five settings: MAI1, BAI3,
PXO86 (defaults), BAI3-1-1 raw (`cterm_min_score = 300`) and BAI3-1-1
corrected (`cterm_min_score = 300, correct_array = TRUE,
max_comparisons = 50`), plus the four segments of
`termini_profile_cases.fa`. Every N-/C-terminus segment from AnnoTALE's
`TALE_Protein_parts.fasta` was searched with `hmmsearch --tblout
--domtblout`; the span of the domains with i-E <= 1e-5 gave the profile
gap on the repeat side. The current rule recomputed this way equals
`array_report.tsv`'s `*_aa_hit` on all 104 segments.
- 108 segments, 99 matching under the current rule, 95 under the
  proposed one. The 4 changes are all N-termini of BAI3-1-1 raw:
  `ROI_00001` (247 aa,
  profile 103-150), `ROI_00002`, `ROI_00007` (287 aa) and `ROI_00008`
  (288 aa), each matching profile positions 1-150 only, E 1e-88 to
  1e-90. In each, the last 138 residues of the segment and of the
  profile are left unaligned. No other terminus changes: MAI1, BAI3,
  PXO86 (truncTALEs included), BAI3-1-1 corrected and the fixture's
  verdicts are identical.
- The same four ROIs in the corrected run match profile 1-288 over the
  whole segment. Their raw DNA carries `GTACGCGCAG` where the corrected
  sequence has `GTACCGCAG` (N-terminus nt ~447, codon ~150): one inserted
  base, and the protein reads `ASPVRAGGSTHARL...` instead of
  `ASPVPQVDLRTLG...` from there on. Whether these are sequencing errors
  is the maintainer's call (BAI3-1-1 is the flagged, non-gold assembly).
- Closest genuine cases to the tolerance: BAI3-1-1 raw `ROI_00003` and
  `ROI_00005`, 285-aa N-termini ending at profile position 286 (gap 2).
  Everything else that matches reaches the end (gap 0). So a tolerance
  of 10 sits between 2 and 138.

**Maintainer, 2026-10-02 (Q56): adopted, all three parts.** Done
2026-10-03:
- `.tale_termini_hmmsearch()` (`R/telltale.R`) adds `--domtblout`; the
  domains with i-E <= `terminus_max_evalue` give
  `nterm_aa_profile_gap`/`cterm_aa_profile_gap`, and a hit needs the gap
  <= `.terminus_max_profile_gap` (10L, internal, argument
  `max_profile_gap` of the helper only). Used by `tell_tales()` and
  `tales_from_annotale()`. `tales_from_telltales()` still reads the hits
  from `array_report.tsv`, so an old output directory keeps its old codes.
- `array_report.tsv` gains the two gap columns, between the E-values and
  the hits (Q56b); docs of `tell_tales()`, `tales_from_annotale()` and
  `tales_anchor_codes()` explain the rule in biological terms.
- Fixture: BAI3-1-1 raw `ROI_00002`'s N-terminus added to
  `termini_profile_cases.fa`; `test_tell_tales.R` asserts gap 138 and no
  hit, gap 0 for the complete and truncated cases (Q56c).
- Golden re-baselined: in both `tell_tales()` runs, `array_report.tsv` and
  `all_ranges.gff` changed (the GFF carries the array metadata as
  attributes: the four array lines gain the two gap attributes, line count
  unchanged at 199); the column fingerprint gains the two columns (all 0
  on the BAI3 sample). Every hit column is unchanged.
- `data-raw/make_telltale_test_fixtures.R` rerun: besides the new columns,
  only dates and temporary paths moved.
- Checked end to end on BAI3-1-1 raw: the four N-termini now `XXXXX`
  (N-terminus codes 2 `NTERM`, 6 `XXXXX`, were 6 and 2).
- Not yet done: `tale_mining.qmd`'s raw BAI3-1-1 section, which the
  maintainer is revising, will show the new codes once re-rendered.

## 43. Pre-1.0.0 decisions (Q52-Q55) **[V]**

Maintainer, 2026-10-02:
- Q52: `tell_tales()` keeps flat arguments, no grouping into lists.
- Q53a: sequences are treated as linear. `?tell_tales` now says so, with
  the consequence (a *tal* gene spanning the junction of a circular
  molecule is cut in two) and the workaround (rotate the sequence). A
  `circular` argument can come later without breaking a call; the
  `TODO` block became a comment pointing here.
- Q53b: `putative_tal_orf.fasta` and `pseudo_tal_cds.fasta` are kept,
  their docs now say exactly what they hold, and the `TODO` block is
  gone. Nothing in the package, tests or articles reads them.
- Q54: ARLEM's duplication/insertion costs (exposable as arguments of
  `tales_tale_distances()`) and §21 option (c) (`tales_consensus()`/
  `tales_consensus_match()` accepting a `tales_msa`, as S3 methods) move
  after 1.0.0: neither breaks an existing call.
- Q55: a 0.99.0 release candidate, yes. The maintainer still has work to
  do on the website and public-facing material; that is compatible, since
  0.99.x only freezes the API in intent (hard renames remain allowed)
  and the site and docs can change between 0.99.x versions.
- §2 and §20 were found already settled by §37's D1 (their functions are
  in `inst/legacy/conversion_retired.R`); START HERE still listed them as
  reserved. Corrected.

## 44. pkgcheck findings (rOpenSci) and their triage (Q57) **[V]**

`pkgcheck` 0.3.2 (from `ropensci-review-tools.r-universe.dev`; its
`pkgstats` needs universal-ctags and GNU global, installed for the run in
a scratch micromamba prefix, not system-wide), run 2026-10-02 on e5022fd
in a detached worktree, 24 min. Verdict "not quite there yet". Coverage
90.7%. Package name available, roxygen2, URL/BugReports, HTML vignette,
website: all fine. The findings and what became of them (maintainer
followed the recommendations, 2026-10-03):
- **P1** no contributing file: `.github/CONTRIBUTING.md` written
  (reporting, development environment and pins, tests that fail rather
  than skip, golden snapshots, cli conditions, rOpenSci naming, NEWS,
  `devtools::check()`). pkgdown links it from the home page sidebar.
- **P2** six exports without examples (`is_tales()`, `is_tales_msa()`,
  `is_pairwise_distances()`, `validate_tales_msa()`,
  `validate_pairwise_distances()`, `tales_msa_width()`): examples added;
  the `tales_msa_width()` one shows the width surviving a subset that
  empties the last columns (an array without a C-terminus).
- **P3** PDF manual failed (LaTeX cannot take the GIF logo of
  `?tantale`): the figure is now inside `\if{html}{}`; `R CMD Rd2pdf`
  builds the manual.
- **P4** R CMD check ERROR, two golden `tell_tales()` tests: pkgcheck
  installs the package under the test session's `tempdir()`, so the
  tempdir rule of `helper-golden.R` dropped the five log lines naming
  the package's files (`n_dropped` 2 -> 7, digest `4dfb575...`,
  reproduced exactly). `.telltale_file_digest()` now rewrites the
  installation directory before the drop; baseline unchanged; a test
  moves the log under `tempdir()` and expects the same digest.
- **P5** 42 Imports, above the 99th percentile: reviewers will ask. Audited
  2026-10-03: every package has a live call site. Three could move to
  Suggests with minimal surgery (`universalmotif`, `gplots`, `biovizBase`);
  the plot cluster (`ggtree`/`tidytree`/`aplot`/`viridis`/`ggnewscale`)
  would need more work; the rest are embedded in core functions and should
  stay. **Decision (maintainer, 2026-10-03): keep all 42 in Imports.** The
  answer to a reviewer is "all are genuinely used". The real size problem
  is R7 (the jars in `inst/tools/`), not the dependency count. Since §46,
  40: the colour style left `biovizBase` and `viridis` without a call site.
- **P6** goodpractice lints (40 kinds). Fixed the bug-prone ones:
  `class(x) %in% ...` -> `inherits()` (`.split_list()`, which errored
  on a multi-class object); every `1:length()`/`1:nrow()`/`1:ncol()`
  in `R/` -> `seq_along()`/`seq_len()` (19 sites; `1:0` is `c(1, 0)`,
  e.g. the hit ids of a run without hits); `sapply()` -> `vapply()` or
  `lengths()` where one value of a known type is expected per element
  (`tell_tales()`'s array metadata, `.split_list()`,
  `.tale_parts_finish()`'s check, `tales_consensus()`'s counts, the
  MAFFT hex encoding, three in `classification.R`,
  `plot_target_preds()`'s RVD labels). `plot_target_preds()` also
  applied `rev()` to each RVD one at a time, a no-op; the minus-strand
  order comes from `xPos` running from end to start, now said in a
  comment. Left as they are: `sapply()` calls that build matrices or
  arrays on purpose, or whose element type varies (`tales_consensus()`
  itself, which takes character or numeric codes); `inst/legacy/`; the
  pure style lints (1325 long lines, `=` assignment, `paste(sep =
  "")`, leading zeros...). False positives: the four "never called"
  functions (`[.tales` and `[.pairwise_distances` are registered S3
  methods, `%||%` is an infix used three times,
  `.repeat_to_rvd_align()` builds a test fixture) and the 20
  "duplicate arguments" (repeated `"i"` names in cli bullet vectors,
  which is how cli takes several bullets).
- **P7** "no continuous integration": pkgcheck looks for a README badge
  or asks GitHub, which it cannot for a private repository. Add the
  badge when the repository goes public; rOpenSci will also expect CI on
  push, which §41 left manual on purpose.
- **P8** BaoVi TramVi's ORCID iD (0000-0002-4319-5544) added to `Authors@R`, 2026-10-03.
- The NOTE about a hidden `.git` is an artefact of checking a worktree
  (`.git` is a file there).
- Not flagged by pkgcheck, still the main obstacle: the 57.8 MB source
  tarball against the 5 MB limit (R7, §34).

---

## 45. `tales_names()` (2026-10-03) **[V]**

The maintainer asked for a `names.tales` method returning the array ids.
Rejected: a `tales` object is a tibble, and dplyr, tibble, `$` and
`print()` all read column names through `names()`, so the method would
break every verb. Added the exported accessor `tales_names(x)`
(`unique(x$array_id)`, row order; family "tales objects"), with tests on
a `tales` and a `tales_msa` in `test_tales_class.R`. The name echoes
`tales_namespace()` and the names of `tales_rvd_strings()`'s output. The
internal `unique(x$array_id)` calls were left as they are.

---

## 46. A colour style for the package's plots (2026-10-03) **[V]**

Raised by the maintainer: the labels of `plot.tales()` are hard to read.
Rule (maintainer, all projects): every palette must be colour-blind safe.
Candidates compared on real data, with deuteranopia, protanopia and
tritanopia simulated by `colorspace` (current colours, Paul Tol's muted
and light schemes, Okabe-Ito, Crameri's scico; viridis excluded by the
maintainer). Under deuteranopia the current `aa_length` fill shows the
20, 33 and 34 aa repeats as one olive.

Decided:
- **Paul Tol's muted scheme** (its wine matches the logo), hex codes kept
  in the package with credit, no new dependency.
- **Colours by biological role** in `plot.tales()`: canonical 34 aa
  repeat sand, final 20 aa half-repeat pale blue, other repeat lengths
  the strong colours; N-termini shades of wine and C-termini shades of
  teal, lighter when shorter. Legend "Part and length" ("repeat, 34 aa").
  Label text black or white by the fill's luminance.
- **Alignment labels, kept simple:** the text colour still says whether a
  residue matches the consensus, black for a match and red `#CC3311` for a
  mismatch; the fills are restricted to pale shades so that both read on
  every fill. The diverging `rvd_sim` scale runs from pale orange (-1) to
  pale blue (+1).

Done (plan agreed 2026-10-03):
- `R/palette.R` holds the colours (`.tol_muted`, `.tol_light`,
  `.tantale_colours`) and `.text_colour_on()` (WCAG luminance, threshold
  0.3). The comparison scripts and sheets lived in the session scratchpad
  only.
- `plot.tales()`, `plot.tales_msa()` (pale teal-sand ramp for
  `domain_clust`, pale YlOrBr for `domain_sim`, pale sunset for
  `rvd_sim`, wider colour bars), `plot_target_preds()` (DNA bases A green,
  C indigo, G wine, T rose: the Tol muted set dark enough for text on white
  that stays most distinct under all three simulations, by CIEDE2000;
  match score 1/2/3 from strong orange to pale yellow; EBE box sand),
  `talomes_heatmap()` (Tol muted ranks, recycled; Tol light for the
  `extra_col` bar; truncation "T" black or white by fill), the hclust
  dendrogram (Tol muted, trunk grey) and the k-medoids silhouette plot
  (cyan, chosen k in wine).
- `biovizBase` and `viridis` left Imports (42 -> 40; P5 above).
- Tests: `test_plot_tales_composition.R` (role colours, text colour),
  `test_plot_tales_msa.R` (label colours, fills pale enough for black).
- The 18 tidyselect "external vector" warnings in
  `test_target_predictions.R` predate this change (seen with it stashed)
  and come from no call in `R/`.
- Second round (maintainer's review, 2026-10-03): `plot.tales()` strips
  on a pale band (`#F2F2F2`) with bold dark text. Dendrogram planned for
  25+ groups: Tol muted cycled in leaf order so neighbours differ, group
  number under each clade, no legend, cut height in the subtitle (the old
  in-plot label overlapped the branches). `plot_target_preds()`: the
  maintainer preferred the original DNA colours (ColorBrewer RdYlBu, as
  biovizBase; they are also more distinct under the simulations than the
  Tol set tried, CIEDE2000 min 23 vs 12), kept as hex codes; match score in
  purple shades (chosen over grey, teal and indigo; orange rejected).
  `talomes_heatmap()`: Tol muted ranks rejected; `colors` now gives the
  ends of a ramp interpolated over the ranks present (a fixed-length
  palette left 2-3 ranks all dark, viridis included), default dark to pale
  wine (chosen over indigo-cyan, teal, indigo-teal-sand).
- Third round: talome ramp inverted (most common pale, rare dark, the
  maintainer's call), pale end deepened to `#E3A9B8`; empty cells white
  (grey and the pale end were CIEDE2000 1.2 apart under deuteranopia,
  white vs `#E3A9B8` 15.4). Dendrograms took `margins` units (5 each)
  against 1 unit per cell, so they dominated a small talome: `margins =
  NULL` (default) now sizes them to a quarter of the heatmap, 1.5-5 cells;
  `plot_type = "single"` needs at least 3.5 cells on top, since
  `heatmap.2()` draws the title in that panel ("figure margins too large"
  in `test_talomes_heatmap.R` otherwise).
- For the site rebuild: `tale_classification.qmd` says "the darkest is
  the most common" and "A grey cell means the strain has no member"; both
  are now the other way round (palest; white). **Done 2026-10-03 (§47)**,
  along with two other colour names the change made wrong: the silhouette
  plot's pick (red, now wine) and `plot.tales_msa()`'s text colours in
  `tales_msa_class.qmd` (cyan/pink, now black/red).
- For the site rebuild: `tale_mining.qmd` line 439 uses `ggplot2::scale_fill_viridis_d()` in
  article code (left to the maintainer, §38).

---

## 47. A structure figure on the home page (2026-10-03) **[V]**

Asked for by the maintainer (2026-10-03, notes from a read of the site): a
figure of the repeat array wrapped around DNA on `pkgdown/index.md`. The
Wikipedia images were not used (licence to check, and the home page says
the figures come from the authors' scripts).

Done: `pkgdown/crd_figure.py` draws PthXo1 bound to its target (PDB 3UGM,
Mak et al. 2012, doi 10.1126/science.1216211, checked on Crossref) with
PyMOL, side view and end-on view from the N-terminal end, into
`pkgdown/assets/crd_dna.png` (pkgdown copies `pkgdown/assets/` to the
site root; `pkgdown/` is in `.Rbuildignore`, so the tarball does not grow).
Colours from `R/palette.R`, as in `plot.tales()`: N-terminal region wine,
repeats alternately sand and olive, RVD residues indigo spheres, DNA grey.
Run from the package root with PyMOL and Pillow (the script's header has
the micromamba line); the output is byte-identical between runs. 256
colours, 220 KB.

Checked against the coordinates before writing the caption: 23 repeats
found by their L-T/P-x-x-Q-V-V-A-I-A-S start, RVDs at 12-13 (N* where 13
is missing); the structure stops inside repeat 23, so neither the final
half repeat nor the C-terminal region is drawn. Repeat centroids rise
3.41 A and turn 32.2 degrees per repeat (one repeat per base pair, 11.2
per turn), right-handed. Chain B is the strand the RVDs read, 5' end on
the N-terminal side. The end view is from the N-terminal end (checked in
camera coordinates).

Left for the maintainer: the caption could name the target, the
*OsSWEET11* (Os8N3) promoter; the crystal's DNA contains
`TGCATCTCCCCCTACTGTACACCAC`, which is that EBE as far as I know, but it
needs the maintainer's confirmation (identifying a specific instance).

Site rebuilt in full the same day (wiped `docs/`, `build_article()` per
file, 14.7 min, no errors). Checked against the render: the home page
figure and caption; the silhouette pick in wine; the talome heatmap
(palest = most common, white = no member, MAI1's other variant in five
groups); `plot.tales_msa()` text black on match, grey on a gap
consensus, termini grey in `rvd_sim`; the dendrogram still nine groups,
26 arrays. The rebuild removed the stale `reference/figures/pipeline.*`
(nothing links to them) and added `CONTRIBUTING.html` (linked from the
news page) and the `tales_names()` page.

---

## 48. `talomes_heatmap()` takes a `tales` object (2026-10-03) **[V]**

The maintainer's question above, surveyed over the 51 exports: three
functions consume RVD strings. `talvez()` and `preditale()` already have
`tales_predict_targets()` in front of them, which renders a `tales` with
`tales_rvd_strings()`; the gap was `talomes_heatmap()`, which
`tale_classification.qmd` fed from a five-line table built by hand.
Maintainer's answer (2026-10-03): option (a), and leave `talvez()` and
`preditale()` exported as they are.

Done: `talomes_heatmap()` accepts a `tales` carrying the group and strain
columns (plus `trunc_tales_col`/`extra_col` if asked for), reduces it to
one row per array, and computes the RVD strings itself; `rvd_col` is then
unused. A missing column, or a column with more than one value in an
array, is a `tantale_error_talome_column` error. Arrays whose RVD string
is empty are dropped (`tales_rvd_strings()` keeps repeats only).

The data-frame path is untouched. Two tests in
`test_talomes_heatmap.R`: the `tales` and table inputs write
byte-identical PNGs in both `plot_type`s, and both error paths. The
article now calls `talomes_heatmap(grouped, group_col = "group",
strain_col = "strain")`; its rendered figure is byte-identical to the
one built from the hand-made table.

---

## 49. `max_comparisons` on BAI3-1-1, and the array that disappears (2026-10-03) **[V]**

Raised in `dev/notes_for_claude.md` (maintainer, 2026-10-03): with the
§42 terminus codes, `max_comparisons = 5` seemed to correct every
BAI3-1-1 array, but one array went missing; find out why, find the value
that gives the best result fastest, and compare with `correct_tales()`.

Sweep (scripts in the session scratchpad only): `tell_tales(BAI3-1-1.fa,
cterm_min_score = 300, correct_array = TRUE, max_comparisons = m)` for m
in 1, 2, 3, 5, 10, 20, 50, 100 and all (~500 shipped references), one run
each on this 8-core machine; installed tantale at 55ff763. Read-outs: the
`rvd_string` of `array_report.tsv` (standard = `NTERM-...-CTERM`) and the
array count of the `tales` object, before and after `sanitize = TRUE`.
The report has 9 candidate regions in every run; `ROI_00004` has no
terminus hit and no ORF in any of them.

| run | seconds | standard | arrays | after sanitize | flagged |
|---|---|---|---|---|---|
| uncorrected | 19 | 0 | 8 | 0 | all 8 |
| m = 1 | 17 | 5 | 6 | 5 | `ROI_00001` |
| m = 2, 3, 5 | 19-22 | 7 | 8 | 7 | `ROI_00001` |
| m = 10, 20 | 26, 34 | 7 | 7 | 7 | none |
| m = 50 | 61 | 8 | 8 | 8 | none |
| m = 100 | 101 | 8 | 8 | 8 | none |
| m = all | 455 | 8 | 8 | 8 | none |
| `correct_tales()` + `tell_tales()` | 29 + 13 | 6 | 8 | 6 | `ROI_00006`, `ROI_00009` |

- Every array except `ROI_00001` gets the same RVD string from m = 2 on.
  `ROI_00001` is the hard one: at m = 2-5 its N-terminus is `XXXXX` and
  it carries 18-19 of its 26 RVDs (ORF coverage 63-67%), so it is
  flagged. At m = 10 and 20 the ORF reaches 72% but AnnoTALE fails to
  parse it; `tell_tales()` drops it, noting it only in `tell_tales.log`
  ("Annotale failed to parse TALE domains for ROI_00001"), and
  `tales_from_telltales()` builds a 7-array object with no warning and no
  anomaly. From m = 50 on, `ROI_00001` is `NTERM`, 26 RVDs, `CTERM`
  (coverage 93%), and all eight RVD strings equal those of the full search.
- This is the missing array of the maintainer's note: `tale_mining.qmd`
  still runs `max_comparisons = 20` (its prose says 5).
- Speed: m = 50 takes 61 s against 42 s for the Java route, which leaves
  the N-termini of `ROI_00006` and `ROI_00009` unmatched. m = 50 gives the
  full search's result 7.5 times faster. The earlier measurement in
  `?tell_tales` (four arrays, 1057 references) also found 50 identical to
  the full search.
- A repeated m = 5 run gave the same result as the sweep's.

Maintainer, 2026-10-03: Q58 (after a check on the other genomes), Q59,
Q60, Q61, Q63 and Q64 agreed; Q62 (the `plot.tales()` legend) to be
discussed.

- **Q58 done.** Check: `max_comparisons = 50` and `NULL` on MAI1, BAI3 and
  PXO86 (defaults, `correct_array = TRUE`; the six runs in parallel, so
  their times are inflated): identical `putative_tal_orf.fasta`, RVD
  strings, ORF coverage and indel counts for all 10, 10 and 19 candidate
  arrays; 210 s against 828 s (MAI1), 212 against 814 (BAI3), 314 against
  1311 (PXO86). BAI3-1-1 at 50 and `NULL` (sequential sweep above) also
  give identical ORF files. Default of `tell_tales()` and
  `.telltale_array_orfs()` now 50; the `@param` rewritten around the
  measurements (the old 1057-reference and 20-reference tables dropped:
  the reference set has changed since). `test_tell_tales_guards.R` asserts
  the default; the golden correction run passes `max_comparisons = NULL`
  (its 20-sequence reference makes 50 and all the same search). Golden
  re-baselined: one row, `tell_tales.log` of the uncorrected run, whose
  parameter echo reads `max_comparisons: 50` instead of `all` (checked by
  diffing the two logs; only the date and path lines also differ, and the
  fingerprint strips those).
- **Q59 done.** `.tale_parts()` warns
  (`tantale_warning_annotale_unparsed`) about each `array_report.tsv` row
  with a terminus DNA hit and no `TALE_Protein_parts.fasta` under
  `annotale/<array_id>/`, with a hint about `max_comparisons` when the
  run corrected. Tests: `test_tales_class.R` (parts files removed from a
  copy of the example output), `test_tell_tales.R` (PXO86 `ROI_00001`,
  which AnnoTALE cannot translate, now expects the warning).
- **Q61 done.** `.tales_anomalies()` sorts by `array_id` then `check`,
  radix order. Golden unchanged.
- **Q64 done.** `plot.tales()` gives `array_id` the levels of
  `.array_factor()` reversed, so the first array is at the top. `array_id`
  is always character in a `tales` (validated), so there is no
  user-supplied factor to respect.
- **Q60 done.** `tale_mining.qmd`: correction inside `tell_tales()` at the
  default; the `@CLAUDE` paragraph on the best method, with the times
  measured during the render (`system.time()`): 69 s against 50 s on a
  first render, 61 s against 40 s on the published one (2026-10-03); "A clean correction without an eight-minute wait"
  replaced by `#sec-max-comparisons` (the sweep table, plus a live
  `max_comparisons = 20` run showing the new warning); "How much did
  either correction actually help?" removed (with it the
  `scale_fill_viridis_d()` of §46); `bai311_best` replaced by
  `bai311_corr`; the `#sec-best-correction` anchor moved to "Correcting
  inside `tell_tales()`", which `tale_classification.qmd` links to. Two
  typos of the maintainer's fixed. In the raw-vs-corrected patchwork, the
  two "Part and length" legends do not merge (different breaks) and the
  combined legend is cut at `fig-height: 6`: left for Q62. `fig-keep: last` on the two patchwork
  figures: `plot.tales()` prints as a side effect, so each chunk produced
  three figures.
- Found while rendering: `.telltale_write_correction_alignments()`
  translated sequences whose length is not a multiple of 3, five
  Biostrings warnings per corrected BAI3-1-1 run; now trimmed to whole
  codons (the translation feeds only the diagnostic HTML alignment). The
  N-substitution warning lacked a space ("orderto").
- Site: `tale_mining`, `trunctale_correction`, `tale_classification`,
  `tales_msa_class` and the "Get started" page re-rendered (the last two
  for the new array order of `plot.tales()`; their prose names no order),
  reference, news, llm docs and search rebuilt; `check_built_site()` and
  `check_pkgdown()` clean. Stale `tale_mining` figures removed by hand.
- `trunctale_correction.qmd`: its `max_comparisons` table named "all 494"
  the default; now "50 (default)".
- **Q63 not started**: the design needs the maintainer (see below).

Q63 data (2026-10-03). `hmmsearch --domtblout` of every AnnoTALE terminus
in the runs above, domains with i-E <= 1e-5; "outer gap" = profile
positions left unmatched at the end away from the repeats (N-terminus:
`hmm_from - 1`; C-terminus: `qlen - hmm_to`). Over MAI1, BAI3, PXO86 (at
50), BAI3-1-1 corrected (at 50) and after `correct_tales()`:
- N-termini: every matched one has outer gap 0, at lengths 230, 264, 283,
  287 and 288 aa (the short ones match the whole profile with internal
  deletions). The two 24-aa N-termini of the Java route match nothing.
- C-termini: outer gap 0 at 286 and 297 aa, 1 at 278 aa; PXO86 `ROI_00001`
  (216 aa) 159 and `ROI_00019` (42 aa, the truncTALE) 242.
- So a `terminus_truncated` check with a 10-position tolerance flags those
  two PXO86 arrays and nothing else on the shipped genomes. Open: where the
  measure lives in a `tales` object (a new column on the terminus rows?),
  its name next to the existing `*_aa_profile_gap`, and the consequence
  that `sanitize = TRUE` would drop truncTALEs.

Maintainer, 2026-10-04: Q62 option c, Q63 "nothing to build, explain it",
Q66 yes. Done the same night:
- **Q62c.** `plot.tales()` fills by length alone ("Length (aa)", levels in
  numeric order), one colour per length whatever the part: 34 aa sand,
  20 aa cyan, the other lengths in increasing order from the remaining Tol
  muted colours (rose, indigo, purple, green, olive, wine, teal) then Tol
  light, recycled past 17. Outline legend in the order N-terminus,
  repeat, C-terminus. The §46 role colours (wine/teal shades for termini)
  are gone. The maintainer's terminus labels (`N-`, `-C`, `??`) kept,
  computed with `match()` instead of `rowwise()`. Two patchworked plots
  still show two fill legends when their length sets differ.
- **Q66.** `plot.tales()` returns the ggplot visibly and no longer prints
  it; `fig-keep: last` removed from `tale_mining.qmd`. `plot.tales_msa()`
  left as it is: it prints a patchwork copy so that the legend of an
  aplot composition goes to the bottom, and returns the aplot for its
  `$plotlist` API; changing it changes the class returned.
- **Q63.** Raw PXO86 (the truncTALE article's run), `hmmsearch` domains:
  both 230-aa N-termini match profile 1-105 and 151-288 over residues
  1-105 and 108-230 (an internal deletion of ~45 positions);
  `ROI_00001`'s 183-aa C-terminus matches profile 1-183 over its whole
  length; `ROI_00019`'s 42-aa C-terminus matches 1-37 over residues 1-37.
  The 216-aa C-terminus the maintainer mentioned is `ROI_00001` after
  `correct_array = TRUE`: there only its first 120 residues match
  (profile 1-120). Explained in `?tales_anchor_codes` (codes record
  relatedness; a short matching terminus has probably lost functional
  regions: T3S signal and degenerate repeats in the N-terminal region,
  NLSs and the activation domain in the C-terminal one), a pointer in
  `?tell_tales` (`terminus_max_evalue`), `trunctale_correction.qmd` §2
  (its "impossible structure" paragraph, stale since §38, rewritten
  around these numbers) and one sentence in `tale_mining.qmd` linking to
  it.

Maintainer, 2026-10-04, after seeing mock-ups of both: **Q62a after all.**
The §46 role colours are back; the fill legend is titled "Length (aa)",
its labels give the length only, and `override.aes` outlines each key in
its domain type's colour (keys ordered N-termini, repeats, C-termini, by
length within each). Test: `test_plot_tales_composition.R` (role colours
restored, legend labels and key outlines). Articles with a `plot.tales()`
figure re-rendered: `tale_mining`, `tale_classification`,
`tales_msa_class`, `tantale`.
- `tale_mining.qmd` (maintainer's suggestion, 2026-10-04): the two
  patchwork comparisons replaced by `tales_bind()` of the tagged objects
  (`array_id` prefixed `raw_`/`java_`/`corr_`, a `method` factor) drawn
  with `facet_by = "method"`: one legend, and the same length gets the
  same colour in both panels (with patchwork each plot coloured its own
  lengths, so 33 aa was rose in one panel and indigo in the other). The
  `library(ggplot2)`/`library(patchwork)` calls, now unused, removed.
  Visible in the new figure: `correct_tales()` gives 288-aa N-termini
  where `correct_array = TRUE` gives 287 aa on six arrays (not
  investigated).

---

## 50. Road map R7: package size (2026-10-04) **[V]**

Measured 2026-10-04 at 43e56cb: `R CMD build --no-build-vignettes` gives
a **57.1 MB** tarball. Tracked files outside `docs/`, `dev/`, `pkgdown/`:
`inst/tools/talecorrect` 34.4 MB (`TALEcorrection.jar` 27 MB, `HMMs/`
6.7 MB, upstream `.java`/`.xml`), `AnnoTALEcli-1.5.jar` 15.8 MB,
`PrediTALE.jar` 14.6 MB, the four genomes 19.3 MB, `tests/testthat`
1.8 MB, `tellTaleExampleOutput` 1.0 MB, `hmmProfile` 0.8 MB, TALVEZ
0.1 MB. A trial build of a copy without the jars, `talecorrect/`, the
four genomes and `inst/legacy` gives **1.0 MB**. `.git` is 85 MB.

References: jars and `talecorrect/` in 3 files of `R/` (defaults of the
wrappers) and `test_correct_tales.R`; genomes in `R/annotale.R`
(examples), `test_correct_tales.R` and five articles; `inst/legacy` only
in comments.

Proposal (builds on §34's parked one-archive plan), put to the
maintainer 2026-10-04; nothing to code before it is agreed:
- R7-1 hosting. The repository is private, so its release assets need
  authentication. Options: a Zenodo record (public, permanent, DOI);
  a public companion repository's releases; this repository's releases
  once public. Recommended: Zenodo.
- R7-2 tools archive: the three jars, TALEcorrection's HMMs and upstream
  sources, README (upstream URL, sha256, version, licence, citation),
  `LICENSES/`, `MANIFEST`. TALVEZ (0.1 MB, redistribution permitted)
  stays in the package.
- R7-3 genomes: a second archive in the same record, fetched on first
  use by an exported helper (name to agree) into
  `tools::R_user_dir("tantale", "data")`, checked by sha256. BAI3-1-1 is
  probably not public, so self-hosting all four keeps one mechanism.
- R7-4 package side: `tantale_setup()` fetches and verifies the tools
  archive (or takes a local copy), one internal resolver gives each
  tool's path, the wrappers' defaults move to it. `test_correct_tales.R`
  gets a BAI3-1-1 excerpt around the TALE loci; tests that need a jar
  fail without it (standing rule); CI caches the archive.
- R7-5 `inst/legacy` into `.Rbuildignore` (stays in the repository).
- R7-6 archives versioned on their own (`tools-1`, `genomes-1`), pinned
  in the package by URL and sha256.
- Not proposed: rewriting git history (85 MB); no review requires it.

Maintainer, 2026-10-04:
- R7-1: the archives go to releases of this repository; it goes public
  again once the website and R7 are done. Until then a download needs
  authentication, so development tests use a local copy of the archive.
- R7-3: both archives are fetched by `tantale_setup()`. BAI3-1-1 is
  private but distributed with the package, so it stays in the genomes
  archive. Accessions found 2026-10-04 (NCBI E-utilities): MAI1 =
  `CP025609.1` (RefSeq `NZ_CP025609.1`), BAI3 = `CP025610.1`
  (`NZ_CP025610.1`), PXO86 = `NZ_CP007166.1`; for all three the shipped
  sequence equals the NCBI record (md5 of the upper-cased sequence),
  lengths 4 735 208, 4 723 880 and that of PXO86. The shipped headers are
  `>MAI1`, `>BAI3`, `>NZ_CP007166` (articles show these as `seqnames`).
- R7-5: the maintainer is considering deleting `inst/legacy` outright.
  Open.
- R7-5 settled: **`inst/legacy/` deleted** (maintainer, 2026-10-04, Q69).
  The last commit holding it is **20755db**. FuncTAL's `Info/2014mat18`
  moved to `data-raw/FuncTAL_2014mat18` (the source of
  `rvd_dna_specificity`; `data-raw/rvd_dna_specificity.R` reads it there,
  same values). Comments in `R/arlem.R` and `R/conversion.R` and
  `?rvd_dna_specificity` updated; `inst/COPYRIGHTS` §3 keeps the dataset.
  The standing rule in `dev/CLAUDE.md` now says: park in
  `R/unused_pending_review.R`, delete after review, name the commit in the
  ledger.

Done 2026-10-04 (maintainer: R7-1 this repository's releases; Q67 a
shorter helper name, taken as `tantale_genome()`; Q68 create the releases):
- `dev/make_archives.R` builds `tantale-tools-1.tar.gz` (50.6 MB: the
  three jars, all of `talecorrect/`, `LICENSES/COPYING.GPL-3`, README,
  MANIFEST) and `tantale-genomes-1.tar.gz` (5.6 MB). Releases `tools-1`
  and `genomes-1` created with `gh` on scunnac/tantale (GitHub labels
  `tools-1` "Latest", the other releases being pre-releases); sha256
  `6ebdfbfe...` and `f7eb3570...`, pinned in `R/tantale_data.R`.
  Anonymous download answers 404 while the repository is private.
- `R/tantale_data.R`: `.tantale_data_dir()` (`TANTALE_DATA_DIR`, else
  `tools::R_user_dir("tantale", "data")`), `.tantale_install_archive()`
  (download or `archive_dir`, sha256, untar, MANIFEST check),
  `.tantale_archive_ok()`, `.tantale_tool()` (default of every wrapper's
  jar/HMM argument), exported `tantale_genome()`. `tantale_setup()` checks
  and with `install = TRUE` installs both archives first, whatever the
  conda state; new argument `archive_dir`; result gains `archives`.
  DESCRIPTION: R >= 4.0.0 (`R_user_dir`).
- Removed from `inst/`: the jars, `talecorrect/`, `COPYING.GPL-3`, the
  four genomes. TALVEZ stays. A trial build gave 1.0 MB.
- `test_correct_tales.R` runs on `data_for_tests/bai311_tale_loci.fa`
  (66 kB, six windows of the TALE loci +-3 kb, built by
  `data-raw/make_bai311_tale_loci.R`): 69 corrections, identical in
  position and edit to 69 of the whole genome's 70; the 70th, an insertion
  at 585387, lies outside every array `tell_tales()` reports. New
  `test_tantale_data.R` (a fake archive in a temporary
  `TANTALE_DATA_DIR`: install, tampering, wrong sha256, missing file,
  failed download, missing tool/genome; plus the real installation).
- CI: the workflow fetches both assets with `gh release download` (works
  while private) and passes them as `archive_dir`.
- Articles, examples: `system.file("extdata", "<genome>.fa", ...)` ->
  `tantale_genome()`. README, home page, `inst/COPYRIGHTS` (§1 rewritten,
  §4 genomes added) updated.
- **Machine quirk:** this Claude session runs inside the VSCodium snap,
  which sets `XDG_DATA_HOME=~/snap/codium/495/.local/share`, so
  `R_user_dir()` resolves there. An R session started elsewhere (terminal,
  RStudio) resolves to `~/.local/share/R/tantale` and will not see this
  installation. Fix: set `TANTALE_DATA_DIR` in `~/.Renviron`, or run
  `tantale_setup(install = TRUE, archive_dir = ...)` in each context.

---

## 51. Optional arguments for the jar wrappers (2026-10-04) **[V]**

Maintainer, 2026-10-04: the wrappers of the Java tools (AnnoTALE predict,
analyze and build, PrediTALE, TALEcorrection) do not let a user pass the
tools' own optional parameters. Add an argument (`opt_arg` suggested): a
string of additional arguments passed to the jar, whose default holds the
values used today. "Soon"; plan to agree before coding. Also to settle
then: which of today's hard-coded values become that default, and how
`preditale()`'s existing `opt_param` fits (it already passes a string).
- Maintainer, 2026-10-04: also show the `position` argument. The
  `correct_tales()` vs `correct_array = TRUE` figure now binds the two,
  aligns them with `tales_align()` (on RVDs; no `dom_code` needed) and
  plots `as_tales(msa)` with `position = "alignment"`,
  `facet_by = "method"`. It shows the one difference: `ROI_00006` and
  `ROI_00009` after `correct_tales()` start with a 24-aa N-terminus and a
  32-aa `nv` repeat in place of the `NN`, `HD` repeats; the prose says so,
  checked against the render. The raw-vs-corrected figure stays
  unaligned.
- Maintainer, 2026-10-04, on seeing it: no alignment there ("an alignment
  makes sense only if the TALEs are related"); back to the plain bound
  plot, prose corrected (five arrays differ by one N-terminal residue,
  `ROI_00008` is 288 aa in both). Added at the maintainer's request: a
  dodged bar chart of part lengths, `facet_grid(method ~ domain_type)`
  (free x per domain type), bars coloured by candidate region (`roi`,
  kept by `tag_method()`; Okabe-Ito). Jittered points and stacked bars
  were tried and set aside.

Maintainer, 2026-10-04 (second session): **Q70** one argument name across
the wrappers, `opt_param` (already used by `preditale()` and `talvez()`;
no rename). **Q71** the Java heap flags (`-Xms`/`-Xmx`) get an argument of
their own, separate from the tool's options. (Asked in the session as Q1a
and Q1b; renumbered here to follow the ledger's Q sequence.)

Optional parameters each jar accepts, read from its own help (installed
`tools-1` archive, 2026-10-04). Every wrapper today leaves them at the
jar's default; the only hard-coded values are heap flags:

| tool, stage | options not set by the wrapper (jar default) | heap flags today |
|---|---|---|
| AnnoTALE 1.5 `predict` (v1.4.2) | `Sensitive` (false) | none |
| AnnoTALE `analyze` (v1.4.1) | `n` run name (cosmetic) | none |
| AnnoTALE `build` (v1.4.1) | `c` class cutoff (5.0), `s` significance level (0.01) | `-Xms512M -Xmx6G` |
| PrediTALE `preditale` | `b`/`bs` background, `t`/`sl`/`n` threshold (sl 1e-4), `Strand` (both), `r` reverse penalty (0.01) | `-Xms512M -Xmx2G` |
| TALEcorrection `correct` | none (s, n, r, c, outdir are all set by the wrapper) | none |

The same letter means different things per tool: `s=` is the strain
(predict), the significance level (build) and the sequences (preditale,
correct).

**Q72** (maintainer, unsure, "I guess not"): `tell_tales()` does not
forward options to its internal AnnoTALE `analyze` calls. Its only
option, `n=`, is a run name, so nothing is lost.

Maintainer, 2026-10-04: **Q73** no `opt_param` on `correct_tales()`
(TALEcorrection has no options). **Q74** the heap flags go in `java_args`,
a string placed before `-jar`; defaults are today's flags (`""` for
predict and `correct_tales()`). **Q75** `opt_param` defaults to the
values used today, `""` where there are none; read as the jar defaults
written out (predict `"Sensitive=false"`, build `"c=5 s=0.01"`),
`preditale()` keeping `""` (its options hang off selections, see Q81).
**Q76** asked how hard a guard against keys the wrapper sets would be.
Tested: PrediTALE given `outdir=` twice silently takes the last one (the
first directory stays empty). The wrapper puts `opt_param` before its own
keys, so a user's `s=` or `outdir=` would be silently ignored.

Maintainer, 2026-10-04: **Q76** yes, the guard; **Q81** yes, the reading
above. **Done** the same day:
- `run_annotale_predict(opt_param = "Sensitive=false", java_args = "")`
  (`opt_param` goes to predict only, `java_args` to both stages);
  `run_annotale_build(opt_param = "c=5 s=0.01", java_args = "-Xms512M
  -Xmx6G")`; `preditale(java_args = "-Xms512M -Xmx2G")`, placed after
  `output_dir` so positional calls keep working; `correct_tales(java_args
  = "")`. Each `@param` lists the tool's options with their defaults.
  The analyze stage's jar path was not `shQuote()`d; it is now.
- `.check_jar_args()` (in `R/annotale.R`, shared by four wrappers): both
  strings must be single strings; keys of `opt_param` (a word at the
  start or after a space, followed by `=`) that the wrapper sets itself
  abort with `tantale_error_jar_args`, naming the R argument to use.
- Real runs (toy regions, sweet promoters): AnnoTALE's
  `protocol_predict.txt` echoes `Sensitive = true/false`; the sensitive
  scan ran out of heap under `-Xmx1G` (no limit: 22 s, same two TALEs as
  the default on the toy file); build `c=0`, `0.5` and `1` give two
  classes where the default gives one; PrediTALE `Strand="forward
  strand"` gives 20 sites, all `+` and all among the default's 24 `+`
  (the default also gives 7 `-`); `-XX:+NoSuchOption` stops each wrapper
  with its error class.
- Tests: new `test_jar_args.R` (the guard, no Java needed); additions to
  `test_annotale.R` (protocol echo, `c=0`), `test_target_predictions.R`
  (forward strand), `test_external_exit_status.R` (bad `java_args`, four
  wrappers). The five files: 32 tests, no failure. NEWS entry.

---

## 52. Readiness for rOpenSci review (2026-10-04) **[P]**

Maintainer, 2026-10-04: the JOSS paper (R9) waits; first make sure the
package would pass an rOpenSci review. §44 (pkgcheck, 2026-10-02) predates
R7, which changed how tools and genomes are installed, so its verdict is
stale.

Scope (devguide "Aims and Scope", read 2026-10-04): tantale fits
"scientific software wrappers" (non-trivial wrappers of field-specific
programs) and "data munging" (tool outputs into tidy objects). "Data
visualization packages are no longer in scope"; plotting as a secondary
feature is not the same thing, but a pre-submission inquiry settles it.

Candidates found before any audit (2026-10-04, by reading the code):
- `tell_tales()`, `run_annotale_predict()`, `run_annotale_build()` default
  `output_dir = getwd()`; `correct_tales()` defaults `corrected_path` to
  the working directory. Writing to the working directory by default is a
  CRAN policy breach that reviewers flag.
- Tests fail rather than skip without Java, conda or the downloaded tools.
  The rOpenSci bot runs pkgcheck in its own container, without them.
- No `inst/CITATION`, no `codemeta.json`; README badges: lifecycle
  "stable" at 0.9.x, coverage as a manual snapshot; CI is manual (§44 P7).

Maintainer, 2026-10-04: **Q77** §51 first, then a fresh pkgcheck on HEAD
and a pass through the devguide's author and reviewer checklists, every
finding recorded here with a recommendation before any fix. **Q78** a
pre-submission inquiry, posted by the maintainer: draft in
`dev/ropensci-presubmission.md` (the issue template of
`ropensci/software-review`, fetched 2026-10-04). It asks about the
plotting and distance/clustering functions, and about how the review
bot's pkgcheck should handle tests that need Java and conda. No R
package for TALE analysis was found (a web search, 2026-10-04).

Maintainer, 2026-10-04: **Q82** tantale is not to be advertised as a
wrapper package only: `tell_tales()`, `tales_align()`, the comparisons,
the grouping methods and the classes with their methods are its own
analysis code. The draft's scope answer now leads with them and asks the
editors directly whether that code fits; a DESCRIPTION rewrite in the
same spirit is proposed (Q85). **Q83** co-author handle `@vibaotram`.
**Q84** the generative-AI paragraph kept as written.
**Q82, Q85** (maintainer, 2026-10-04): DESCRIPTION's Title and
Description rewritten to present the whole package (discovery, classes,
alignment, distances, grouping, then the wrappers), program names in
single quotes as CRAN asks; `tools:::.check_package_description()` clean.

### pkgcheck on ac8c534 (2026-10-04)

pkgcheck 0.3.2, ctags/global from a scratch micromamba prefix, detached
worktree, 22 min. Passed: name available, contributing file, roxygen2,
URL/BugReports, HTML vignette, examples for all functions, website,
ORCIDs, no R CMD check warning, coverage 91.1% (90.7% in §44). Findings:
- **F1** (error) The golden `tell_tales()` correction test fails: the
  last fingerprint row (`tell_tales.log`) drops 3 lines instead of 2.
  Cause, by elimination: `correction_ref` is `test_path(...)`, and
  pkgcheck's R CMD check runs the tests from a directory under the test
  session's `tempdir()`, so the tempdir rule drops that log line. §44's
  P4 fix rewrites only the installation directory (`system.file()`), and
  was verified by a simulation, not a pkgcheck rerun; the uncorrected
  golden test, whose inputs all come from `system.file()`, now passes.
  Proposal: `.telltale_file_digest()` rewrites the tests directory the
  same way, plus a helper test; verified by a second pkgcheck run.
- **F2** (fail) No continuous integration: the workflow runs on manual
  dispatch only (§41) and pkgcheck cannot query a private repository
  (§44 P7). rOpenSci expects CI on push and pull requests, with a badge.
- **F3** (note) Hidden `.git`: an artefact of checking a worktree (§44).
- **F4** (goodpractice) Packages imported whole: `import(Biostrings)`,
  `import(cli)`, `import(fs)`, `import(magrittr)` in NAMESPACE. Measured
  with `codetools::findGlobals()` over every package function: no
  unqualified call to `fs` or `Biostrings` (its `start`/`end` hits are
  dplyr column names), one bare `cli_inform()`, and `magrittr`'s `%>%`
  and `%<>%`.
- **F5** (lintr, 1396 lints) Mostly style: 1100 long lines, 88
  `expect_equal()` where `expect_identical()` is preferred, 40 implicit
  assignments, 26 `:::` (tests reaching internals). Checked for bugs:
  the two "missing argument" lints are trailing commas inside `mutate()`
  and `geom_tile()` (harmless); the two "unreachable code" lints are
  `next()` as the last statement of an `else` (false positive); the 20
  "duplicate arguments" are cli bullet vectors (§44). Small real
  cleanups: 10 `%<>%`, `ifelse(x, TRUE, FALSE)` (classification.R 464),
  `any(duplicated())` (2), `any(!x)` (3), two `<<-` in `R/` (closures in
  `.tales_anomalies()` and the `tryCatch()` handler of
  `.telltale_run_annotale()`).
- **F6** (goodpractice) Long or complex functions: cyclomatic complexity
  `plot.tales_msa()` 43, `talomes_heatmap()` 34, `tantale_setup()` 29,
  `plot.tales()` 21, `tales_group_kmedoids()` 16, `summary.tales()` 15;
  over 150 lines: `plot.tales_msa()` 370, `talomes_heatmap()` 263,
  `.build_repeat_msa()` 212, `tell_tales()` 189, `plot_target_preds()`
  178. `plot.tales_msa()`'s matrix helpers are §21, reserved for the
  maintainer.
- **F7** 40 Imports (> 99th percentile): kept, §44 P5.
- **F8** "Unused" internals: the same four false positives as §44 P6.
- **F9** Package size: 4636 lines of R (97th percentile), 217 functions
  (98th). Informational.

Maintainer, 2026-10-04: all recommendations agreed. **Q87** fix F1 and
confirm with a second pkgcheck run; **Q88** CI on push to `main` and on
pull requests (Linux; macOS stays manual), badge in the README; **Q89**
the manual coverage badge removed; **Q90** F4 as proposed; **Q91** the
small F5 cleanups only; **Q92** F6 left until a reviewer asks; **Q93**
the devguide checklist pass next, findings recorded here.

**Done 2026-10-04** (Q87-Q91):
- F1: `.telltale_file_digest()` rewrites the test data directory to
  `/tests` before the tempdir rule; new test in `test_golden.R` (fails
  without the rewrite: the line is dropped).
- F2: `R-CMD-check.yaml` also runs on push to `main` and on pull
  requests (Linux; `paths-ignore` for `docs/`, `dev/`, `pkgdown/`,
  `.claude/`), a newer push cancels the older run; README CI badge.
- Q89: the manual coverage badge and its note removed from the README.
- F4: `import()` of `fs`, `cli`, `magrittr` and `Biostrings` replaced:
  `importFrom(magrittr, "%>%")`; `cli::cli_inform()` in `startup.R`;
  `importFrom(Biostrings, ...)` with the 19 base-named generics the code
  calls unqualified, so name resolution is unchanged (the comment in
  `R/tantale.R` says to extend the list for a new unqualified call). The
  full suite caught one dependency codetools cannot see: `{qty()}`
  inside two cli strings of `.tale_parts()` (cli evaluates them in the
  caller's frame); now `{cli::qty()}`. A scan of every glue/cli string in
  `R/` for `fs`, `Biostrings` or `magrittr` calls found nothing else.
- F5: the 18 `%<>%` written as `x <- x %>% ...`; `grepl()` in place of
  `ifelse(grepl(), TRUE, FALSE)`; `anyDuplicated() > 0L` (2); `!all()` (3).
  The two `<<-` in `R/` left as they are (closures over a local list).
- Full suite after F4/F5: 462 tests, 0 failures, 3 errors, all the `qty`
  one; the two files rerun clean after the fix.
- URLs (`urlchecker`, run because pkgcheck skipped it): redirects fixed
  (site URL without its trailing slash in DESCRIPTION, README, `?tantale`
  and the "Get started" page; the old PubMed URL in the AnnoTALE wrappers,
  now the DOI; jstacs.de), and the README's link to the prototype bundle,
  which lacked the file name (404). Two doi.org links answer 403 to the
  checker (publisher bot blocking), left.
- **Second pkgcheck run, on 06cb505** (21 min): R CMD check no error,
  no warning, one note (the worktree's `.git`); coverage 91.15%; the
  whole-package import check passes; 1388 lints (1396 before). Still
  flagged: no CI (pkgcheck cannot see the workflow of a private
  repository; the README badge is in place), 40 Imports, the four
  "unused" internals (F7, F8).

### Devguide checklist pass (2026-10-04, Q93)

Read: the author checklist of the submission template, the reviewer
template, the packaging guide (`pkg_building`) and the CI chapter.
Already met: Title Case title without "in R", program names quoted
(Q85), roxygen2, examples and `@return` for every export (two internal
Rd topics, `.tales_check_crd_unique` and `.tales_report_anomalies`,
have no `\value`), `?tantale`, an HTML vignette, pkgdown, `@family`,
NEWS, CONTRIBUTING, URL/BugReports/Maintainer, ORCIDs, `main` branch,
cli conditions, object_verb naming, coverage 91%. Findings:
- **C1** "MUST NOT use startup messages except when necessary": the
  `.onAttach()` greeting ("Attaching the tantale package", "Email ... for
  comments"). Recommend removing it.
- **C2** "MUST indicate system dependencies in DESCRIPTION": no
  `SystemRequirements`. Java (the jars are class version 52, Java 8),
  Perl, conda/mamba/micromamba. Recommend adding the field.
- **C3** "MUST run on Windows, macOS and Linux", exceptions for documented
  reasons: Windows impossible (no bioconda builds of the pins; `OS_type:
  unix`), macOS never run, README says "written with only Linux in mind".
  Recommend a macOS CI run once public (free for public repositories) and
  a README sentence giving the Windows reason; the inquiry already asks.
- **C4** "CI must also be linked to a code coverage service", coverage
  badge: Codecov (preferred). Needs the maintainer to enable the Codecov
  GitHub app for the repository; then a coverage workflow (with the
  tools installed, as in R-CMD-check) and a badge. Q89's removal of the
  manual badge stands.
- **C5** CI on R release, previous and devel. Recommend a Linux matrix.
- **C6** Examples run in CI with `--run-donttest` (recommended): the CI
  has the tools, so add it to the check arguments.
- **C7** CITATION (should, bibentry): none. Recommend `inst/CITATION` for
  the package, with a header asking to cite the wrapped tools' papers.
- **C8** "MUST acknowledge bundled code authors with cph role": TALVEZ
  3.2 and the FuncTAL table are by Alvaro L. Pérez-Quintero; TALVEZ's
  `simplescancode/` is Java code credited to "Matzieu and Hatzigeorgiu
  2010" in TALVEZ's header. Recommend `cph` entries in `Authors@R`.
- **C9** README "MUST include ... usage": no usage example; also
  recommended: repostatus badge, comparison with similar software,
  citation guidance, code of conduct. Items reviewers will query: the
  "Connector with daTALbase (to be done)" bullet, "contributing: send me
  an email" next to the CONTRIBUTING file, "only Linux in mind". The README
  is the maintainer's; the site home is `pkgdown/index.md`.
- **C10** "MUST identify data sources (organisation, access URL)":
  BAI3-1-1 is "an unpublished assembly of a BAI3 derivative" in
  `inst/COPYRIGHTS`; who produced it is not stated. The maintainer has
  the facts.
- **C11** Console messages "SHOULD provide an opt-out (option or
  environment variable)" and verbosity levels: none beyond
  `suppressMessages()`. Recommend deferring until a reviewer asks.
- **C12** "MUST NOT use print() or cat() except in print methods":
  `.tales_group_hclust` prints a tree when `plot_tree = TRUE`;
  `plot.tales_msa()` prints a patchwork copy (Q66 kept it). Both print
  plots, not text; recommend leaving and explaining if asked.
- **C13** Writing to the working directory by default (CRAN policy):
  `tell_tales()`, `run_annotale_predict()`, `run_annotale_build()`
  (`output_dir = getwd()`), `correct_tales()` (`corrected_path` in
  `getwd()`). Removing the defaults breaks calls that rely on them.
- **C14** Tests "SHOULD skip_on_cran() for external resources" against
  the standing fail-not-skip rule: the inquiry asks the editors.
- **C15** Code style: Air or styler recommended, not required (1100 long
  lines, F5). Leave.

Maintainer, 2026-10-04 (Q94-Q103; the repository is public since the
night of 2026-10-04): all yes. Q97: Codecov enabled by the maintainer
(`CODECOV_TOKEN` secret present). Q99: no tantale paper to cite yet.
Q100: the maintainer believed `simplescancode/` had been recoded in R;
it has not (ARLEM was, §33): `talvez()` runs TALVEZ_3.2.pl, which calls
`java simplescancode.Background` and `simplescancode.Scan`, so the cph
question for its authors stays open (Q104). Q102: BAI3-1-1 was produced
by the package authors and is available only in tantale. **Done:**
- C1 `R/startup.R` removed (the greeting); dev/CLAUDE.md's cli rule
  updated.
- C2 `SystemRequirements: Java (>= 8), Perl, and conda, mamba or
  micromamba` (the jars are class version 52).
- C4 `.github/workflows/test-coverage.yaml` (covr with the tools
  installed as in R-CMD-check, upload with `CODECOV_TOKEN`); Codecov
  badge in the README draft.
- C5/C6 R-CMD-check matrix over R release, oldrel-1 and devel;
  `--run-donttest`.
- C7 `inst/CITATION`: a Manual entry for the package, footer with the
  DOIs of AnnoTALE, TALEcorrection, PrediTALE, TALVEZ, DisTAL/FuncTAL.
- C8 Alvaro L. Pérez-Quintero as `cph` in `Authors@R`.
- C10 BAI3-1-1's origin in `inst/COPYRIGHTS` and `?tantale_genome`.
- C13 `output_dir` required in `tell_tales()`, `run_annotale_predict()`,
  `run_annotale_build()` (and the internal `.run_annotale_analyze()`),
  `corrected_path` in `correct_tales()`, which checks it before running
  anything (`tantale_error_missing_output`, tested). No call in tests,
  examples or articles relied on the defaults. NEWS entries.
- Checks: targeted tests (error_conditions, correct_tales, annotale,
  jar_args, golden: 40 tests, 0 failures); `rcmdcheck` without tests or
  vignettes: 0 errors, 0 warnings, 0 notes.
- C3/C9 README draft left uncommitted for the maintainer: repostatus and
  Codecov badges, platforms paragraph (Linux and Intel macOS; why not
  Windows), a Usage section (run for real: 31 s, 9 MAI1 arrays, TALVEZ
  predictions on the SWEET promoters), Related software, Citation, the
  daTALbase "to be done" bullet removed, contributing points to
  CONTRIBUTING and the issues.
- Maintainer, 2026-10-04: asks whether README should be generated from a
  README.qmd ("discussed before and concluded not"; no record of that
  decision in the ledger or the history). Q105.

Maintainer, 2026-10-04: **Q104** could not find the reference. Found it:
"Matzieu" is Megraw. PlantTFBS (Molly Megraw, Artemis G. Hatzigeorgiou;
chapter "MicroRNA Promoter Analysis", Methods Mol Biol, Plant MicroRNAs;
page https://megraw.cgrb.oregonstate.edu/node/714) is "freely available
under the GNU Public License". Both added as `cph`; `inst/COPYRIGHTS` and
the README licence section say so. Open: the upstream zip now returns an
HTML page, so the bundled classes were not compared with it, and tantale
ships GPL class files without their source (Q111). **Q105** README.qmd
adopted: `format: gfm` with `wrap: preserve`, install blocks not run,
the Usage example run at each render (32 s; figure in
`man/figures/README-usage-plot-1.png`); README.qmd, README_files and
.quarto build-ignored. **Q106** draft approved. **Q107** push.
C13 revised (maintainer's proposal, **Q108** go ahead, **Q109** no
message: the docs warn, **Q110** yes): `output_dir`/`corrected_path`
default to `tempfile(...)` under `tempdir()`; `run_annotale_predict()`
and `run_annotale_build()` return `output_dir` invisibly (was `0`), as
`tell_tales()` and `correct_tales()` already returned their output, so
`tales_from_telltales(tell_tales(genome))` works; the explicit
`corrected_path` check and its test removed; new test of the default in
`test_annotale.R`; README example shortened; NEWS rewritten. 57 tests
in the affected files, 0 failures.
**Q111** dropped (maintainer: enough has been done on PlantTFBS). **Q112**
version 0.99.0 (the §43 release candidate: the repository is public),
full site rebuild, push. **Q113** an archive with a DOI will come from
DataSud, not Zenodo.

---

## 53. Correction strategies inside `tell_tales()` (future release) **[P]**

Maintainer, 2026-10-04, for a release after this one: prototype new
correction methods in `tell_tales()`:
1. `correct_tales()` on the whole input up front, then the usual search.
2. `correct_tales()` on each candidate ORF separately (not known whether
   TALEcorrection works on a single ORF).
3. In either, an ORF that is still broken after `correct_tales()` goes to
   the DECIPHER correction (today's `correct_array = TRUE`).

Evidence that motivates the hybrid (§49, truncTALE article): on BAI3-1-1
`correct_tales()` left the N-termini of `ROI_00006` and `ROI_00009`
unmatched (24 aa) where DECIPHER repaired them; on PXO86 DECIPHER
over-corrected the truncTALE `ROI_00001` (C-terminus 183 -> 216 aa) and
`correct_tales()` left it alone. Each method fixes what the other gets
wrong on these genomes.

Open before any code: a definition of "broken" for step 3. Candidates
already computed per array: a missing or unmatched terminus
(`*_aa_profile_gap`, §42), ORF coverage, AnnoTALE failing to parse the
ORF (§49 Q59), the non-standard structures `tales_anomalies()` reports.
A rule that sends a genuine truncTALE to DECIPHER reproduces the
over-correction, so the definition has to be checked against
`ROI_00001`/`ROI_00019` of PXO86. Also to settle: how the output records
which correction each array received.

**Constraint (maintainer, 2026-10-04): the change must not break the
current API.** Every existing call keeps working with the same result:
`correct_array = FALSE` (no correction) and `correct_array = TRUE`
(DECIPHER) keep their meaning, and the new methods are reached only
through something added. Two additive shapes: `correct_array` also
accepting method names (`TRUE` staying an alias of the DECIPHER one), or
a new argument whose default reproduces today's behaviour. The
`tell_tales()` arguments are flat (§43), so either fits. After 1.0.0 any
later rename follows the lifecycle rules (dev/CLAUDE.md, API
conventions).

---

## 54. TALVEZ in R (future release) **[P]**

Maintainer, 2026-10-04: plan a port of TALVEZ to R for a future release.
universalmotif, already imported (`tales_to_universalmotif()`), makes it
more feasible than before. It would remove the bundled Perl script, the
Perl dependency of the conda environment, and the PlantTFBS Java classes
whose source tantale does not ship (§52 Q104, Q111). Starting points: the
two RVD-DNA matrices `mat1`/`mat2`, TALVEZ's position correction (`-l`)
and scoring threshold (`-t`, `-a`), and its background model from
`simplescancode.Background`; validate against `talvez()` on the SWEET
promoters, as ARLEM was validated against its executable (§33).

**Q111 settled (maintainer, 2026-10-05):** no action on the PlantTFBS
classes; if someone raises it, the answer is this port, which removes
them.

---

## 55. Updating the RVD-to-nucleotide tables (2026-10-04) **[P]**

Maintainer, 2026-10-04: look for published data, especially from the
TALEN literature, that could update the RVD/base association tables.
Search done (Q116); nothing changed, the maintainer judges the sources.

**What tantale ships today.** Two tables, both from Pérez-Quintero's work
and unchanged since 2015:
- `inst/tools/TALVEZ_3.2/mat1` and `mat2`, used by `talvez()`: 17 rows
  (HA NA HD ND HG IG NG HI NI NK HN NN NS H* N* OO XX), four weights each.
- `rvd_dna_specificity` (QueTAL FuncTAL, the `rvd_sim` fill of
  `plot.tales_msa()` and `tales_align()`): 404 rows, all amino-acid pairs
  plus `H*`, `N*`, `OO`, `XX`. 86 rows are uniform (no information), and
  only 58 have a value >= 30; the median row maximum is 3, so most of the
  404 carry no real preference.

**Gaps measured, not assumed** (RVD counts over the shipped RVD sets:
`Sample_TALEs_RVDSeqs_AnnoTALE.fasta`, the SWEET control TALEs,
`TalA_RVDSeqs_AnnoTALE.fasta`; 10 distinct RVDs):
- **`NV` occurs in real TALEs (4 of 957 RVDs) and is absent from
  `mat1`/`mat2`.** TALVEZ_3.2.pl lines 143 and 147 fall back to
  `$RVDs{XX}`, which is `1 1 1 1`, so `talvez()` scores `NV` as entirely
  non-specific. `rvd_dna_specificity` does give `NV` a maximum of 50.
- **`N*` is the fifth most frequent RVD (4.7%)** and its
  `rvd_dna_specificity` row is `1 4 1 4` (maximum 4), among the weakest
  in the table, while `mat1` gives it `1 4 1 3`.

**Candidate sources found** (provenance only; none assessed as correct,
which is the maintainer's call, see
`feedback_tantale_biology_ground_truth`):
- **PrediTALE's own model, already shipped in `PrediTALE.jar`**:
  `projects/tals/prediction/preditale_quantitative_PBM.xml` (132 KB).
  Its RVD alphabet has 420 symbols (every amino-acid pair) and an
  `indexMap` assigning all 420 to **6 parameter groups**, with 6 RVDs
  given separate parameters: the rare RVDs inherit the specificity of
  common ones sharing the 13th amino acid, which is TALVEZ's idea fitted
  to quantitative data. Learned in the *Xanthomonas*/rice system (Erkes
  A. et al. 2019, PLoS Comput Biol 15, e1007206,
  doi:10.1371/journal.pcbi.1007206; CC-BY; source in
  github.com/Jstacs/Jstacs). Closest provenance to tantale's own use, and
  already redistributable. Reading the parameters out needs the Jstacs
  encoding.
- **Streubel J. et al. (2012), Nat Biotechnol 30, 593-595**
  (doi:10.1038/nbt.2304), "TAL effector RVD specificities and
  efficiencies": plant reporter assays, introduces `NH` as
  guanine-specific and reports that `NK` performs worse than `NN`. Plant
  system, so the closest assay context.
- **Miller J.C. et al. (2015), Nat Methods 12, 465-471**
  (doi:10.1038/nmeth.3330), expanded RVD repertoire for TALENs: many
  non-canonical RVDs, measured in human cells. Broadest RVD coverage
  found; the system differs from a TALE activating a plant promoter.
- **Rogers J.M. et al. (2015), Nat Commun 6, 7440**
  (doi:10.1038/ncomms8440), protein-binding microarrays, CC-BY, GEO
  GSE56978: quantitative, but only the four common RVDs (NI, HD, NN,
  NG); its value is neighbour context, not new RVDs. It is the PBM data
  PrediTALE was trained on.

**Open.** Which table to improve first (the target-prediction weights
used by `talvez()`, or the `rvd_sim` similarity used for alignment and
distances: the second is a different question from base preference);
whether a TALEN assay in human cells transfers to a TALE activating a
rice promoter; and whether adding `NV` to `mat1`/`mat2`, rather than
leaving it at `XX`, is worth doing on its own. Note that `preditale()`
already applies a data-derived model, so a user wanting current
specificities has that route today.

### How the two matrices were built, and what PrediTALE's model says (Q117, Q118)

**`rvd_sim` is derived from TALVEZ, not independent** (maintainer,
2026-10-04; `data-raw/sysdata.R`): `rvdSimDf` is the Spearman correlation
between the rows of `mat1`. So the alignment/distance similarity and the
target-prediction weights are one table seen two ways, and `rvdSimDf`
covers the same 17 RVDs. A change to `mat1` propagates to both.

**TALVEZ (Pérez-Quintero A.L. et al. 2013, PLoS ONE 8, e68464).** `mat1`
is Moscou and Bogdanove's counts matrix modified by hand, assigning the
same vector to RVDs sharing the 13th amino acid (the residue that
contacts the base), with position 0 adjusted to allow cytosine. `mat2` is
a scaled-down version "adopted after preliminary tests against other
possible options". There is no separate training set: the same 72
curated TALE-EBE pairs (35 positive, 37 negative) were used to choose the
position-correction value and to measure performance (ROC, plus
whole-genome ranking of known targets and two rice microarray series).
The paper does not discuss `XX` or rare RVDs; it notes that five control
pairs rank very low because the EBE "deviates from the canonical
RVD-nucleotide association model", "not currently well understood".
This supports the maintainer's caution (Q117): the matrices are tuned,
not fitted.
- Checked against the file: the "same 13th residue, same vector" rule
  holds only for D (HD/ND) and G (HG/IG/NG). `HA`/`NA`, `HI`/`NI`,
  `HN`/`NN` and `H*`/`N*` differ, the N-prefixed ones being much
  stronger, so `mat1` also carries repeat *efficiency*, not base
  preference alone.
- **`NV`'s `XX` fallback follows TALVEZ's own scheme rather than being an
  oversight**: the 13th residues covered by `mat1` are `* A D G I K N O S
  X`; `V` is absent, so `NV` has no sibling to inherit from.

**PrediTALE (Erkes A. et al. 2019).** Trained on pooled quantitative
experiments, including the PBM data of Rogers et al. (21 TALEs built from
NI, HD, NN, NG). Validated on *independent* data: RNA-seq of rice after
infection by three Xoo and ten Xoc strains, taking up-regulated genes as
ground truth, benchmarked against Target Finder, Talvez and TALgetter.
Rare RVDs again inherit from the 13th amino acid. Stated limit: position-0
dependencies are fitted only for R0 = {HD, NN, NG, NI, NS}. So PrediTALE
is the one of the two with validation on data it was not tuned on.

**Q118 done.** `dev/preditale-specificities.py` reads the model out of
the jar (`projects/tals/prediction/preditale_quantitative_PBM.xml`):
`specsThirteen` is 21 x 4 keyed on the 13th residue, `separateSpecs` is
6 x 4 for HD, NN, NG, HG, NI, NK, both unnormalised log-scores. (The
sibling `a`/`b` vectors, length 905, are optimiser state; reading them as
specificities gives HD -> G and other values that contradict the biology,
which is how the error announces itself.) Validated: D->C, G->T, I->A,
K->G, N->G; HD->C, NG->T, NI->A, NN->G, NK->G. The CSV is regenerated by
the script rather than committed: it derives from a GPL-3 work, and
tantale is MIT.

Comparison of the preferred base, TALVEZ against PrediTALE's 13th-residue
table, over the 17 RVDs of `mat1`: 11 agree, 6 differ.
- `HA`, `NA`: TALVEZ C, PrediTALE T.
- `HN`, `NN`: TALVEZ A, PrediTALE G. For `NN` this is barely a
  disagreement: `mat1` gives `30 10 30 1`, A and G tied, and `which.max`
  takes the first.
- `NS`: TALVEZ A (`20 7 10 7`), PrediTALE G; `NS` is promiscuous in both.
- `H*`: TALVEZ T, PrediTALE C.
- **`NV`**: PrediTALE's 13th-residue row for `V` is nearly uniform
  (A 0.309, C 0.252, G 0.243, T 0.196 against a 0.25 baseline), so its
  data gives little ground for a confident `NV` vector either. Leaving
  `NV` at the flat `XX` is defensible on the evidence found.

---

## 56. The R version matrix finds two failures (2026-10-04) **[P]**

The CI matrix added in §52 (C5) ran for the first time on 0.99.0
(bab07a7, run 37213179171). **R release passed; R oldrel-1 and R devel
failed.** Neither is a defect in tantale's logic, and both need a
decision before the rOpenSci submission.

**R devel: `merge()` on an S4 `DataFrame`.** 19 test failures and
`checking examples ... ERROR`, all one cause:

    Error in `.local(x, row.names, optional, ...)`: unused argument (validRN = FALSE)
    .telltale_add_array_measures() -> base::merge() -> merge.default ->
      base::merge(as.data.frame(x), as.data.frame(y), ...)

R devel's `merge.data.frame` now calls `as.data.frame(y, validRN =
FALSE)`; the S4Vectors method for `DataFrame` does not accept that
argument. Triggered by `R/telltale.R:1047` and `:1051`, which call
`merge()` directly on `S4Vectors::mcols(by_array)`. Upstream
incompatibility between R devel and S4Vectors, which tantale meets by
merging a `DataFrame` rather than a `data.frame`.
Options: (a) convert with `as.data.frame()` before merging and put the
result back into `mcols()` -- contained, and removes the dependency on
how S4Vectors tracks base's signature; (b) use a dplyr join on a plain
data frame; (c) wait for S4Vectors to catch up and keep devel red.
Prefer (a).

**Devel fixed 2026-10-04 (Q138), option (a).** Three call sites, not two:
`.telltale_add_array_measures()` (two `merge()` calls) and the
correction path in `.telltale_array_orfs()` (the indel counts), all
passing `S4Vectors::mcols(by_array)` straight to `merge()`. Each now
wraps it in `as.data.frame()`; `mcols<-` coerces back. Behaviour-neutral
on R release: golden 18/18 unchanged, `test_tell_tales.R`,
`test_tell_tales_correction.R`, `test_tell_tales_guards.R` 31/31.
Confirmation on devel needs a CI run.

**Oldrel solved and fixed 2026-10-04 (Q139, Q145).** A throwaway
workflow on branch `ci-oldrel-gff` printed the file under both versions:
of 199 lines **exactly one differed**, line 2, `##source-version
rtracklayer 1.72.0` against `1.70.1`. rtracklayer stamps its own version
into the GFF header, which is a record of how the run happened, like
`##date` and HMMER's `# Version:` already dropped beside it. Added
`^##source-version` to `.run_specific_pattern()`.
Re-baselined, every changed row explained first: **two files in each of
the two fingerprints**, `all_ranges.gff` and `hits_report.gff` (both
rtracklayer GFFs), each `n_dropped` 1 -> 2 with a new digest; all 173
snapshot vectors were compared position by position, and nothing else
moved. Checking positionally rather than reading the diff is what caught
the second file: the first account of this change named only
`all_ranges.gff`. Full suite after: 470 tests, 0 failures.

**R oldrel-1: the golden `all_ranges.gff` digest.** 2 failures
(`test_golden.R:142` and `:249`), both row 1 of the fingerprint,
`all_ranges.gff`: the digest differs while `n_lines` (199) and
`n_dropped` (1) are unchanged, so the file has the same shape and
different content, most likely rtracklayer/GenomeInfoDb writing the GFF
differently on the older stack. The baseline is therefore tied to the
package versions that recorded it, not only to the machine.
Options: (a) find the differing line (download the oldrel check
artefact, or run an oldrel container) and extend the normalisation if it
is a version or formatting artefact; (b) restrict the golden tests to
the release version; (c) drop `all_ranges.gff` from the fingerprint.
Prefer (a): the point of the baseline is that a content change is
explained before it is accepted (dev/CLAUDE.md), and the artefact is
attached to the run.

### Third pkgcheck run, on 9e82f88 / 0.99.0 (2026-10-04)

21 min. **R CMD check: no errors, no warnings** (the §52 F1 fix holds).
Coverage 91.1%. Everything else passes except the two already settled:
40 Imports (§44 P5) and the goodpractice lints.

**The "no continuous integration" flag is an artefact of running
pkgcheck locally, not a finding.** `pkgcheck:::pkgchk_ci_badges()`
downloads the README from the *remote* default branch and returns NULL
early unless `gh::gh_token()` is non-empty; the background runs here had
no token, so it reported none. Re-run with `GITHUB_PAT=$(gh auth token)`
it finds the badge:
`[![R-CMD-check.yaml](.../R-CMD-check.yaml/badge.svg)](...)`. rOpenSci's
own bot runs authenticated, so it will see CI. Set `GITHUB_PAT` for any
future local pkgcheck run, or discount that one line.

### Why DESCRIPTION keeps `biocViews: Software` (2026-10-04)

Maintainer asked whether to remove it. Keep it: `remotes` (and so
`devtools::install_github()`) decides whether to add the Bioconductor
repositories by `is_bioconductor <- function(x) !is.null(x$biocviews)`
(r-lib/remotes `R/utils.R`; used in `R/deps.R`, lines 137-144). Without
the field, the ~15 Bioconductor Imports would not resolve for a user
without Bioconductor repositories configured. `pak`, which the README
uses, resolves them either way. R CMD check, CRAN and rOpenSci ignore the
field. A Bioconductor submission would need more specific biocViews terms
than `Software`.

### Pre-submission inquiry posted (2026-10-04)

The maintainer posted it as ropensci/software-review#813 (16:55 UTC,
label `0/presubmission`). The bot acknowledged it and listed the "most
similar" rOpenSci and CRAN packages from `pkgmatch`; none is TALE-related
(geojson, rixpress, stats19, tepr...), which is consistent with the
inquiry's "no other R package for TALE analysis". Editor reply expected
within 5 business days.

---

## 57. `tale_annotations`, a curated reference talome (2026-10-04) **[V]**

Maintainer, 2026-10-04: Bao Tram Vi's curated table of the TALEs of ten
published *X. oryzae* genomes, wanted in the package "if only to have a
safe place for it". Source file
`.../Tram/former_xanthopore-scripts/talomes_analysis/VXO_genomes/tal_analysis/ref_tale_annotation.tsv`,
copied unchanged to `data-raw/ref_tale_annotation.tsv`. Note the sibling
`all_tale_annotation.tsv` (305 rows, 20 strains) **must not** be used: it
includes ten unpublished VXO genomes.

Decisions (Q125-Q128): **Q125** first public appearance, Tram's
unpublished work, so `@source` credits her and points at the GenBank
assemblies each row names. **Q126** ship the gaps as they are (19 rows
without `tal_name`, 54 without `annotale_group`, 19 without `pubmed`;
the maintainer filled some PubMed ids first, BAI3 among them). **Q127** a
plain tibble, not a `tales` object. **Q128** `distal_group` dropped: empty
in every row here, and how the sibling file's values were computed is not
recorded, so they cannot be vouched for; recomputing is
`tales_compare_distal()` + `tales_group_hclust()`.

Maintainer, same day: `telltale2_id` and `annotale_id` dropped too. Both
name an array within one run (`ROI_00001`, `tempTALE1`) and neither
survives a rerun. `annotale_group` is kept, being AnnoTALE's own class
names (`TalAH30`). `strain` + `label` is the key: unique over all 128
rows, never missing. `truncTALE` **keeps its spelling** (maintainer,
2026-10-04: "no rename of truncTALE outside of R code; this is the
standard way to call them, with iTALEs"). The snake_case rule governs R
code; a term the literature writes `truncTALE` is not R code.
`talomes_heatmap()`'s argument `trunc_tales_col` stays snake_case, being
an argument name.

Result: 128 rows x 10 columns, 2.8 KB. `data-raw/tale_annotations.R`
builds it with the invariants asserted; `?tale_annotations` documents
every column; `test_tale_annotations.R` has 6 tests; NEWS entry;
`@family TALE discovery` puts it on the site with no `_pkgdown.yml`
change. `check_pkgdown()` clean, `rcmdcheck` (no tests, no vignettes)
0/0/0.

**The table feeds `talomes_heatmap()` with no reshaping**, verified by
drawing it: `talomes_heatmap(subset(tale_annotations,
!is.na(annotale_group)), group_col = "annotale_group", strain_col =
"strain", rvd_col = "rvd_seq", trunc_tales_col = "truncTALE")` gives the
four Asian strains that carry AnnoTALE classes (PXO61, PXO71, PXO86,
PXO99A; 74 rows, 31 groups) clustered group by group, truncTALEs marked
"T". That call is now an example in `?tale_annotations`, and the Rd
examples were run to confirm. The African strains have no
`annotale_group`, which is where the 54 gaps are.

**Found while testing: the RVD strings carry lowercase RVDs** (`hd`,
`ng`, `ni`, `nn`, `ns`, `n*`; 17 occurrences over 16 of the 128 arrays).
Lowercase is the convention for a repeat of non-standard length, which
tantale already reads that way (`grepl("[a-z]", rvdseq)` in
`classification.R`, and §49's 32-aa `nv` repeat). A first test asserted
`^[A-Z][A-Z*]$` and failed on real data, not the other way round. Both
the test and `?tale_annotations` now say case carries meaning and must
not be normalised away. Only 4 of those 16 arrays also carry an
`unusual_feature` note, so the two records of "atypical" are independent:
the note is the curator's, the lowercase is the annotation pipeline's.

### Filling the missing AnnoTALE classes (Q130-Q133, 2026-10-04) **[A]**

The 54 rows without `annotale_group` are the six African strains
(AXO1947, BAI3, MAI1, MAI68, MAI129, MAI134), nine TALEs each; the four
Asian ones (PXO61, PXO71, PXO86, PXO99A) have all 74. Maintainer asked
whether a script could fetch those genomes and compute the missing
classes.

**It can, with a tool tantale does not yet wrap.** `AnnoTALEcli-1.5.jar`
holds eleven tools, not the three tantale calls. Two matter here, and
both are already in the installed jar (`ClassAssignmentTool.class`,
`LoadAndViewClassesTool.class`), so wrapping them adds no download, no
licence question and no change to the `tools-1` release:
- **`loadAndView`** downloads AnnoTALE's curated catalogue of TALE
  classes (the shared `TalAA`, `TalAB`... naming) and writes it as a
  Class builder XML. Slow: several minutes, since it rebuilds the class
  definitions after downloading. A wrapper should cache the XML.
- **`assign`** places a set of TALEs into those existing classes, opening
  new ones for TALEs that fit none, and proposes systematic names from
  the strain name given. This is what produced `TalAH30` and the rest.
  Note `build`, which tantale does wrap, builds classes *de novo*
  (`Class_1`, `Class_2`) and cannot produce catalogue names.

Pipeline: `loadAndView` once -> per strain, fetch the genome by accession
-> `run_annotale_predict()` -> `assign` with the strain name.
Genome fetching is available: `biomartr` is already in Suggests and
NCBI's datasets API answers for these accessions.

Decisions: **Q130** recompute all 128, not just the 54, so the column has
one class-definition vintage; the catalogue grows over time and Tram's
assignments date from 2021, so mixing vintages would make the column
mean two different things. **Q131** the script reports where a recomputed
value disagrees with Tram's and changes nothing without the maintainer's
sign-off. **Q133** match rows on `rvd_seq`, verifying the match rate first
on the four strains that already have classes. Maintainer: that is how
this was done with Tram; a duplicated TALE gene matches twice but both
copies fall in the same class, so either serves. The case to watch is a
near-identical pair differing by a repeat, which could land in different
classes -- the table already carries such a note (an allele of
`tal7b`/`tal8b` with a deletion in the repeat region).

**Q132 open:** whether `assign` and `loadAndView` become exported
wrappers (`@family external TALE tools`, `opt_param`/`java_args` as §51
settled) or stay inside the `data-raw/` script. Recommended: wrap them,
since "which published class is this TALE" is a question tantale cannot
answer today. The interface is declared stable, so this is an addition,
not a change.

### The catalogue answers it without computing anything (2026-10-04)

`loadAndView` took **17 minutes** and wrote **440 MB**: a 382 MB
`Class_builder_download.xml`, 429 per-class PDFs, a tree of classes, and
plain-text lists. The XML is readable in principle (`xml2`) but is Jstacs
serialisation holding 8,866 `<TALE>`, 168,050 `<Repeat>` and **262,960
cached pairwise `<StringAlignment>` with costs** -- the saved state
`assign` compares against, which is why it is that size and that slow. It
is not where a human reads class membership.

`Lists_of_classes,_strains_and_TALEs/List_of_classes.txt` (7 MB, plain
text) is: per class, an aligned RVD matrix, one row per member, each
tagged with its TALE id and strain. Parsed: **12,101 entries over 433
classes**. So the classes come out of a text file, with no genome fetch,
no `predict` and no `assign` -- the whole Q130 pipeline collapses.

**Matching (Q133) works.** On `strain` + RVD string, case-insensitive:
**126 of 128 rows match.** The two that do not are explainable: PXO99A
`Tal7b` is the 5-repeat allele whose `unusual_feature` records "a
deletion in the central repeat region", too degenerate to be catalogued;
MAI1 `TalH` (12 repeats) matches an entry the catalogue attributes to
another strain. Matching on RVD alone instead gives 127 but picks the
wrong strain's copy where a TALE is shared, so the strain constraint
matters.

**What `annotale_group` actually encodes, and why it cannot be stored.**
Of the 73 rows carrying a value, the **class letters agree 73/73** and
the **trailing number differs 73/73**. The number is a member index
within the class, running 1..n: class TalAH has 116 members numbered
1-116, and 281 of 432 classes number theirs exactly 1..n. So `TalAH30` is
"the 30th member of TalAH **in the catalogue of that day**". The
catalogue grows, members renumber, and Tram's 2021 numbers no longer
hold: `TalAH30` is `TalAH25` today, `TalBH12` is `TalBH1`. Storing the
number in a published dataset stores something already wrong and
guaranteed to drift; the class letters are the stable part and they
agree perfectly.

**Done (Q140-Q142).** `annotale_group` replaced by `annotale_class`,
holding the class alone, derived for all 128 rows by matching `strain` +
upper-cased `rvd_seq` against the catalogue. 126 carry one; the two that
do not (MAI1 `TalH`, PXO99A `Tal7b`) are asserted by name in
`data-raw/tale_annotations.R` and documented, so a silent change in the
match would fail the build. The member index is not kept in a second
column (Q142): it invites the mistake this section just found.
`List_of_classes.txt` is kept gzipped in `data-raw/` (677 KB), so the
dataset rebuilds without the 17-minute download, and `?tale_annotations`
and `@source` give the catalogue date. Tests 7/7.

The completed column is what makes the dataset useful: the
`talomes_heatmap()` example now draws **all ten strains over 38 classes**
where only four could be plotted before, and the clustering separates the
Asian strains (PXO61, PXO71, PXO86, PXO99A) from the African ones (MAI1,
MAI68, MAI129, MAI134, BAI3, AXO1947), which share almost no classes
with them.

### What the completed table shows, and an argument for tantale's grouping

Maintainer, 2026-10-04, on the `talomes_heatmap()` figure: the biology is
as expected. Checked against the data:
- The Asian strains (PXO61, PXO71, PXO86, PXO99A) and the African ones
  (MAI1, MAI68, MAI129, MAI134, BAI3, AXO1947) share **no AnnoTALE class
  at all**: 27 classes against 11, zero in common.
- The six African strains share a **core of nine classes** (TalBS, TalDG,
  TalDH, TalDI, TalDJ, TalDK, TalDL, TalDM, TalDN), which is Doucouré et
  al.'s "all sequenced strains shared nine evolutionary related tal
  effector genes" (Frontiers in Microbiology 9:1657,
  \doi{10.3389/fmicb.2018.01657}, PMID 30127769 -- already the `pubmed`
  value on those rows).
- Two classes sit outside that core: **TalFR in MAI129 and MAI68**, and
  **TalFU in MAI68 alone**. (The maintainer recalled the pair as MAI68 and
  MAI134; the table says MAI129.)
- All nine truncTALEs are in Asian strains.

**The maintainer's point, worth keeping for an article:** AnnoTALE gives
those divergent MAI68/MAI129 TALEs *classes of their own*, so they show
up as isolated columns saying nothing about how far out they are -- a
class is a hard partition. `tales_group_hclust()` on DisTAL distances
would keep them inside a consistent African group while showing their
divergence, which is information the class label discards. That contrast,
drawn on this very dataset, would make a concrete case for the package's
grouping as a complementary view. Vignettes come last, so this is noted,
not started.

### The two class wrappers (Q132, Q134-Q137) **[V]**

`run_annotale_load_classes()` and `run_annotale_assign()` added to
`R/annotale.R`, `@family external TALE tools`, so `_pkgdown.yml` needed
no change (`check_pkgdown()` clean). Both take `opt_param`/`java_args` as
§51 settled, refuse the keys they set themselves, and check their input
files before Java starts. `java_args` defaults to `-Xmx8G` for both: the
catalogue is large and rebuilding it was seen holding 1.8 GB resident.
`run_annotale_load_classes()` returns the class builder path invisibly so
it feeds straight into `run_annotale_assign()`, which returns its
`output_dir`; `class_builder = ` reads a kept copy instead of
downloading. Its `@return` points at `List_of_classes.txt` as the plain-text
way to answer "which class is this TALE in" for a TALE already
catalogued, which is how §57 filled `annotale_class` without running
`assign` at all. Guard tests in `test_jar_args.R` (5/5, no Java needed);
the download itself is not tested, being a quarter of an hour of network
(Q137).

**Verified by running it (2026-10-04).** `run_annotale_assign()` on
`Sample_TALEs_RVDSeqs_AnnoTALE.fasta` against the downloaded catalogue:
**9.4 minutes**, wrote `TALE_names_(<strain>).tsv`, the per-class reports
and an augmented class builder. It assigned exactly the nine core African
classes §57 had found from `List_of_classes.txt` (TalBS, TalDG, TalDH,
TalDI, TalDJ, TalDK, TalDL, TalDM, TalDN), so the text-file route and
`assign()` agree on the vocabulary -- a useful cross-check of how
`annotale_class` was filled.

For MAI1 the two routes agree on eight of nine, and `assign()` puts
`TalH`, the row `tale_annotations` leaves `NA`, in **TalDN**. The `NA`
was **not** filled from this: the sequence `assign()` classified is not
the one in the table. The shipped sample file's `TalH_MAI1` has 16
repeats (`NN-HD-NV-HD-NI-NG-NI-NN-NS-HD-HD-NI-NG-NI-NG-NI`), the table's
MAI1 `TalH` has 12 (the same first twelve, without `NG-NI-NG-NI`). So
Tram's curation and AnnoTALE's catalogue record different repeat counts
for what looks like the same gene, which is a question for the maintainer
(**Q146**) and not something to resolve by assuming. The `NA` stands.

**Q146 answered from the catalogue (2026-10-04).** Tram's 12-repeat MAI1
`TalH` string matches **four catalogue entries exactly, all class
`TalDN`**: TalDN2 (Xoo BAI3), TalDN8 (CFBP7320), TalDN9 (CFBP7321),
TalDN12 (CFBP7324). So that sequence is a real catalogued TALE, recorded
in BAI3 and three CFBP strains rather than in MAI1, whose own entry
`TalDN23` has 16 repeats. The class is unambiguous from three
directions: the exact string is TalDN elsewhere, MAI1's own entry is
TalDN, and `assign()` put the 16-repeat version in TalDN. Only the
*repeat count* stays open -- whether MAI1 carries a 12-repeat allele like
BAI3's or the 16-repeat one. **Q147 put to the maintainer:** fill
`annotale_class` as `TalDN` for that row, noting it matched other
strains' entries rather than MAI1's own, which would make the column 127
of 128. Not done without sign-off.

**The length settled by running it (maintainer's suggestion, 2026-10-04).**
A fresh `tell_tales()` on the shipped MAI1 genome finds that locus as
`ROI_00007` with **16 repeats**,
`NN-HD-NV-HD-NI-NG-NI-NN-NS-HD-HD-NI-NG-NI-NG-NI`, matching AnnoTALE's
`TalDN23` for MAI1. Tram's 12-repeat string is an exact prefix of it: the
final `NG-NI-NG-NI` is absent from her record, so the table's `rvd_seq`
for that row is truncated against the genome it came from.

**Q147 done.** `annotale_class` for MAI1 `TalH` set to `TalDN`, which both
routes agree on (the 12-repeat string is TalDN in four other strains; the
16-repeat array at that locus is MAI1's own TalDN). The fill is explicit
in `data-raw/tale_annotations.R`, guarded by `stopifnot()` so it fails if
the row ever stops being the unmatched one. Column now **127 of 128**;
only PXO99A `Tal7b` has none. Docs, test and NEWS updated; tests 7/7.

**Q148 open:** whether to correct that `rvd_seq` to the 16 repeats the
genome shows. Arguments for: it is demonstrably what the sequence says,
and the row would then match MAI1's own catalogue entry with no special
case. Against: the table is a record of Tram's curation, and editing a
curated value is a different act from completing a missing one. Left for
the maintainer; `?tale_annotations` documents the discrepancy either way.

### §56 devel, second round (2026-10-05)

The first fix (as.data.frame() before merge()) was right but not
sufficient: CI on f1626ac had **oldrel-1 passing** (the
`##source-version` fix worked) and **devel still failing** at the same
19 tests and the `tell_tales` example. The backtrace went one level
deeper: `data.frame(..., longest_orf_seq = full_orf)` in
`.telltale_add_array_measures()` puts an `XStringSet` in a base
`data.frame()`, and R devel's `data.frame()` calls
`as.data.frame(x[[i]], optional = TRUE, validRN = FALSE)` on each column.
Now `unname(as.character(full_orf))`; output-neutral (golden 18/18).

**Upstream cause and status.** CI's R devel job installs Bioconductor
*release* packages (Biostrings 2.80.2, S4Vectors 0.50.3, rtracklayer
1.72.0) next to R devel (2026-10-02). Biostrings fixed exactly this in
devel: **2.81.9, 2026-09-05, "Fix issue with as.data.frame.XStringSet()"**,
whose method now takes `validRN`. So the failure was R devel against
an unpatched release Biostrings; the change here makes tantale work
whichever Biostrings is installed.

**Searched for the rest rather than waiting for the next CI failure.**
Since R devel cannot be run locally, the full suite was run with a
tracer on `S4Vectors:::as.data.frame.Vector`, recording each call with
base `data.frame()` on the stack. After the fix, the only remaining path
(30 calls) is inside rtracklayer: `export()` builds a `data.frame()`
from an `Rle`. That one is safe: S4Vectors 0.50.3's `as.data.frame` for
`Rle` decodes to a plain vector and passes `...` (so `validRN`) on to
base R's own `as.data.frame()`. So the XStringSet was the only path
that breaks; CI is the confirmation. Suite under the tracer: 472 tests,
0 failures.

### §57 Q148 done (2026-10-05)

MAI1 `TalH` corrected to the 16 repeats of the genome
(`NN-HD-NV-HD-NI-NG-NI-NN-NS-HD-HD-NI-NG-NI-NG-NI`), in
`data-raw/tale_annotations.R`, guarded by `stopifnot()` on the old value;
the source TSV stays as Tram left it. The special-case class fill of
Q147 is gone: the corrected string matches MAI1's own catalogue entry
(TalDN23) by the normal route, asserted. 127 of 128 classed; only PXO99A
`Tal7b` has none. `?tale_annotations` says one curated value was
corrected against the genome, and why.

### Closing the session (2026-10-05)

**CI green on all three R versions** on f2e402b (run 37239246389): the
§56 fixes are confirmed. Site reference pages built for the two class
wrappers and `tale_annotations` (162b9f1). **`v0.99.0` released** as a
GitHub pre-release on 162b9f1, which differs from the CI-verified commit
only in `docs/` (Q113, Q150; the maintainer asked what a release adds
over the pushed branch: a tag pinning one commit, since `main` had moved
24 commits under the same version number). macOS check dispatched,
`macos-15-intel`, run 37242907805; result not yet read.


## 58. macOS check, first run (2026-10-05) **[V]**

Run 37242907805 (`macos-15-intel`, on b090fd5) failed on all three R
versions. Three separate causes, read from the job logs:
- **gdtools binary needs XQuartz.** `plot.tales_msa()` with a tree panel
  and `tales_group_hclust(plot_tree = TRUE)` load ggtree, whose chain
  loads `gdtools.so`, linked against `/opt/X11/lib/libXrender.1.dylib`;
  the runner has no XQuartz. Three tests (R release, oldrel-1). R devel
  passes them. A runner setup matter: install XQuartz in the macOS job.
- **R oldrel-1 segfaults** later, in a `dyn.load()` called by pillar while
  formatting a tibble, after "stack imbalance" warnings that follow the
  failed loads. Probably a consequence of the failed `dyn.load()`;
  recheck once XQuartz is installed.
- **Three golden failures on every version.** (i) `tales_compare_distal()`
  on four arrays: `domain_distances`, column `dissim`, 43 distinct values
  on Linux, 38 on macOS; `tales` identical. (ii) and (iii) both
  `tell_tales()` runs: `n_terminus_aa_alignment.html` differs in content,
  same line count. Both outputs come from `DECIPHER::AlignSeqs()` (the
  default `aln_method = "DECIPHER"` of the domain distances, and
  `.telltale_align_termini()`, which writes the HTML with `BrowseSeqs()`);
  MAFFT and HMMER are not involved. Suspects: a different DECIPHER version
  in the macOS binaries (3.8.0 here), or platform-dependent arithmetic
  inside DECIPHER. A first guess, `ceiling()` in `.arlem_cost_matrix()`,
  was rejected before the failing table was identified correctly.
- **Diagnostic run (Q154, agreed):** branch `ci-macos-diag` (temporary,
  delete once the cause is known) carries the XQuartz step in
  `R-CMD-check.yaml` and a workflow `diag-macos-golden.yaml` that dumps,
  on Linux and macOS, package versions, the conda list, the four-array
  `domain_distances`/`tale_distances`/ARLEM costs at 17 digits and the
  golden `n_terminus_aa_alignment.html`, uploaded as artifacts. Runs
  37290907358 (diagnostic) and 37290907077 (R CMD check, macOS, on the
  branch). Local Linux reference produced with the same script.

- **Diagnostic result (run 37290907358).** Same versions on both runners
  (R 4.6.1, DECIPHER 3.8.1, Biostrings 2.80.2; MAFFT 7.453 and HMMER
  3.3.2, platform builds). Linux CI reproduces the local values exactly.
  On macOS, 328 of 2304 `domain_distances` pairs differ, by up to 0.53
  (one alignment column on ~264), and the N-terminus protein alignment
  differs on three lines (`ROI_00001` and two consensus lines). So
  `DECIPHER::AlignSeqs()` places some gaps differently on macOS with the
  same version: platform arithmetic inside DECIPHER (compiler, floating
  point tie-breaks). The ARLEM costs differ as a consequence, yet
  `tale_distances` is identical on both. Both alignments are valid; the
  golden baseline is what cannot be shared. Proposed: per-platform
  snapshots (`expect_snapshot_value(variant = )`) for the DECIPHER-based
  golden expectations, the macOS variant taken from the check's uploaded
  snapshots (Q159).

## 59. `ROI_00019` of PXO86: C-terminus evidence at the DNA level (2026-10-05) **[P]**

Maintainer: the truncTALE article shows `ROI_00019` with a C-terminus
protein hit (`cterm_aa_evalue` 5.9e-18, code `CTERM`) and
`cterm_dna_hit = FALSE`; how can both hold? Reopens the cut-off for
`NTERM`/`CTERM` against `XXXXX` (§42).

Checked by running `tell_tales()` on PXO86 with `cterm_min_score` 200
(default) and 0. nhmmer does find the C-terminus profile on DNA downstream
of the last repeat (array on the minus strand, 2812633-2815116), in two
adjacent pieces:
- profile 1-97 at 2812653-2812557, 94.3 bits, E 3.4e-27;
- profile 783-861 at 2812559-2812482, 56.1 bits, E 1.1e-15.

Both are below `cterm_min_score = 200`, so neither reaches
`hits_report.tsv` and `cterm_dna_hit` is `FALSE`. All other C-terminus
hits in the genome score 228-1115 bits (complete ones ~1100). The
profile positions skipped between the two pieces (98-782, 685 nt) are not
a multiple of three. Translation agrees: AnnoTALE's 43-residue segment is
`SIVAQLSRRDPALAALTNDQLVALACLGGRPA` (canonical start of the C-terminal
region) followed by `PHSRKRKSHD*`; the second piece, read in another
frame, gives `PAFKEEEIA*LMELLPQ`, the canonical end of the TALE
C-terminus (`...PAFNEEELAWLMELLPQ`) with a stop in place of the Trp.
So the locus carries the two ends of a C-terminus joined by a
frame-shifting deletion of ~685 nt. Whether this describes the genuine
gene is the maintainer's call.

Consequences:
- The article's "nothing downstream of its last repeat resembles a TALE
  C-terminus to any of the three profiles" is wrong (Q151).
- The two flags measure different things at different stringency:
  `cterm_dna_hit` is a bit-score cut-off of nhmmer on the genome, tuned
  for discovery (200 bits, about a fifth of a complete C-terminus);
  `cterm_aa_hit` is an E-value on hmmsearch over a handful of protein
  segments, where a 37-residue match passes easily (Q152, Q153).

**The maintainer's proposal (2026-10-05).** Maintainer: (1) report more about the terminus hits, DNA and protein, in
`array_report.tsv` and/or `hits_report.tsv` (to be designed together);
(2) code a terminus `NTERM`/`CTERM` only when it has both the protein
match and a DNA hit, documented in the function docs. Claude's reading:
(1) yes; (2) to be measured first, because the DNA cut-offs are
discovery thresholds (C 200 bits, about a fifth of a complete
C-terminus), the §42 rule deliberately codes a terminus truncated at its
far end as `CTERM`, and `tales_from_annotale()` has no DNA hits to apply
the condition to. Alternative put to the maintainer: keep the protein
code, report the discordance as a `tales_anomalies()` check. Q170-Q172.

**Maintainer, 2026-10-05: a new paradigm for terminus calls.** Q171: yes,
step (1) for both termini. Q172: `NTERM`/`CTERM` are to mean a canonical
terminus, one the maintainer can be reasonably confident performs its
usual function; anything else is `XXXXX`, and the TSV files say what it
may be. So the protein match must be full-length or close to it, with
"close" learned from the match features of a set of canonical TALEs and
applied to new sequences. This replaces §42's rule (a terminus truncated
at its far end is `CTERM`). Proposed calibration set: `tale_annotations`
(10 published genomes, 128 TALEs; 9 `truncTALE`, 9 with an
`unusual_feature`, some of which concern the termini). Q173-Q176.

**Q173-Q176 (maintainer, 2026-10-05).** Q173 yes; Q174 exclude both
termini of every TALE with an `unusual_feature`; Q175 keep `NTERM`/`CTERM`
unless better names come up; Q176 no joint design needed. Maintainer's
remark: what `sanitize = TRUE` drops should follow the new definition.

Step (1) done: `array_report.tsv` gains, per terminus, `*_dna_score`,
`*_dna_evalue`, `*_dna_cover`, `*_dna_pieces` (every nhmmer hit of the
terminus profile next to the repeats, whatever its score, within twice the
profile length) and `*_aa_score`, `*_aa_far_gap`, `*_aa_cover`,
`*_aa_domains`; `hits_report.tsv` (and the two GFFs) gain `score`,
`evalue`, `hmm_from`, `hmm_to`. On PXO86 `ROI_00019`: `cterm_dna_score`
94.3, `cterm_dna_cover` 0.204, 2 pieces, as found by hand above. Golden:
only `all_ranges.gff`, `array_report.tsv`, `hits_report.gff`,
`hits_report.tsv` changed, in both runs; every existing column's digest
is unchanged; the table fingerprints only gain the new columns.

Calibration (`data-raw/terminus_calibration.R`, data in
`../tantale_calibration/`): all 128 curated TALEs matched an array by RVD
string; arrays: 112 canonical, 9 truncTALE, 7 unusual, 11 unannotated
(10 without any terminus segment). `*_aa_cover` of the canonical termini:
- C-terminus: 111 of 112 at 0.99 or more; one at 0.452, PXO61
  tal6c/PthXo3 (133 aa, DNA cover 1.0). truncTALE C-termini 0.133 (42 aa)
  and 0.656 (183 aa); PXO71 `ROI_00005` (truncTALE) has a 14-aa segment
  and no match, DNA cover 1.0.
- N-terminus: 103 at 0.99 or more; 7 at 0.927, TalC/TalE of the African
  strains (264-265 aa, two domains, an internal deletion); one at 0.833,
  PXO99A tal6b (240 aa, far gap 48, DNA cover 1.0); PXO71 `ROI_00008`
  (curated canonical) has a 22-aa segment and no match, DNA cover 0.999.
  truncTALE N-termini 0.844 (230 aa, two domains).
- `unusual` TALEs: every terminus at 1.0 except PXO61 tal1b's C-terminus
  (133 aa, 0.430).
Proposed rule: `*_aa_cover` at least 0.9 for both termini (between the
lowest accepted canonical value, 0.927, and the highest truncTALE
N-terminus, 0.844). It rejects PthXo3 (PXO61), tal6b (PXO99A) and tal1b
(PXO61), each with a full-length DNA match and a short protein segment;
whether those genes are genuine is the maintainer's call. Q177-Q180.

**Q177-Q180 (maintainer, 2026-10-05).** Q177 yes: the rule is adopted,
`terminus_min_cover = 0.9` in `tell_tales()` and `tales_from_annotale()`
(`.terminus_min_cover`), on top of the E-value and the repeat-side gap.
Q178 no: the three outliers (PXO61 PthXo3 and tal1b, PXO99A tal6b) are a
separate job. Q179: the anomaly check is renamed `terminus_noncanonical`;
`sanitize = TRUE` must not drop truncTALEs, so a more permissive sanitize
is to be defined (Q181). Q180: TalC/TalE are canonical, an unquestionable
transcription factor (maintainer).
Done: rule, argument, docs (`?tales_anchor_codes` rewritten for the new
meaning), tests (the PXO86 truncTALE C-terminus is now `XXXXX`; the
frameshift test needs both tolerances loosened). Golden: the reference
anomaly renamed, and `tell_tales.log` one line longer (the new parameter)
in both runs; the BAI3 sample's four TALEs keep their codes. On the
calibration set the rule recodes 18 termini to `XXXXX`: both termini of
the 8 truncTALEs that had matched, PthXo3 and tal1b (PXO61, C), tal6b
(PXO99A, N). Articles not yet re-rendered (truncTALE, mining).

**Q181 (maintainer, 2026-10-05): a permissive sanitize.** Q181a yes, Q181b
yes (`terminus_absent` is kept), Q181c the filter on `tales_anomalies()`
is enough for canonical-only, to be seen once implemented. Done:
`tales_anomalies()` has a `kind` column (`noncanonical` for
`terminus_noncanonical` and `terminus_absent`, `integrity` for every other
check); `sanitize = TRUE` drops only arrays with an `integrity` anomaly
(warning `tantale_warning_tales_sanitized`) and keeps the others with a
message (`tantale_message_tales_noncanonical`); the example of
`?tales_anomalies` shows the canonical-only filter. Golden: the reference
anomaly table gains the column.
Q182 (maintainer): `sanitize = "canonical"` rather than a filter, easier
for moderately proficient R users. Done: `sanitize` takes `FALSE`, `TRUE`
or `"canonical"` (anything else is `tantale_error_bad_argument`);
`"canonical"` drops every array with an anomaly, with the
`tantale_warning_tales_sanitized` warning; the construction warning names
both levels; the `?tales_anomalies` example uses it.

**Site revised (2026-10-06).** Articles: `trunctale_correction.qmd`
(both truncTALEs now flagged, coverage evidence, `sanitize` levels, and
`ROI_00019`'s C-terminus fragments, correcting the Q151 sentence),
`tale_mining.qmd` (meaning of the codes, report columns, both `sanitize`
levels; raw BAI3-1-1 now has six arrays with both termini `XXXXX`, and
`ROI_00002` is flagged too), `tales_class.qmd` (anomaly kinds, `sanitize`
levels). Article cache rebuilt from scratch (backup in the session
scratchpad); `tale_classification`, `tale_msa`, `tales_msa_class` and
`tale_target_prediction` render identical text, so the grouping is
unchanged. Home, all reference pages (now stamped 0.99.0.9000), news,
llm docs and search rebuilt without wiping `docs/`; `check_built_site()`
clean.

## 60. rOpenSci #813: first editor reply and the bot's pkgcheck (2026-10-05) **[P]**

Adam Sparks (editor, 2026-10-05 06:08 UTC) ran `@ropensci-review-bot
check package` (report 2026-10-04 23:57), is discussing the fit with the
editorial team, and asks to run goodpractice locally, singling out the
"duplicate arguments" lints as severe. Items the bot marks blocking:
- **R CMD check error, coverage failed:** 76 test failures, all from the
  absent external tools (the case the inquiry asked about). 64 of them
  surface as reticulate's base `stop()` "Unable to find conda binary",
  wrapped in a `%in%` evaluation error, from `.create_tantale_env()`
  (`R/tantale_conda_env.R`); 7 are tantale's own
  `tantale_error_tool_missing`. The 64 deserve the same cli condition.
- **`\dontrun{}`** in `run_annotale_load_classes()` and
  `run_annotale_assign()`, whose examples also write to the working
  directory (`output_dir = "annotale_classes"`).
- **Unused internal functions, unreachable URLs:** urlchecker finds only
  `doi.org/10.1142/S0219720009004060` (World Scientific answers 403 to
  scripts; the link works in a browser). The unused functions were
  false positives in §44 (S3 methods, `%||%`, a fixture builder).
- **The 20 duplicate arguments** (rechecked with lintr): all repeated
  `"i"`/`"x"` names in cli bullet vectors, which is how cli takes several
  bullets. The two "missing argument" lints are trailing commas
  (`R/distalr.R:465`, `R/target_predictions.R:465`).

**Maintainer, 2026-10-05: Q155-Q158 yes.** Done the same day, local
commits on `main` (version 0.99.0.9000, NEWS entry for the conda error):
- Q155: `.tantale_env_prefix()` checks `reticulate::conda_binary()` first
  and aborts with `tantale_error_tool_missing`, pointing to
  `tantale_setup(install = TRUE, conda = TRUE)`; test in
  `test_tantale_setup.R`.
- Q156: both examples write under `tempfile()`; the trailing commas are
  gone. `\dontrun{}` kept for now: under `--run-donttest` (CI, and
  `--as-cran`) the two examples would add the catalogue download twice
  (~15 min each, plus the assignment, ~9 min) to every check job, and R
  devel already takes ~57 of its 120 minutes (Q160).
- Q157: lint pass with lintr 3 on the goodpractice categories, 66 edits
  plus 77 `expect_equal()` -> `expect_identical()`. The 11
  `expect_equal()` left compare integers with doubles (`nrow()` against
  `n^2`, `choose()`, ARLEM ids). Left on purpose: the 21 duplicate-argument
  lints (cli bullets), `length(unique()) == nrow()` in `summary.tales`
  (a count, not a duplicate test), the URL built with `paste(sep = "/")`,
  one `<<-` in a calling handler, `%in% TRUE` (NA-safe on purpose),
  `:::` in tests, line length. Full suite after the pass: only the 11
  type-strict expectations failed, since reverted; golden 46/46.
- Q158: draft reply in the session, for the maintainer to post.

**Follow-up (2026-10-05).** macOS check on `ci-macos-diag` (run
37290907077): XQuartz installed, the three plot tests pass, the oldrel-1
segfault is gone; only the three DECIPHER golden expectations fail
(release, oldrel-1). The reply draft and the bot's report are kept
outside the repository, in `../tantale_ropensci/` (maintainer: not in the
git tree until things are settled). Q159: maintainer accepts per-platform
snapshots if nothing lighter exists. Q161: push deferred until the
batch is complete.

**Q159b, Q160 (maintainer, 2026-10-05), done.**
- Q159b with mmseqs2 for the distance golden (Biostrings is too slow,
  maintainer). `test_golden.R` runs `tales_compare_distal(aln_method =
  "mmseq2")`; `helper-golden.R` keeps only the line counts of
  `*_terminus_aa_alignment.html` (`.PLATFORM_DEPENDENT`). Re-baselined,
  every changed row explained: `domain_distances` gains mmseqs2's output
  columns (`evalue`, `raw`, `pident`...; 664 of 2304 pairs without an
  mmseqs2 hit, so `NA` there) and `dissim` has 39 distinct values (43
  with DECIPHER); `tale_distances` changes in `dissim` and `arlem_score`
  only, following from those distances; in both `tell_tales()`
  fingerprints the two protein terminus alignments now carry `NA` digests,
  nothing else moved. Golden 46/46 twice. The DECIPHER backend keeps its
  non-golden tests. Whether mmseqs2's osx-64 build gives the same values
  is checked by the macOS run on `ci-macos-diag`.
- Q160: `@examplesIf interactive()` for both class wrappers; no
  `\dontrun{}` left in `R/`.

**Local pkgcheck on 0f33e9c (2026-10-05).** pkgcheck 0.3.2, ctags/global
from a scratch micromamba prefix, detached worktree, 21 min. Report
rendered with `checks_to_markdown()` to
`../tantale_ropensci/pkgcheck-local-2026-10-05.html` (outside the repo),
network page beside it. With the tools installed: R CMD check no errors,
no warnings (one NOTE, `.git` in a worktree); coverage 90.4%; URLs all
reachable from here. Still failing: "no continuous integration" (an
artefact of the detached worktree, where pkgcheck finds no default
branch to query; the bot saw the badge and the runs) and "all internal
functions are used", the same four as §44: `.repeat_to_rvd_align()`
(used only by `test_plot_tales_msa.R`), `[.tales` and
`[.pairwise_distances` (registered S3 methods), `%||%` (an infix, called
four times). Dependency table: fs, gplots, cluster, ggnewscale and
rtracklayer all have call sites; pkgstats misses them, and attributes
calls to mgcv and Matrix, which tantale does not use. macOS check on
`ci-macos-diag` (run 37301244765): R release and oldrel-1 pass with the
mmseqs2 golden; R devel still running.
Final: macOS run 37301244765 green on R release, oldrel-1 and devel;
Linux run 37301468672 on `main` (0f33e9c) green on all three. Branch
`ci-macos-diag` deleted (the diagnostic workflow survives in commit
9fe8e4b). §58 closed.

**Reply rewritten (Q158, 2026-10-05).** Blocking items first, then the
lints by severity. goodpractice 1.2.0.1 was already installed; its lint
check is `lintr::lint_package(linters = goodpractice:::linters_to_lint())`,
rerun on 0f33e9c. The earlier draft claimed every lint but line length
was fixed; untrue. Left in `R/`: 10 `sapply()`, 4 implicit assignments
(among them `ifelse(..., refName <- ...)` used for control flow,
`R/tales_plot.R:860`), 2 `<<-` (`.tales_anomalies()`'s `add()`, the
calling handler in `R/telltale.R:765`), 1 nested pipe, 3 `%in% TRUE`.
Tests: 35 implicit assignments (`expect_warning(x <- f())`), 26 `:::`,
4 `<<-` (caches). Open for the maintainer: Q162 fix the `sapply()` calls,
the `ifelse` assignment and the nested pipe before posting; Q163 move
`.repeat_to_rvd_align()` into a test helper, which clears the "unused
internal functions" item of its one real case; Q164 run the macOS check
on `main` before the reply claims macOS (pushes run Linux only).

**Q162, Q163, Q164 yes (maintainer, 2026-10-05).**
- Q162 done: the 10 `sapply()` replaced (`vapply()` with a declared type,
  `lapply()` for the XStringSet lists in `.telltale_align_termini()` and
  `.telltale_array_gaps()`, `do.call(rbind, lapply())` in
  `.build_repeat_msa()`, `match()` for `talomes_heatmap()`'s side bar,
  vectorised `ifelse()` for its truncation labels); the four implicit
  assignments and the nested pipe in `R/` rewritten. A test of the
  suspicion that the outer `sapply()` of `.telltale_align_termini()` could
  simplify two equal-sized XStringSets into a list-matrix: it does not
  (S4 elements stay a named list). `talomes_heatmap()`'s side bar now
  takes a strain's first `extra_col` value, where several values used to
  make `sapply()` return a list. Left in `R/`: 21 duplicate arguments,
  2 `<<-`, 3 `%in% TRUE`, 1 `paste(sep = "/")`, 1 `length(unique())`.
- Q163 done: `.repeat_to_rvd_align()` moved from `R/conversion.R` into
  `test_plot_tales_msa.R`, its only user.
- Q164: needs the commit on `main`, so the push comes first.
- Q168 (maintainer: yes): `.github/workflows/pkgcheck.yaml`, run by hand.
  Job container `ghcr.io/ropensci-review-tools/pkgcheck-action:latest`
  (the official action's image, built on the bot's), so every step runs
  in rOpenSci's environment; micromamba, the pak install `/check.R` does,
  `tantale_setup(install = TRUE)`, then the image's own `/check.R`. The
  report goes to the run's summary page and an artifact; the job fails
  when the summary has a blocking item, as the action does. The official
  action itself cannot take a setup step: it is a Docker action, so
  earlier steps run on the host, outside its container.
  First run (37360975924): setup, `tantale_setup()` and goodpractice's
  steps (covr, rcmdcheck) went through; pkgcheck then aborted in its
  ORCID check. The image carries pkgcheck 0.2.0.044 (c4743df), whose
  `paste(x$given, x$family)` gives two names for "Bao Tram"; fixed
  upstream 2026-09-29, and the bot runs 0.3.2. The workflow now updates
  pkgcheck from GitHub first. Java is not in the image either: pak
  installed `default-jdk` as a system requirement of tantale's
  dependencies. Also, the hand-started macOS run had cancelled the Linux
  run of the same push (one concurrency group per branch);
  `R-CMD-check.yaml` now adds the runner to the group.

**End of session (2026-10-06).** CI on 11f93d8 still running at session
end: pkgcheck 37380949134, Linux 37380940975, macOS 37380952515;
test-coverage 37380940798 passed. e162f49 (README) committed, not pushed.
Seen in passing, not acted on: `tell_tales.log`'s "Number of gaps of size
below 500nt" divides by 2 and can print 1.5 (pre-existing, in the shipped
example log too).

**CI on 11f93d8, read 2026-10-06.** R-CMD-check green on Linux and macOS
(`macos-15-intel`), R release, oldrel-1 and devel. pkgcheck 37380949134
failed, but after every check had run (goodpractice's covr, rcmdcheck,
lintr etc. all prepared): the image's `/check.R` (pkgcheck-action
`R/check.R`, line 43) died in `writeLines(md[1:(s_break - 1)],
"summary.md")` with "can only write character objects". So
`checks_to_markdown()` of the current pkgcheck (taken by the "Update
pkgcheck" step) returns something that is not a character vector, which
the image's script, written for the older pkgcheck, does not expect. No
report was written, hence no artifact and an empty summary page. Not a
verdict on tantale either way.
Q186 (maintainer: a): the workflow no longer calls `/check.R`; its "Run
pkgcheck" step runs the same steps inline (`Rscript {0}`), flattens the
report to character if `checks_to_markdown()` returns anything else (and
logs the non-character elements, material for an upstream issue, Q186c),
writes `full.md`/`summary.md` and sets `status`. Locally (pkgcheck 0.3.2,
cached goodpractice report for 3570f42, with or without `GITHUB_PAT`)
`checks_to_markdown()` returns a character vector, so the list arises
only on the runner. Q187 (maintainer: yes): the #813 reply waits for a
green run.
