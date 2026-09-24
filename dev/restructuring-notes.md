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

- **`tales_rvd_strings(rvd_only =)`** means "repeats only"; its sibling
  says `repeats_only` (§8.2b).
- **§2** retire `repeat_to_rvd_map()`, and **§20** rename/rewrite
  `tale_parts_to_rvd()`. Reserved for the maintainer.
- **`plot.tales_msa()`'s x-axis title** still reads "Position in array";
  the axis is the alignment position (§22, `R/tales_plot.R:406`).
- **§21 option (c)**: `tales_consensus()`/`tales_consensus_match()` taking
  a `tales_msa` directly (`dev/class-design.md` §4.6 has it as `[A]`).
- **`tell_tales()`'s 17 arguments** were never regrouped, and two `TODO`
  blocks remain in its body (circular molecules; what two output files
  should contain) (§5.3).
- **ARLEM's duplication and insertion costs** are fixed at 10
  (`.arlem_dup_cost`, `.arlem_indel_cost`), not exposed. Their ratio to
  substitution costs (up to 99) is a modelling choice about how TALE
  arrays evolve (§6).
- **Distribution channel** and the one-archive plan (§34).
- **A beta of 1.0.0**, if wanted: number it **0.99.0** (then 0.99.1, ...).
  R versions are numeric only (`1.0.0-beta` is invalid, `1.0.0-1` sorts
  after 1.0.0); 0.9.9010 < 0.99.0 < 1.0.0; Bioconductor also requires
  0.99.z for new submissions. Mark it with a `# tantale 0.99.0` NEWS
  heading, a `v0.99.0` tag and a GitHub *pre-release*. Hard renames are
  still allowed at 0.99.x (lifecycle rules start at 1.0.0).

### Housekeeping pending (2026-09-24)

- **`docs/` is behind 0.9.9010**: reference pages (the `run_annotale_*()`
  return values, `talomes_heatmap()`), home page (README's stability
  note) and news. Partial rebuild per `dev/CLAUDE.md` (reinstall first;
  no page added or removed, so no wipe). The site may be offline anyway
  while the repository is private.
- **Re-run `dev/function-graph-dataflow.R`** (now ~10 min, the suite is
  longer) so the §29.2 data-flow view records the functions tested since:
  `talomes_heatmap()`, `preditale()`, `plot_target_preds()`,
  `run_annotale_*()`.

### Worth investigating, no deadline

- `correct_array = TRUE` fabricated a ~489 nt region on a two-array PXO86
  excerpt that does not exist in the genome; the fixture
  `dev/fixtures/pxo86_roi18_19_excerpt.fa` reproduces it in under a
  minute (§25). Did not reproduce on the full genome.
- Which corrected BAI3-1-1 object the article cache should build on (§25b).
- Validate the 136-sequence correction reference against the 494 one
  (maintainer's own, §8.1b).
- Each toy region yields a spurious single-hit array at its 3' end,
  because `min_domain_hits` filters subject sequences (§8.1).
- `.rvds_from_annotale_file()` (`R/tales_ingest.R:55`) has no caller in
  `R/` or `tests/`: a parking candidate (§29.1).
- Consumers of `tales_rvd_strings()` other than
  `tales_to_universalmotif()` were never checked for silently dropping
  an array with no repeats (§12b).
- Seven condition sites still carry only the generic `tantale_error`
  class (`conversion.R` 2, `talecorrection_java.R` 3, `tales_ingest.R` 1,
  `tales_plot.R` 1) (§9.5).
- Rendered error messages from an installed package show a source path
  ("at tantale/R/tales_class.R:818:3"); cosmetic (§30).
- Commented-out developer snippets still carry `/home/cunnac/...` paths
  (§34).
- `man/figures/pipeline.svg`/`.png` are stale and referenced from nowhere
  outside `dev/` (§7.6). The data-flow view of §29.2 could replace them.

### Parked, reserved or deferred

- **Reserved for the maintainer, do not start unasked:** §21 items 1-4
  (the matrix helpers of `plot.tales_msa()` and the tests pinning their
  shape), §20, §2.
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

## 2. Conversion functions **[P]**

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
- `man/figures/pipeline.svg`/`.png`: see §7.6.

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
- **`man/figures/pipeline.svg`/`.png`** show pre-restructuring function
  names and are referenced from nowhere outside `dev/`. Redraw, remove,
  or replace with the §29.2 data-flow view. Re-exporting the PNG needs
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

## 20. `tale_parts_to_rvd()` -- candidate for a rename and a refactor/rewrite **[P]**

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
position; renamed `alignment_position`. **Held for the maintainer:** the
visible axis title still reads "Position in array"
(`R/tales_plot.R:406`), where `plot.tales()` says "Position in
alignment" for the same coordinate.

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

## 29. Function dependency diagram, to rebuild the maintainer's mental map -- phase 1 (29.1) and phase 2 (29.2) built **[V]**

Maintainer's request (2026-09-22): a map of the package covering
internals, and the class objects if legible. Decisions: a dev-only `.qmd`
under `dev/`, interactive (`visNetwork`), generated from the code so it
cannot go stale, opening grouped by file. Distinct from
`man/figures/pipeline.svg` (§7.6), a curated workflow figure.

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
