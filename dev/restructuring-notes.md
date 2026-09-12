# tantale — restructuring notes and action ledger

Working document for the pre-publication overhaul. Records findings, agreed
actions and deferred questions so they don't live only in conversation.

Branch: `dev`. Last updated: 2026-09-12.

Status markers used below:

- **[V]** verified empirically against the code/data in this repo
- **[A]** agreed direction, not yet executed
- **[P]** parked — needs a judgement call or a dedicated review pass

---

## 0. Framing

`tale_parts` was introduced late in development, in response to an interface
that had grown too many conversion functions. It emerged while `distalr()` was
being written. Before that, computations went through `runDistal()` plus a
family of conversion functions whose job was to reshape its output into
whatever MSA or plotting needed.

The current code is therefore transitional: several exported conversion
functions are archaeology from the pre-`tale_parts` era, and can be retired or
unexported once `tale_parts` is genuinely central.

Two cautions that apply throughout, and that earlier analysis got wrong:

- **Absence from `pipeline.svg` says nothing about a function's value.** The
  package also ships standalone utilities useful outside that workflow.
- **Zero internal call sites says nothing either.** Exported functions are
  meant to be called by users; unused internals may be dormant rather than
  dead. Judge by what a function *does*, not by how often the package calls it.

---

## 1. Shrink `distalr()`'s returned list

`distalr()` currently returns six elements. Three are derivable from a fourth.

| element | finding |
|---|---|
| `coded.repeats.str` | **[V]** pure projection of `tale_parts`. Rebuilt from the stored `tale_parts` and compared: identical names, byte-identical sequences, identical widths, `all.equal()` TRUE. `identical()` differs only in the `XStringSet` `pool` slot (an allocation detail, not content). |
| `repeats.code` | **[V]** pure projection. `identical()` TRUE outright. |
| `repeats.cluster` | **[V]** zero consumers anywhere in `R/`, `tests/` or `vignettes/`. The only code needing cluster IDs — `.repeat_to_cluster_align()`, for plot colouring — recomputes them from `repeat.similarity`, and at a different default cut height (`distalr()` uses `h_cut = 10`; both plotting functions default to `90`). |
| `tale_parts`, `tal.similarity`, `repeat.similarity` | **[V]** necessary and sufficient. The two similarity tables are the irreducible expensive products (pairwise protein alignment, then ARLEM on coded arrays); `tale_parts` is the substrate the rest project from. |

**Actions**

- **[A]** Drop `coded.repeats.str`, `repeats.code`, `repeats.cluster` from the
  returned list.
- **[A]** Re-expose the two projections as *methods* rather than stored fields,
  so the reconstruction lives in exactly one place.
- **[P]** `repeats.cluster` embeds a cut-height choice that consumers override
  anyway. Cut height looks like a display/analysis parameter, not an intrinsic
  property — confirm before deciding whether any clustering belongs in the
  object at all.

Note on the reconstruction: it is safe today only because `distalr()` hard-stops
when any `aaSeq` is `NA`, so `domCode` is never `NA`. If that guard were relaxed,
`paste(domCode, collapse = " ")` would silently emit the literal string `"NA"`
as a repeat code. Another reason for the logic to exist once, not at each call
site.

---

## 2. Conversion functions

### Retire

- **[V]** `repeat_to_rvd_map()` is redundant given `tale_parts`. Both paths were
  computed from the same fixture and compared: 251 rows each, identical once
  sorted, zero rows unique to either. The function exists to *re-derive* by
  melt-and-rejoin a correspondence that `tale_parts` already holds on a single
  row (`domCode` and `rvd` side by side). It is exactly the kind of function
  that was necessary when `runDistal()` handed back `coded.repeats.str` and RVD
  fasta files as unrelated artifacts.
- **[A] Migrate its assertion first.** It `stopifnot`s that each `repeatID` maps
  to exactly one RVD. That invariant is currently enforced *nowhere else* — move
  it onto `tale_parts` validation rather than losing it.
- **[A]** `repeat_to_rvd_map_distalr()` survives, but its name is misleading: it
  depends on `domCode` being present, not on `distalr()` having been run.

### Keep internal

- **[V]** `.repeat_to_sim_align()`, `.repeat_to_cluster_align()` — genuinely
  plotting-only internals, correctly unexported.

### Dormant but valuable — repair, then decide (NOT removal)

An earlier pass recommended deleting these two on call-count alone. That was
wrong; both implement capability available nowhere else in the package.

- **`.rvd_to_match_align()`** is the structural twin of `.repeat_to_sim_align()`
  at a different level of biological abstraction. The repeat-level version asks
  *how similar are these two repeats' protein sequences*; this one asks *how
  similar are these two RVDs' DNA-binding preferences*, using `rvdSimDf`
  (correlation of base-preference profiles: `NI`/`NI` = 1.0, `NI`/`NN` = 0.63,
  `NI`/`NG` = −0.4). Those come apart — repeats can be sequence-divergent yet
  share an RVD (same specificity), or near-identical yet differ at positions
  12–13 (different target base). **[V]** Every fill mode in both plotting
  functions is repeat-level (`fill_type` ∈ {`repeat_sim`, `repeat_clust`};
  `plot_type` ∈ {`repeat.similarity`, `repeat.clusters`,
  `repeat.clusters.with.rvd`}) — there is no RVD-level mode. This is an unwired
  feature, not a duplicate. It is also the sole consumer of `rvdSimDf`; removing
  it orphans that dataset entirely.
  - **[A]** Fix `tantale::rvdSimDf` → `rvdSimDf` (internal data accessed with
    `::`, which errors — the function cannot ever have run as written).
  - **[P]** Consider exposing as `fill_type = "rvd_sim"`. Integration wrinkle:
    `Cor` is signed on [−1, 1] whereas repeat `Sim` is [0, 100], so the colour
    scale can't be reused as-is.
  - **[P]** Confirm the provenance of `rvdSimDf$Cor` — some values look
    surprising relative to a plain correlation of `rvdToNtAssocMat` rows.

- **`.rvd_to_repeat_align()`** is the inverse of the exported
  `repeat_to_rvd_align()`. The package currently supports only one direction
  (align on repeat codes, view as RVDs). This provides the other (align on RVDs,
  view as repeats), which matters because the two alphabets differ greatly — in
  the sample data, ~17 RVD symbols versus 251 repeat codes, mostly singletons —
  so the two alignments are genuinely different problems. `plot_tales_msa()`
  consumes *both* matrices simultaneously, so an RVD-based alignment currently
  has no way to obtain its matching repeat matrix.
  - **[A]** Add a guard: it assumes each row's non-gap count equals
    `length(repeat_vecs[[r]])` and does not check. A mismatch silently
    misregisters codes against positions.

---

## 3. Legacy cemetery (`inst/legacy/`)

**[A]** The AnnoTALE↔QueTAL format shims:
`.annotale_to_quetal_rvd()`, `.quetal_to_annotale_rvd()`, `.reformat_array_report()`.

"Shim" here means a thin file-format adapter between two external tools'
incompatible on-disk conventions for identical content — no biological or
computational transformation. AnnoTALE writes RVDs as standard FASTA (name and
sequence on separate lines); QueTAL/FuncTAL wants one line per TALE with a tab
between ID and RVD string. `.annotale_to_quetal_rvd()` is literally
`paste0(">", name, "\t", seq)`.

Reasons they go:

- **[V]** Dead as a *set* — they feed the FuncTAL branch, which is itself
  dormant (`functal()` is blocked on uninstalled Perl deps).
- **[V]** They carry lab-specific hardcodings presented as general converters:
  `(MAI\d{1,3}).*TALE(\d{1,3})` name rewriting, and a `Sebra/` directory regex
  for deriving strain names.
- **[V]** `.reformat_array_report()` is stale against its own upstream. It
  strips extremity flags using `^(BBB-)*([^(ZZZ)]*)(-ZZZ)*$`, but `BBB`/`ZZZ`
  appears nowhere else in the package and current `tell_tales()` output uses
  `NTERM`/`CTERM` (0 occurrences of `BBB`/`ZZZ` in `arrayReport.tsv`; `NTERM`/
  `CTERM` throughout `rvdSequences.fas`). It would strip nothing today.

If the FuncTAL branch is ever revived, the shim should not be: with
`tale_parts` + `tale_parts_to_rvd()`, producing QueTAL input is a small writer
off a `tale_parts` object, not a file-to-file translation between two foreign
formats.

---

## 4. Plotting: `msa_heatmap()` is superseded

**[A]** Label `msa_heatmap()` obsolete. It was written first, in base graphics
(`gplots::heatmap.2`), with a clumsy style and cryptic comments.
`plot_tales_msa()` came later using ggplot grammar and returns a composable
`aplot` object (vignette p2 demonstrates disassembling `p1$plotlist`, restyling
a panel and reassembling).

**[V]** Both call the *same* internals — `.repeat_to_sim_align()` and
`.repeat_to_cluster_align()` — so the underlying data transformations are
already shared. The difference is purely rendering.

First-pass feature delta (what `msa_heatmap()` has that `plot_tales_msa()` does
not). **[P]** — needs confirmation before deprecating:

| feature | assessment |
|---|---|
| `save_path` (direct-to-file, dimensions computed, format from extension) | trivially replaced by `ggsave()` on the returned object |
| `note_colors` (configurable matched/mismatched colour pairs; viridis ramp in similarity mode) | genuine gap — `plot_tales_msa()` has no colour-scheme argument |
| `...` passthrough to `heatmap.2` | backend-specific escape hatch; not meaningful to port |
| bespoke legend via `extrafun = extra.key()` | ggplot handles legends natively |
| `plot_type = "repeat.clusters.with.rvd"` (combined cluster fill + RVD text) | appears **already covered**: `plot_tales_msa()` expresses the combination by *which arguments are supplied* (both `repeat_align` and `rvd_align`, with `fill_type = "repeat_clust"`) rather than by a mode string — arguably the better design |
| row dendrogram | present in both (`Rowv` vs `ggtree` + `aplot::insert_left`) |

Still to check: whether `consensus = TRUE` and the `ref_pattern` reference-row
mechanism behave identically in both.

---

## 5. OOP restructuring

> Detailed class definitions — identity, invariants, constructors, method
> policy — now live in **`dev/class-design.md`**. This section stays the
> ledger: what is settled, what is open, and why.

**[A]** Settled so far:

- Additive S3 — classes tag the native type (tibble/matrix stays a
  tibble/matrix), preserving dplyr compatibility.
- Per-stage classes; no top-level session/project object.
- Long table as the canonical form for alignments.
- Data column names to be canonicalised (`TAL1`/`RepU1`/`RepID`/`Rep_clust` →
  a shared schema). This is the data-level counterpart to the code-level
  rename already done.
- `annout` migrated to S3 — **[V]** with the caveat that its
  `contains = "AAStringSet"` inheritance is load-bearing: `telltale.R:651` does
  `unlist(Biostrings::AAStringSetList(annoTaleOut))`, which works only because
  each `annout` *is-a* `AAStringSet`. An S3 class cannot inherit that way, so
  the collection step needs rewriting, not relabelling. Note `annout` is purely
  a private bundling device inside `tell_tales()`'s per-array loop (an
  `AAStringSet` plus a `domainsReport`), unpacked immediately afterwards — it
  never reaches the public pipeline despite being `exportClasses`-tagged.

**[P]** Open design questions:

- Primary justification. **[V]** the dominant structural problem is that every
  core shape exists in two representations (wide matrix ↔ long table), with the
  conversion re-implemented ad hoc at each call site: **17** `melt`/`acast`/
  `dcast` calls and **16** positional `colnames(x) <- c(...)` assignments across
  `R/`. `acast(tal_sim, TAL1 ~ TAL2, value.var = "Sim")` appears at four sites.
  Semantics are carried by convention, not structure — if `melt()` ever changed
  its column order, the package would silently mislabel data rather than error.
  This argues the classes need canonical representations and conversion
  *methods*, with validation as the second benefit rather than the first.
- ~~Whether `repeat_align`/`rvd_align` become parent + subclass (shared methods
  written once, since both are character matrices with `arrayID` rownames
  differing only in cell meaning) rather than two peers.~~ **Resolved — neither.**
  They are two *value layers over one alignment geometry*, which a long
  `tales_msa` carries simultaneously; the matrices become `as.matrix()` views.
  See `class-design.md` §4. This is why `plot_tales_msa()` currently needs both
  `repeat_align` and `rvd_align` as separate arguments — a matrix can only hold
  one layer.
- ~~Whether `tal.similarity` and `repeat.similarity` unify into one class.~~
  **Resolved — yes**, parent `pairwise_sim` with `tale_sim`/`repeat_sim`
  subclasses carrying entity semantics only, and canonical id columns. See
  `class-design.md` §3. Their divergence was indeed accidental: different ID
  column names (`TAL1`/`TAL2` vs `RepU1`/`RepU2`), different column *order*
  (`RepU2` precedes `RepU1`), different extras
  (`arlemScore`/`maxLength`/`normArlemScore` vs `Dissim`) — yet downstream code
  treats them interchangeably.
- **[V]** `domCode` is a whole-set-dependent surrogate key (`cur_group_id()`
  over `aaSeq`). Recomputing it on a subset renumbers everything and silently
  breaks the join to both similarity tables. It must be carried, never
  recomputed — a real invariant for a class to protect, and an argument that a
  `tale_parts` and its companion similarity tables must be subset coherently.
  **Resolved** — "carried, never recomputed" is now invariant 4 of `tales`
  (`class-design.md` §2.4), and cross-object coherence is *enforced* by a
  `dom_code` namespace tag: a content hash stamped on the `tales` and on each
  `dom_code`-keyed companion, compared by methods that consume two of them
  (`class-design.md` §3.5). **No container object is needed** — the
  "no top-level session/project object" decision above stands. **[V]** The real
  hazard was never subsetting (which fails loudly, or not at all) but *mixing
  runs*: `cur_group_id()` mints `1..N` every run, so a cross-run join succeeds
  and silently maps repeats to the wrong sequences.

Downstream consequence: a good part of the conversion functions can then be
unexported.

---

## 6. Correctness review backlog **[P]**

Deferred to a dedicated pass on "computations that may not match intent":

- `hclust(as.dist(Sim))` in both `.cluster_repeats()` and
  `.repeat_to_cluster_align()` — a similarity (`Sim = 100 - Dissim`) is fed
  where a distance is expected, so high similarity reads as far apart.
- `rvdSimDf` is orphaned while `build_repeat_msa()` forces identity scoring for
  RVD alignments (`if (is.null(repeat_sims) || repeatType == "rvds")
  maffMatOpt <- ""`), even though `rvdSimDf` is exactly the biologically
  informed RVD substitution matrix that case would want.
- `aaSeq` ↔ `rvd` is 1:1 in the sample data but not enforced. The two come from
  independent sources (protein-parts file vs `rvdSequences.fas`) and
  `tale_parts()`'s own comments say the disagreement is deliberately kept
  visible. If they ever diverged, `repeats.code` would gain duplicate `code`
  rows and quietly stop being a key-per-code table.
- `build_repeat_msa()` infers whether its input is RVDs or repeat codes by
  testing against a hardcoded list of six frequent RVDs
  (`NN NG HD NI N* NS`). A small alignment of unusual TALEs containing none of
  them would be silently misclassified. A typed object removes the guess.
- `functal()` is blocked on uninstalled Perl deps (`List::MoreUtils`,
  `Bio::Perl` — the latter absent from the conda BioPerl package). Currently
  documented by its test.
- **[V]** `tales_consensus_match(long = TRUE)` mislabels the alignment
  coordinate as `positionInArray`
  ([msa.R:74-76](../R/msa.R#L74-L76)). It melts an MSA matrix whose columns are
  *alignment* positions — `build_repeat_msa()` sets
  `colnames(...) <- 1:ncol(...)` on the **gapped** matrix
  ([msa.R:233](../R/msa.R#L233)) — then asserts the name positionally. The two
  coordinates diverge as soon as a gap is inserted: in
  `sampleRepeatMsaByGroup.rds`, array `BAI3-1-1_ROI_00006` has 14 parts spread
  over 18 columns, agreeing up to part 12 and then jumping — part 13 sits at
  alignment position 17, part 14 at 18. Nothing errors, because the name is
  asserted rather than derived. Another instance of the 16 positional
  `colnames(x) <- c(...)` assignments noted in §5. **Expected to resolve
  itself** when the alignment gains a real `alignment_position` coordinate
  distinct from `position_in_array` (`class-design.md` §4) — listed here so it
  is not lost if that design changes, not as separate work.

### ARLEM score semantics

**[V]** `arlemScore` is a **cost (a distance), not a similarity**. ARLEM takes
`-cfile cost_file`; its error strings refer to "the cost matrix" and "distance
matrix"; `distalr()` feeds it the `Dissim` matrix as that cost file
([distalr.R:571](../R/distalr.R#L571)); and every self-comparison yields
`arlemScore = 0`. Three consequences:

- **`tal.similarity$Sim` is crushed into the top ~8% of its nominal scale.**
  [distalr.R:666-667](../R/distalr.R#L666-L667) computes
  `normArlemScore = arlemScore/maxLength` (mean per-position alignment cost)
  then `Sim = 100 - normArlemScore`. The subtraction presumes the cost is on a
  0–100 scale. It could be in principle — the cost matrix spans 0–99 — but
  aligned TALEs mostly match, so the observed mean per-position cost is ~6.
  Measured on the fixture: `arlemScore` 0–231, `normArlemScore` 0–8.25,
  therefore **`Sim` 91.75–100**. All discriminating signal sits in 8 units of
  100. Note this affects the *array-level* table only; the repeat-level
  `Sim = 100 - Dissim` ([distalr.R:565](../R/distalr.R#L565)) is a genuine
  full-range similarity.
- **Every consumer immediately undoes it.** All three users of
  `tal.similarity$Sim` convert straight back to a distance —
  [classification.R:24](../R/classification.R#L24),
  [msa.R:332](../R/msa.R#L332), [msa.R:760](../R/msa.R#L760) — recovering
  exactly `normArlemScore`. The `Sim` column is a round trip that costs
  interpretability and buys nothing; exposing the distance directly would be
  more honest. *(This also narrows the `as.dist(Sim)` question above: at array
  level everyone does invert correctly. Only the repeat-level clustering fails
  to, and there the values really are 1–100 similarities — so that one looks
  like a genuine inversion bug rather than a scale quibble.)*
- **Three magic numbers in the cost model.**
  [distalr.R:612-613](../R/distalr.R#L612-L613) hardcodes
  `"# Indel align 10"`, `"# Indel hist 10"`, `"# Dup 10"` into the cost file.
  So an indel or duplication costs 10 while a substitution between dissimilar
  repeats costs up to 99. Since TALE arrays evolve largely by duplication and
  loss, the relative weighting of duplication against substitution is a real
  biological modelling choice — currently fixed, undocumented and not exposed
  as a parameter.

---

## 7. Documentation and artifacts

- **[V]** `man/figures/pipeline.svg` has been updated on `dev` (names refreshed;
  `runDistal()`, `tree` and six orphaned arrows removed).
- **[A] TODO — re-export `man/figures/pipeline.png` from the current
  `pipeline.svg` in Inkscape.** It must be done from Inkscape rather than
  scripted: librsvg's text metrics differ and would silently change the
  figure's typography. The SVG carries the export settings already
  (`export-filename`, 300 dpi, 3555x1621). *(A re-export was done once and then
  destroyed by an erroneous `git checkout --` on my part; it needs redoing.)*
- **[V]** `p3_tantale_objects.Rmd` embeds that figure by absolute local path
  (`/home/cunnac/...`), so it renders for nobody else. The vignette body is
  otherwise just `TODO!!!`.
- **[V]** The committed pkgdown site (`docs/`) was built 2023-09-17. 25 of its
  28 reference pages document functions that no longer exist under those names,
  including pages for deleted functions (`runDistal`, `buildDisTalGroups`,
  `read_distal_aligns`, `tellTaleLegacy`, `tellTale2`). `_pkgdown.yml` itself
  needs no change. A rebuild is the fix — but see below.
- **[V]** A clean rebuild is **blocked on vignette reproducibility**, not on the
  rename: vignettes 2, 3 and 4 all begin with
  `load(file.path(outdir, "mining.RData"))` from `~/TEMP/test_tantale/`, which
  does not exist on a fresh machine. The vignettes are currently
  personal-workstation notebooks rather than portable documents.
- **[A]** Self-deprecatory / quality-shadowing language to revise:
  - `README.md:54` — "work in progress ... not necessarily fully and properly
    implemented!!!"
  - `R/tantale.R:6` — "(IDEALLY)" in the package description, which renders on
    `?tantale`
  - `README.md:10` — "nightmarish experience" (about the ecosystem, not tantale,
    but prominent and early)
  - `R/AnnoTALE_QueTAL_functions_library.R` — names a collaborator's unpublished
    data ("Hinda's Malian strains") and a tool crash in shipped source
  - `R/tellTale_utilities.R:136` — `## !! THIS SHOULD BE MADE OBSOLETE ...`
  - `R/telltale.R:317` — "I do not know why but it fails to work ..."
  - `inst/legacy/tellTaleLegacy.R` — "far from optimal", "curiosity only"
    (lower priority; already unexported with no man page)
- **[V]** Local roxygen2 is 8.0.0 and rewrote `RoxygenNote` →
  `Config/roxygen2/version`. Worth checking against the usual toolchain.
- **[P]** `talomes_heatmap()`'s roxygen describes `group_col` as displayed "as
  rows" and `strain_col` "as columns", but the code does
  `dcast(tale_annotation, strain ~ group)` — rows are strains, columns are
  groups. The doc text looks reversed (pre-existing, not introduced by the
  rename).

---

## 8. Tests

- **[V]** `test_plot_tales_msa.R` contains no `expect_*` calls — it registers as
  an empty/skipped test. It runs code inside `try()` but asserts nothing.
- **[A]** Coverage gaps, deliberately deferred until class shapes settle:
  `msa.R` 21%, `conversion.R` 22%, `target_predictions.R` 17%. (Baseline
  measured before `test_group_tales.R` and `test_functal.R` were added; overall
  was 37%.)
- **[A]** Conventions established: run targeted test files rather than the full
  suite unless a change ripples broadly; tests depending on the `tantale` conda
  environment fail loudly with install instructions rather than skipping
  silently (`skip_on_cran()` retained, as it is inert off CRAN infrastructure).

---

## 9. Long-term systematic passes **[A]**

Whole-codebase sweeps, to be done deliberately rather than opportunistically.
Deferred until the class design settles, since it will dictate several of the
names.

### 9.0 Governing convention — rOpenSci package API guidelines **[A]**

Source: <https://devguide.ropensci.org/pkg_building.html#package-api>. Adopted
as the reference standard for 9.1, 9.2, and — importantly — for the generics
design in §5, which has to be decided *before* any further renaming, since the
class design dictates several of the names.

The five rules, each audited against the current `dev` state:

1. **`object_verb()` naming scheme** for functions sharing a data type or API.
   *Status: inconsistent, but do not act yet.* The package currently splits:
   - object-first: `talomes_heatmap()`, `msa_heatmap()`, `tales_consensus()`,
     `tales_consensus_match()`, `tale_parts_to_rvd()`, `repeat_to_rvd_align()`,
     `repeat_to_rvd_map()`
   - verb-first: `plot_tales_msa()`, `plot_tale_composition()`,
     `plot_target_preds()`, `build_repeat_msa()`, `group_tales()`,
     `correct_tales()`, `diagnose_tale_parts()`, `tell_tales()`,
     `run_annotale_*()`, `split_list()`

   **The OOP transition dissolves most of this split rather than renaming
   through it.** Every verb-first name above is a verb applied to one of our
   prospective classes; under S3/S4/S7 those become *methods on a generic*
   (`plot()`, `autoplot()`, `summary()`) dispatching on the object, so the
   object moves out of the function name and into the signature. Renaming
   `plot_tales_msa()` → `msa_plot()` now would be churn we then undo. Decide
   §5 first; re-audit this rule against whatever survives as a plain function.

2. **Data/object as the first argument**, for pipe compatibility.
   *Status: broadly satisfied, two real violations* — both cases of the same
   pair of objects taken in opposite order, which also breaks rule 5:

   | | first arg | second arg |
   |---|---|---|
   | [`msa_heatmap()`](../R/msa.R#L297) | `tal_sim` | `repeat_align` |
   | [`plot_tales_msa()`](../R/msa.R#L630) | `repeat_align` | `tal_sim` |
   | [`.repeat_to_sim_align()`](../R/conversion.R#L237) | `repeat_align` | `repeat_sim` |
   | [`.repeat_to_cluster_align()`](../R/conversion.R#L270) | `repeat_sim` | `repeat_align` |

   The second pair is the worse of the two: the name reads
   `repeat_ -> _cluster_align`, so the alignment is the subject, yet the
   similarity matrix is passed first. §4 already marks `msa_heatmap()` as
   superseded, so the first pair may resolve by deletion rather than by
   reordering.

   Path-taking entry points (`tale_parts(telltale_dir)`,
   `run_annotale_*(fasta_file)`, `tell_tales(subject_file)`,
   `correct_tales(uncorrected_path)`, `functal(TALfile)`) are **not**
   violations — they are constructors/readers, and the file *is* the input.

3. **snake_case throughout.** Confirms the direction already taken for
   functions; 9.1 and 9.2 are the unfinished remainder.

4. **No name conflicts with base or popular packages.** *Status: clean.*
   Verified by set-intersecting all 24 exports against the exports of `base`,
   `stats`, `utils`, `ggplot2`, `dplyr`, `magrittr`, `data.table`, `tidyr`,
   `purrr` and `Biostrings` — zero collisions. Worth re-running after any
   rename pass, and worth keeping in mind for generic names specifically,
   where the point is to *deliberately* collide (i.e. register a method on an
   existing generic) rather than to shadow.

5. **Consistent argument naming and order across functions with similar
   inputs.** *Status: naming largely unified by the rename pass* (`tal_sim`,
   `repeat_align`, `repeat_sim`, `h_cut` are now used uniformly — except
   `h.cut`, see 9.1); *ordering is not*, per the table in rule 2. This rule is
   the one that makes 9.1 and 9.2 a single coordinated pass rather than two
   independent ones: the argument vocabulary and the column vocabulary have to
   agree, so that `repeat_align` the argument and `repeat_align` the column
   mean the same thing.

### 9.1 Argument names — finish the snake_case conversion

The rename pass converted the public API, but six functions still carry
non-conforming argument names (verified by `formals()` introspection, so this
is the complete list):

| function | argument(s) |
|---|---|
| `functal()` | `TALfile` |
| `.compute_match_string()` | `RVDSeq`, `EBESeq` |
| `.extract_seqs_from_hits()` | `DNAsequences` |
| `.annotale_to_quetal_rvd()` | `inputFile` |
| `.quetal_to_annotale_rvd()` | `inputFile` |
| `.repeat_to_cluster_align()` | `h.cut` — **dot.case**, and inconsistent with `h_cut` used everywhere else |

`h.cut` is the one with real bite: the same parameter is `h_cut` in `distalr()`,
`msa_heatmap()` and `plot_tales_msa()`, but `h.cut` in the internal that
actually performs the clustering.

### 9.2 Column names — adopt snake_case across all tables

Currently **no** column in `distalr()`'s output is snake_case. Present state:

| table | columns |
|---|---|
| `tale_parts` | `arrayID`, `domainType`, `positionInCrd`, `dnaSeq`, `sourceDirectory`, `positionInArray`, `aaSeq`, `rvd`, `seqnames`, `domCode` |
| `repeats.code` | `code`, `AA Seq`, `rvd` |
| `repeat.similarity` | `RepU2`, `RepU1`, `Dissim`, `Sim` |
| `tal.similarity` | `TAL1`, `TAL2`, `arlemScore`, `maxLength`, `normArlemScore`, `Sim` |
| `repeats.cluster` | `RepID`, `Rep_clust`, `Rep_order` |

Four different conventions coexist: camelCase (14 names), uppercase-acronym
(`TAL1`, `RepU1`, `RepID`), capitalised-snake (`Rep_clust`, `Rep_order`), and
one name **containing a literal space** (`AA Seq`, which forces backtick
quoting everywhere it is used).

This is the data-level counterpart to the code-level rename already completed,
and it is a prerequisite for unifying the two similarity tables into one class
(they currently differ in id-column naming *and* column order). It is also a
breaking change for anyone with scripts reading these columns, so it belongs
in the same release as the class work.

Note the interaction with 9.1: several *argument* names deliberately mirror
*column* names (`tale_parts`, `rvd_map`), so the two sweeps should agree on a
single vocabulary rather than be done independently.

### 9.3 Decide `@internal` vs `@noRd` per function

Currently near-absent as a policy: **2** `@noRd` tags in the whole package
([conversion.R:269](../R/conversion.R#L269),
[target_predictions.R:460](../R/target_predictions.R#L460) — and the latter is
written `##' @noRd`, a typo that stops roxygen seeing it at all), **0**
`@internal`, **0** `@keywords internal`.

The distinction to apply case by case across the ~31 non-exported functions:

- **`@noRd`** — no help page generated at all. Right for trivial plumbing whose
  roxygen block is really just a source comment.
- **`@keywords internal`** — a help page *is* generated and `R CMD check`
  validates it, but it is kept out of the package index and the pkgdown
  reference. Right for internals that are non-obvious enough to deserve real
  documentation (the three `.pairwise_align_*` backends, `.compute_match_string()`,
  the alignment converters), and it means examples and `@param` coverage get
  checked rather than silently rotting.

Related, and worth settling at the same time: which of the currently-exported
functions should stop being exported once `tale_parts` is central (see §2), and
whether the S4 class `annout` should remain `exportClasses`-tagged given it
never reaches the public pipeline.

### 9.4 Use `@family` wherever justified **[A]**

Current state: **0** `@family`, **0** `@seealso`, and exactly **two**
`\code{\link{}}` cross-references in the entire package
([AnnoTALE_QueTAL_functions_library.R:65](../R/AnnoTALE_QueTAL_functions_library.R#L65)
and [:125](../R/AnnoTALE_QueTAL_functions_library.R#L125)) — a hand-maintained
reciprocal pair between `run_annotale_predict()` and `run_annotale_build()`,
which is precisely the case `@family` automates.

`_pkgdown.yml` also has no `reference:` section, so the website emits one flat
alphabetical list of all 24 exports. **There is currently no grouping of the
API anywhere** — not in R help, not on the site.

One `@family` tag does double duty (verified against roxygen2):

```
#' @family plotting functions
```
emits into every member's Rd both
```
\seealso{Other plotting functions: \code{\link[=beta2]{beta2()}}}
\concept{plotting functions}
```

- the `\seealso{}` block is **bidirectional and auto-maintained** — adding a
  member updates every sibling, which hand-written `@seealso` cannot do;
- the `\concept{}` is what pkgdown's `has_concept("plotting functions")`
  selector reads, so the same tag can drive a grouped `reference:` index on the
  website instead of maintaining that list separately in `_pkgdown.yml`.

Provisional families (to be confirmed — **do not apply before §5**, since
functions that become methods on generics are documented differently and may
not warrant their own topic at all):

| family | members |
|---|---|
| TALE discovery | `tell_tales`, `correct_tales`, `tale_parts`, `diagnose_tale_parts` |
| external TALE tools | `run_annotale_predict`, `run_annotale_build`, `functal` |
| similarity and grouping | `distalr`, `group_tales` |
| repeat alignment | `build_repeat_msa`, `tales_consensus`, `tales_consensus_match` |
| format conversion | `repeat_to_rvd_align`, `repeat_to_rvd_map`, `repeat_to_rvd_map_distalr`, `tale_parts_to_rvd`, `split_list` |
| TALE plots | `plot_tales_msa`, `plot_tale_composition`, `talomes_heatmap`, `msa_heatmap` |
| target prediction | `talvez`, `preditale`, `plot_target_preds` |

Points to settle when applying:

- A function may carry **several** `@family` tags; `plot_target_preds()` is the
  obvious dual member (plots *and* target prediction).
- The family name is interpolated verbatim into "Other <name>:", so it must
  read as a plural noun phrase — "TALE plots" works, "plotting" does not.
- Interacts with 9.3: `@family` only has an effect on topics that generate an
  Rd, so anything marked `@noRd` is out of scope, while `@keywords internal`
  topics *can* carry a family.
- §4 marks `msa_heatmap()` as superseded. Either leave it out of the family so
  the grouping does not imply parity with `plot_tales_msa()`, or keep it in and
  lean on an explicit deprecation note.

### 9.5 Unify the user-messaging system **[A]**

Current state — four idioms coexisting, ~190 call sites:

| idiom | count | channel |
|---|---|---|
| `logger::log_*()` | 86 (info 24, debug 12, warn 16, error 34) | logger appender |
| `stop()` | 57 | condition |
| `stopifnot()` | 16 | condition |
| `warning()` | 10 | condition |
| `cat()` | 10 | **stdout** |
| `message()` | 7 | stderr |
| `print()` | 6 | **stdout** |
| `cli_*()` | 1 | condition |

`logger` is the *dominant* idiom, not a minority — concentrated in
[distalr.R](../R/distalr.R) (30) and [msa.R](../R/msa.R) (22). `cli` is
already in `Imports` but used in exactly one place,
[startup.R:3](../R/startup.R#L3).

#### Correctness problems, not just style **[V]**

These need fixing regardless of which system wins:

1. **22 bare `stop()` and 1 bare `warning()`** produce a condition whose
   message is the **empty string** (verified). The diagnostic text was sent to
   the logger on the preceding line, so the R condition itself carries nothing.
   Consequences: `tryCatch`/`conditionMessage` see nothing, `expect_error(regexp=)`
   cannot assert on them, and if the user's logger threshold or appender sends
   ERROR elsewhere the failure is **completely silent with a blank message**.
   The `log_error(...)` + `stop()` pairing at
   [distalr.R:66](../R/distalr.R#L66), [:149](../R/distalr.R#L149) and the
   `log_warn(...)` + `warning()` at [:144](../R/distalr.R#L144) are the pattern
   to look for.

2. **[distalr.R:563](../R/distalr.R#L563)** —
   `logger::log_errors() && stop("'aln_method' parameter must be either ...")`.
   `log_errors()` (plural) is not a predicate; it installs a global error
   handler via `globalCallingHandlers()`. The `&&` short-circuits on its return
   value, so **the intended message is unreachable** — the user passing a bad
   `aln_method` gets `should not be called with handlers on the stack`
   (verified). Also a hidden global side effect fired from inside a function.

3. **`logger` is never configured by the package** — no `log_threshold()`,
   `log_appender()`, `log_layout()`, and crucially no **namespace**. tantale
   therefore reads and writes the *user's global* logger settings: someone who
   sets `log_threshold(WARN)` for their own code silently loses tantale's
   progress output, and tantale cannot adjust verbosity without stomping their
   configuration. The `namespace =` argument exists for exactly this and is
   unused.

4. **`cat()`/`print()` write to stdout**, so they cannot be silenced with
   `suppressMessages()` and they contaminate captured output. In a package,
   user-facing narration belongs on stderr.

5. **`import(cli)` and `import(logger)`** in NAMESPACE are whole-package
   imports (~200 symbols from cli alone). Call sites are also inconsistently
   qualified — `logger::log_error()` at [distalr.R:66](../R/distalr.R#L66)
   versus bare `log_error()` at [:107](../R/distalr.R#L107). Should be
   `importFrom` or fully qualified.

6. [startup.R:9](../R/startup.R#L9) wraps `cli_inform(class = "packageStartupMessage")`
   inside `packageStartupMessage()`. `cli_inform()` already signals the
   condition and returns `NULL`, so the outer call re-emits an empty message.

#### The framing to apply

"logger vs cli" is not quite the right axis — they solve different problems:

- **logger** is a *logging framework*: severity thresholds, appenders
  (destinations, incl. files), namespaces, layouts. For diagnostics that are
  filtered by level and may persist.
- **cli** is *console UI + condition signalling*: semantic bullets, glue
  interpolation, pluralisation, progress bars, and `cli_abort()` / `cli_warn()`
  / `cli_inform()`, which are classed-condition wrappers.

So split the ~190 calls by **what they are**, not by which package is fashionable:

| kind | count | target |
|---|---|---|
| conditions the caller might catch | 57 `stop` + 10 `warning` + 16 `stopifnot` | `cli::cli_abort()` / `cli::cli_warn()` |
| narration of what is happening | `cat`, `print`, `message`, most `log_info` | `cli::cli_inform()` / `cli_alert_*()` |
| long-running external tool progress | parts of distalr/telltale | `cli::cli_progress_bar()` |
| level-filtered developer diagnostics | `log_debug`, table dumps via `skip_formatter(kable(...))` | the open question below |

Converting the condition-signalling third is unambiguous and carries the most
value: it makes the 23 empty messages impossible by construction, and classed
conditions make errors assertable in the test suite (§8) via
`expect_error(class = )` instead of brittle regex matching.

#### The one real decision

Whether to keep `logger` at all, for the level-filtered debug channel.

- **Drop it** — one system, one mental model; removes a dependency; `cli` covers
  narration and progress better. Cost: lose threshold filtering and the
  `skip_formatter(kable(...))` table dumps, which have no direct cli
  equivalent (`cli_verbatim()` is the closest).
- **Keep it, namespaced** — `logger::log_threshold(..., namespace = "tantale")`
  and debug-level only, with **every** user-facing message and every condition
  moved to cli. Justified if a persistent log file of a long pipeline run is
  genuinely wanted.

Either is defensible; what is not defensible is the current state, where the
two overlap with no boundary. Recommendation: **cli for everything
user-facing, and drop logger unless the log-file use case is real** — the
package's own usage is dominated by narration and errors, not by diagnostics
anyone filters by level.

---

## 10. Explicitly ruled out

- Deleting dormant internals such as the unused HMMER wrappers
  (`.write_hmm_file()`, `.run_hmmer_search()`, `.run_hmmalign()`,
  `.extract_seqs_from_hits()`). Obsolete now, plausibly useful later.
- Treating `run_annotale_predict()` / `run_annotale_build()` as dead. They have
  no internal call sites but are legitimate standalone user utilities — the type
  case for "useful outside the `pipeline.svg` workflow".
