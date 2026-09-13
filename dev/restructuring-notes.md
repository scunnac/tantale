# tantale — restructuring notes and action ledger

Working document for the pre-publication overhaul. Records findings, agreed
actions and deferred questions so they don't live only in conversation.

Branch: `dev`. Last updated: 2026-09-13.

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
  it onto `tales` validation rather than losing it.
  **[V]** The assertion is now pinned by a test
  (`test_untested_exports.R`), which feeds it a repeat code carrying two
  different RVDs and confirms it errors. So retiring the function can no longer
  drop the check silently: the test will fail until the invariant has a new
  home. Note the `tales` validator already has invariant 8, the `aa_seq` ↔
  `dom_code` bijection, which is the same *shape* of constraint one level down —
  that is the natural place for it.
- **[A]** `repeat_to_rvd_map_distalr()` survives, but its name is misleading: it
  depends on `domCode` being present, not on `distalr()` having been run.

### Keep internal

- **[V]** `.repeat_to_sim_align()`, `.repeat_to_cluster_align()` — genuinely
  plotting-only internals, correctly unexported.

### Dormant but valuable — REPAIRED **[V]**

Both repairs are done; the "decide whether to wire them up" question stays open.

**`.rvd_to_match_align()`** — the `tantale::rvdSimDf` default is fixed to
`rvdSimDf`. **[V]** The bug was real and total: `rvdSimDf` lives in
`R/sysdata.rda`, so `tantale::rvdSimDf` raises *"not an exported object from
'namespace:tantale'"*. The function could never have run as written, in any
version. It now executes (checked on a small hand-built matrix). Note that only
*execution* is verified -- its semantics are still unexercised by any caller.

**`.rvd_to_repeat_align()`** — the missing length guard is added. The
back-mapping is positional (k-th non-gap cell = k-th repeat code), which is
well defined only if the counts agree; a mismatch previously misregistered
codes against positions silently. Two guards now fire, both tested: a row of
the alignment with no entry in `repeat_vecs`, and a row whose non-gap count
disagrees with its repeat-vector length.

Both also gained roxygen with `@keywords internal`, so `R CMD check` validates
their documentation rather than letting it rot.

#### Original notes

### Dormant but valuable — original notes

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

## 3. Legacy cemetery (`inst/legacy/`) — DONE **[V]**

The three AnnoTALE<->QueTAL shims are moved to
`inst/legacy/annotale_quetal_shims.R`, joining the existing
`inst/legacy/tellTaleLegacy.R`.

**[V]** Verified uncalled before moving: zero references to
`.reformat_array_report()`, `.annotale_to_quetal_rvd()` or
`.quetal_to_annotale_rvd()` anywhere in `R/`, `tests/` or `vignettes/` outside
their own definitions.

Moved rather than deleted. The file conventions they encode -- how AnnoTALE and
QueTAL disagree about writing the same RVD content -- are not recorded anywhere
else, and that is worth keeping even though nothing calls them. `inst/legacy/`
ships with the package but is not sourced, so they cost nothing at load time.

#### Original notes

## 3-original Legacy cemetery

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

## 4. Plotting: `msa_heatmap()` is superseded — GAP IDENTIFIED **[V]**

The ledger asked what `msa_heatmap()` provides that `plot_tales_msa()` does
not. Answered.

**Most of the apparent difference is not a difference.** `msa_heatmap()`'s six
`plot_type` values are combinations of orthogonal features that
`plot_tales_msa()` exposes as separate arguments -- which is the better design:

| `plot_type` | `plot_tales_msa()` equivalent |
|---|---|
| `repeat.clusters` | `fill_type = "repeat_clust"` |
| `repeat.similarity` | `fill_type = "repeat_sim"` |
| `with.rvd` | pass `rvd_align` |
| `repeat.clusters.with.rvd` | both of the above |
| `reference` | pass `ref_pattern` |
| `consensus` | `consensus = TRUE` -- **but see below** |

**The one real blocker [V]: `consensus` does not work in `plot_tales_msa()`.**
The argument is accepted and documented as *"NOT IMPLEMENTED YET"*; the body
carries `#### TODO: Bind a 'consensus' tibble or a consensus plot if requested`.
`msa_heatmap()` does render a consensus row.

Note the groundwork is already there: `plot_tales_msa()` computes
`tales_consensus(rvd_align)` and builds a tibble of it, using it to colour
matches. What is missing is only *displaying* it as a row.

**Smaller gaps**, all arguably out of scope for a ggplot function:

- `save_path` -- writes the plot to a file. A ggplot is returned as an object,
  so `ggsave()` covers this; not a real gap.
- `note_colors` -- customises the matched/mismatched colours. A ggplot caller
  adds a scale instead; the `p2` example in `p2_multiple_alignments.Rmd` shows
  exactly that.
- `...` passed to `gplots::heatmap.2`. Nothing equivalent, and nothing should be.

**Conclusion.** `msa_heatmap()` cannot be retired until `consensus` display
lands in `plot_tales_msa()`. That is the single prerequisite, and it is a
contained piece of work. Everything else it offers is either already present or
better served by returning a ggplot.

I did not implement it unattended: it is a visual feature whose result needs a
human eye, not a passing test.

#### Original notes

## 4-original Plotting: `msa_heatmap()` is superseded

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

### 5.1 Does `diagnose_tale_parts()` survive the `tales` class? **[P]**

**[V]** It checks three things — rows with `NA` in `aaSeq`, `dnaSeq` or `rvd` —
and has two modes: report the offending arrays (default), or *remove* them
(`sanitize = TRUE`).

Against the `tales` contract (`class-design.md` §2.4):

| check | status under the class |
|---|---|
| `rvd` is `NA` | **redundant** — hard invariant 3 errors on an `NA` residue column, so a valid `tales` cannot reach the check |
| `aa_seq` is `NA` | **not covered** — deliberately a *precondition* of `tales_relatedness()`, not an invariant, since a `tales` built from RVD strings has no `aa_seq` at all |
| `dna_seq` is `NA` | **partly** — soft invariant 10 warns, does not error |

So it is not made redundant, but its remit has narrowed to two things no
invariant provides:

1. **Array-level triage.** Validation is per row; this reports the *whole
   array* when any one of its parts lacks a sequence — which is the right
   granularity, since a partial array cannot be aligned.
2. **Repair.** `sanitize = TRUE` drops the bad arrays so the rest can proceed.
   The class errors instead; it has no "carry on without the broken ones" mode.

**To decide:** whether to keep it as an explicitly-named triage/repair tool
(dropping the now-unreachable `rvd` branch, and renamed — it is not a
validator, and sharing vocabulary with `validate_tales()` would mislead), or
to fold the `sanitize` behaviour into a `tales` helper and retire the rest.
Note both call sites use it as a *guard* (`conversion.R:205`, `:373`), which
is a third use again — and that guard is stricter than the class, since it
rejects the whole input if any array is affected.

It also carries a bare `warning()` with an empty message, one of the §9.5
cases.

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
- **[V] FIXED** — `plot_tales_msa()` aborted on **every** call under ggplot2
  4.0.3: `msa.R:830` passed `palette =` to `ggplot2::scale_fill_manual()`,
  which has no such argument (it takes `values`), so the name collided with the
  `palette` that `discrete_scale()` supplies internally —
  *"formal argument 'palette' matched by multiple actual arguments"*. The scale
  was built before the `fill_type` branch, so both branches died; verified by
  calling the function directly on the package's own fixture matrix. Replaced
  with `discrete_scale()`, which is the scale that actually accepts a palette
  *function*. This is ggplot2 API drift of the same kind as the Bioconductor
  drift found in Phase 1. It also explains the assertion-free
  `test_plot_tales_msa.R`: nothing there *could* have asserted, since every
  call aborted — that file now has real assertions, including that the plot
  renders to a file.
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

### RESOLVED: the as.dist() inversion bug **[V]**

Closed. Of the sites that looked suspect, only one was a live bug:

| site | verdict |
|---|---|
| `msa.R` (both dendrogram sites) | correct -- they compute `100 - tal_sim` first |
| `classification.R` | correct -- its matrix is now a distance |
| `conversion.R` `.repeat_to_cluster_align()` | **was the bug**; fixed |
| `distalr.R` `.cluster_repeats()` | same defect, but unreachable after `distalr()` was removed; deleted |

Measured effect of the fix on the reference output: 45 clusters instead of 59,
93.8% pair-agreement. `h_cut` defaults moved 90 -> 10 in both plotting
functions, since the height is now read on a distance scale, and the dot.case
`h.cut` argument was renamed `h_cut` at the same time (part of 9.1).

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

### 7.1 Vignettes 1-4 cannot be built **[V]**

Found while removing the deprecated functions, and it predates that work.

`R CMD build` (with vignettes) **fails today**, at `2_tale_classification.Rmd`.
The cause is not an API problem: vignettes 1-4 are chained through session
state written to a hardcoded user path.

- `1_tale_mining.Rmd` ends with `save.image(file.path(outdir, "mining.RData"))`
  where `outdir <- fs::dir_create("~/TEMP/test_tantale")`.
- `2`, `3` and `4` each begin with `load(file.path(outdir, "mining.RData"))`.
- That file is absent on a clean machine, so 2-4 error immediately.

Consequences:
- `R CMD build` must be run with `--no-build-vignettes` to succeed at all.
- This, not anything about pkgdown itself, is the real blocker behind the
  "pkgdown rebuild blocked on vignette reproducibility" note.
- The four vignettes could not be verified after the API migration, since they
  cannot execute. Their call sites were updated mechanically and are unchecked.

By contrast `p1_tales_compare.Rmd` and `p2_multiple_alignments.Rmd` are
self-contained -- they read fixtures from `inst/extdata` -- and both were
re-knitted successfully after migration.

The fix is to make each vignette stand alone: read its inputs from
`inst/extdata` (or build them in a setup chunk) rather than inheriting a
`save.image()` from the previous one. Worth doing before any pkgdown rebuild.

---

### 7.2 `R CMD check` results **[V]**

First full check of the package (with `_R_CHECK_FORCE_SUGGESTS_=false`, since
`ggcorrplot` and `corrr` are not installed here -- an environment gap, not a
package defect).

**Fixed tonight:**

| finding | fix |
|---|---|
| NOTE: `exportClasses(annout)` requires `methods` | added `methods` to Imports |
| WARNING: `::` imports not declared -- `BiocGenerics`, `BiocParallel`, `GenomeInfoDb`, `S4Vectors`, `rtracklayer` | all five added to Imports (all were already installed, so this was purely a declaration gap) |
| WARNING: malformed cross-reference `\link[Biostrings::BStringSet]{BStringSet}` | corrected to `\link[Biostrings]{BStringSet}` |
| WARNING: `correct_tales.Rd` documents `telltale_dir`, which is not an argument | stale `@param` removed |
| WARNING: undocumented `rvd_vecs` in `repeat_to_rvd_map.Rd` | documented |
| WARNING: undocumented `plot_tree`, `k`, `k_range`, `method` in `tales_group.Rd` | restored from the deleted `group_tales()`, which had them |
| WARNING: undocumented `plot_type` in `talomes_heatmap.Rd` | documented |

Note the `tales_group.Rd` gap was **pre-existing**, not caused by the removal:
`tales_group()` never carried those `@param` tags, while the `group_tales()`
alias it superseded did. Deleting the alias merely made the omission visible.

**Second pass** brought it to 5 WARNINGs / 3 NOTEs, then fixed two more:

- `@param plot_type` had been inserted into the wrong roxygen block. Both
  `tales_group()` and `talomes_heatmap()` live in `classification.R`, and
  `talomes_heatmap()`'s block has no `@return`, so an "insert before the nearest
  preceding `@return`" heuristic landed it in `tales_group()`. Anchoring on
  `@export` instead fixed it. Third time tonight that a positional heuristic
  found a plausible-but-wrong target silently.
- `methods` was declared (needed for `exportClasses`) but never imported from.
  Added `@importFrom methods setClass` on the `annout` class definition, which
  is the only `methods` machinery the package uses.

**Left open:**

- ~~*Imports declared but not imported from*~~ **RESOLVED [V]**. All six were
  checked individually and all six were genuinely unused, so all six were
  removed:

  | package | evidence |
  |---|---|
  | `GenomicFeatures` | zero occurrences in `R/`, `inst/`, `vignettes/`, `tests/` |
  | `RColorBrewer` | zero occurrences anywhere |
  | `dichromat` | zero occurrences anywhere |
  | `optparse` | zero occurrences anywhere |
  | `scales` | one occurrence, inside a **commented-out** line (`msa.R:334`) |
  | `msa` | 38 textual hits, but **all of them our own identifiers** -- `tales_msa`, `plot_tales_msa`, `rvd_msa_by_group`, `.tidy_biostrings_msa`. No `msa::`, no bare call to any of its functions. Alignment shells out to MAFFT, not to this package. |

  Also confirmed for each: no `library()`/`require()`/`requireNamespace()` call,
  and no `@import`/`@importFrom` roxygen tag. `msa` and `GenomicFeatures` are
  substantial Bioconductor packages that users were being made to install for
  nothing.
- The vignette WARNINGs all trace to 7.1 (no `inst/doc`, because the vignettes
  cannot build).
- WARNINGs on executable files and non-portable file names -- these are the
  long test-fixture paths and bundled tool binaries, both known. The path
  lengths are the 60 over-long tar entries already recorded; fixing them means
  renaming fixture directories such as
  `tellTaleErrorMissingAnnotaleDnaDomain/`.
- ~~NOTE on `R code for possible problems`~~ **TRIAGED AND LARGELY FIXED [V]**.
  It had two halves and they were nothing alike:

  | half | count | verdict |
  |---|---|---|
  | "no visible global function definition" | ~40 | **every one real.** Not a single false positive. |
  | "no visible binding for global variable" | 126 | all NSE column names -- genuine false positives |

  The first half is how `plot_tale_composition()`'s breakage was found. It also
  turned up `methods::Quote`, used ten times in `telltale.R` and never
  imported. All are fixed by declaring what the code calls.

  The second half is now declared in `R/globals.R` (81 names). That is not
  cosmetic: 126 lines of noise are exactly what let ~40 real unresolved calls
  sit unread. A check output nobody can read is a check nobody runs.
- NOTE on `package subdirectories` -- relates to `inst/`.

---

### 7.3 A systemic habit: unqualified calls to non-imported packages **[V]**

Three separate defects tonight had the identical shape, which makes it a habit
rather than three accidents:

| where | call | consequence |
|---|---|---|
| `plot_tale_composition()` | `mutate()`, `ggplot()` | **failed outright** unless the user had dplyr attached |
| `telltale.R` (x10) | `Quote()` | `methods::Quote`, never imported |
| `msa_heatmap()` (`msa.R:356`) | `countMatches()` | `S4Vectors::countMatches` -- and `msa.R:43` gets it *right*, 300 lines earlier |

All three lived in code no test exercised. None was visible by reading -- the
calls look perfectly ordinary; only the namespace resolution is wrong.

**Root cause.** Listing a package in `Imports` makes it *installable*, not
*visible*. Without an `@import` or `@importFrom`, a bare call to it resolves
through the caller's search path, so the function works in an interactive
session where the user has done `library(dplyr)` and fails everywhere else.

**Review heuristic worth keeping: in this package, a bare call to anything
outside base is suspect.** `R CMD check`'s "no visible global function
definition" is the tool that finds them, and it is only readable once the NSE
false positives are declared away (see 7.2) -- which is the real argument for
`R/globals.R`.

---

## 8. Tests — error conditions now covered **[V]**

`tests/testthat/test_error_conditions.R` added (18 assertions). It exists
because 11 of the package's `tantale_error_*` classes had **zero** test
coverage, including six created during the cli conversion. A classed condition
that nothing asserts on buys nothing.

**It immediately earned its keep: it found a live bug in two error handlers.**

`.pairwise_distances_rename_legacy()` and `.tales_rename_legacy()` both build a
message whose `"i"` bullet carries a `{?s}` plural marker with no quantity to
count. cli processes each bullet separately, so it raised
*"Cannot pluralize without a quantity"* — as a plain `simpleError`. The
intended `tantale_error_*_name_clash` class was never signalled, and the user
saw a cli internals complaint instead of the real problem.

Fixed with an explicit `{cli::qty(clash)}`. Both now raise their proper class.

Worth remembering when writing cli messages: an inline style span such as
`{.fn tales}` is **not** a quantity. A first scan for this bug missed the
`tales_class.R` instance for exactly that reason.

#### Original notes

## 8-original Tests

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

### 9.1 Argument names — DONE **[V]**

Complete. Verified by `formals()` introspection over every function in the
namespace: no argument anywhere in the package now contains a capital letter or
a dot. The last six were `functal(TALfile)`, `.compute_match_string(RVDSeq,
EBESeq)`, `.extract_seqs_from_hits(DNAsequences)`, the two
`inputFile`s, and `.repeat_to_cluster_align(h.cut)` -- the one with real bite,
since it was dot.case *and* disagreed with `h_cut` everywhere else. Fixed as
part of the clustering bug fix.

#### Original notes

### 9.1-original Argument names — finish the snake_case conversion

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

#### Second category: inconsistent vocabulary, not casing **[P]**

These are all *already* snake_case, so the list above does not catch them, but
they violate rule 5 (consistent argument naming across functions taking similar
inputs) — and since `tale_sim` and `repeat_sim` are now **class names**
(`class-design.md` §3), the arguments and the classes disagree:

| name | where | exported? |
|---|---|---|
| `repeat_sim` | `msa_heatmap()`, `plot_tales_msa()` | yes — collides exactly with the `repeat_sim()` constructor |
| `repeat_sim` | `.repeat_to_sim_align()`, `.repeat_to_cluster_align()` | internal |
| `repeat_sims` | `build_repeat_msa()`, `tales_align()` | yes — singular/plural split for the same thing |
| `tal_sim` | `group_tales()`, `msa_heatmap()`, `plot_tales_msa()` | yes — one letter from the `tale_sim` class, same meaning |

**[V]** No correctness risk, verified: a parameter bound to a data frame does
not shadow a same-named function, because R skips non-function bindings when
resolving a symbol used in call position. This is a readability problem.

**[D]** `tales_align()`'s `repeat_sims` is *new* code that inherited the plural
from `build_repeat_msa()`, rather than inherited debt. Deliberately left
unrenamed so the whole vocabulary is settled in one pass here, per this
section's own warning that the argument and column sweeps should agree rather
than be done piecemeal. It has no users yet, so it is free to change.

The classes give the sweep a fixed point to converge on: whatever the arguments
become, they should agree with `tale_sim` / `repeat_sim` / `tales_msa` rather
than each other.

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

### 9.3 Decide `@internal` vs `@noRd` per function — PARTLY DONE **[V]**

Current state over 55 non-exported functions: 21 `@noRd`, 12
`@keywords internal`, and roughly two dozen with no roxygen at all.

**A finding that changes what this item means [V].** `@keywords internal` on
its own does *nothing*. roxygen generates no Rd for a block with no title, so
of the pre-existing `@keywords internal` tags on dot-prefixed helpers
(`.tales_check_key()` and its siblings), none produces a help page. They read
as a policy decision but have no effect: those functions are documented exactly
as if they carried `@noRd`.

So the real distinction is not `@noRd` vs `@keywords internal` — it is
**whether the block has a title at all**:

| block | Rd generated? | checked by `R CMD check`? |
|---|---|---|
| `@noRd`, with or without title | no | no |
| `@keywords internal`, **no title** | no | no |
| `@keywords internal`, **with title** | yes, as `man/dot-<name>.Rd`, hidden from the index | yes |

Only the third row buys anything. The two functions repaired in §2 use it
deliberately and are the first in the package to generate `man/dot-*.Rd`.

Note this does not conflict with the class-implementation session's §3.7
decision to use `@noRd` for its new helpers: that decision explicitly said it
did not pre-empt this item.

**Remaining work**, and why I left it: deciding which of the two dozen
undocumented internals deserve real, check-validated documentation is a
per-function judgement. Bulk-adding a bare `@noRd` would only restate what
already happens.

#### Original notes

### 9.3-original Decide `@internal` vs `@noRd` per function

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

### 9.4 Use `@family` wherever justified -- DONE **[V]**

42 `@family` tags added across eight families:

| family | n |
|---|---|
| tales objects | 8 |
| TALE alignment | 7 |
| tales projections | 6 |
| TALE plots | 5 |
| pairwise distances | 5 |
| target prediction | 4 |
| TALE discovery | 4 |
| external TALE tools | 3 |

Two implementation notes:

- Blocks carrying `@rdname` were deliberately skipped. They share a topic with
  their parent, so a tag on each would emit duplicate `\concept{}` entries into
  one Rd.
- `plot_target_preds()` carries two tags (target prediction *and* TALE plots),
  the dual membership anticipated when this item was written.

Each tag emits both a bidirectional `\seealso{Other <family>: ...}` and a
`\concept{<family>}`, so a grouped `reference:` section can now be added to
`_pkgdown.yml` with `has_concept("<family>")` rather than maintaining the list
separately. That pkgdown change is **not** done and is the natural follow-up.

#### Original notes

### 9.4-original Use `@family` wherever justified

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

### 9.5 Unify the user-messaging system -- DONE **[V]**

Executed. `logger` is removed from `R/` and from `DESCRIPTION`; `cli` is the
single messaging system.

| before | after |
|---|---|
| 87 `logger::log_*()` | 0 |
| 22 bare `stop()` + 1 bare `warning()` (empty messages) | 0 |
| `log_errors() && stop(...)` (message unreachable) | `cli_abort()` naming the bad value |
| `cat()` / `print()` narration on stdout | `cli_inform()` on stderr |

Conversion rules applied:

- `log_error(msg)` immediately followed by `stop()` -> `cli_abort(msg, class = "tantale_error")`.
  This is what fixes the empty-message class of bug: the text was going to the
  logger while the condition itself carried nothing.
- `log_error(msg)` *not* followed by `stop()` -> `cli_warn()`, since it never aborted.
- `log_warn(msg)` + `warning()` -> `cli_warn(msg)`.
- `log_info()` -> `cli_inform()`.
- `log_debug()` -> **deleted**. logger's default threshold is INFO, so these
  were already invisible; removing them changes nothing a user could see. The
  `skip_formatter(kable(...))` table dumps went with them, their information
  folded into the neighbouring `cli_warn()` bullets where it mattered.

Open follow-up: converted legacy sites carry the generic `tantale_error` class
only. Giving them specific subclasses (as the class-system code already does,
e.g. `tantale_error_tales_type`) would make them individually assertable in
tests. Worth a pass, but each needs a judgement call about the right name.

### 9.5-superseded Unify the user-messaging system **[A]**

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

### 9.6 "similarity" and "repeat" are both wrong names -- DONE **[V]**

Both halves resolved and implemented.

**(a)** confirmed by measurement (71 of 251 ids, 28%, are terminus domains) and
fixed by renaming to the `distances` vocabulary:
`pairwise_sim`/`tale_sim`/`repeat_sim` -> `pairwise_distances`/`tale_distances`/`domain_distances`,
`tales_relatedness()` -> `tales_compare()`, helpers `sim_*` -> `distances_*`.

**(b)** resolved in favour of storing the distance. `dissim` is the required
column; `sim` and `norm_arlem_score` are folded in on ingest and dropped, since
all three were exact restatements of one quantity. Legacy tables are still
accepted and converted.

Naming decisions taken with the user, in order: the function keeps a verb
(`tales_compare`), the slots mirror the class labels, the constructors mirror
them too, and the entity word stays singular while the head noun is plural
(`tale_distances`, not `tales_distances`) -- both because English puts the
modifier in the singular and because `tales_*` already means "operates on a
tales object" in this package.

#### Original notes

### 9.6-original "similarity" and "repeat" are both wrong names

Two independent problems with the same family of names — the class
`repeat_sim`, the slot `repeat.similarity`, the column `sim`, and the
functions and arguments built on them. They want deciding together, because
any fix touches class names, function names *and* column names at once.

#### (a) "repeat" understates what the object covers **[V]**

`repeat.similarity` compares **every domain type, not just repeats**. Measured
on the fixture: of its 251 ids, **71 (28%) are terminus domains** and 180 are
repeats — and the two sets are disjoint, so no id is both. Termini are
`dom_code`-ed like any other part and participate in the alignment
(`class-design.md` §4.5), so this is by design, not an accident.

`domain` is the accurate word and is already the vocabulary elsewhere —
`domain_type`, `dom_code`. So: `repeat_sim` → `domain_sim`,
`repeat.similarity` → `domain_sim`, and the `repeat_*` arguments of §9.1
follow. Note `build_repeat_msa()`/`tales_align()` are *also* misnamed on the
same grounds, since they align termini too.

#### (b) similarity may be the wrong quantity to store **[P]**

The ARLEM analysis in §6 already established that `arlemScore` is a cost, that
`tal.similarity$Sim` is crushed into 91.75–100, and that **all three consumers
immediately convert it back to a distance** — "a round trip that costs
interpretability and buys nothing".

That argues for storing dissimilarity. But §6 also drew a distinction worth
keeping: at the **array** level `Sim` is near-useless, while at the **domain**
level `Sim = 100 - Dissim` is a genuine full-range similarity. So this is a
strong case for `tale_sim` and a real choice for `domain_sim`, not one verdict.

**Consequence for the class, if dissimilarity wins.** `pairwise_sim` currently
requires `sim` and treats `dissim` as optional (`class-design.md` §3.3). That
inversion is the substantive part of this change — the renames are mechanical,
but which quantity is *required* is a contract decision. Whichever is chosen,
only one should be stored: keeping both invites them drifting out of step.

**Also to decide:** whether the class name should say the quantity at all.
`domain_sim`/`tale_sim` presume similarity; `domain_dist`, or a neutral
`domain_relatedness`, would not. A neutral name would survive changing the
stored quantity later.

---

## 10. Explicitly ruled out

- Deleting dormant internals such as the unused HMMER wrappers
  (`.write_hmm_file()`, `.run_hmmer_search()`, `.run_hmmalign()`,
  `.extract_seqs_from_hits()`). Obsolete now, plausibly useful later.
- Treating `run_annotale_predict()` / `run_annotale_build()` as dead. They have
  no internal call sites but are legitimate standalone user utilities — the type
  case for "useful outside the `pipeline.svg` workflow".
