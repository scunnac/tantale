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

### Dormant but valuable — REPAIRED AND WIRED UP **[V]**

Both are now reachable, and `rvdSimDf` has two independent consumers rather
than none.

**`.rvd_to_match_align()` is `fill_type = "rvd_sim"`.** It colours each cell by
how alike that RVD's DNA-binding preference is to the reference TALE's RVD at
the same position -- the RVD-level counterpart of `repeat_sim`, which scores
protein sequence. The two genuinely differ: `HD` and `ND` are distinct repeats
with identical specificity, while repeats differing only at 12-13 are
near-identical proteins targeting different bases.

**[V]** Working example on the fixture: at position 23 the reference carries
`NG` (T-binder, 5/10/1/50) and MAI1 carries `NN` (A/G-binder, 30/10/30/1),
scoring **-0.95**. A repeat-level fill renders that as merely "a different
repeat"; the RVD fill shows the specificities are opposed.

It needs `rvd_align` but **not** `repeat_sim`, so it works without having run
`tales_compare()`. It gets a diverging scale centred on zero, since the score
is signed on [-1, 1] -- the sequential 0-100 palette would flatten "opposite"
and "somewhat different" together. Cells with no score render grey. On the reference fixture those are the
termini (`NTERM`/`CTERM`), which are not RVDs and so have no specificity to
compare -- not `XX`, which is a genuine RVD with a uniform profile. The cell
keeps its own label either way, so nothing reads as `NA`.

**`rvdSimDf` also scores RVD alignments** -- see 7.5. That use came out of
noticing that RVD alignments had no scoring matrix at all.

The earlier instinct to delete these on call-count alone was wrong in a
specific way worth remembering: they were not dead, they were **unwired**, and
the capability existed nowhere else.

#### Original notes

### Dormant but valuable — original notes

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

## 4. Plotting: `msa_heatmap()` — RETIRED **[V]**

Done. Moved to `inst/legacy/msa_heatmap.R` (345 lines) once the consensus
prerequisite was implemented in `plot_tales_msa()`. Vignette 3's three calls
were migrated: `plot_type = "repeat.similarity"` -> `fill_type = "repeat_sim"`,
`"repeat.clusters"` -> `fill_type = "repeat_clust"`, and
`"repeat.clusters.with.rvd"` -> the same plus `rvd_align` and `consensus = TRUE`.

The internals it used (`.repeat_to_sim_align()`, `.repeat_to_cluster_align()`,
`.pick_ref_name()`) are shared with `plot_tales_msa()` and stay.

#### The analysis that justified it


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

**Conclusion — PREREQUISITE NOW MET [V].** `consensus` display is implemented
in `plot_tales_msa()`, so nothing blocks retiring `msa_heatmap()` any more.

Implementation note worth keeping, because the obvious approach cannot work:
the consensus is a **separate `aplot` panel**, not an extra row of the
alignment. `aplot::insert_left()` reorders the main plot's y axis onto the
tree's leaves, and a y level with no matching leaf is *silently dropped* --
measured, the composed y levels came back `NA | 1 | 2 | 3 | NA` with the
consensus row simply absent. There is no variant of "just add a row" that
survives the tree.

Two further details:

- The consensus follows whatever the cells are labelled with: taken from
  `rvd_align` when supplied, `repeat_align` otherwise, with the same 3-character
  padding the cells use for `dom_code`.
- `aplot`'s `height` is a *ratio* of the main plot, so a fixed value grows with
  the array count -- several rows tall for a large group. It is
  `1.0 / countOfTales`, clamped, which holds the consensus at about one row
  whatever the count. Verified at 3 and 12 arrays.

**[V]** `tales_consensus()` was cross-checked against an independent
`table()`-based mode calculation over all 28 positions of the fixture: identical.

What remains before `msa_heatmap()` can actually go is only the decision, plus
`save_path`/`note_colors` having no ggplot equivalent -- and §4 already argues
`ggsave()` and an added scale cover those.

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

### 5.1 Does `diagnose_tale_parts()` survive the `tales` class? — RETIRED **[V]**

**Resolved.** The function is gone; `tales_anomalies()` reports the same
conditions, `tales(sanitize = TRUE)` drops the offending arrays, and the
checks are part of the class rather than a separate diagnostic. Nothing in
`R/` or `NAMESPACE` mentions it.

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

### 5.2 Talome-wide MSA summary plot **[P]**

The idea: show every group's alignment in one figure, faceted by group, using
`plot_tales_composition(position = "alignment")`. That layout puts domain type
and amino-acid length onto the alignment coordinate, which neither existing
plot does -- aberrant repeats line up as a column instead of scattering.

**Prerequisite, already done:** `group` is now a recognised `tales` column,
validated as constant within an array (it is an array-level property, like
`seqnames`). It was always usable as a free column; what is new is the
invariant. The name matches what `tales_group()` already returns.

#### The dilemma: what does the summary function take?

**[V]** Each group aligns independently, so each alignment has its own width
*and its own coordinate system*. Measured on the fixture: group 1 is 28 columns
numbered 1-28, group 2 is 22 columns numbered 1-22.

A single concatenated `tales_msa` is therefore **structurally valid but
semantically incoherent**: `alignment_position = 5` means unrelated things in
different groups, `alignment_width` is one scalar that cannot describe two
alignments, and `as.matrix()` would build a single grid spanning both.

Two honest shapes, undecided:

| | shape | cost |
|---|---|---|
| **A** | a list of `tales_msa`, one per group | honest about there being N alignments; needs a `tales_msa_list` type or just a plain list |
| **B** | concatenate to a plain `tales`, carrying `group` and `alignment_position` as ordinary columns, and facet with `scales = "free"` | closest to the original idea; the object stops claiming to be one alignment, which is the accurate claim |

B is reachable today: demoting a `tales_msa` with `as_tales()` keeps
`alignment_position`, and `position = "alignment"` consumes it. Note the facet
currently uses `scales = "free_y", space = "free"`, which **shares** the x
axis -- per-group alignments need `free_x` too, or the narrower group is padded
out to the wider one's width.

#### Decided: how `group` gets populated **[V]**

Maintainer's call, and the right one: `tales_group()` takes the `tales`
object whose comparison produced `tal_sim` and returns it with `group`
filled, rather than a separate `add_group()` combining a bare mapping with a
`tales` after the fact.

Implemented as `tales_group(x, tal_sim, ...)`, returning `x` with `group`
added. Previously it took `tal_sim` alone and returned a
`data.frame(name, group)`.

The argument order is a judgement made when implementing, not something the
maintainer specified: `x` first, following the rOpenSci data-first
convention (§9.0) and the rest of the `tales_*` API. Trivially flipped.

**Why the combining variant is worse, concretely.** Taking `x` is the only
point at which the correspondence between the distances and the object can
be checked. `tales_group()` now errors (`tantale_error_group_mismatch`) if
any array in `x` is ungrouped or any grouped name is absent from `x`, which
means "these distances did not come from this object". An `add_group()`
called later could do the same check, but nothing would *oblige* the caller
to route through it, and the failure mode it prevents -- a partly-grouped
object -- surfaces far downstream and confusingly.

The bare mapping is not lost: `unique(out[c("array_id", "group")])`.

**Still open:** this settles how `group` is populated, not the A-vs-B
question below, which is about what a multi-group *alignment* is.

#### Connected: a group-aware `tales_align()`

Rather than making the user loop, `tales_align()` could notice a `group` column
and align each group separately, returning either a reassembled object or a
list. That is the natural home for the loop.

**Blocked on a prerequisite that does not exist yet:** there is no `c()` or
`bind_rows()` method for `tales` or `tales_msa`. Reassembly needs one, and it
is not trivial:

- **[V]** `array_id` uniqueness is a hard invariant, so a bind must reject
  colliding ids rather than silently fanning out.
- The `dom_code` namespace attribute must agree across the parts being bound --
  §3.5 exists precisely because mixing runs joins wrongly and silently.
- For `tales_msa`, `alignment_width` cannot survive a bind of two alignments,
  which is the same incoherence as above; a bind would have to demote to
  `tales`, or be refused.

So the ordering is: decide A vs B, then add the bind method the chosen shape
needs, then make `tales_align()` group-aware. Not before.

---

### 5.3 `tell_tales()` refactoring — MECHANICAL PASS DONE **[V]**

| | before | after |
|---|---|---|
| lines in `tell_tales()` | 745 | **160** |
| deepest indentation | 65 | 48 |
| `if`/`else` branches | 18 | 4 |
| `for` loops | 5 | 0 |
| arguments | 17 | 17 (untouched) |

Seventeen internals, each named for what it does, and the body now reads as
a pipeline: prepare the subject, read the profiles, find the hits, put them
on the genome, merge, group into arrays, find the ORFs, run AnnoTALE, finish
the RVD strings, align the termini, measure, report, log.

**Verification.** Every extraction was checked the same way, not assumed: a
baseline captured on the refactored code was re-run against the
pre-extraction commit (`git stash push R/telltale.R`). Passing in both
directions means the two produce identical output. Every step also ran the
full suite, and each is its own commit, so any one of them reverts alone.

**Two defects found, neither of which any test would have caught**

1. The GFF export decided whether to include the unmerged hits by testing
   `exists("reducedOlapGr")` -- an intermediate variable of the merge branch.
   Lifting that branch into a function removed the variable from the frame,
   `exists()` silently became `FALSE`, and `allRanges.gff` lost half its
   records: 199 lines to 103. It asks `merge_hits` now.
2. Removing the `annout` S4 class removed an `@import Biostrings` that its
   roxygen block had been carrying for the whole package. Four unqualified
   calls depended on it; the worst was `nchar()` on an `XStringSet`, which
   reads as base R and only differs for an S4 argument. The import is now
   declared deliberately (7.3 still wants it narrowed).

**Also fixed on the way:** the unguarded `min_domain_hits` filter (8.0); the
three terminus anchor codes, hardcoded here as literals, now read from
`tales_anchor_codes()`; the nested `AnnoTALEanalyze()` promoted to
`.run_annotale_analyze()` (8.0b); the DNA and protein terminus alignments,
written twice, unified; `methods::Quote()` dropped from the package imports.

**Not done, deliberately.** The 17 arguments are untouched: grouping them is
a judgement call, it interacts with 9.1, and it changes the user-facing
interface, which the rest of this pass did not. Two `TODO` blocks remain in
the body -- circular molecules, and what two output files should contain --
both of which are questions for the maintainer rather than cleanups.

#### Original notes

### 5.3-original `tell_tales()` needs refactoring **[A]**

**[V]** Measured, so the scale is on record rather than impressionistic:

| | |
|---|---|
| lines in the one function | **745** |
| lines in the whole of `R/telltale.R` | 903 |
| arguments | **17** |
| `if`/`else` branches | 18 |
| `for` loops | 5 |
| deepest indentation | **65 spaces** |

So a single function is 82% of its file, and at 65 spaces of indent the tail of
it is unreadable at normal width. It is the entry point of the whole pipeline
and the least tractable code in the package.

It leans on only four package internals -- `.run_nhmmer_search()`,
`.hits_report_to_gff()`, `.correction_tibble()`, plus one `system()` call --
which suggests the bulk is inline orchestration that could be named and lifted
out rather than genuinely irreducible logic.

Not attempted yet. Worth noting the ordering against other items: it produces
the input to `tales_from_telltale()`, and §5.1's `sanitize` work exists
precisely because its output is sometimes malformed. Straightening it might
reduce how much sanitising is needed, but the two are independent -- sanitising
guards against bad input whatever its provenance, so neither blocks the other.

A sensible first pass would be purely mechanical: extract the named stages into
internals, without changing behaviour, and get the indentation down. The 17
arguments are a separate question and interact with 9.1.

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

### 7.1 Vignettes 1-4 cannot be built **[V]** — DEFERRED BY POLICY

**Vignettes do not constrain the code.** They will be rebuilt from the
finished functionality and the sharpened interface, not the other way round.
A broken vignette is not a regression to chase, and no API decision should be
made to keep one knitting. This section and 7.5 are downstream of everything
else here.

The diagnosis below stands, for whenever that rebuild happens.


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

### 7.4 Shrink the payload — DONE **[V]**

`inst/tools` is **122 MB -> 63 MB**. MAFFT and HMMER now come from the
`tantale` conda environment, which already existed and already declared both.

**The versions were checked, not assumed.**

- **HMMER 3.3 -> 3.3.2.** Whole pipeline run both ways: of 36 output files, 2
  differed, and within those exactly one line -- the version banner. Every
  hit identical.
- **MAFFT 7.450 -> 7.520 changed results**, reproducibly and in both
  directions (28 vs 29 columns on `dom_code`, 29 vs 28 on `rvd`). Worse,
  7.520 left the termini unanchored: with 7.450 every array starts at column
  1 with its N-terminus and ends at the last column with its C-terminus,
  while 7.520 staggered them. Gap penalties do not explain it -- `--op` 0-5
  and `--ep` 1-10 all gave the same wrong answer. `--globalpair --op 1 --ep 1`
  recovers the column count and the anchoring, and is identical for
  well-populated arrays, but still differs on sparse ones.
- **MAFFT 7.453 from bioconda is byte-identical to the bundled 7.450**, on
  both layers. That is what the yaml pins, with a comment saying why.

**What this changes for users.** The core pipeline no longer works offline
out of the box: the first `tales_align()` or `tell_tales()` creates the conda
environment. Conda was already a stated requirement, and TALVEZ and functal
already depended on it, but they were optional paths and these are not.

`mafft_path` and `hmmer_path` default to `NULL`, meaning "use the conda
environment". Passing a directory still works for a standalone MAFFT, and
the docs warn that the version matters.

**A trap found on the way.** conda and micromamba keep separate roots, so
`tantale` can exist twice with different contents -- it did here, and
`reticulate::conda_list()` returned both. `.tantale_env_prefix()` prefers the
one belonging to the binary in use rather than picking arbitrarily.

**What is left, and why.** `arlem` (144 KB) has no conda package, so the
`R CMD check` executable-files WARNING remains -- but for one small file
rather than roughly 1300. The jars stay: AnnoTALE (16 MB), PrediTALE (14 MB)
and TALEcorrection (27 MB) have no conda packages, and `correct_tales()` is
wanted for the vignettes.

#### Original notes

### 7.4-original Shrink the payload: get MAFFT and HMMER from conda **[A]**

`inst/` is 163 MB, and `inst/tools` is 122 MB of it -- the bulk of what a user
downloads.

| | size | already in the conda env? |
|---|---|---|
| `mafft-linux64` | 34 MB | **yes** -- `mafft=7.520` |
| `talecorrect` | 33 MB | no |
| `hmmer-3.3` | 27 MB | **yes** -- `hmmer=3.3.2` |
| `AnnoTALEcli-1.5.jar` | 16 MB | no |
| `PrediTALE.jar` | 14 MB | no |
| `arlem`, `TALVEZ_3.2`, `QueTAL_v1.1` | < 150 KB each | no |

**[V] The duplication is already paid for.** `inst/tools/tantale_conda_env.yaml`
declares `mafft=7.520` and `hmmer=3.3.2`, so any user who has run
`.create_tantale_env()` has both installed -- and the package still ships its
own copies and calls those instead:

- `msa.R:95` -- `mafft_path = system.file("tools", "mafft-linux64", ...)`, with
  `mustWork = TRUE`
- `tellTale_utilities.R:4` -- `.get_hmmer()` returns
  `system.file("tools", "hmmer-3.3", "bin", ...)`, also `mustWork = TRUE`

So **61 MB of the 122 MB is redundant with a dependency the package already
creates**. Dropping both would roughly halve `inst/tools`.

Both are reached through a single indirection (`mafft_path`, `.get_hmmer()`),
so the change is contained: resolve to the conda env first and fall back to a
bundled copy only if present. The `mustWork = TRUE` has to go either way.

**Caveats worth checking before doing it:**

- MAFFT is used in `--text` mode via `mafft.bat` plus `hex2maffttext` and
  `maffttext2hex` from `mafftdir/libexec`. Confirm the conda build ships those
  helpers and the `.bat` wrapper -- the text mode is unusual and may not be in
  every distribution.
- Versions differ (bundled HMMER 3.3 vs conda 3.3.2; MAFFT unstated vs 7.520).
  Worth confirming output is unchanged before switching.
- This makes the conda env a hard requirement for alignment and TALE mining,
  where today it is only needed for mmseq2 and the Perl tools. That is a real
  change in the install story, not just a size saving.

**The jar files are a separate question.** `AnnoTALEcli` and `PrediTALE` are
30 MB together; bioconda has no package for either as far as I know, so they
would have to stay, be downloaded on demand, or move to a data package.

---

### 7.4a `tantale_setup()` -- DONE **[V]**

A single entry point that checks, and optionally builds, everything the
package needs outside R. Decided after 7.4 made the conda environment a
prerequisite of the core pipeline.

**The reason to build it is correctness, not convenience.**
`.create_tantale_env()` today tests only whether an environment *named*
`tantale` exists. If one does, it prints "can be used for analysis" and
returns success **without looking inside it**. An environment built by an
older version of this package holds MAFFT 7.520, which silently produces
different alignments -- unanchored termini, a different column count (7.4).
Nothing would report it. This happened three times during 7.4 and was caught
only because a golden baseline existed to compare against; a user has no such
thing.

So the pins in `tantale_conda_env.yaml` are currently aspirational. Verifying
them is the point of this function.

**Shape**

```r
tantale_setup(install = FALSE, conda = FALSE, conda_bin = "auto")
```

Diagnostic by default -- called bare it reports and changes nothing:

```
✔ conda binary     /home/cunnac/bin/micromamba       # what reticulate drives
ℹ default root     /home/cunnac/micromamba           # where `-n` would create
✔ tantale env      /home/cunnac/mamba/envs/tantale   # what is actually used
✔ mafft            7.453   (required 7.453)
✔ hmmer            3.3.2   (required 3.3.2)
✔ mmseqs2          14.7e284
✖ java             not on PATH -- needed by AnnoTALE, PrediTALE, TALEcorrection
✔ perl             5.36.0
ℹ Run tantale_setup(install = TRUE) to build or repair the environment.
```

**Requirements**

1. **Check versions against the yaml pins**, not merely presence. Parse the
   pins out of `inst/tools/tantale_conda_env.yaml` so there is one source of
   truth; a hardcoded second list would drift.
2. **Repair, not just create.** `install = TRUE` on an environment with the
   wrong MAFFT must fix it. Note that `micromamba create` on an existing
   environment does *not* downgrade a package -- 7.4 learned this the hard
   way; an explicit `install` of the pinned version does.
3. **Check Java and Perl too.** They are hard requirements of the AnnoTALE,
   PrediTALE and TALEcorrection wrappers, they are not conda's business, and
   they currently fail deep inside a `system()` call with nothing useful said.
   This is the only place they would ever be checked.
4. **Report three separate paths, not one.** The binary, the default root,
   and the environment actually in use are different things, and on a machine
   with any history they diverge. Measured on the development machine:

   | | |
   |---|---|
   | binary `reticulate` drives | `/home/cunnac/bin/micromamba` |
   | micromamba's own root (`MAMBA_ROOT_PREFIX`) | `/home/cunnac/micromamba` |
   | the `tantale` env reticulate resolves to | `/home/cunnac/mamba/envs/tantale` |

   `reticulate::conda_list()` scans several known locations, so it returns
   environments from every root it finds -- two rows named `tantale` here.
   This is not hypothetical: during 7.4 three rebuilds appeared to succeed,
   with honest logs saying MAFFT 7.453 was linked, while the package kept
   using an env in the other root that still held 7.520.

   **Corollary for the implementation: always operate on `-p <prefix>`, never
   `-n <name>`.** `-n` creates under the binary's default root, which is not
   necessarily the root the environment lives in.
5. **Installing conda itself stays opt-in** behind its own argument. Putting a
   package manager on someone's machine is a larger side effect than building
   an environment, and should be asked for. Note that
   `reticulate::install_miniconda()` installs **miniconda**, not mamba -- do
   not let the docs promise otherwise.
6. **The lazy path stays.** Users cannot be made to call this, so
   `.create_tantale_env()` must still work on demand. It gains the version
   check, and its failure message should point at `tantale_setup()`.

Precedent for the idiom: `keras::install_keras()`,
`tensorflow::install_tensorflow()`, `spacyr::spacy_install()`. All of them
require the user to ask before touching the system.

**Knock-on:** this shrinks 7.4b a lot. The README stops needing to explain
conda roots and lazy environment creation, and says instead: install tantale,
run `tantale_setup(install = TRUE)`.

**Built**, in `R/tantale_setup.R`, to the shape specified above. All six
requirements met:

1. *Versions, not presence.* `.tantale_pins()` parses
   `inst/tools/tantale_conda_env.yaml` so there is one source of truth;
   `.tantale_installed()` reads `conda-meta/`, whose filenames are
   `name-version-build.json`. No subprocess, and it does not need conda to be
   working in order to report that conda is not working. Both hard cases are
   covered and tested: hyphenated names (`perl-statistics-r`) and
   non-numeric versions (`14.7e284`).
2. *Repair, not just create.* `.tantale_repair()` runs an explicit
   `conda_install` of the pinned specs, because `create` against an existing
   environment will not downgrade.
3. *Java and Perl checked*, with what needs them named in the failure line.
4. *Three paths reported.* `.tantale_conda_root()` exists because
   `dirname(dirname(bin))` is **not** the root: micromamba's binary sits in
   `~/bin` while its root is `MAMBA_ROOT_PREFIX`. Getting this wrong is
   precisely what made 7.4's rebuilds land in a different root from the one
   in use. Everything operates on `-p <prefix>`.
5. *Conda install is opt-in* behind its own `conda` argument, and the docs
   say `install_miniconda()` installs miniconda, not mamba.
6. *The lazy path still works.* `.tantale_warn_if_unpinned()` rides along
   with `.tantale_env_prefix()`, warning once per session and pointing at
   `tantale_setup(install = TRUE)`.

**Removed while here:** `.create_tantale_env()` announced on every run that
an environment "has been found on your system and can be used for analysis"
-- without having looked inside it. Noise when true and a false assurance
when false, which is the exact failure 7.4a was written to catch. Now silent.

**Not exercised by the test suite:** the `install = TRUE` and `conda = TRUE`
branches, which need network and would modify the machine. The parsing and
comparison they depend on are tested against fixtures; the install call
itself is one `reticulate::conda_install()`.

### 7.4b Tell users how to get conda, and that they now need it -- DONE **[V]**

7.4 changed what a user must have before the package works at all. Before,
MAFFT and HMMER shipped inside it and the core pipeline ran on a bare Linux
box; now the first `tales_align()` or `tell_tales()` builds the `tantale`
conda environment, so **conda or mamba, plus a working network connection, is
a hard prerequisite of the main workflow** rather than of optional extras.

The README and the pkgdown site both need updating, and neither currently
says enough:

- `README.md` line 59 says only "**Conda and Mamba** must be installed as
  well.", in a list of caveats, with no instructions.
- The package-level doc in `R/tantale.R` links to `install_miniconda()`'s help
  page, which is better but still only a link, and it is buried under
  "CAUTIONARY NOTES" alongside remarks about Java and Perl.

**What to write.** The point users need is that they do not have to install
conda by hand or know anything about it -- `reticulate` will do it from
inside R:

```r
install.packages("reticulate")
reticulate::install_miniconda()      # or point at an existing installation
```

and that `reticulate::conda_binary()` is how tantale finds it afterwards, so
an existing conda/mamba/micromamba is used if there is one. Worth stating
explicitly:

- it happens **once**, and the environment is built on first use, not at
  install time, so the first call is slow and needs the network;
- which tools come from it -- MAFFT, HMMER, mmseqs2 and the Perl
  dependencies -- so a failure to build it has an understandable consequence
  rather than an opaque one;
- that conda and micromamba keep **separate roots**, and an environment named
  `tantale` in one is not the one in the other. This bit us during 7.4 and
  will bite a user who has both.

Also worth revisiting while there: the README's last bullet still explains
that Perl libraries are bundled and "cause tantale to occupy quite some disk
space". After 7.4 the size story has changed and that sentence should be
re-checked against what is actually shipped.

**Written.** `README.md` gains a proper three-step Installation section
(install the package, make sure conda is available, run `tantale_setup()`),
replacing the single caveat bullet. The package-level doc in `R/tantale.R`
gains a `@section Setting up:` that says the same thing, and the two remaining
cautionary bullets now point at `tantale_setup()` rather than at a
`reticulate` help page.

Everything 7.4b asked to be stated explicitly is stated: that it happens
once and on first use rather than at install time; which tools come from the
environment; that miniconda is not mamba; that an existing conda/mamba is
found automatically; and that conda and micromamba keep separate roots, with
`tantale_setup()`'s three-path report as the way to see it.

**Two stale things corrected while there.**

- The README claimed bundled Perl libraries "cause tantale to occupy quite
  some disk space". False since 7.4 -- there are no bundled Perl libraries.
  The ~60 MB is three Java programs (AnnoTALE 16 MB, PrediTALE 14 MB, TALE
  correction 27 MB), none of which has a conda package. Corrected.
- `_pkgdown.yml` still listed `annout-class`, whose Rd disappeared when the
  S4 class was retired to `inst/legacy`. That is a **hard error** in
  `pkgdown::build_site()`, so the site could not have been rebuilt. Removed,
  and a "Setting up" section added for `tantale_setup()`.
  `pkgdown::check_pkgdown()` is now clean.

Worth keeping as a habit: **run `pkgdown::check_pkgdown()` after retiring or
adding an exported topic.** Nothing else catches a dangling reference entry,
and `R CMD check` does not look at `_pkgdown.yml`.

### 7.5a An article on the `tales` class, and what a `dom_code` is **[A]**

A full pkgdown article on the `tales` class, written for an audience that is
biologists first. The `tales_msa` class gets the same treatment later, with
the alignment material.

**The section that matters most: what a `dom_code` is.** Nothing currently
explains it to someone who is not already reading the source, and it is the
concept the whole comparison machinery rests on. What it has to say:

- **A `dom_code` names a distinct *domain* sequence -- not a repeat.** This
  is the point the name is making and it must be said first. A TALE part is
  an N-terminus, a repeat, or a C-terminus, and all three get codes on the
  same footing; "domain" is the word chosen precisely to cover them
  indiscriminately. Two parts with the same amino acid sequence get the same
  code, whatever kind of part they are.

  The distinction is not pedantic. On the reference fixture, 251 distinct
  codes cover **180 repeats and 71 termini** -- describing the total as a
  repeat count overstates it by nearly a third. (This exact error was made
  and caught while writing `summary.tales()`.)

- **It is what makes a TALE alignable.** Aligning TALEs residue by residue is
  meaningless -- the repeats are near-identical, so everything matches
  everything. Giving each distinct repeat a symbol turns an array into a
  *sequence of repeat units*, and that can be aligned the way a protein
  sequence is, with insertions and deletions of whole repeats. This is why
  `tales_align()` works on `dom_code` (or `rvd`) rather than on `aa_seq`.

- **It is finer than an RVD, and defined where an RVD is not.** The RVD is
  residues 12-13 and says what base the repeat binds. Two repeats can carry
  the same RVD -- the same specificity -- while differing elsewhere in the
  repeat, and they get different `dom_code`s. So `rvd` is the functional
  layer and `dom_code` the identity layer, which is why the class carries
  both and why the plots let you choose.

  And a terminus has no RVD at all: the `rvd` column holds `NTERM`, `CTERM`
  or `XXXXX` there, which are placeholders standing in for "not a repeat",
  whereas its `dom_code` is a real identifier of a real sequence. Another
  reason the two layers are not interchangeable.

- **How they are computed, plainly**: group the parts by `aa_seq`, number the
  groups. `dplyr::cur_group_id()`, nothing cleverer.

- **And the consequence that bites**: the numbers depend on which arrays were
  in the table when they were assigned. Code 42 from one run is not code 42
  from another. This is not a wart to apologise for but a fact to state
  early, because it is why `tales_compare()` stamps a
  `dom_code_namespace` and why mixing a similarity table from one run with
  codes from another is an error the package tries to catch. The redundancy
  is worth a number: the reference fixture has **251 distinct domains across
  955 parts**, which is also what makes the pairwise comparison affordable --
  it runs over distinct domains, not over parts.

**The rest of the article** should cover: what a `tales` is (one row per
part, not per TALE, and why); the column contract and which columns are
optional; the three ways to build one (`tales_from_telltale()`,
`as_tales()`, `tales_compare()`); that it is a tibble and dplyr verbs work on
it; `tales_anomalies()` and `sanitize`; and the projections
(`tales_rvd_strings()`, `tales_coded_strings()`, `tales_domain_codes()`).

Relates to 8.5b: if `tales_compare()` is broken into three exported steps,
the code-assignment step becomes the natural place to link this explanation
from.

### 7.5 Worked examples: vignettes and `@examples` **[A]**

**Found during the 9.2 sweep:** the "Overview of TALE composition by genome"
chunk in vignette 2 hand-rolls, in about fifteen lines of `ggplot()` calls,
exactly the figure `plot_tales_composition()` now produces in one. It also
passed `color = isNaAaSeq`, a variable defined nowhere in the vignette or the
package -- dead since it was written, and invisible only because these
vignettes do not build (7.1). The stray argument is removed; replacing the
chunk with a `plot_tales_composition()` call belongs to this rewrite.


`plot_tales_msa()` is the most capable function in the package and the hardest
to use: three independent things determine the rendering (cell text, text
colour, block fill), each with its own inputs. Its `@details` now explains the
mechanism, but explanation is not the same as demonstration.

Two gaps, both out of scope for now:

**The pkgdown MSA section is obsolete.** It was written against
`msa_heatmap()`, which is retired, and against the legacy slot names. Vignette
3's calls were migrated mechanically but the surrounding prose still describes
the old workflow, and none of it can be verified while 7.1 stands. It needs
rewriting as a set of commented examples covering the combinations a user
actually reaches for -- each `fill_type`, with and without a tree, with and
without a consensus panel, RVD versus repeat-code labels.

**No exported function has `@examples`.** Nothing in `man/` carries a runnable
example, so `R CMD check` exercises none of the documented API and a reader has
nothing to copy. This matters most for the plotting and class constructors,
where the argument combinations are the hard part.

Note the dependency: useful `@examples` need small, fast, self-contained
fixtures. `inst/extdata` has some, but the plotting examples would want a tiny
alignment that does not require running MAFFT. Worth building that fixture
first; it would serve the vignettes too.

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

## 8.5 Internals audit — function census **[V]**

A census of every top-level definition in `R/` (111: 51 exported, 60
internal), counting call sites within the package.

**Single-caller internals: 27.** But 18 of those already sit in the same file
as their caller, so the maintenance cost is concentrated in the 9 that do
not:

| internal | defined in | only caller |
|---|---|---|
| `.repeat_to_sim_align()` | `conversion.R` | `msa.R` |
| `.repeat_to_cluster_align()` | `conversion.R` | `msa.R` |
| `.rvd_to_match_align()` | `conversion.R` | `msa.R` |
| `.tale_parts()` | `distalr.R` | `tales_class.R` |
| `.build_repeat_msa()` | `msa.R` | `tales_msa_class.R` |
| `.tales_dom_code_namespace()` | `tales_class.R` | `distalr.R` |
| `.tales_msa_contract_holds()` | `tales_msa_class.R` | `tales_class.R` |
| `.run_nhmmer_search()` | `tellTale_utilities.R` | `telltale.R` |
| `.hits_report_to_gff()` | `tellTale_utilities.R` | `telltale.R` |

**The test is not "how many callers".** Some single-caller helpers earn their
name: `.pairwise_align_biostrings/mmseq2/decipher()` are three siblings behind
a `switch` and their symmetry is the readable part; the `.tales_check_*()`
validators surface in error provenance. The useful question is whether the
helper has a name the reader needs. If it does, it should live *next to* its
caller; if it does not, it is a paragraph of the caller that was given a name
for no reason.

**Callers with no caller: 5.** Moved to `R/unused_pending_review.R`, not
deleted -- see the header of that file for what is known about each.

### 8.3 Regression baseline — DONE **[V]**

`tests/testthat/test_golden.R` plus `helper-golden.R`. Twenty snapshots
covering the tales column contract, the anomaly report, the requirements
table, the five projections, `tales_compare()`, both `tales_align()` layers,
`tales_group()`, the data behind `plot()` on a `tales_msa`, and the consensus
of the reference alignment. Runs in about 19 seconds; MAFFT and arlem are
exercised for real.

These assert nothing about what the values *should* be. They record what the
pipeline produces, so a refactor meant to change nothing can be shown to have
changed nothing, and one that does change something says where.

**Whole tables are not snapshotted.** Each is reduced to a fingerprint: one
row per column carrying type, length, distinct count, missingness and an md5
of the values. Snapshots stay readable (23 KB in total) and a diff names the
artefact *and* the column that moved. Doubles are rounded before digesting so
that last-bit differences between machines do not register. Small artefacts
worth reading -- the column contract, the requirements table, the consensus --
are snapshotted whole.

`expect_golden()` forces `cran = TRUE`: `expect_snapshot_value()` skips on CRAN
by default, and a baseline that quietly does not run is worse than none.

Verified to work by reintroducing the `tales_consensus()` tie-break bug: two
snapshots failed, naming the consensus and the plot data. This replaces an
ad-hoc baseline kept in a session scratchpad, which was lost when the session
restarted -- the reason it now lives in the repository.

To accept an intended change: inspect the diff, then
`testthat::snapshot_accept("golden")`.

### 8.0b Quoting paths in shell commands — DONE **[V]**

`tell_tales()` had a nested function definition that shelled out to
AnnoTALE's "analyze" stage. It is now `.run_annotale_analyze()` at top level.

It overlaps the exported `run_annotale_predict()`, which runs predict *and*
analyze starting from a genome. They are not duplicates -- `tell_tales()` has
already found the ORF by the time it calls AnnoTALE, so it needs analyze on
its own -- but the two build their `java -jar` command lines separately, and
they disagree: the exported one wraps paths in `shQuote()` and the internal
one does not. A path with a space in it works through one and not the other.

**Resolved by quoting, not by sharing.** No common helper: the maintainer's
call was to add `shQuote()` where it was missing rather than build an
abstraction over command construction.

The audit found it missing well beyond AnnoTALE. Every shell command the
package builds now quotes its interpolated paths:

| file | command |
|---|---|
| `telltale.R` | `.run_annotale_analyze()`, `.run_nhmmer_search()` |
| `tales_msa_class.R` | the MAFFT pipeline and its `--textmatrix` (2 sites) |
| `distalr.R` | arlem, and the four mmseqs calls |
| `talecorrection_java.R` | nhmmer, and the TALEcorrection jar |
| `target_predictions.R` | PrediTALE, TALVEZ |
| `AnnoTALE_QueTAL_functions_library.R` | functal |

`run_annotale_predict()` and `run_annotale_build()` already quoted theirs.
The parked HMMER wrappers in `unused_pending_review.R` were left alone.

Two of these needed restructuring rather than a wrapped variable, because
they pasted a directory and a filename into one string: the MAFFT binaries
(`{mafft_path}/mafft.bat`) and the TALEcorrection nhmmer outputs
(`{outputFolder}/out_nhmmer.{domains}.txt`). A path is built with
`file.path()` first and quoted whole; quoting only the directory would have
left the separator outside the quotes.

### 8.0 `tell_tales()` has an unguarded filter — FIXED **[V]**

Found while checking that the new baseline is actually sensitive.

`min_domain_hits` filters at `telltale.R:319`:

```r
nhmmerTabularOutput <- subset(nhmmerTabularOutput,
  target_name %in% temp_df[temp_df$V1 > min_domain_hits, "target_name"])
```

If it removes everything, the run does not stop. It carries on and dies
several stages later inside Bioconductor:

```
Error: Rle of type 'NULL' is not supported
  call: new_Rle(values, lengths)
```

A bare `simpleError`, no `tantale` class, nothing naming the argument that
caused it. The two filters immediately before this one -- "no TALE cds hit"
and "no record remains after filtering on score" -- are both guarded and warn
properly, so the pattern to follow is already in the function a few lines up.

Two further things worth noticing about this argument:

- It filters on `target_name`, the **subject sequence**, not on the array. The
  name reads as "minimum hits per TALE array"; it is actually "minimum hits per
  contig". `min_domain_hits = 12` on a fixture with four TALEs changes nothing
  at all, because the two contigs carry far more than twelve hits between them.
  Worth asking whether the documented meaning and the implemented meaning are
  the same thing.
- The comparison is `>` where the name says "minimum", so `min_domain_hits = 4`
  keeps sequences with **five** hits or more.

**Done**, during 5.3. The guard lives in `.telltale_find_domain_hits()`, the
internal that now owns this stage, and says which argument caused it and that
it counts per subject sequence. All three of `tell_tales()`'s give-up points
have tests now (`test_tell_tales_guards.R`); none had any before.

**The other two observations, resolved by the maintainer:**

- **The off-by-one is fixed.** The filter compares with `>=` now, so a
  subject sequence carrying exactly `min_domain_hits` hits is kept, which is
  what the documentation always said. Nothing changes at the default: real
  TALE contigs carry dozens of hits, and only a sequence sitting exactly on
  the threshold behaves differently.
- **Per contig vs per array was not a bug** -- the documentation and the code
  agreed, and my first reading of it here was wrong. The author's own note at
  that line asked whether short *arrays* should also be dropped. They should,
  optionally: `min_array_length` is a second argument, defaulting to `0`,
  which drops arrays with too few repeats after grouping. The two filters are
  complements, not alternatives -- one is a cheap pre-filter on input
  sequences, the other a quality filter on the arrays found in them, and the
  per-contig filter can never catch a 3-repeat fragment sitting on a contig
  that also holds two real TALEs.

  It counts **repeat units**, not all hits, so an array is not penalised for
  having had a terminus missed, and because the repeat count is what "array
  length" means for a TALE -- it determines how long a target box it
  recognises. Default 0 because whether a short array is noise or a truncated
  TALE is a judgement about the biology: a pseudogene with three surviving
  repeats is real, and may be exactly what someone is looking for.

### 8.1 A purpose-built fixture for `tell_tales()` -- DONE **[V]**

`tell_tales()` is slow, and the slowness is not where it looks.

| run on the current 116 kb fixture (2 regions, 4 TALEs) | elapsed |
|---|---|
| `correct_array = FALSE` | 8.8 s |
| `correct_array = TRUE` | **254.6 s** |

The correction is 29x the rest of the pipeline. But it does **not** scale with
the size of the subject sequence. `CorrectFrameshifts()` is called with
`maxComparisons = length(AAref)` against `decipher_ref_tales_aa.fa.gz`, which
holds **1057 reference TALEs**, so every array is aligned against all of them.
The cost is (number of arrays) x (size of the reference set):

| reference TALEs | elapsed |
|---|---|
| 1057 (default) | 254.6 s |
| 100 | 60.9 s |
| 20 | 20.3 s |

So there are two independent levers, and they fix different things:

1. **A trimmed reference set** is what makes the correction path testable at
   all -- 12x, and `correction_ref` is already an argument, so a test can pass
   its own without touching the package default. **DONE**:
   `tests/testthat/data_for_tests/correction_ref_20.fa.gz`, 20 sequences taken
   at even intervals through the shipped file rather than the first 20, so the
   subset is not biased by whatever ordering that file happens to have. The
   correction branch now has a baseline and runs in ~19 s. It pins the code
   path, not the biology: a 20-sequence reference is not claimed to correct as
   well as the full one.
2. **A toy subject sequence** -- two TALEs, one carrying a single-nucleotide
   insertion in its ORF -- shortens the HMMER stage and, more importantly,
   gives the correction a **known right answer** to assert against. Today
   nothing checks that the correction corrects anything; it is only checked
   that the call does not error.

Suggested composition for the toy, to be built from the existing BAI3 regions
rather than synthesised, so the sequences stay biologically real:

- two complete TALE ORFs with short flanks, enough for `extend_len = 300`
  to have something to extend into;
- one left intact, as the negative control -- correction must not change it;
- one with a single nucleotide inserted in a known repeat, far enough from
  the ends that the frameshift truncates the ORF and is detectable. Record
  the insertion point in the fixture's name or a companion file so the test
  can assert the correction restores exactly that position;
- optionally a third region with no TALE at all, which would pin the
  "no hits" branch on real sequence rather than the random DNA the current
  test generates.

Relates to 5.3: the refactor needs a characterisation baseline, and the
baseline needs a fixture that can run in seconds.

**Item 2 built.** `data-raw/make_toy_tale_regions.R` cuts three regions from
the shipped BAI3 sequences into
`tests/testthat/data_for_tests/toy_tal_regions.fasta` (14.6 kb against the
old 116 kb), with a companion `toy_tal_regions_truth.tsv` recording the
answer:

| region | what it is |
|---|---|
| `toy_intact` | one complete TALE, untouched |
| `toy_frameshift` | the same TALE, one `A` inserted mid-array |
| `toy_no_tale` | 4 kb of the same genome with no TALE in it |

**Two copies of the same TALE is the design.** The intact one is the control,
so any difference between the two is attributable to the inserted base rather
than to the arrays being different TALEs. That is what lets the test assert
something real without needing the correction to be perfect in absolute
terms -- which it is not, against the deliberately small 20-sequence test
reference.

Measured, and now asserted in `test_tell_tales_correction.R`:

| | intact | frameshifted |
|---|---|---|
| longest ORF, correction **off** | 4305 nt (93%) | **2631 nt (57%)** |
| longest ORF, correction **on** | 4305 nt | **4305 nt** |
| RVD string, correction on | — | **identical to intact** |
| predicted insertions | 2 | **3** |

So the inserted base truncates the ORF to 57% coverage, and correction
recovers the original TALE exactly: same ORF length, same RVD string, and
precisely one more insertion charged.

**What this replaces.** The correction branch was covered only by a golden
digest. A digest pins the code path but cannot notice correction silently
ceasing to work -- once the changed digest is accepted it simply records the
new wrong answer. These tests fail instead.

The `toy_no_tale` region also pins the no-hits branch on real genomic
sequence, where the previous test used randomly generated DNA, which is a
much easier negative than the real thing.

*Incidental finding:* each toy region also yields a spurious single-hit
array at its 3' end, because `min_domain_hits` filters **subject sequences**,
not arrays -- a contig with enough hits overall keeps all of its arrays,
however small. Pre-existing and not touched here; the tests select the real
arrays explicitly. Worth a look if short spurious arrays ever become a
nuisance.

### 8.1b Curating the shipped correction reference **[A]**

Separate from the test fixture above, and a biological question rather than
an engineering one.

`inst/extdata/decipher_ref_tales_aa.fa.gz` holds **1057 reference TALEs**,
1.25 M amino acids, median length 1198. `CorrectFrameshifts()` is called with
`maxComparisons = length(AAref)`, so every candidate array is aligned against
every one of them. That is the entire reason correction costs 255 s where the
rest of the pipeline costs 9 s, and the cost is borne by every user on every
run.

Questions worth answering before touching it:

- **Where did the 1057 come from?** **Answered by the maintainer:** it is the
  raw output of `tell_tales()` run over a large set of *Xanthomonas oryzae*
  genomes, assembled long ago, **with no curation applied**. The file's own
  contents corroborate that exactly -- 70 genome accessions, AnnoTALE-style
  `<accession>-tempTALE<n>` names, `(Pseudo)` markers left in place, and
  fragments far too short to be TALEs.

  So the question is not whether to *re*-curate but whether to curate at all,
  for the first time.
- **How redundant is it?** TALEs are highly similar by construction. If the
  set collapses to a few dozen clusters at high identity, most of those 1057
  alignments are re-deriving the same answer.
- **Does a smaller set correct as well?** This is the measurable one: correct
  a set of known-frameshifted arrays against the full reference and against
  candidate subsets, and compare the corrected sequences. If a curated 100
  reproduces the full set's output, the default could change and every user's
  run gets ~10x faster.
- **Is `maxComparisons = length(AAref)` the right call at all?** DECIPHER's
  own default is lower. Capping it is a one-line change that does not require
  touching the reference file, though it makes which references get compared
  depend on ordering.

Do not trim the shipped file on speed grounds alone: a reference that is fast
but corrects worse is a bad trade, and correction rewrites the user's
sequences.

**Resolved.** Answers to the four questions above, and what was done.

*Provenance.* Confirmed by the maintainer and by the file: uncurated
`tell_tales()` output over 70 *X. oryzae* genomes. 1057 sequences, of which
555 are byte-identical duplicates and 8 are too short to be TALEs (shortest
23 aa against a median of 1198).

*Redundancy.* Very high, as expected. After dedup and a 300 aa floor, 494
remain; clustering the intact ones at 0.04 collapses them to 65 clusters.

*Does a smaller set correct as well?* **Reference size is not the lever.** 494
and 1057 give the same corrections at almost the same cost, because the cost
is dominated by alignment, not by the cheap pre-screen every reference goes
through.

*Is `maxComparisons = length(AAref)` right?* **This is the lever**, and the
finding that mattered. Reading the DECIPHER source settled what the docs do
not say: `CorrectFrameshifts()` scores *all* references with a cheap distance,
`order(d, widths, decreasing = TRUE)`, truncates to `maxComparisons`, and only
then aligns, stopping early at `acceptDistance`. So the cap does not pick
arbitrary references -- it bounds how deep a ranked search goes. Measured on
four arrays against the 1057 set, all byte-identical: 252 s uncapped, 10 s at
20.

(Method note, worth keeping: I first asserted this ranking behaviour, then
wrongly retracted it when the documentation was silent. The maintainer's
correction -- *"reading the source code of the function may teach you far more
that you ask for with the tests"* -- was right, and reading it confirmed the
original claim. Prefer the source to inference when the docs are silent.)

**What shipped.**

- `data-raw/make_correction_references.R` builds both sets from
  `data-raw/tale_correction_ref_source.fa.gz` (the renamed original).
- `inst/extdata/tale_correction_ref.fa.gz` -- 494, the new default.
- `inst/extdata/tale_correction_ref_representative.fa.gz` -- 136, a
  diversity-sampled subset (all pseudogenes + one longest member per cluster).
- Pseudogenes are kept in both, deliberately: the genomes were high quality,
  so the frameshifts are real biology, and a reference of only intact TALEs
  risks "repairing" a genuine pseudogene into an ORF no strain carries.
- `Clusterize()` settings: `includeTerminalGaps = TRUE`,
  `penalizeGapLetterMatches = NA`, `method = "overlap"` (inert when
  `includeTerminalGaps` is TRUE, stated for clarity).
- `max_comparisons` promoted to a real `tell_tales()` argument, defaulting to
  `NULL` = all references. It had to be a real argument, not passed through
  `...`, which collided ("formal argument matched by multiple actual
  arguments").
- It is echoed into `tell_tales.log`, since it changes results.

**The trade-off, measured.** The cap is not a free speedup, and it fails in
the worse direction -- not by leaving an array uncorrected but by correcting
it against a poor reference, which still looks like a corrected ORF. Against
a deliberately small 20-sequence reference:

| `max_comparisons` | indels called per array |
|---|---|
| all (20), 20, 10 | 2, 2, 0, 1 |
| 5 | 2, 2, 0, **2** |
| 2 | **9, 11, 0, 15** |

So what matters is not the ratio to the reference set but whether the closest
`max_comparisons` are genuinely close: 20 of 1057 is ample, 5 of 20 is not.
This is in `@param max_comparisons` and pinned by a test in
`test_tell_tales_guards.R`.

*`processors`.* Tried, per the maintainer's recollection that it once broke
things. It no longer breaks, and gives ~10% -- it does not scale. Not worth
exposing; `max_comparisons` is the real lever.

**Deferred to the maintainer:** validating the 136-sequence set against the
494 on real frameshifted arrays. Their call -- *"I will do the test myself
later on"*.

### 8.1d Golden baseline records machine-specific paths -- DONE **[V]**

`tell_tales.log` echoes absolute paths -- the three HMM files, and
`correction_ref`. Under `load_all()` these are the source tree; from an
installed package they are the library path. `.RUN_SPECIFIC` in
`helper-golden.R` drops `/tmp/` and `Rtmp` lines but not these, so the golden
snapshot for `tell_tales.log` only reproduces on the machine that recorded it.

Pre-existing, not introduced by 8.1b -- but it surfaced there, because
renaming the reference file changed that line and nothing else.

A baseline that fails for everyone but its author is worth much less than one
that travels. Options: log basenames instead of paths (loses provenance),
filter those lines (loses the ability to catch a default change -- which is
exactly what caught the rename), or teach the fingerprint to rewrite absolute
paths to a placeholder rather than drop whole lines. The third keeps both
properties and is probably right.

**Fixed**, by the third option: `helper-golden.R` now rewrites absolute
directories to `<path>/` before digesting, keeping the basename. Both
properties are preserved -- the baseline travels between machines, and a
change of *which* reference file is used still shows up.

Also added `^# Current dir:` to the drop list; that is HMMER echoing the
working directory, which is where the run happened rather than what it found.

**The interesting part was getting the pattern narrow enough.** A first
attempt matched "anything between two slashes", which also rewrote
`</title></head>`, HMMER's `//` record separators and the `//` in
`http://hmmer.org/`. It would have produced a perfectly stable digest while
quietly destroying content -- a weaker baseline wearing the appearance of a
more portable one. The pattern now requires a slash that starts a token and
at least one non-empty `segment/` group, and every case the broad version got
wrong is a test.

Audited after the change: exactly one file (`tell_tales.log`) and exactly the
four lines this section identified are touched.

### 8.1c `...` must not hide arguments behind an internal **[V]** — audited

**The rule.** When an exported function forwards `...` to something the user
cannot see, the `@param ...` has to name the arguments themselves, not point
at the callee. "Passed to `.build_repeat_msa()`" is useless advice: the reader
cannot call that function, cannot read its help, and has no way to discover
what it accepts.

**The case that prompted it.** `tales_align()` took `x`, `residue_col`,
`repeat_sims` and `...`, and its `@param ...` read "Passed to
`tales_align()` (e.g. `mafft_opts`)" -- circular, and naming an argument
without saying what it does or what its default is. MAFFT's options were
therefore reachable but undiscoverable, which matters because the default
sets gap penalties (`--op 0 --ep 5`) that are unusual on purpose and that a
user may well want to change.

Fixed by promoting `mafft_opts` and `mafft_path` to real arguments of
`tales_align()`, documented in terms of what they do to an alignment of
TALE repeats rather than as a pass-through.

**The audit.** Every exported function that forwards `...`:

| function | `...` reaches | verdict |
|---|---|---|
| `tales_align()` | `.build_repeat_msa()` (internal) | **was the problem; fixed** |
| `as_tales()`, `as_tales.data.frame()` | `tales()` | fine -- exported and documented |
| `tales_predict_targets()` | `talvez()`, `preditale()` | fine -- both exported, and `@param ...` links to them |
| `tell_tales()` | `DECIPHER::CorrectFrameshifts()` | fine -- names the external function and links to its help |

So this was one occurrence, not a pattern. The rule stands for anything added
later: **if `...` lands somewhere the reader cannot open, the arguments belong
in the signature or spelled out in the docs.** Worth re-running the audit
(`scratchpad/dots.R` in the session notes, trivially rebuilt) whenever a new
exported wrapper appears.

### 8.5b Break `tales_compare()` into three composable steps **[A]**

Maintainer's proposal, and I agree with it. `.tales_compare_core()` does three
things that are separable and each independently useful:

1. assign the `dom_code`s;
2. compute the `domain_distances`;
3. compute the `tale_distances`.

**One correction to the ordering, which improves the design rather than
complicating it.** Steps 2 and 3 are not parallel: **the TALE distances are
built *from* the domain distances.** At `distalr.R:498-510` the pairwise
repeat dissimilarities are cast to a matrix, passed through
`stats::dist(method = "minkowski", p = 3.5)` to force the triangle
inequality, rescaled to 0-100, and written as ARLEM's cost matrix. ARLEM then
aligns the repeat-code strings *using that matrix* as its substitution cost.

So the real chain is **1 -> 2 -> 3**, and that is worth exposing rather than
hiding, because it states something biological the current monolith conceals:
two TALEs are compared by aligning their repeat arrays, where the cost of
substituting one repeat for another is how different those repeats are as
proteins. The repeat-level comparison is not a by-product of the TALE-level
one; it is its input.

Composed, the three would read:

```r
x  <- tales_assign_domain_codes(x)          # 1
dd <- tales_domain_distances(x, aln_method) # 2
td <- tales_tale_distances(x, dd)           # 3, consumes dd
```

and `tales_compare()` stays as the convenience wrapper that runs all three.

**What already exists, and what does not.** `tales_domain_codes()` is taken
but does something else -- it *reads back* the `dom_code`/`aa_seq`
correspondence from an object that already has codes. Step 1 needs a
different name.

**The hazard to think about before exporting step 1.** Codes are assigned
with `dplyr::cur_group_id()` over `aa_seq`, so they depend on which arrays
were in the table at the time. That is exactly why the `dom_code_namespace`
stamp exists -- to catch similarity tables from one run being used with codes
from another. Exporting the assignment makes that run-dependence part of the
public API, so the function must stamp a namespace and its documentation must
be blunt: **these codes are meaningful only within one call, and comparing
them across calls is an error the namespace is there to catch.**

Worth doing. Every step is separately useful -- someone may want the
repeat-level distances without paying for ARLEM at all -- and the
decomposition documents the model.

### 8.2b `tales_coded_strings()` needs a `sep` argument -- DONE **[V]**

The two projections are siblings and should take the same arguments, but do
not:

```r
tales_rvd_strings(x, sep = "-", rvd_only = TRUE)
tales_coded_strings(x)
```

`tales_coded_strings()` hardcodes `collapse = " "`
(`tales_projections.R:25`). Add `sep = "-"`... but **check the default before
changing it**, because the separator is not cosmetic here:

- `.build_repeat_msa()` is handed repeat-code strings built with `sep = " "`,
  and splits them back on the same character. `tales_align()` constructs its
  own strings rather than calling `tales_coded_strings()`, so it is probably
  insulated -- confirm that before assuming it.
- A repeat code is a bare integer rendered as text, so `"1 2 3"` and
  `"1-2-3"` are both unambiguous; unlike RVDs, no code contains either
  character. The choice is therefore free, which is exactly why it should be
  the caller's.

Whether the default should match `tales_rvd_strings()`'s `"-"` (consistency
between siblings) or stay `" "` (not changing output for existing callers) is
a judgement about which matters more. Matching the sibling reads better but
changes what the function returns today.

Same pass should check `rvd_only`: `tales_rvd_strings()` has it and
`tales_coded_strings()` does not, and it is not obvious whether dropping the
terminus codes makes sense for repeat codes.

**Resolved.** `tales_coded_strings(x, sep = " ", repeats_only = FALSE)`.

*Was the separator safe to change?* Checked, and yes -- but it was not
changed. `tales_coded_strings()` has **no callers inside `R/` at all**, only
tests; `tales_align()` builds its own `" "`-joined strings at
`tales_msa_class.R:288` rather than going through the projection, so it was
never coupled. The default stays `" "` because the function's documented job
is the encoding ARLEM and MAFFT `--text` consume, and those split on spaces.
Changing it would have silently altered output for existing users with no
error anywhere.

So the siblings now take the same two arguments with *different defaults*,
which is the honest answer rather than a compromise -- they feed different
consumers, and the docs say so in a small table.

| | `tales_rvd_strings()` | `tales_coded_strings()` |
|---|---|---|
| `sep` | `"-"` (AnnoTALE) | `" "` (ARLEM, MAFFT) |
| filter | `rvd_only = TRUE` | `repeats_only = FALSE` |

*And `rvd_only`?* Added as `repeats_only`, defaulting to `FALSE`. Target
prediction concerns repeats only, so dropping termini is right for RVDs;
alignment is the consumer here and the termini are the most reliable anchors
an alignment of TALE arrays has. It filters on `domain_type == "repeat"` and
errors if that column is absent.

**Naming debt noted:** `tales_rvd_strings()`'s `rvd_only` means "repeats
only" and would be better named `repeats_only` to match. Not renamed -- it is
an exported argument and §9.1 is closed. Worth folding into any future
breaking pass.

**A bug this uncovered.** The shared `coded_tales()` test fixture carried
`rvd = NTERM, NI, NTERM` for one array -- two N-termini, the second
mid-array. Not a TALE. It had gone unnoticed because the fixture had no
`domain_type` column, and without one the anomaly checks cannot run. Adding
the column made `tales()` flag it at once (`terminus_duplicated`,
`terminus_misplaced`). Fixture rebuilt as one complete and one incomplete
array, preserving every property the tests relied on (three distinct codes,
code 7 -> `MDP`, a recurring code within an array).

Worth generalising: **a fixture that omits the columns the validators key on
is not exercising the validators.** Other minimal fixtures in the suite are
likely in the same position.

### 8.2 Silence MAFFT by default — DONE **[V]**

`.build_repeat_msa()` runs MAFFT through `system()` with
`ignore.stderr = FALSE` (`tales_msa_class.R`, the `res <- system(...)` call).
MAFFT writes its banner, strategy notice and per-sequence progress to stderr,
so every alignment floods the console with dozens of lines the user did not
ask for. The alignment itself is already redirected to a file with `>`, so
stdout carries nothing of interest either.

Wanted: a `mafft_verbose = FALSE` argument on `.build_repeat_msa()`, surfaced
through `tales_align()`. Note the spelling — the package converted every
argument to snake_case in 9.1, so `mafft_verbose`, not `mafftVerbose`.

**One thing to get right.** Simply setting `ignore.stderr = TRUE` also
discards MAFFT's error messages, and the current failure path is already thin:
when the output file comes back empty the code aborts with nothing but the
exit status, so a silenced run would report *that* it failed and never *why*.
Better to redirect stderr to a temporary file (`2> {logfile}` in the command,
or `stderr = TRUE` on a captured call) and replay its contents only when the
run fails. That gives silence in the normal case and more diagnostics than
today in the failing one.

**Done.** `mafft_verbose = FALSE` on `.build_repeat_msa()`, surfaced as a
real argument of `tales_align()`. Measured: **61 lines of stderr per
alignment, down to 0.**

`--quiet` was not used. Redirecting stderr to a temporary file does more: it
covers the banner and the strategy notice as well as the progress, and the
captured text is replayed when the run fails. Silence therefore costs nothing
diagnostically -- the failure path is *better* than before, which reported an
exit status and nothing else. Verified against a deliberately bad option:
the abort carries MAFFT's own usage output and is classed
`tantale_error_mafft_failed`.

### 8.4 `print()` methods for `tales` and `tales_msa` — DONE **[V]**

Both classes currently fall through to the tibble print method, so the screen
says `# A tibble: 955 x 10` and nothing about what the object *is*. Everything
the class knows that a tibble does not is invisible:

- that it is a `tales` at all, rather than a data frame that happens to have
  these columns;
- how many arrays it holds, as opposed to how many parts (rows);
- which residue layers are present (`rvd`, `dom_code`), which is what decides
  what `tales_align()` and `plot()` can do with it;
- the `dom_code` namespace stamp, whose whole purpose is to catch tables from
  different runs being mixed, and which is invisible until something fails;
- for a `tales_msa`, the alignment width, and how gappy it is.

A print method is the cheapest place to surface all of that, and it is the
first thing a user sees. Worth doing before the vignettes are rebuilt, since
the printed object will appear throughout them.

**Done.** `R/tales_print.R`, 25 tests.

```
<tales> 44 arrays, 955 parts
  layers: rvd, dom_code   |   namespace: 4a3059c6   |   6 other columns
                      dom_code
  BAI3_ROI_00001      194 68 152 68 94 154 151 157 152 94 60 34 153 154 15 ...
  BAI3_ROI_00002      199 64 50 64 94 149 8 54 5 94 127 115 158 39 56 50 5 ...
  ...                 ...
  PXO86_ROI_00018     213 131 120 128 140 79 145 84 10 40 96 116 18 166 96 ...
  PXO86_ROI_00019     214 139 116 116 114 50 21 50 39 17 18 65 93 115 94 1 ...
```

A `tales_msa` prints the same way but draws its gaps and pads every cell to a
common width, so the columns line up down the page -- an alignment is
something you recognise by looking at it.

The preview follows Biostrings: first two, `...`, last two, eliding nothing
when the object holds four arrays or fewer. It previews `dom_code` over
`rvd` because repeat codes discriminate better -- two arrays can share an RVD
sequence while being built from different repeats.

**Two bugs caught while writing the tests**, both worth remembering:

- The header was built with `cli::cli_text()`, which writes to the **message
  connection**. A `print()` method must write to stdout: as it was, the header
  would interleave wrongly under redirection and `capture.output()` could not
  see it at all. `cli::format_inline()` plus `cat()` keeps the styling and the
  pluralisation while going to the right place.
- I briefly "fixed" a non-existent portability problem with `%||%`, thinking
  it was base-only since R 4.4 while `DESCRIPTION` allows 3.6.3. The package
  defines its own at `tales_msa_class.R:206`, so it was never at risk.

**`format()` methods too**, for both classes, and not as polish: both
inherit a `format()` from tibble, so adding `print()` alone left the two
halves of one operation disagreeing. `print(x)` showed the view above while
`format(x)` still returned `# A tibble: 955 x 10`. Anyone writing
`cat(format(x), sep = "\n")` -- the idiomatic way to get a printed form as
text -- got the wrong one.

`format()` now builds the lines and `print()` only emits them. The tibble
rendering stays reachable as `format(tibble::as_tibble(x))`.

A third bug fell out of that: the `tales_msa` header used a `\\` line
continuation inside the `cli` template, which left an **embedded newline in
one element**. `format()` reported six lines while printing seven, so the
vector lied about its own length -- exactly what breaks a caller indexing
lines to put in a log. Tests now assert that no element contains a newline
and that `format(x)` equals `capture.output(print(x))` for both classes.

**`summary()` methods too**, `R/tales_summary.R`, 26 tests. They return an
object that a `print` method renders, so the numbers are usable and not
merely visible.

```
<tales> summary
  arrays / parts            44 / 955
  distinct domains          251 of 955 parts  (C-terminus 37, N-terminus 34, repeat 180)
  distinct RVDs             17
  repeats per array         min 12   median 19   max 27
  arrays with both termini  44 of 44
  source sequences          4
  anomalies                 none
```

Contents were cross-checked against what `tell_tales()` already thinks worth
logging, which independently confirmed "complete arrays" and array-length
min/median/max. It also turned up a commented-out line in that log --
"Total number of distinct types of RVD" -- which the author wanted and lost
when the table it needed was disabled. It is back, here.

For a `tales_msa` the measure that earns its place is **columns with no
consensus, per layer**. On a gappy four-array alignment: 17 of 28 columns
have no `dom_code` consensus but only 10 have no `rvd` one. That gap is the
biology -- repeats that are distinct proteins can share a base preference --
and it is a one-line answer to "is this alignment telling me something, or is
it disagreement all the way down".

**Left out deliberately:** the RVD frequency table (composition analysis,
belongs in its own function returning data), and per-array breakdowns (they
scale with the object; a summary should not).

Still open: nothing. `format()`, `print()` and `summary()` are all in place
for both classes.

### 8.6b Co-locating single-caller internals — PARTLY DONE **[V]**

Done, and the file responsibilities now line up:

- `msa.R` is **drawing an alignment**: `plot.tales_msa()`, the three
  fill-layer builders moved in from `conversion.R`, the consensus functions
  and the reference picker.
- `tales_msa_class.R` is **the class and how to build one**: constructor,
  validators, `as.matrix()`, `tales_align()`, and the MAFFT runner moved in
  from `msa.R`.
- `conversion.R` is down to five functions, about projecting a `tales` onto
  strings and maps, which is what its name suggests.

The remaining six need a decision about which file owns what, so they are
left alone. My reading of each:

| internal | situation | suggestion |
|---|---|---|
| `.tale_parts()` + its two helpers | **DONE** — the trio and `tales_from_telltale()` are now `R/tales_ingest.R`. `distalr.R` no longer reads anything off disk | — |
| `.build_repeat_msa()` | **DONE** — moved to `tales_msa_class.R` beside `tales_align()`, taking `.as_mafft_score_table()` and `.rvd_score_table()` with it, since it is their only caller. Co-located, not inlined: 157 lines of MAFFT plumbing inside a 57-line `tales_align()` would have made the caller harder to read, and "run MAFFT in text mode and return a matrix" is a name the reader needs. Six tests also call it directly with fasta paths and bare sequence lists, exercising edge cases `tales_align()` cannot reach | — |
| `.tales_dom_code_namespace()` | in `tales_class.R`, called from `distalr.R` | leave. It is a property of the class read by another module, which is normal |
| `.tales_msa_contract_holds()` | in `tales_msa_class.R`, called from `tales_class.R` | leave, same reason |
| `.run_nhmmer_search()`, `.hits_report_to_gff()` | in `tellTale_utilities.R`, called only from `telltale.R` | `tellTale_utilities.R` is down to five members after the parking. Either fold what is left into `telltale.R` and drop the file, or leave it |

### 8.6 Legacy preconditions leaking through class methods — DONE **[V]**

`plot.tales_msa()` decomposes its object and hands the pieces to
`plot_tales_msa()`, whose argument checks are written for a caller assembling
matrices by hand. Measured against what the class guarantees:

| check in `plot_tales_msa()` | reachable via the method? |
|---|---|
| neither `repeat_align` nor `rvd_align` given | no -- the method always builds `repeat_align` |
| `repeat_align` was coerced to a vector | no -- `as.matrix.tales_msa()` uses `matrix()`; a one-array subset still returns a `1 x n` matrix |
| `rvd_align` was coerced to a vector | no -- same |
| fewer than one sequence | **yes** -- a zero-row `tales_msa` is valid |

So three of four are unstateable, and the one that fires reports a problem
with `repeat_align`, an argument a `plot(x)` caller never supplied and cannot
inspect.

Decided: `plot_tales_msa()` is folded into `plot.tales_msa()` and unexported.
Done.

**The other entry points, audited the same way:**

- `plot.tales()` -> `plot_tales_composition()`: **clean.** A one-line
  pass-through to a function that already takes the object and checks through
  `.tales_require()`; every message names a column, not an argument the caller
  did not supply. The only question here is a different one -- two public
  names (`plot(x)` and `plot_tales_composition(x)`) for one operation, both
  taking the same object. Unlike the msa case there is no legacy matrix
  interface to remove, so this is an API-surface choice, not a defect.

- `tales_compare()` -> `.tales_compare_core()`: **one dead check.** The core
  aborts with "Your tale arrays identifers are probably not unique. Make sure
  that there is only one part per position per array_id." Tested: `tales()`
  already rejects a duplicated `array_id`/`position_in_array` pair, so this
  cannot fire through `tales_compare()`. (It also misspells "identifiers".)
  **Deleted.** The class is where that invariant belongs, and duplicating it
  in a private function only created a second place for it to go stale.

  Its two *other* checks are live and must stay: `tales()` accepts `NA` and
  `""` in `aa_seq` (tested), so "Some of the provided TALE parts have no amino
  acid sequence" is reachable and doing real work.

- `tales_align()` -> `.build_repeat_msa()`: **fixed.** Its messages named
  `input_seqs`, an internal argument a `tales_align()` caller has no way to
  inspect.

**Incidental finding, resolved:** `tales()` accepts `NA` in a residue column,
but not silently -- `.tales_anomalies()` reports it as `missing_rvd` /
`missing_dom_code`, `tales()` warns, and `sanitize = TRUE` drops the array.
Working as designed; empty strings are covered by the same check.

### 8.7 `tales_consensus()` depended on row order **[V]** — FIXED

Found while checking that folding `plot_tales_msa()` into `plot.tales_msa()`
preserved behaviour: the plot data differed, and the difference was real.

`tales_consensus()` scored candidates with
`unique(allElements)[which.max(freq)]`. `unique()` returns values in order of
first appearance and `which.max()` takes the first maximum, so a tie was won
by whichever array happened to be the top row. Permuting the rows of an
alignment changed its consensus. Fixed by sorting the candidates first; the
counting is untouched.

**The deeper question is left open.** On the three-array fixture, two of 28
columns carry three *distinct* repeats — no majority exists at all, yet a
value is still reported and the figure then colours cells by whether they
"match the consensus" at a position that has none. A deterministic arbitrary
pick is better than a non-deterministic one, but it is still arbitrary.
Options, for the maintainer:

- return `NA` where no element is strictly more frequent than the rest, so
  the figure shows no consensus rather than a fictitious one;
- keep a value but mark weak columns (the usual sequence-logo convention);
- leave as is, treating it as "modal element" rather than "consensus", and
  say so in the name.

This is a biological question about what the figure should claim, not a
coding one, which is why it is parked here.

## 9. Long-term systematic passes **[A]**

Whole-codebase sweeps, to be done deliberately rather than opportunistically.
Deferred until the class design settles, since it will dictate several of the
names.

### 9.0 Governing convention — rOpenSci package API guidelines **[reference]**

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

### 9.2 Column names — DONE for the classes, PARTLY for the report files **[V]**

The `tales` class and the distance tables use snake_case throughout, and the
two egress bridges that let legacy-named internals survive are deleted.

**Correction (found later):** the sentence that stood here said "every table
the package produces". That overstated it -- in the three TSVs `tell_tales()`
writes, only `array_id` was renamed. See §9.2b for the inventory.

**What the vocabulary is now**

| table | columns |
|---|---|
| `tales` | `array_id`, `domain_type`, `position_in_crd`, `dna_seq`, `source_directory`, `position_in_array`, `aa_seq`, `rvd`, `seqnames`, `dom_code` |
| `tales_domain_codes()` | `dom_code`, `aa_seq`, `rvd` |
| `domain_distances` / `tale_distances` | `id1`, `id2`, `dissim` (+ `arlem_score`, `max_length`) |
| `tales_group()` | `name`, `group` |

The two distance tables now agree on their id columns, which was the stated
prerequisite for unifying them into one class -- that unification landed in
9.6, and this sweep removes the last places that still spoke the old
vocabulary behind it.

**Bridges: one kept, two deleted**

- `.tales_rename_legacy()` — **kept**. Ingest only. A `tell_tales` output
  directory written by an older version still has camelCase headers, so
  `tales()` must keep accepting them. `.as_mafft_score_table()` and
  `.pairwise_distances_rename_legacy()` are the same courtesy for the
  distance tables.
- `.tales_to_legacy()` — **deleted**. Existed only because
  `.tales_compare_core()` was written against camelCase.
- `.distances_to_legacy()` — **deleted**. Existed only because
  `plot_tales_msa()` was written against `TAL1`/`RepU1`/`Sim`. That function
  now normalises both of its similarity arguments through
  `pairwise_distances()` at entry, so it accepts either spelling and its
  internals speak one.

**On-disk formats changed too.** `arrayReport.tsv`, `domainsReport.tsv` and
`hitsReport.tsv` now write `array_id`, and the derived GFFs carry an
`array_id` attribute. Authorised explicitly -- this release breaks things by
design. The fixtures under `tests/testthat/data_for_tests/` and
`inst/extdata/` were rewritten to match.

**Defects this uncovered**

1. `.tale_parts_from_file()` named the column `arrayIDs` in its empty-input
   branch and `arrayID` in the populated one, so the two returns had
   incompatible schemas.
2. `repeat_to_rvd_map_distalr()` and `tale_parts_to_rvd()` are exported and
   documented as taking a `tales_compare()` result, but read camelCase -- so
   both had been broken against that result since the class work landed.
   Neither had a test that would notice.
3. Two test assertions went silently vacuous when the fixture moved:
   `d$tale_parts$arrayID` returns `NULL`, and `expect_identical(NULL, NULL)`
   passes. Both now index with `[[ ]]`, which errors on a missing column.
   **This is the third distinct way `$` has hidden a defect in this project;
   prefer `[[ ]]` in tests.**
4. The tree panel of `plot_tales_msa()` had no test at all, so the dendrogram
   rewrite was flying blind until one was added.

**Not renamed, deliberately:** `repeatClusterId`, `repeatSimVsRef`,
`rvdSimVsRef`, `matchConsensusRepeat`, `matchConsensusRvd`. These are columns
of the intermediate plot-data tibble inside `plot_tales_msa()`, not of any
table the package returns. They are reachable as `p$data`, so they are worth
a later pass, but renaming them changes nothing a documented API promises.

**Follow-up needed from the maintainer:** `man/figures/pipeline.svg` (and the
PNG exported from it) label the `tale_parts` box with the old column names —
`arrayID`, `domainType`, `positionInCrd`, `dnaSeq`, `aaSeq`,
`positionInArray`, `domCode`, `sourceDirectory`. The figure needs re-exporting
once the labels are updated. Not touched here.

#### Original notes

### 9.2-original Column names — adopt snake_case across all tables

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

### 9.2b The sweep stopped at `array_id` in the report files **[A]** — needs your call

§9.2 above says "Every table the package produces now uses snake_case". That
is true of the `tales` class and the distance tables, and **not** true of the
three TSVs `tell_tales()` writes. What 9.2 actually changed there was
`array_id`; the rest of the columns were left alone.

Measured on a corrected run:

| file | snake_case | still legacy |
|---|---|---|
| `domainsReport.tsv` | all 4 | -- |
| `hitsReport.tsv` | 10 of 12 | `nhmmerHitID`, `hitID` |
| `arrayReport.tsv` | 3 of 17 | the other 14 |

`arrayReport.tsv` in full: `OriginalSubjectName`, `Start`, `End`, `Strand`,
`NumberOfHits`, `ArraySeq`, `AllDomains`, `SeqOfRVD`, `aberrantRepeat`,
`N.terminusAAlength`, `C.terminusAAlength`, `LongestOrfLength`,
`OrfCovOverArrayLength`, `LongestORFSeq`. (The three that are already right
are `array_id` and the two `predicted_*_count` columns, which correction
adds.)

Note `N.terminusAAlength` is not merely camelCase -- the dots are what
`data.frame()` does to `N-terminus`, so that name is an accident rather than
a choice.

**Proposed mapping**, if you want it finished:

| now | proposed |
|---|---|
| `OriginalSubjectName` | `seqnames` (matches every other table) |
| `Start`, `End`, `Strand` | `start`, `end`, `strand` |
| `NumberOfHits` | `n_domain_hits` |
| `ArraySeq` | `array_seq` |
| `AllDomains` | `has_all_domains` |
| `SeqOfRVD` | `rvd_string` |
| `aberrantRepeat` | `has_aberrant_repeat` |
| `N.terminusAAlength` / `C.terminusAAlength` | `nterm_aa_length` / `cterm_aa_length` |
| `LongestOrfLength` | `longest_orf_length` |
| `OrfCovOverArrayLength` | `orf_coverage` |
| `LongestORFSeq` | `longest_orf_seq` |
| `nhmmerHitID`, `hitID` | `nhmmer_hit_id`, `hit_id` |

**Why this is not done unattended.** These are the column names of the
package's primary output files, and the person who knows what reads them
downstream is the maintainer, not me. §9.2 broke this format once already
"by design", so doing it again is defensible -- but it is a decision, not a
chore. The work itself is mechanical, has golden coverage, and would take
one pass.

Two things to decide: whether to do it at all, and whether
`.tales_rename_legacy()` should learn the old spellings so that directories
written by the current version still load after the change.

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

### 9.7 Roxygen markdown enabled -- DONE **[V]**

`DESCRIPTION` had no `Roxygen:` field, so markdown was off -- while the
roxygen comments had been written for years as if it were on. Backticks,
`*emph*` and `**strong**` were reaching the rendered help as literal
characters: `**precondition**` showed up in `?tales_assert_complete` with the
asterisks visible.

Turned on with `Roxygen: list(markdown = TRUE)`.

**How it was verified**, since this reparses every block in the package and
can silently change any of them: copy `man/` aside, flip the flag,
re-document, and diff. 31 of 64 Rd files changed. Classified:

- Most of the diff is whitespace re-wrapping -- no content change.
- The rest is the fix: `` `dom_code` `` and friends became `\code{}` (7
  occurrences of `dom_code` alone), `*valid*`/`*unreadable*` became
  `\emph{}`, `**not**`/`**precondition**` became `\strong{}`.
- **Two regressions, caught and fixed.** Square brackets in prose are link
  syntax under markdown, so `[eg PacBio, ONT]` and `[see the min_gap
  parameter]` in `tell_tales`'s description became `\link{}` to targets that
  do not exist -- an R CMD check WARNING. Rewritten as plain prose.

Checks that it is clean: the set of `\link` targets is now identical to
before the switch, all 64 Rd files pass `tools::parse_Rd()`, and no literal
`**` remains anywhere in `man/`.

Snake_case survived unharmed: CommonMark does not treat intraword
underscores as emphasis, so `dom_code` and `array_id` render as written.

**The rule for anything written later:** square brackets in roxygen prose are
a link. Use parentheses, or escape them.

## 10. Explicitly ruled out

- Deleting dormant internals such as the unused HMMER wrappers
  (`.write_hmm_file()`, `.run_hmmer_search()`, `.run_hmmalign()`,
  `.extract_seqs_from_hits()`). Obsolete now, plausibly useful later.
- Treating `run_annotale_predict()` / `run_annotale_build()` as dead. They have
  no internal call sites but are legitimate standalone user utilities — the type
  case for "useful outside the `pipeline.svg` workflow".
