# tantale — class system design

Detailed design for the S3 class overhaul. Companion to
`restructuring-notes.md` §5 (OOP restructuring), which stays the action ledger;
this document holds the class definitions themselves — identity, invariants,
constructors and method policy.

Branch: `dev`. Last updated: 2026-09-13.

Status markers (same as the ledger):

- **[V]** verified empirically against the code/data in this repo
- **[A]** agreed direction, not yet executed
- **[P]** parked — open question

---

## 0. Governing principles **[A]**

**P1 — Additive S3.** Classes tag the native type; a `tales` *is* a tibble, a
matrix view *is* a matrix. dplyr compatibility is preserved rather than
replaced.

**P2 — Invariant vs precondition.** The distinction that resolved most of the
design questions below:

> **Invariants** are properties no legitimate operation breaks; they are
> validated. Properties that hold only of a *complete* object — and that
> ordinary subsetting rightly breaks — are **method preconditions**, checked
> by the methods that need them.

Applied: key uniqueness is an invariant; array completeness, matrix
squareness and coordinate contiguity are preconditions. The test is whether a
*reasonable* user operation breaks the property. `filter(domain_type ==
"repeat")` is reasonable, so completeness cannot be an invariant. Two rows
claiming the same slot of the same array is never meaningful, so uniqueness
can be.

**P3 — Long form is canonical.** Matrices are views produced by `as.matrix()`,
never alternative return values. This kills the dual-return branches and is
the structural answer to the 17 ad-hoc `melt`/`acast`/`dcast` sites recorded
in ledger §5.

**P4 — Preconditions are checked where they matter, with specific messages.**
Not on every dplyr verb. The codebase already works this way — `msa.R:189`
asserts squareness by hand exactly where the mafft matrix is written.

---

## 1. Naming **[A]**

Per ledger §9.0 (rOpenSci API guidelines), `object_verb()` with the object
first:

| current | becomes | note |
|---|---|---|
| `split_list()` | `as_tales()` | ingests fasta / `XStringSet` / list |
| `tale_parts(dir)` | `tales_from_telltale(dir)` | may be called inside `tell_tales()` |
| `tell_tales()` | `tell_tales()` | unchanged |
| `distalr()` | `tales_relatedness()` | see 1.2 |
| `group_tales()` | `tales_group()` | |
| `build_repeat_msa()` | `tales_align()` | |
| `plot_tales_msa()` | `plot.tales_msa()` | becomes a method, not a function |
| `talvez()`, `preditale()` | `tales_predict_targets()` | generic; both survive as backends, see 1.3 |

Base generics stay verb-first by necessity: `print()`, `format()`, `plot()`,
`as.matrix()`, `[`.

### 1.2 `distalr()` → `tales_relatedness()` **[A]**

**[V]** The current name refers to a tool the function no longer runs. Its one
external call is to ARLEM (`distalr.R:627`); repeat dissimilarities now come
from DECIPHER / mmseq2 / Biostrings via `aln_method`. DISTAL survives only as
an *input format* — `.format_repeat_dist_mat()` reads a `*_Repeatmatrix.mat`
file produced by Distal-1.2, used by `build_repeat_msa()` (`msa.R:178`). That
compatibility is kept; the name is not.

`tales_relatedness()` over `tales_similarity()` deliberately: the latter would
collide confusingly with the `pairwise_sim` class it *returns*.

**Documentation requirement.** The rename must not obscure the provenance —
the function is largely an R rewrite of the DisTAL Perl code. The existing
description already states this well (`distalr.R:455-459`: "An R
re-implementation of the original DisTal Perl program: it still uses the Arlem
binary for the same repeat-array alignment step…"); what it lacks is a
**citation**. Add the QueTAL suite paper, already used elsewhere in the
package for exactly this tool — <https://doi.org/10.3389/fpls.2015.00545>
(`README.md:28`, `AnnoTALE_QueTAL_functions_library.R:156`).

**ARLEM's attribution exists but is in the wrong place** **[V]**. It is emitted
as three `logger::log_info()` lines immediately before the system call
(`distalr.R:624-626`), relaying ARLEM's own console banner:

```
Running ARLEM version 1.0 :
Copyright by Mohamed I. Abouelhoda
Plz. cite Abouelhoda, Giegerich, Behzadi, and Steyaert
```

A runtime log line is not a citation: it is invisible to anyone reading
`?tales_relatedness`, and only appears at all if the logger threshold admits
INFO. Move it into `@references` (keeping the log line is fine — it faithfully
relays the tool's own banner).

**[V]** The publication, confirmed against the Crossref registry from that
author list:

> Abouelhoda M.I., Giegerich R., Behzadi B., Steyaert J.-M. (2009). *Alignment
> of minisatellite maps based on run-length encoding scheme.* Journal of
> Bioinformatics and Computational Biology **7**(2), 287–308.
> <https://doi.org/10.1142/S0219720009004060>

(An earlier version of the same work by the same four authors appeared at APBC
2007, doi:10.1142/9781848161092_0028.) The fit is not incidental: ARLEM aligns
*minisatellite maps* — tandem repeat arrays — which is structurally what a TALE
central repeat domain is, and presumably why DisTAL reached for it.

### 1.3 `talvez()` / `preditale()` → `tales_predict_targets()` **[A]**

**[V]** The two outputs are *already* deliberately unified. Both roxygen blocks
carry the same sentence — "column names have been modified … **in order to
homogenize column names across TALE target prediction programs in tantale**" —
and both rename into the same vocabulary (`subjSeqId`, `score`, `ebeSeq`,
`taleId`, `start`, `end`, `strand`; `target_predictions.R:67`, `:186`). The
package already treats these as one operation with interchangeable backends;
the API just never said so.

```r
tales_predict_targets(x, subject, method = c("talvez", "preditale"), ...)
```

Tool-specific arguments (`talvez_dir`/`conda_bin` vs `predictor_path`, and the
differing `opt_param` defaults `"-t 0 -l 19"` vs `""`) pass through `...`, and
**must be documented per backend** rather than left to `...`'s usual vagueness.

**`talvez()` and `preditale()` remain exported, with their own help pages.**
This is deliberate, not a transitional courtesy: their documentation carries
the citations to the original tools — Talvez
<https://doi.org/10.1371/journal.pone.0068464> (`target_predictions.R:92`) and
PrediTALE <https://www.jstacs.de/index.php/PrediTALE>
(`target_predictions.R:10`) — and folding them into one page would orphan
those attributions. They are also the natural place to document each backend's
own options. Naming them after the tools they wrap is a legitimate exception to
`object_verb` (rule 1): they are proper nouns.

**[P]** PrediTALE's reference is a tool page, not a paper; Talvez's is a DOI.
Worth levelling up if PrediTALE has a citable publication.

### 1.4 Provenance is a design constraint **[A]**

Generalising from 1.2 and 1.3: `tantale` is largely a wrapper and
reimplementation layer over published tools, so **no rename or merge may
orphan a tool citation**. Any function that disappears into a generic must
leave its references somewhere a user can still find them — normally by
surviving as a documented backend.

### 1.1 Column names **[A]**

Ledger §9.2 settles snake_case for all table columns, and notes it is a
prerequisite for the similarity-table unification. The target schema for
`tales`:

| current | target |
|---|---|
| `arrayID` | `array_id` |
| `positionInArray` | `position_in_array` |
| `positionInCrd` | `position_in_crd` |
| `domainType` | `domain_type` |
| `aaSeq` | `aa_seq` |
| `dnaSeq` | `dna_seq` |
| `domCode` | `dom_code` |
| `sourceDirectory` | `source_directory` |
| `rvd` | `rvd` |
| `seqnames` | `seqnames` — **kept**, deliberately, for Bioconductor familiarity (`GenomicRanges`/`Biostrings` spell it this way). The one sanctioned exception to snake_case in this schema. |

The rest of this document uses the target names.

---

## 2. Class `tales`

Class vector: `c("tales", "tbl_df", "tbl", "data.frame")`.

### 2.1 Identity **[A]**

> A `tales` object is a **long table of TALE array parts**: one row per
> (array, slot), describing the ordered decomposition of one or more TALE
> arrays into their constituent domains.

**Primary key: `array_id` × `position_in_array`.** **[V]** — `distalr()`
already enforces exactly this and stops if violated (`distalr.R:501-509`);
that check moves into the validator.

Row order is not part of identity — `position_in_array` carries the order — so
`arrange()` is always safe.

### 2.2 Two coordinate systems **[V]**

Deliberate and both needed:

- `position_in_array` — position among **all** parts, including termini. Starts
  at 1 (the N-terminus, when present).
- `position_in_crd` — position within the **central repeat domain**: numbers
  the repeats only, `NA` on termini, starting at 1 at the first repeat.

They are related by **[V]** (holds for all 44 arrays of the test fixture):

```
position_in_crd == position_in_array - (number of non-repeat parts
                                        at lower position_in_array in that array)
```

collapsing to `position_in_array - 1` only when an N-terminus is present.

This is why a repeats-only subset is a *legitimate* `tales` and not a broken
one: `filter(domain_type == "repeat")` leaves `position_in_crd` contiguous
from 1 while `position_in_array` no longer starts at 1 **[V]**. Any invariant
demanding contiguity of `position_in_array` would wrongly outlaw it.

### 2.3 Column contract **[A]**

| tier | columns | rule |
|---|---|---|
| required | `array_id` (chr), `position_in_array` (int) | error if absent or wrong type |
| required, ≥1 of | `rvd`, `dom_code` (chr) | error if neither present |
| optional, validated if present | `domain_type`, `position_in_crd`, `aa_seq`, `dna_seq`, `seq_name`, `source_directory` | error only if present *and* malformed |
| free | anything else (strain, host, clade, …) | silently preserved, never warned about |

Deliberately **not** required:

- `seq_name` — produced only by the telltale path (`distalr.R:161-165`); a
  fasta-derived `tales` has no source sequence.
- `rvd` specifically — the `dom_code` workflow needs the other one.

Deliberately **no warning on extra columns**: users will `mutate()` on
metadata, and a class that complains each time trains them to ignore warnings.
Errors are reserved for a broken contract.

Fix on ingest: `position_in_array` is currently double **[V]**; coerce to
integer.

### 2.4 Invariants **[A]**

**Hard — always validated:**

1. `array_id` × `position_in_array` unique **[V]**.
2. `position_in_array` is a positive integer, never `NA`.
3. `array_id` never `NA`; the required residue column never `NA`.
4. **`dom_code` is carried, never recomputed** **[V]** (ledger §5). It is a
   whole-set-dependent surrogate key — `cur_group_id()` over `aa_seq`
   (`distalr.R:512-516`) — so recomputing it on a subset renumbers everything
   and silently breaks the join to both similarity tables. No constructor or
   method may regenerate it from a subset.

An empty `tales` (zero rows) is **valid** — `filter()` legitimately returns
nothing.

**Conditional — validated only when the column is present:**

5. `domain_type` ⊆ {`"N-terminus"`, `"repeat"`, `"C-terminus"`}, with **at
   most one** N-terminus and **at most one** C-terminus per array. If present,
   the N-terminus has `position_in_array == 1`; the C-terminus has
   `position_in_array` greater than every other row of its array. Both are
   absolute claims, so both survive row subsetting.
6. `position_in_crd` is `NA` exactly on non-repeat rows **[V]**, and unique
   within `array_id` where non-`NA`.
7. **Order agreement** between the two coordinate systems: among the repeats of
   an array, ranking by `position_in_crd` gives the same order as ranking by
   `position_in_array`.

   *Corrected during implementation.* This invariant originally stated the
   **exact offset relation** of §2.2. That is wrong: the relation is not closed
   under subsetting. After `filter(domain_type == "repeat")` the N-terminus is
   gone, so "non-repeat parts before it" is 0 on the subset while
   `position_in_crd` still carries the offset from the complete array — the
   relation fails on an object we agreed is valid. Order agreement is the part
   that survives, and the exact relation moves to a precondition
   (`tales_assert_complete()`). Caught by a test, not by inspection.
8. `aa_seq` ↔ `dom_code` is **bijective** — same `aa_seq` ⟺ same `dom_code`
   **[V]**. This is what makes `repeat_to_rvd_map_distalr()` well-defined
   instead of a guess.
9. `seq_name` constant within `array_id` **[V]**.

**Soft — warn:**

10. `dna_seq` present but `NA`/empty — matches `distalr()`'s current behaviour
    (`distalr.R:494-496`).

**Explicitly not invariants (preconditions instead):**

- **Completeness** — "array holds all its parts, numbered from 1", *and* the
  exact `position_in_crd` offset relation of §2.2. Required by
  `tales_align()`, since it is what makes the mafft back-mapping
  well-defined. Broken legitimately by `filter(domain_type == "repeat")`.
  Implemented as `tales_assert_complete()`.
- **Non-empty `aa_seq`** — required by `tales_relatedness()`, which errors
  today (`distalr.R:487-493`). A `tales` built from RVD strings has no
  `aa_seq` at all.
- Uniform array length — arrays legitimately differ in repeat count.

**Closure property [V]:** with completeness demoted, every invariant above is
closed under row subsetting — now *tested*, not just reasoned about, for
uniqueness, integrality, at-most-one-terminus, the bijection, constant
`seqnames` and order agreement. This is what makes the dplyr policy in §2.6
cheap. The one property that failed this test is recorded in invariant 7
above.

### 2.5 The anchor sentinels **[V]**

Termini carry sentinel values in `rvd`, and there are **three**, not two:
`"NTERM"`, `"CTERM"`, and `"XXXXX"` — the last meaning *terminus detected in
the CDS but no HMMer hit, identity unknown* (`telltale.R:141`,
`telltale.R:182`). `tale_parts()` already treats all three as anchor codes
(`distalr.R:132`).

They live in the same column as real RVDs, so the set belongs in one exported
constant — `tales_anchor_codes()` — rather than being retyped across four
files. Since termini participate in the alignment (§4), `"XXXXX"` is an
ordinary alignable symbol.

### 2.6 dplyr policy **[A]**

- `dplyr_reconstruct()` — checks **only the column contract** (§2.3). Keeps the
  class if it holds, degrades silently to tibble otherwise. Never re-checks
  rows: by the closure property, it doesn't need to.
- `[.tales` — the same contract check. **[V]** Necessary, and not redundant
  with the above: dplyr's `dplyr_col_select()` calls `dplyr_reconstruct()`
  *only* for plain `data.frame`/`data.table` (verified, dplyr 1.2.1); for a
  tibble subclass it relies entirely on the class's own `[`. Without
  `[.tales`, `select(x, -array_id)` silently keeps the class on an object that
  no longer satisfies the contract. Pinned by a test, since it depends on
  dplyr internals.
- `dplyr_col_modify()` — additionally checks key uniqueness, one
  `anyDuplicated()`. `mutate()` overwrites values in place and *can* collide
  the key (`mutate(array_id = seq_name)` merges arrays from one contig); row
  operations cannot. Cheap, and never a false positive.
- Degradation is **silent**, not warned: the only code that trips it is doing
  so deliberately (see §3.4).

`[` keeps tibble semantics (row/column indexing). Array-level selection is
`filter(array_id %in% ids)`, not an overloaded `[` — overriding `[` to mean
array-selection would break data-frame expectations for a class that *is* a
data frame.

### 2.7 Constructors **[A]**

Standard three-layer, all returning `tales`:

- `new_tales(x)` — bare, no checks, internal.
- `validate_tales(x)` — all of §2.4, errors naming the offending `array_id`s.
- `tales(x, ...)` — user-facing: data frame in, coerces types, validates.
- `as_tales(x, ...)` — fasta / `XStringSet` / list path (today's
  `split_list()`, `conversion.R:18-30`). Assigns `position_in_array` from
  vector order; no `domain_type`, no `aa_seq`.
- `tales_from_telltale(dir)` — today's `tale_parts()`, behaviour unchanged.

---

## 3. Class `pairwise_sim` (+ `tale_sim`, `repeat_sim`)

### 3.1 Identity **[A]**

> A `pairwise_sim` is the **long form of a square pairwise similarity** over
> one entity set — one row per ordered pair.

Parent `pairwise_sim`; subclasses `tale_sim` and `repeat_sim` add nothing
structural, only entity semantics for validation and printing, and a hook if
methods ever need to diverge.

> **Names provisional** — see ledger §9.6 **[P]**. Two open problems:
> `repeat_sim` understates its content (**[V]** 28% of the ids in the real
> table are *terminus* domains, not repeats, so `domain_sim` is the accurate
> name), and the stored quantity may become dissimilarity rather than
> similarity, which would invert the required column of §3.3. The structure
> below is unaffected either way; only the names and which value column is
> required would change.

### 3.2 Why a class at all **[V]**

The same long→square-matrix operation is hand-written **four times**, differing
only in hardcoded id column names: `classification.R:24`, `msa.R:331`,
`msa.R:759`, `conversion.R:271`. "Filter to a subset of entities" appears three
more times: `msa.R:329-330`, `msa.R:757-758`, `msa.R:188`. Plus the hand-rolled
squareness assertion at `msa.R:189`. One `as.matrix()` method and one
subsetting method collapse eight sites.

### 3.3 Canonical columns **[A]**

Current state **[V]**:

| | ids | extras |
|---|---|---|
| `repeat.similarity` | `RepU2`, `RepU1` (that order) | `Dissim`, `Sim` |
| `tal.similarity` | `TAL1`, `TAL2` | `arlemScore`, `maxLength`, `normArlemScore`, `Sim` |

Decision: **canonical id columns** `id1` and `id2` **[A]**, with entity type
carried by the subclass, not by column naming. Rejected the
alternative of keeping `TAL1`/`RepU1` and recording the names in an attribute:
an attribute naming columns goes stale the moment one is renamed, and would
have to be policed in `dplyr_col_modify()`. Canonical names also fix the
rule-5 wart that the two tables order their id columns differently.

Column contract — **minimal**:

| tier | columns |
|---|---|
| required | `id1`, `id2`, `sim` |
| optional | `arlem_score`, `max_length` |

Minimal on purpose: `msa.R:176` legitimately selects down to three columns,
dropping `dissim`. Requiring it would make correct existing code invalid.

*Corrected before implementation.* This table originally read `Sim`,
`Dissim`, `arlemScore`, `maxLength`, `normArlemScore` — the existing spellings,
carried over unchanged because only the *id* columns were being canonicalised
at the time. That contradicted §1.1 and ledger §9.2, which settled snake_case
for every table column. The value columns are renamed on the same terms as the
ids; legacy spellings are accepted on input and normalised, exactly as for
`tales`.

### 3.4 Invariants and preconditions **[A]**

**Invariant:** the column contract only.

**Preconditions, not invariants:**

- **Squareness** — n² rows over the n ids present **[V]** (63001 = 251²,
  1936 = 44²). Checked by `as.matrix()`/`sim_matrix()`, with a message naming
  the unbalanced ids — replacing `msa.R:189`'s opaque `all.equal` assertion.
  Cannot be an invariant: `conversion.R:248` filters asymmetrically
  (`RepU1 == refState`) *inside an `apply()` over every alignment column*, and
  `msa.R:329-330` filters one id column then the other, passing through an
  invalid intermediate on purpose.
- Diagonal `Sim` = 100, symmetry of `Sim` **[V]** — same treatment.

Phrasing squareness as "complete over the ids it contains" (not complete in
some absolute sense) is what lets a symmetric subset stay square.

**Degradation is expected and correct here:** `msa.R:190-191` overwrites both
id columns with raw hex bytes and then sets `colnames(...) <- NULL`,
deliberately demolishing the tibble into a mafft input file. Silent degradation
to a plain data frame is the right behaviour; a warning would be noise.

### 3.5 Cross-object coherence — the `dom_code` namespace **[A]**

**[V]** Three companion tables are keyed by `dom_code`, not one:
`repeat.similarity`, `repeats.cluster` (`RepID`) and `repeats.code` (`code`);
`tal.similarity` is keyed by `array_id`.

**[V]** "Subset coherently" turned out to be the wrong framing — neither
direction of subsetting fails silently:

- Shrinking a `tales` is safe. `tales_align()` requires the similarity table to
  *cover* the residues present, then narrows it itself
  (`msa.R:187-189`): extra rows are fine.
- Shrinking the similarity table fails **loudly**, at that same `stopifnot()`.

The genuinely silent failure is **mixing objects from two different runs**.
`dom_code` is `cur_group_id()`, so every run mints `1..N`; the ranges overlap,
a cross-run join *succeeds*, and repeats are mapped to the wrong sequences with
nothing to signal it. Recomputation is just the special case where a run is
remade.

**Decision — a namespace tag, not a container.** Each `tales_relatedness()`
run is stamped with a `dom_code` namespace identifier: a content hash of the
sorted unique `aa_seq` set at construction, carried as an attribute on the
`tales` and on every `dom_code`-keyed companion. Methods consuming two of them
compare tags and error on mismatch.

Properties:

- Enforces coherence without a container, so §5's "no top-level
  session/project object" stands.
- Exactly parallel to `dom_code` itself — **stamped once, carried, never
  recomputed** — and therefore stable under subsetting, because it is stamped
  rather than derived.
- Content-hashed rather than random, so two runs over identical data are
  correctly recognised as compatible instead of falsely flagged.
- **[V]** Costs no machinery: custom attributes already survive `filter()`,
  `mutate()`, `arrange()`, `[` and `select()` on a tibble, so the tag rides
  along for free and `dplyr_reconstruct()` preserves it.

Limits, on the record: the tag is only as good as the methods that check it,
and stripping attributes defeats it. It catches accidents, not adversaries —
the right target.

---

## 4. Class `tales_msa`

Class vector: `c("tales_msa", "tales", "tbl_df", "tbl", "data.frame")`.

### 4.1 Identity **[A]**

> A `tales_msa` is a `tales` that additionally carries an **alignment
> coordinate**: one row per aligned part, with gaps represented by the
> *absence* of a row.

### 4.2 Gaps are implicit **[A]**

This supersedes an earlier decision that gaps would be explicit rows with
`alignment_position` set and the residue columns / `position_in_array` /
`position_in_crd` set to `NA`. That form is **incompatible with being a
subclass**: `tales` invariant 2 forbids `NA` `position_in_array`, and an array
has many gaps, so every gap row would carry the same `NA` and break invariant
1's key uniqueness too. A `tales_msa` would not have been a valid `tales`, and
every inherited method would have needed a special case.

Rejected alternatives: making `tales_msa` a sibling class (honest, but
duplicates the whole contract and inherits nothing), and relaxing `tales`'
invariants to permit `NA` positions (dismantles the parent's key to
accommodate the child).

With gaps implicit, **every parent invariant holds unchanged** and the
subclass relation is real. Both keys coexist: `array_id` × `position_in_array`
(inherited) and `array_id` × `alignment_position` (new) are each unique.

The `NA`-filled rectangular form does not disappear — it is what `as.matrix()`
produces, and what `build_repeat_msa(gap_symbol = NA)` already returns. It
simply stops being the canonical storage.

### 4.3 Column contract **[A]**

Everything in §2.3, plus:

| tier | column | rule |
|---|---|---|
| required | `alignment_position` (int) | error if absent or wrong type |

Attributes: the `dom_code` namespace tag (§3.5), inherited, plus the
**alignment width** `L` (see 4.4).

### 4.4 Invariants **[A]**

All of §2.4, unchanged — that is the point of 4.2. Plus:

11. `alignment_position` is a positive integer, never `NA`.
12. `alignment_position` is unique within `array_id`.
13. **Order agreement**: within an array, ranking rows by `alignment_position`
    gives the same order as ranking by `position_in_array`. An alignment may
    insert gaps but may never reorder parts. Nothing in the current code
    checks this.
14. `alignment_position` ⊆ `1:L`.

`L` (alignment width) is stamped as an **attribute**, not derived. Deriving it
as `max(alignment_position)` is correct only while some array still occupies
the last column; subsetting arrays can silently shrink it. Carried like the
namespace tag, for the same reason. Read with `tales_width()`.

**Degradation is graded** **[A]**, a refinement found while implementing:
dropping `alignment_position` while keeping the `tales` contract demotes the
object to a plain `tales` rather than all the way to a tibble. Both `[` and
`dplyr_reconstruct()` go through one shared regrade step, so `tales_msa`
needs no `[` method of its own — it inherits `[.tales`.

**Preconditions, not invariants:**

- ~~**Grid completeness**~~ — **struck during implementation.** With gaps
  implicit, `as.matrix()` needs no completeness precondition at all: a missing
  row *is* a gap, so the grid is always well defined, and the method simply
  fills `arrays × 1:L` and places the rows it has. This is an unlooked-for
  benefit of §4.2 — the implicit-gap representation removed a precondition
  rather than adding one.
- **No all-gap column** — true of mafft output, but subsetting arrays can
  empty a column. Honest under the implicit-gap representation: such a column
  simply has no rows.

### 4.5 Construction — `tales_align()` **[A]**

Takes a `tales`, returns a `tales_msa`: **one input type, one return type**.

The mafft back-mapping is **positional**: the k-th non-gap cell of an aligned
row is the k-th part fed in (`msa.R:218-232`), recoverable as
`cumsum(!is.na(row))`. That is well-defined only because `tales_align()`
builds the mafft input from the `tales` itself rather than accepting a
pre-made fasta — which is why the input type is not a union.

- **Precondition**: array completeness (§2.4) — each array must hold all its
  parts, contiguously from 1 — since that is what makes the k-th cell
  correspond to the k-th part.
- A residue-column argument selects what is aligned (`rvd` or `dom_code`),
  replacing `build_repeat_msa()`'s inference from a hardcoded list of six
  frequent RVDs (ledger §6).
- **Termini participate** — current behaviour, since `coded.repeats.str` is
  built from all parts including termini (`distalr.R:526`).
- A fasta reaches `tales_msa` by composing `as_tales()` first; no separate code
  path, and no matrix-only return.

### 4.6 Views and methods **[A]**

- `as.matrix(x, value = "rvd")` — materialises the grid; gaps become `NA` (or a
  `gap_symbol`). This is the one place the rectangular form is built, replacing
  the ad-hoc `melt`/`acast` sites of ledger §5.
- `plot(x, fill = , label = )` — replaces
  `plot_tales_msa(repeat_align, rvd_align, tal_sim, repeat_sim, …)`. The
  four-argument signature exists only because a matrix holds one value layer;
  a long `tales_msa` carries them all. Requesting a layer the object lacks
  (e.g. `dom_code` on an msa built from an RVD fasta) is a method precondition,
  not a structural difference.
- `tales_consensus_match()` becomes a method returning a `tales_msa` with a
  logical layer, collapsing its `long` flag — and removing the mislabelled
  `positionInArray` column recorded in ledger §6.

---

## 5. Open questions

All structural questions are closed. What remains is a small documentation
debt, carried from §1.2 and §1.3.

| # | question | ref |
|---|---|---|
| ~~1~~ | ~~`seqnames` → `seq_name`?~~ **Resolved** — kept as `seqnames` | §1.1 |
| ~~2~~ | ~~Canonical similarity id column names~~ **Resolved** — `id1`/`id2` | §3.3 |
| ~~3~~ | ~~Is cross-object `dom_code` coherence enforced or documented?~~ **Resolved** — enforced by a namespace tag; no container needed | §3.5 |
| ~~4~~ | ~~`tales_msa` key~~ **Resolved** — `array_id` × `alignment_position`; the inherited key remains valid too, since gaps are implicit | §4.2 |
| ~~5a~~ | ~~`distalr()`'s new name~~ **Resolved** — `tales_relatedness()`, with the DisTAL provenance cited | §1.2 |
| ~~5b~~ | ~~One generic over both prediction backends?~~ **Resolved** — yes; `talvez()`/`preditale()` stay exported and documented | §1.3 |
| ~~6~~ | ~~ARLEM has no citation anywhere in the package~~ **Resolved** — it was in runtime log lines, not the docs; publication identified and verified, move to `@references` | §1.2 |
| 7 | **[P]** PrediTALE is cited by tool page, not paper — level up if a publication exists | §1.3 |
