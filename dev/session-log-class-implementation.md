# Implementation session log — `tales` / `tales_msa`

Written while you were away, for review. Covers commits `2e92e47`, `16d8458`,
`1d06ca1`, `f0e0b0a`. Everything below is on branch `dev` and committed; nothing
is pushed.

Status markers as elsewhere: **[V]** verified empirically, **[A]** agreed
direction, **[D]** a decision I took unilaterally and that you may want to
revisit.

---

## 1. What was built

| file | contents |
|---|---|
| `R/tales_class.R` | `tales()`, `new_tales()`, `validate_tales()`, `as_tales()`, `tales_from_telltale()`, `is_tales()`, `tales_namespace()`, `tales_anchor_codes()`, `tales_assert_complete()`, `[.tales`, `dplyr_reconstruct.tales`, `dplyr_col_modify.tales` |
| `R/tales_msa_class.R` | `tales_msa()`, `new_tales_msa()`, `validate_tales_msa()`, `is_tales_msa()`, `tales_width()`, `as.matrix.tales_msa()`, `tales_align()` |
| `tests/testthat/test_tales_class.R` | 71 passing |
| `tests/testthat/test_tales_msa_class.R` | 42 passing, including live MAFFT |

Adjacent suites re-run and unaffected: `test_tale_parts.R` (4),
`test_split_list.R` (6), `test_build_repeat_msa.R` (5), `test_group_tales.R`
(18). All pass.

---

## 2. Design errors implementation caught

These are corrections to `class-design.md`, already applied there.

**2.1 Invariant 7 was wrong — the offset relation is not closed under
subsetting** **[V]**

§2.4 stated the exact `position_in_crd = position_in_array - (non-repeat parts
before it)` relation as an invariant, and claimed all invariants had been
"checked individually" for closure. That claim was reasoning, not testing, and
it was false for this one: after `filter(domain_type == "repeat")` the
N-terminus is gone, so the offset computed *on the subset* is 0 while
`position_in_crd` still carries the offset from the complete array. The
relation fails on an object we had explicitly agreed is valid.

Replaced by **order agreement** (ranking by `position_in_crd` matches ranking
by `position_in_array`), which does survive. The exact relation moved into
`tales_assert_complete()`. A test caught this, not inspection — the doc now
says "tested, not just reasoned about".

**2.2 `dplyr_reconstruct()` is not the hook for column-dropping** **[V]**

Verified in dplyr 1.2.1's source: `dplyr_col_select()` calls
`dplyr_reconstruct()` *only* for plain `data.frame`/`data.table`; for a tibble
subclass it relies entirely on the class's own `[`. So `select(x, -array_id)`
kept the class until `[.tales` was added. Pinned by a test, since it depends on
dplyr internals.

**2.3 Grid completeness is not needed after all** **[V]**

§4.4 listed it as a precondition of `as.matrix()`. With gaps implicit, a
missing row *is* a gap, so the grid is always well defined — the method fills
`arrays × 1:L` and places whatever rows exist. The implicit-gap choice removed
a precondition rather than adding one. Struck from the doc.

---

## 3. Decisions I made on my own **[D]**

Listed most consequential first.

**3.1 Constructors accept legacy camelCase and normalise on the way in.**
`tales()` renames `arrayID` → `array_id`, `positionInArray` →
`position_in_array`, etc. via a lookup table. Without this the class would
validate nothing the current pipeline produces until the §9.2 rename sweep
lands, so the class would be untestable against real data. `seqnames` is
excluded from the mapping per your call. Easy to delete once §9.2 is done.
*Reversible, low risk.*

**3.2 `as_tales()` takes an explicit `residue_col` argument rather than
guessing.** Values land in `rvd` (default) or `dom_code`, chosen by the caller.
This is the design's §4.5 intent — replacing `build_repeat_msa()`'s inference
from a hardcoded list of six frequent RVDs (ledger §6) — applied to ingestion
as well as alignment. *Worth a look: the default of `"rvd"` may or may not be
what you want.*

**3.3 `as_tales()` is a generic with a `data.frame` method that forwards to
`tales()`.** So `as_tales()` works as a single entry point for every input
type. Not specified in the design either way.

**3.4 `tales_from_telltale()` wraps the existing `tale_parts()` rather than
reimplementing it.** It calls `tale_parts(dir)` and passes the result through
`tales()`. Keeps the parsing logic in one place; `tale_parts()` can be
deprecated later without touching this. *No behaviour change to `tale_parts()`
itself.*

**3.5 Degradation is graded.** Dropping `alignment_position` demotes
`tales_msa` → `tales` rather than straight to a tibble. Both `[` and
`dplyr_reconstruct()` share one regrade step, so `tales_msa` needs no `[`
method of its own. Felt clearly right; recorded in §4.4.

**3.6 `tales_align()` joins arrays into space-separated strings internally.**
MAFFT's input is built by pasting each array's residues with `" "` regardless
of layer, since neither RVDs nor repeat codes contain a space. The `sep`
distinction (`"-"` for RVDs, `" "` for codes) is thereby an ingestion concern
only, invisible in alignment. *If any RVD vocabulary can contain a space this
is wrong — I believe it cannot, but you would know better.*

**3.7 Private helpers use `@noRd`, not `@keywords internal`.** Only for the new
`.`-prefixed helpers, to avoid generating `man/dot-*.Rd` pages. Does not
pre-empt ledger §9.3, which is about auditing existing functions.

**3.8 New code uses `cli::cli_abort()` with classed conditions.** Classes are
`tantale_error_<specific>` plus a shared `tantale_error` parent, so tests assert
with `expect_error(class = )` rather than regex. This is §9.5's recommended
direction; new code seemed the right place to start exemplifying it rather than
adding to the four-idiom mess. No existing call sites were converted.

**3.9 Kept a local `%||%`.** Base R gained it in 4.4.0 but `DESCRIPTION`
declares `R (>= 3.6.3)`, so relying on base would break the stated minimum.

---

## 4. Known limitations, not addressed

- **[V] `tales_align(residue_col = "dom_code")` will fail on large sets.**
  `build_repeat_msa()` maps each distinct residue to one of ~247 usable ASCII
  bytes. The `sampleDistalrOutput` fixture already has **251** distinct
  `dom_code` values, which exceeds that. RVD alignment is unaffected (a few
  dozen distinct values). Pre-existing, not introduced here, but the typed API
  makes it reachable more directly — worth a ledger entry.
- `plot.tales_msa()` is not implemented; `plot_tales_msa()` is untouched.
- `pairwise_sim` is not implemented, so `tales_align()` passes `repeat_sims`
  through to `build_repeat_msa()` unchanged and unvalidated.
- The `dom_code` namespace tag is carried and checked for *presence*, but
  nothing yet *stamps* it — that belongs to `tales_relatedness()`, and the
  content hash is unimplemented (it needs a hashing dependency; see §5).
- No `print()`/`tbl_sum()` method, so a `tales` prints as a plain tibble.

---

## 5. Questions for you

1. **Hashing dependency.** §3.5's namespace tag is specified as a content hash
   of the sorted unique `aa_seq` set. Base R has no general hash function;
   `rlang::hash()` or `digest::digest()` would each be a new `Imports` entry
   (both are already installed transitively). Which, or something else?
2. **`as_tales()` default `residue_col = "rvd"`** — right default?
3. **Deprecation timing.** `tale_parts()`, `split_list()` and `build_repeat_msa()`
   are now each wrapped by a newer function but remain exported and undeprecated.
   Deprecate now, or once the whole class layer is in place?
4. Whether the `dom_code` ASCII-budget limit in §4 deserves a ledger §6 entry
   (I did not add one unprompted).
