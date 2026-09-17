# Unattended session, 2026-09-17

Five commits on `dev`, none pushed. Full suite **642 passing, 0 failing**
(was 555). `R CMD check`: **0 errors**, 1 warning, 2 notes — all three
pre-existing and unrelated (bundled executables, long paths under
`data_for_tests`, empty `NEWS.md`).

Everything below is written up at length in `restructuring-notes.md`; this is
the map.

| commit | ledger |
|---|---|
| `8f88aca` Curate the correction reference; expose `max_comparisons` | 8.1b, 9.7 |
| `c63c109` `tantale_setup()`; `tales_group()`; `sep` for coded strings | 7.4a, 5.2, 8.2b, 8.1d |
| `68f5d52` Document the conda prerequisite | 7.4b |
| `72c6e71` Toy fixture with a known frameshift | 8.1 |
| `5185722` Break `tales_compare()` into three steps | 8.5b |

---

## What needs you

**§9.2b — the snake_case sweep is not finished, and §9.2 said it was.**

This is the one thing I stopped on rather than doing. §9.2 claims "every
table the package produces now uses snake_case". In the three TSVs
`tell_tales()` writes, what it actually changed was `array_id`:

| file | snake_case | still legacy |
|---|---|---|
| `domainsReport.tsv` | all 4 | — |
| `hitsReport.tsv` | 10 of 12 | `nhmmerHitID`, `hitID` |
| `arrayReport.tsv` | 3 of 17 | the other 14 |

I corrected the overstated sentence in §9.2 and put the full inventory and a
proposed mapping in §9.2b. I did not do the rename: these are the column
names of the package's primary output files, and you are the one who knows
what reads them downstream. The work itself is mechanical and has golden
coverage — one pass, once you say so.

**§5.2 A-vs-B** is still open. Your message settled how `group` is
*populated* (done, below); it did not settle what a multi-group *alignment*
is, which is the other half.

---

## Your `tales_group()` call, implemented

`tales_group(x, tal_sim, ...)` now returns the `tales` object with `group`
filled, instead of a bare `data.frame(name, group)`.

The argument order is a judgement I made, not something you specified: `x`
first, following the data-first convention in §9.0 and the rest of the
`tales_*` API. Trivially flipped if you want it the other way.

Your option 1 turns out to buy more than symmetry. Taking `x` is the only
point where the correspondence between the distances and the object can be
checked, so `tales_group()` now errors if any array in `x` is ungrouped or
any grouped name is absent from `x` — i.e. "these distances did not come
from this object". A separate `add_group()` could do the same check, but
nothing would oblige a caller to route through it, and the failure it
prevents — a partly-grouped object — surfaces far downstream.

The bare mapping is still `unique(out[c("array_id", "group")])`.

---

## Done

**`tantale_setup()`** (§7.4a) — checks and optionally repairs the conda
environment, verifying versions against the yaml pins rather than presence.
All six requirements in the spec are met. Two things worth knowing:

- `.create_tantale_env()` used to announce on every run that an environment
  "has been found on your system and can be used for analysis" — without
  having looked inside it. Noise when true, false assurance when not, which
  is the exact failure 7.4a was written to catch. Now silent.
- `.tantale_conda_root()` exists because `dirname(dirname(bin))` is *not*
  the root — micromamba's binary sits in `~/bin` while its root is
  `MAMBA_ROOT_PREFIX`. That confusion is what made the 7.4 rebuilds land in
  a different root from the one in use.

Not exercised by tests: the `install = TRUE` and `conda = TRUE` branches,
which need network and would modify the machine.

**Conda documentation** (§7.4b) — README gains a real three-step install
section; `R/tantale.R` gains a matching `@section Setting up:`. Two stale
things fixed on the way: the README's claim that bundled Perl libraries
account for the package size (false since 7.4 — it is three Java jars), and
a `_pkgdown.yml` entry for `annout-class`, whose Rd vanished when the S4
class was retired. That last one is a **hard error** in
`pkgdown::build_site()`, so the site could not have been rebuilt.

**Toy fixture with a known answer** (§8.1) — three regions cut from the BAI3
sequences: one intact TALE, the same TALE with one base inserted mid-array,
and 4 kb with no TALE. Two copies of the *same* TALE is the design, so the
intact one is the control.

| | intact | frameshifted |
|---|---|---|
| longest ORF, correction off | 4305 nt (93%) | **2631 nt (57%)** |
| longest ORF, correction on | 4305 nt | **4305 nt** |
| RVD string, correction on | — | **identical to intact** |

Correction had never been asserted to correct anything — the branch was
covered only by a golden digest, which cannot notice correction silently
breaking, since once the changed digest is accepted it just records the new
wrong answer.

**`tales_compare()` decomposed** (§8.5b) — into
`tales_assign_domain_codes()` → `tales_domain_distances()` →
`tales_tale_distances()`, with `tales_compare()` as the wrapper. Golden
passed **44/44 with no snapshot changes**, so it is byte-identical to the
monolith.

Exporting step 1 makes `dom_code`'s run-dependence public, so the namespace
stamp is now enforced rather than advisory: `tales_tale_distances()` refuses
arguments from different runs. That matters more than it looks — step 3
indexes domains by code, so distances from another run would not fail, they
would silently compare the wrong domains.

**`tales_coded_strings(sep, repeats_only)`** (§8.2b) — added, defaults
unchanged. It has no callers inside `R/` at all, so nothing was coupled to
the separator. The siblings now take the same arguments with deliberately
*different* defaults, because they feed different consumers (`"-"` for
AnnoTALE, `" "` for ARLEM and MAFFT), and the docs say so.

**Golden baseline made portable** (§8.1d) — absolute directories are
rewritten to `<path>/` before digesting, keeping the basename.

---

## Three things I got wrong and caught

Recording these because each was caught by something specific, and the
something is reusable.

**1. The path normaliser ate content.** My first pattern matched "anything
between two slashes", which also rewrote `</title></head>`, HMMER's `//`
record separators and the `//` in `http://hmmer.org/`. It produced a
perfectly stable digest while quietly destroying content — a *weaker*
baseline wearing the appearance of a more portable one. Caught by auditing
which lines actually changed instead of trusting the tests to pass. Every
case the broad version got wrong is now a test.

**2. Rich docs written into a `@noRd` block.** The `correction_ref` and
`max_comparisons` prose was in the internal `.telltale_array_orfs()` block,
where no user would ever see it, while the exported `tell_tales()` had
legacy one-liners and no `@param max_comparisons` at all.

**3. A test fixture that hid a bug by being incomplete.** `coded_tales()`
carried `rvd = NTERM, NI, NTERM` — two N-termini, the second mid-array. Not
a TALE. It went unnoticed because the fixture had no `domain_type` column,
and without one the anomaly checks cannot run. Adding the column made
`tales()` flag it immediately.

That third one generalises: **a fixture that omits the columns the
validators key on is not exercising the validators.** Other minimal fixtures
in the suite are likely in the same position — worth a sweep sometime.

---

## Also

**Roxygen markdown is on** (§9.7), as you asked. It had been *written for*
for years without being enabled — `**precondition**` was reaching
`?tales_assert_complete` with the asterisks visible. Enabling it changed 31
of 64 Rd files, mostly rewrapping plus those fixes, and created two real
regressions: `[eg PacBio, ONT]` and `[see the min_gap parameter]` became
`\link{}` to targets that do not exist, because `[...]` is link syntax.
Both fixed. Verified by diffing the full set of link targets against the
pre-switch baseline (now identical) and parsing all 64 Rd files.

The rule for later: **square brackets in roxygen prose are a link.**

**Incidental, not acted on:** each toy region yields a spurious single-hit
array at its 3' end, because `min_domain_hits` filters *subject sequences*,
not arrays — a contig with enough hits overall keeps all of its arrays,
however small. Pre-existing. Worth a look if short spurious arrays ever
become a nuisance.

---

## Queue after this

- **§9.2b** — needs your call (above).
- **§5.2** — A vs B, needs your call.
- **§7.5a** article on the `tales` class and `dom_code`; **§7.5** vignettes
  and `@examples`. Both deliberately last, and §8.5b just changed the API
  they would describe, so waiting was right.
- `pipeline.svg` still shows old column names.
