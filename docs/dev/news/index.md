# Changelog

## tantale (development version)

### Breaking change: the `tales` class system replaces the list-of-tables API

The pipeline now passes typed objects rather than named lists of plain
data frames.
[`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
(formerly `distalr()`) returns three of them:

| element            | class                                                          |
|--------------------|----------------------------------------------------------------|
| `tales`            | `tales` — the parts table, stamped with a `dom_code` namespace |
| `domain_distances` | `domain_distances` / `pairwise_distances`                      |
| `tale_distances`   | `tale_distances` / `pairwise_distances`                        |

#### Removed functions

These were deprecated in an earlier development commit and are now gone.
There is no alias; `master` still carries the old API if you need a
reference.

| removed              | replacement                                                                                       |
|----------------------|---------------------------------------------------------------------------------------------------|
| `distalr()`          | [`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)             |
| `tale_parts()`       | [`tales_from_telltale()`](https://scunnac.github.io/tantale/dev/reference/tales_from_telltale.md) |
| `split_list()`       | [`as_tales()`](https://scunnac.github.io/tantale/dev/reference/as_tales.md)                       |
| `build_repeat_msa()` | [`tales_align()`](https://scunnac.github.io/tantale/dev/reference/tales_align.md)                 |
| `group_tales()`      | [`tales_group()`](https://scunnac.github.io/tantale/dev/reference/tales_group.md)                 |

`distalr()`’s six-element list is reduced to the three typed elements
above. `coded.repeats.str`, `repeats.code` and `repeats.cluster` were
verified to be pure projections of the parts table or to have no
consumer at all, and are replaced by
[`tales_coded_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_coded_strings.md),
[`tales_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_codes.md)
and an explicit clustering call respectively.

#### Bug fix: repeat clustering was computed on an inverted matrix

`.repeat_to_cluster_align()`, which colours the repeat-cluster fill in
`plot_tales_msa()` and `msa_heatmap()`, passed a **similarity** matrix
straight to [`as.dist()`](https://rdrr.io/r/stats/dist.html).
[`as.dist()`](https://rdrr.io/r/stats/dist.html) expects a distance, so
the dendrogram was built upside down and the resulting clusters were
wrong.

Fixed by inverting to a distance first. Because the cut height is now
read on a distance scale, the `h_cut` default changes from **90 to 10**
in both plotting functions; if you pass `h_cut` explicitly, subtract it
from 100.

The clusters change: on the reference output, cutting the corrected tree
yields 45 clusters where the old code gave 59, with 93.8% pair-agreement
between the two partitions.

`.cluster_repeats()` had the same defect but became unreachable when
`distalr()` was removed, and has been deleted.

#### `msa_heatmap()` is retired

Superseded by `plot_tales_msa()` and moved to `inst/legacy/`. Its six
`plot_type` values were combinations of features that `plot_tales_msa()`
exposes as orthogonal arguments; the one thing it did that
\`plot_tales_msa() could not – draw the consensus row – is now
implemented there.

| `msa_heatmap(plot_type =)`   | `plot_tales_msa()`           |
|------------------------------|------------------------------|
| `"repeat.similarity"`        | `fill_type = "repeat_sim"`   |
| `"repeat.clusters"`          | `fill_type = "repeat_clust"` |
| `"with.rvd"`                 | pass `rvd_align`             |
| `"repeat.clusters.with.rvd"` | both of those                |
| `"reference"`                | pass `ref_pattern`           |
| `"consensus"`                | `consensus = TRUE`           |

`save_path` has no equivalent because none is needed: the returned
ggplot is `ggsave()`-able. `note_colors` likewise – add a scale to the
returned plot.

#### `plot_tales_msa()` can draw the consensus

`consensus = TRUE` now adds a consensus row above the alignment, where
it previously did nothing and the argument was documented as “NOT
IMPLEMENTED YET”. The consensus is the most frequent element per column,
taken from `rvd_align` when supplied and `repeat_align` otherwise, so it
always matches what the cells are labelled with.

It is drawn as its own panel above the alignment rather than as an extra
row. That is not cosmetic: `aplot` reorders the alignment’s y axis onto
the tree’s leaves, so a row the tree has no leaf for is silently
dropped.

This was the one feature `msa_heatmap()` still provided that
`plot_tales_msa()` did not.

#### Messaging is now cli throughout

The package mixed four messaging idioms: `logger`,
[`cat()`](https://rdrr.io/r/base/cat.html), base
[`message()`](https://rdrr.io/r/base/message.html) and base
[`stop()`](https://rdrr.io/r/base/stop.html)/[`warning()`](https://rdrr.io/r/base/warning.html).
It now uses `cli` only, and `logger` has been dropped from `Imports`.

This also fixes a class of bug rather than being cosmetic. Twenty-three
call sites wrote their diagnostic to the logger and then raised a bare
[`stop()`](https://rdrr.io/r/base/stop.html) or
[`warning()`](https://rdrr.io/r/base/warning.html) – conditions whose
message was the empty string.
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) saw nothing,
tests could not assert on them, and if your logger threshold excluded
ERROR the failure was silent. Errors now carry their message on the
condition, with a `tantale_error` class.

One guard was simply broken: `distalr()` tested `aln_method` with
`logger::log_errors() && stop("...")`. `log_errors()` installs a global
error handler rather than returning a predicate, so the helpful message
was unreachable and an invalid `aln_method` produced an unrelated
complaint about calling handlers.

#### Similarity became distance

`pairwise_distances` (formerly `pairwise_sim`) stores `dissim`, not
`sim`, and `tale_sim`/`repeat_sim` are now
`tale_distances`/`domain_distances`.

Two reasons. First, `repeat_sim` understated its content: 71 of the 251
ids in the reference output (28%) are terminus domains rather than
repeats, so `domain` is the accurate word. Second, the distance is the
*native* quantity in both tables — the aligners emit a dissimilarity,
and the TALE-level score is an ARLEM cost divided by array length —
while `Sim` was a derived `100 - x` that every consumer immediately
inverted back. Storing one quantity rather than two also removes the
risk of the two drifting apart.

Legacy tables are still accepted on input: a data frame carrying `Sim`,
`TAL1`/`TAL2`, `RepU1`/`RepU2` or `normArlemScore` is converted to the
canonical `id1`/`id2`/`dissim` form at construction.

Note the changed reading of the numbers: what used to display as a TALE
similarity of 91.75–100 is the same information shown as a distance of
0–8.25.

### Breaking change: every column is now `snake_case`

The tables `tantale` returns used four naming conventions at once, one
of them with a literal space in a column name. They now use one.

#### The `tales` table

`arrayID`, `domainType`, `positionInCrd`, `dnaSeq`, `sourceDirectory`,
`positionInArray`, `aaSeq` and `domCode` become `array_id`,
`domain_type`, `position_in_crd`, `dna_seq`, `source_directory`,
`position_in_array`, `aa_seq` and `dom_code`.

[`tales()`](https://scunnac.github.io/tantale/dev/reference/tales.md)
still accepts the old spellings and renames them, so a
[`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)
output directory written by an earlier version still loads.

#### The distance tables

`domain_distances` and `tale_distances` previously disagreed about how
to name the pair being compared – `RepU1`/`RepU2` in one, `TAL1`/`TAL2`
in the other, and in opposite column orders. Both now use `id1`, `id2`
and `dissim`, with `arlem_score` and `max_length` alongside for the TALE
table.

[`pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md)
accepts `TAL1`/`RepU1`/`Sim`/`Dissim`/`arlemScore` and renames them, and
`plot_tales_msa()` puts its `tal_sim` and `repeat_sim` arguments through
it, so passing an old-format table still works.

#### On-disk output

[`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)
writes `array_id` instead of `arrayID` in `arrayReport.tsv`,
`domainsReport.tsv` and `hitsReport.tsv`, and the GFFs derived from them
carry an `array_id` attribute. **Scripts that read these files by column
name need updating.**

#### Fixed along the way

- `.tale_parts_from_file()` named its column `arrayIDs` when the input
  was empty and `arrayID` otherwise, so the two returns had incompatible
  schemas.
- [`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map_distalr.md)
  and
  [`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/dev/reference/tale_parts_to_rvd.md)
  are documented as taking a
  [`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
  result but read the pre-class column names, so both had been broken
  against that result.
- `.repeat_to_cluster_align()` recovered a distance by computing
  `100 - Sim` before clustering, which assumed the scale ran 0-100. It
  now reads the stored distance.
- Six `if (class(x) == "...")` comparisons became
  [`inherits()`](https://rdrr.io/r/base/class.html).

### Breaking change: package-wide naming overhaul

The package had accumulated two incompatible naming styles (camelCase,
snake_case, and some dot.case) across its history, which made the API
hard to predict and, in argument names, occasionally collided with R’s
own S3 dispatch conventions. Every exported function, and most internal
ones, have been renamed to a single consistent `snake_case` style. There
is no backward-compatible alias for any old name — this is a clean
break, not a deprecation cycle. If you have scripts using the old API,
use the `master` branch, which still has the old names, as a reference
while you update.

#### Exported function renames

| Old name                                                                                                                                                            | New name                                                                                                      |
|---------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
| `FuncTAL()`                                                                                                                                                         | [`functal()`](https://scunnac.github.io/tantale/dev/reference/functal.md)                                     |
| `analyzeAnnoTALE()`                                                                                                                                                 | [`run_annotale_predict()`](https://scunnac.github.io/tantale/dev/reference/run_annotale_predict.md)           |
| `buildAnnoTALE()`                                                                                                                                                   | [`run_annotale_build()`](https://scunnac.github.io/tantale/dev/reference/run_annotale_build.md)               |
| `buildRepeatMsa()`                                                                                                                                                  | `build_repeat_msa()`                                                                                          |
| `convertRepeat2RvdAlign()`                                                                                                                                          | `repeat_to_rvd_align()`                                                                                       |
| `correcTales()`                                                                                                                                                     | [`correct_tales()`](https://scunnac.github.io/tantale/dev/reference/correct_tales.md)                         |
| `diagnoseTaleParts()`                                                                                                                                               | `diagnose_tale_parts()`                                                                                       |
| `getRepeat2RvdMapping()`                                                                                                                                            | [`repeat_to_rvd_map()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map.md)                 |
| `getRepeat2RvdMappingFromDistalr()`                                                                                                                                 | [`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map_distalr.md) |
| `getTaleParts()`                                                                                                                                                    | `tale_parts()`                                                                                                |
| `ggplotTalesMsa()`                                                                                                                                                  | `plot_tales_msa()`                                                                                            |
| `groupTales()`                                                                                                                                                      | `group_tales()`                                                                                               |
| `heatmap_msa()`                                                                                                                                                     | `msa_heatmap()`                                                                                               |
| `heatmap_talomes()`                                                                                                                                                 | [`talomes_heatmap()`](https://scunnac.github.io/tantale/dev/reference/talomes_heatmap.md)                     |
| `matchConsensus()`                                                                                                                                                  | [`tales_consensus_match()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus_match.md)         |
| `plotTaleComposition()`                                                                                                                                             | `plot_tale_composition()`                                                                                     |
| `plotTaleTargetPred()`                                                                                                                                              | [`plot_target_preds()`](https://scunnac.github.io/tantale/dev/reference/plot_target_preds.md)                 |
| `taleAlignConsensus()`                                                                                                                                              | [`tales_consensus()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus.md)                     |
| `taleParts2RvdStringSet()`                                                                                                                                          | [`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/dev/reference/tale_parts_to_rvd.md)                 |
| `tellTale()`                                                                                                                                                        | [`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)                               |
| `tellTale2()`                                                                                                                                                       | removed (was a deprecated alias for `tellTale()`)                                                             |
| `toListOfSplitedStr()`                                                                                                                                              | `split_list()`                                                                                                |
| [`preditale()`](https://scunnac.github.io/tantale/dev/reference/preditale.md), [`talvez()`](https://scunnac.github.io/tantale/dev/reference/talvez.md), `distalr()` | unchanged                                                                                                     |

#### Argument renames (package-wide conventions)

- `condaBinPath` -\> `conda_bin`
- `outputDir` / `outDir` -\> `output_dir` (previously inconsistent
  across functions)
- `taleSim` / `talsim` -\> `tal_sim`
- `repeatAlign` / `rvdAlign` / `repeatSim` -\> `repeat_align` /
  `rvd_align` / `repeat_sim`
- `repeat.clust.h.cut` / `repeats.cluster.h.cut` (dot.case) -\> `h_cut`
- `rvdSeqs` / `subjDnaSeqFile` -\> `rvd_seqs` / `subj_file`
- `taleParts` -\> `tale_parts`

Plus function-specific cleanups, notably: -
[`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)’s
verbose scoring arguments (e.g. `TALE_NtermDNAHitMinScore` -\>
`nterm_min_score`, `minDomainHitsPerSubjSeq` -\> `min_domain_hits`). -
[`talomes_heatmap()`](https://scunnac.github.io/tantale/dev/reference/talomes_heatmap.md)’s
cryptic column-selector arguments (`col`/`row`/`value` -\>
`group_col`/`strain_col`/`rvd_col`; `truncTaleLab` -\>
`trunc_tales_col`; `mapcol`/`sepwid`/`sepcol`/`mar.side` -\>
`colors`/`sep_width`/`sep_color`/`margins`).

#### Internal functions

Roughly 35 non-exported helper functions were also renamed to
`snake_case`, and every non-exported function now consistently starts
with a leading `.` (this was previously applied only in some files).
Since these were never part of the public API, this isn’t a breaking
change for package users, but it’s listed here for anyone reading the
source or a diff against an older version. Notable ones, since they
change *behavior*-adjacent naming rather than pure casing: -
`.distalPairwiseAlign()` / `2` / `3` -\> `.pairwise_align_biostrings()`
/ `.pairwise_align_mmseq2()` / `.pairwise_align_decipher()` (named after
the actual backend instead of an arbitrary number) -
`convertRepeat2SimAlign()`, `convertRepeat2ClusterIDAlign()`,
`convertRvd2RepeatAlign()`, `convertRvd2MatchAlign()` -\>
`.repeat_to_sim_align()`, `.repeat_to_cluster_align()`,
[`.rvd_to_repeat_align()`](https://scunnac.github.io/tantale/dev/reference/dot-rvd_to_repeat_align.md),
[`.rvd_to_match_align()`](https://scunnac.github.io/tantale/dev/reference/dot-rvd_to_match_align.md)
(joining the `_to_` naming family used by the exported conversion
functions) - `computeRVDSeqEBESeqMatchQualityString()` -\>
`.compute_match_string()` - `systemInCondaEnv()` / `createTantaleEnv()`
-\> `.run_in_conda()` / `.create_tantale_env()`

#### Other fixes made alongside the rename

- `group_tales()` and
  [`talomes_heatmap()`](https://scunnac.github.io/tantale/dev/reference/talomes_heatmap.md)
  called [`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html),
  [`viridis::scale_color_viridis()`](https://sjmgarnier.github.io/viridis/reference/scale_viridis.html)
  and
  [`tidytree::MRCA()`](https://rdrr.io/pkg/tidytree/man/MRCA.html)/`groupClade()`
  without declaring `cluster`, `viridis`, or `tidytree` in
  `DESCRIPTION`. Fixed.
- Fixed a handful of prose/log messages that referred to the AnnoTALE
  tool by name and would otherwise have been corrupted by the `annoTALE`
  -\> `annotale_jar` argument rename (e.g. “Now running AnnoTALE predict
  for…” was at risk of becoming “Now running annotale_jar predict
  for…”).
