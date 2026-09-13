# tantale (development version)

## Breaking change: the `tales` class system replaces the list-of-tables API

The pipeline now passes typed objects rather than named lists of plain data
frames. `tales_compare()` (formerly `distalr()`) returns three of them:

| element | class |
|---|---|
| `tales` | `tales` — the parts table, stamped with a `dom_code` namespace |
| `domain_distances` | `domain_distances` / `pairwise_distances` |
| `tale_distances` | `tale_distances` / `pairwise_distances` |

### Removed functions

These were deprecated in an earlier development commit and are now gone. There
is no alias; `master` still carries the old API if you need a reference.

| removed | replacement |
|---|---|
| `distalr()` | `tales_compare()` |
| `tale_parts()` | `tales_from_telltale()` |
| `split_list()` | `as_tales()` |
| `build_repeat_msa()` | `tales_align()` |
| `group_tales()` | `tales_group()` |

`distalr()`'s six-element list is reduced to the three typed elements above.
`coded.repeats.str`, `repeats.code` and `repeats.cluster` were verified to be
pure projections of the parts table or to have no consumer at all, and are
replaced by `tales_coded_strings()`, `tales_domain_codes()` and an explicit
clustering call respectively.

### Bug fix: repeat clustering was computed on an inverted matrix

`.repeat_to_cluster_align()`, which colours the repeat-cluster fill in
`plot_tales_msa()` and `msa_heatmap()`, passed a **similarity** matrix straight
to `as.dist()`. `as.dist()` expects a distance, so the dendrogram was built
upside down and the resulting clusters were wrong.

Fixed by inverting to a distance first. Because the cut height is now read on a
distance scale, the `h_cut` default changes from **90 to 10** in both plotting
functions; if you pass `h_cut` explicitly, subtract it from 100.

The clusters change: on the reference output, cutting the corrected tree yields
45 clusters where the old code gave 59, with 93.8% pair-agreement between the
two partitions.

`.cluster_repeats()` had the same defect but became unreachable when
`distalr()` was removed, and has been deleted.

### Messaging is now cli throughout

The package mixed four messaging idioms: `logger`, `cat()`, base `message()`
and base `stop()`/`warning()`. It now uses `cli` only, and `logger` has been
dropped from `Imports`.

This also fixes a class of bug rather than being cosmetic. Twenty-three call
sites wrote their diagnostic to the logger and then raised a bare `stop()` or
`warning()` -- conditions whose message was the empty string. `tryCatch()` saw
nothing, tests could not assert on them, and if your logger threshold excluded
ERROR the failure was silent. Errors now carry their message on the condition,
with a `tantale_error` class.

One guard was simply broken: `distalr()` tested `aln_method` with
`logger::log_errors() && stop("...")`. `log_errors()` installs a global error
handler rather than returning a predicate, so the helpful message was
unreachable and an invalid `aln_method` produced an unrelated complaint about
calling handlers.

### Similarity became distance

`pairwise_distances` (formerly `pairwise_sim`) stores `dissim`, not `sim`, and
`tale_sim`/`repeat_sim` are now `tale_distances`/`domain_distances`.

Two reasons. First, `repeat_sim` understated its content: 71 of the 251 ids in
the reference output (28%) are terminus domains rather than repeats, so
`domain` is the accurate word. Second, the distance is the *native* quantity in
both tables — the aligners emit a dissimilarity, and the TALE-level score is an
ARLEM cost divided by array length — while `Sim` was a derived `100 - x` that
every consumer immediately inverted back. Storing one quantity rather than two
also removes the risk of the two drifting apart.

Legacy tables are still accepted on input: a data frame carrying `Sim`,
`TAL1`/`TAL2`, `RepU1`/`RepU2` or `normArlemScore` is converted to the
canonical `id1`/`id2`/`dissim` form at construction.

Note the changed reading of the numbers: what used to display as a TALE
similarity of 91.75–100 is the same information shown as a distance of
0–8.25.

## Breaking change: package-wide naming overhaul

The package had accumulated two incompatible naming styles (camelCase, snake_case, and some dot.case) across its history, which made the API hard to predict and, in argument names, occasionally collided with R's own S3 dispatch conventions. Every exported function, and most internal ones, have been renamed to a single consistent `snake_case` style. There is no backward-compatible alias for any old name — this is a clean break, not a deprecation cycle. If you have scripts using the old API, use the `master` branch, which still has the old names, as a reference while you update.

### Exported function renames

| Old name | New name |
|---|---|
| `FuncTAL()` | `functal()` |
| `analyzeAnnoTALE()` | `run_annotale_predict()` |
| `buildAnnoTALE()` | `run_annotale_build()` |
| `buildRepeatMsa()` | `build_repeat_msa()` |
| `convertRepeat2RvdAlign()` | `repeat_to_rvd_align()` |
| `correcTales()` | `correct_tales()` |
| `diagnoseTaleParts()` | `diagnose_tale_parts()` |
| `getRepeat2RvdMapping()` | `repeat_to_rvd_map()` |
| `getRepeat2RvdMappingFromDistalr()` | `repeat_to_rvd_map_distalr()` |
| `getTaleParts()` | `tale_parts()` |
| `ggplotTalesMsa()` | `plot_tales_msa()` |
| `groupTales()` | `group_tales()` |
| `heatmap_msa()` | `msa_heatmap()` |
| `heatmap_talomes()` | `talomes_heatmap()` |
| `matchConsensus()` | `tales_consensus_match()` |
| `plotTaleComposition()` | `plot_tale_composition()` |
| `plotTaleTargetPred()` | `plot_target_preds()` |
| `taleAlignConsensus()` | `tales_consensus()` |
| `taleParts2RvdStringSet()` | `tale_parts_to_rvd()` |
| `tellTale()` | `tell_tales()` |
| `tellTale2()` | removed (was a deprecated alias for `tellTale()`) |
| `toListOfSplitedStr()` | `split_list()` |
| `preditale()`, `talvez()`, `distalr()` | unchanged |

### Argument renames (package-wide conventions)

- `condaBinPath` -> `conda_bin`
- `outputDir` / `outDir` -> `output_dir` (previously inconsistent across functions)
- `taleSim` / `talsim` -> `tal_sim`
- `repeatAlign` / `rvdAlign` / `repeatSim` -> `repeat_align` / `rvd_align` / `repeat_sim`
- `repeat.clust.h.cut` / `repeats.cluster.h.cut` (dot.case) -> `h_cut`
- `rvdSeqs` / `subjDnaSeqFile` -> `rvd_seqs` / `subj_file`
- `taleParts` -> `tale_parts`

Plus function-specific cleanups, notably:
- `tell_tales()`'s verbose scoring arguments (e.g. `TALE_NtermDNAHitMinScore` -> `nterm_min_score`, `minDomainHitsPerSubjSeq` -> `min_domain_hits`).
- `talomes_heatmap()`'s cryptic column-selector arguments (`col`/`row`/`value` -> `group_col`/`strain_col`/`rvd_col`; `truncTaleLab` -> `trunc_tales_col`; `mapcol`/`sepwid`/`sepcol`/`mar.side` -> `colors`/`sep_width`/`sep_color`/`margins`).

### Internal functions

Roughly 35 non-exported helper functions were also renamed to `snake_case`, and every non-exported function now consistently starts with a leading `.` (this was previously applied only in some files). Since these were never part of the public API, this isn't a breaking change for package users, but it's listed here for anyone reading the source or a diff against an older version. Notable ones, since they change *behavior*-adjacent naming rather than pure casing:
- `.distalPairwiseAlign()` / `2` / `3` -> `.pairwise_align_biostrings()` / `.pairwise_align_mmseq2()` / `.pairwise_align_decipher()` (named after the actual backend instead of an arbitrary number)
- `convertRepeat2SimAlign()`, `convertRepeat2ClusterIDAlign()`, `convertRvd2RepeatAlign()`, `convertRvd2MatchAlign()` -> `.repeat_to_sim_align()`, `.repeat_to_cluster_align()`, `.rvd_to_repeat_align()`, `.rvd_to_match_align()` (joining the `_to_` naming family used by the exported conversion functions)
- `computeRVDSeqEBESeqMatchQualityString()` -> `.compute_match_string()`
- `systemInCondaEnv()` / `createTantaleEnv()` -> `.run_in_conda()` / `.create_tantale_env()`

### Other fixes made alongside the rename

- `group_tales()` and `talomes_heatmap()` called `cluster::pam()`, `viridis::scale_color_viridis()` and `tidytree::MRCA()`/`groupClade()` without declaring `cluster`, `viridis`, or `tidytree` in `DESCRIPTION`. Fixed.
- Fixed a handful of prose/log messages that referred to the AnnoTALE tool by name and would otherwise have been corrupted by the `annoTALE` -> `annotale_jar` argument rename (e.g. "Now running AnnoTALE predict for..." was at risk of becoming "Now running annotale_jar predict for...").
