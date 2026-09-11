# tantale (development version)

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
