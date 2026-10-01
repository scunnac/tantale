# tantale 0.9.9011

## Terminus codes say whether a terminus resembles a TALE terminal domain

`tell_tales()` now searches the segment AnnoTALE reports on each side of
the repeats with the TALE N- and C-terminal protein profiles shipped in
`inst/extdata/hmmProfile/` (`hmmsearch`). `NTERM` and `CTERM` mark a
segment that matches its profile with an E-value at most
`terminus_max_evalue` (new argument, default 1e-5); `XXXXX` marks a
segment that does not match. The codes used to record whether an nhmmer
hit of that terminus type was found anywhere in the array's DNA. That
labelled as termini some segments of unrelated sequence, where the ORF
starts or ends in a frameshifted region, and coded as `XXXXX` genuine
termini too short for the DNA search, such as a C-terminus truncated to
42 residues.

`tales_from_telltale()` reads the codes from `array_report.tsv`. A
directory written by an earlier `tell_tales()` is an error; run
`tell_tales()` again on the same sequences.

## `array_report.tsv` and `tell_tales()`: renamed and new columns

`n_domain_hits` is now `n_dna_hits`, and `has_all_domains` is replaced by
`nterm_dna_hit` and `cterm_dna_hit`. New columns `nterm_aa_evalue`,
`cterm_aa_evalue`, `nterm_aa_hit` and `cterm_aa_hit` hold the result of
the protein-profile search. `?tell_tales` describes every column. The
argument `min_domain_hits` of `tell_tales()` is now `min_dna_hits`.

## `tales_from_telltale()`: absent termini and inconsistent arrays

When AnnoTALE reports no terminus on one side of the repeats, the array
now has no part on that side, with a warning, and its positions are
counted from its first part. It used to receive an empty part. An array
whose AnnoTALE protein and DNA parts disagree is left out with a warning,
and `tell_tales()` deletes the DNA parts AnnoTALE writes for an ORF it
could not translate. Repeat RVDs are read from AnnoTALE's own RVD file.

## `tell_tales()`: fewer spurious warnings, no index file next to the input

The warning about overlapping nhmmer hits now concerns only hits of the
same domain type, which make `n_dna_hits` count a domain twice; it
appears only with `merge_hits = FALSE`. A terminus hit overlapping the
adjacent repeat hit by a few nucleotides is normal and no longer
reported. Two warnings from the underlying Bioconductor packages are
gone: "invalid seqlevels ... ignored", issued when a subject sequence
carried no TALE, and "GRanges object contains ... out-of-bound ranges",
issued when an array lies within `extend_len` of a sequence end (the
extended range was, and still is, clipped to the sequence).

`tell_tales()` no longer writes a `.fai` index next to the subject file,
and Rsamtools is no longer a dependency. Sequence lengths now come from
the full fasta headers: with a header containing a space, they used to
be lost, so an array near the end of such a sequence could be extended
past it. A subject file with duplicated or empty sequence names is an
error.


## Licences of the bundled programs are now stated

`inst/COPYRIGHTS` lists every program and data file tantale bundles from
other projects, with its authors, licence, upstream download and source
code: AnnoTALE, PrediTALE and TALEcorrection (GNU GPL 3 or later, from
the Jstacs project), TALVEZ 3.2 and QueTAL FuncTAL (redistributed by
permission of their author). The GPL text ships as
`inst/tools/COPYING.GPL-3`. tantale's own code stays under the MIT
licence.

## `rvd_dna_specificity`: the RVD NA has its row back

The row for the RVD `NA` (Asn-Ala) had a missing name, because the table
was read with the string `"NA"` taken as a missing value. As a result,
`tales_to_universalmotif()` and `tales_compare_functal()` gave every `NA`
repeat the flat `XX` profile instead of its own (1/2/1/0). Fixed; no
other row changed.

## `rvd_sim` fill: a rare RVD identical to the reference now scores 1

`plot.tales_msa(fill_type = "rvd_sim")` takes its RVD similarities from
TALVEZ's table, which covers 17 RVDs. Any other RVD (`NV`, for example)
was left grey, even where it was identical to the reference's. It now
scores 1 there, as it already did in `tales_align()`'s RVD scoring
matrix; other pairs involving it stay grey. The legend says what grey
means.

## The interface is now stable

From version 1.0.0 on, changes to exported functions and their arguments
follow the [lifecycle](https://lifecycle.r-lib.org/articles/stages.html)
conventions: deprecation with a warning first, removal in a later release.

## `tales_to_universalmotif()` works without `library(tantale)`

`tales_to_universalmotif()` and `tales_compare_functal()` read the
`rvd_dna_specificity` table by a name that resolved only once tantale was
attached, so calling them as `tantale::tales_to_universalmotif()` failed
with "object 'rvd_dna_specificity' not found". Fixed.

## `talomes_heatmap()` fixes

With the default `plot_type = "all"`, the `title`, `x_lab` and `y_lab`
arguments were ignored (the default text was drawn whatever was passed),
and saving to `save_path` left an extra graphics device open. Both fixed;
the default figure is unchanged. An explicit `save_path = NULL` now draws
on the current device, and an unknown `plot_type` is an error instead of
drawing nothing.

## External programs: a failure now stops the call

nHMMER (in `tell_tales()` and `correct_tales()`), AnnoTALE
(`run_annotale_predict()`, `run_annotale_build()`), PrediTALE
(`preditale()`) and TALEcorrection (`correct_tales()`) were run without
checking their exit status, so a failing program could leave empty or
partial output for later steps to trip over. Each call now stops with an
error naming the program. Inside `tell_tales()`, an AnnoTALE failure on
one array still skips that array with a warning, which now gives
AnnoTALE's own error as its cause. `run_annotale_predict()` and
`run_annotale_build()` return `0` invisibly on success; their paths are
now quoted in the shell command, so paths with spaces work.

## Plots: no more ggplot2 deprecation warnings

`plot.tales_msa()` and `plot_target_preds()` used arguments ggplot2 has
deprecated (`label.size`, and `size` for a line width), so every plot
raised warnings. They now use `linewidth`, and the figures are unchanged.
tantale now requires ggplot2 3.5.0 or later.

## `array_report.tsv`: terminus lengths no longer count the stop codon

`nterm_aa_length` and `cterm_aa_length` counted AnnoTALE's `*` as a residue
whenever the stop codon fell inside the terminal part, which is common for
C-termini. They now count residues only, matching the `tales` object's
`aa_seq`: for example, a 278-residue C-terminus is reported as 278, not 279.

## Biostrings backend: free end gaps, as in DisTAL

`tales_domain_distances(aln_method = "Biostrings")` now aligns with free end
gaps (DisTAL's "sliding ends") instead of a global alignment that charged a
length difference twice. A half-repeat is now about 41 from the full repeat
it matches (it was 64.7), as with the other two backends. Pairs differing
only by substitutions are unchanged.

## Domain distances now count missing residues (DECIPHER backend)

`tales_domain_distances()`'s default `"DECIPHER"` backend ignored gaps, so
a half-repeat was at distance 0 from any full repeat it is a prefix of, and
indels between domains went uncounted. It now follows DisTAL's published
definition (percentage of amino acids that change, normalised by the longer
domain), as the `"mmseq2"` backend already did: a 20-residue half-repeat is
about 41 from the 34-residue repeat it matches. Domain and array distances
change for pairs that differ in length or align with gaps; on the three
genomes in the articles the grouping is unchanged. With a
`domain_distances` scoring matrix, `tales_align()` no longer pulls a
half-repeat away from its identical match.

## Reference pages reviewed

Every exported function's help page was re-read against the code. Among
the corrections: the `pairwise_distances` family is now described as the
distance table it is (it stores `dissim`), `repeat_to_rvd_map_distalr()`
and `tale_parts_to_rvd()` document the input they actually take,
`rvd_dna_specificity` explains its special rows (`N*`/`H*`, `OO`, `XX`),
`talomes_heatmap()` explains what its colours mean, and the package page
and README no longer claim that tantale bundles none of the programs it
drives (the Java tools ship with it).

## Array alignment computed in R; the ARLEM executable is no longer bundled

`tales_tale_distances()`, and so `tales_compare_distal()`, aligns the
coded arrays with an R implementation of ARLEM's minisatellite alignment
model (Abouelhoda et al., [2009](https://doi.org/10.1142/S0219720009004060))
in place of the ARLEM 1.0 executable that shipped in `inst/tools/arlem/`.
The scores are identical to the executable's: checked on about 6000 random
array pairs and on real TALE arrays, and recorded in a test fixture. The
executable ran only on Linux x86-64, and its licence did not clearly allow
redistribution; it has been removed from the package. This also clears the
"undeclared executable file" warning from `R CMD check`.

The R version takes about 4 ms per pair of 20-domain arrays, roughly five
times the executable's time; the domain-level alignment in
`tales_domain_distances()` remains the slow step. `matrixStats` is now
imported. Column names are unchanged (`arlem_score`, `max_length`).

## Bug fix: demoting a `tales_msa` now drops its alignment width

`as_tales()` (and `tales()`) on a `tales_msa` returned a plain `tales`
that still carried the alignment's width, so `tales_width()` kept
answering on an object that no longer claims to be an alignment. The
width is now removed on demotion, as it already was when
`alignment_position` is dropped with `select()`.

## Articles checked against their own output

Every article, the README and the home page were re-read against what
their code actually renders. Statements that disagreed with a table or a
figure are corrected, notably: the truncTALE article (`ROI_00001`'s
N-terminus is reduced like `ROI_00019`'s; `correct_tales()` leaves
`ROI_00001`'s C-terminus length unchanged), the alignment article (a
`domain_distances` scoring matrix does not make this group's alignment
more compact; it moves the final half-repeat away from its identical
match), and the home-page dendrogram description (it shows nine groups). The
classification article's motif tree now has tip labels, and its
silhouette plot is described and drawn on every build.

## New article: genuine truncTALEs and frameshift correction

TALEs are not always broken when they are short. `PXO86` carries two
naturally truncated TAL effectors (truncTALEs), documented in
[Ji et al. 2016](https://doi.org/10.1038/ncomms13435) and
[Read et al. 2016](https://doi.org/10.3389/fpls.2016.01516), and the two
differ enough at the DNA level that they respond differently to
automatic frameshift correction: `tell_tales(correct_array = TRUE)`
extends one of them as if it were an assembly error, while
`correct_tales()` leaves both alone. See the new
[Genuine truncTALEs and frameshift correction](articles/trunctale_correction.html)
article for the full comparison.

## `tale_mining.qmd`'s correction sections updated for the `correct_tales()` fix

The "Correcting frameshifts, two ways" and "How much did either
correction actually help?" sections were written against the pre-fix,
swapped-flags `correct_tales()` (see the bug fix entry below) and
described it as fixing only one of BAI3-1-1's two frameshifted arrays.
With the flags corrected, `correct_tales()` fixes both on its own; the
article's prose, timing table and coverage figure are updated to match,
verified against a real re-run rather than assumed.

## `plot.tales_msa()`: internal column names corrected, consensus computation simplified

Some of `plot.tales_msa()`'s internal column names claimed a
repeat-specificity the underlying computation never had (termini are
scored identically to repeats) -- `matchConsensusRepeat`, `repeatClusterId`
and `repeatSimVsRef` are now `matchConsensusDomain`, `domainClusterId` and
`domainSimVsRef`. This only affects code reading `plot(msa)`'s returned
`$data`/`$plotlist[[i]]$data` by column name to compose further layers, per
its own documented pattern -- the figure itself, and every other argument
and return value, are unchanged. The consensus/match computation behind
the alignment's text colour no longer round-trips through a matrix
internally; not user-visible, part of an ongoing internal simplification
(`dev/restructuring-notes.md` §21, not finished this round).

## Bug fix: `correct_tales()` had its nHMMER inputs swapped

`correct_tales()` built its call to `TALEcorrection.jar` with the repeat-
and C-terminus-domain nHMMER search results assigned to the wrong CLI
flags -- confirmed against the tool's own printed usage, which documents
`r=` as wanting the repeats file and `c=` the C-terminus file, exactly
backwards from what every prior call supplied. Every correction run has
therefore been telling the external tool the wrong domain for those two
inputs. Fixed; corrected sequences from `correct_tales()` may now differ
from previous runs on the same input.

## In-depth documentation review of every exported function

Every exported function and S3 method's documentation was checked for
accuracy against its actual current behaviour, completeness (arguments,
return value, a real runnable example), and consistency of tone --
several stale or copy-pasted descriptions, wrong defaults, and missing
`@examples` fixed along the way. See `dev/restructuring-notes.md` §17
for the full record.

## TALE comparison by predicted binding specificity, reworked

`functal()`, the wrapper around QueTAL's vendored Perl `FuncTAL` script, is
retired: moved to `inst/legacy/` along with the Perl tool itself
(`inst/tools/QueTAL_v1.1/` -> `inst/legacy/QueTAL_v1.1/`), unexported, no
longer built. It could not be fixed in place -- the script needs
`Bio::Perl`, dropped by BioPerl in the 1.7 reorganisation, and no
perl/bioperl combination available restores it.

In its place: `tales_to_universalmotif()` converts a `tales` object into a
list of `universalmotif` PWMs (one per array, from each array's RVD-to-DNA
binding-specificity weights), and `tales_compare_functal()` compares them
via `universalmotif::compare_motifs()`. This is not a port of the old
script's numbers -- `compare_motifs()`'s Pearson correlation over matched
columns is a different statistic from FuncTAL's single correlation over
the whole flattened, padded region, and the two do not agree. Once a
`tales` object is a list of real `universalmotif` motifs, the rest of that
package's toolkit becomes directly usable: `motif_tree()`, `view_motifs()`,
`scan_sequences()`, `merge_motifs()`, and more.

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

### `msa_heatmap()` is retired

Superseded by `plot_tales_msa()` and moved to `inst/legacy/`. Its six
`plot_type` values were combinations of features that `plot_tales_msa()`
exposes as orthogonal arguments; the one thing it did that `plot_tales_msa()
could not -- draw the consensus row -- is now implemented there.

| `msa_heatmap(plot_type =)` | `plot_tales_msa()` |
|---|---|
| `"repeat.similarity"` | `fill_type = "repeat_sim"` |
| `"repeat.clusters"` | `fill_type = "repeat_clust"` |
| `"with.rvd"` | pass `rvd_align` |
| `"repeat.clusters.with.rvd"` | both of those |
| `"reference"` | pass `ref_pattern` |
| `"consensus"` | `consensus = TRUE` |

`save_path` has no equivalent because none is needed: the returned ggplot is
`ggsave()`-able. `note_colors` likewise -- add a scale to the returned plot.

### `plot_tales_msa()` can draw the consensus

`consensus = TRUE` now adds a consensus row above the alignment, where it
previously did nothing and the argument was documented as "NOT IMPLEMENTED
YET". The consensus is the most frequent element per column, taken from
`rvd_align` when supplied and `repeat_align` otherwise, so it always matches
what the cells are labelled with.

It is drawn as its own panel above the alignment rather than as an extra row.
That is not cosmetic: `aplot` reorders the alignment's y axis onto the tree's
leaves, so a row the tree has no leaf for is silently dropped.

This was the one feature `msa_heatmap()` still provided that
`plot_tales_msa()` did not.

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

## Breaking change: every column is now `snake_case`

The tables `tantale` returns used four naming conventions at once, one of them
with a literal space in a column name. They now use one.

### The `tales` table

`arrayID`, `domainType`, `positionInCrd`, `dnaSeq`, `sourceDirectory`,
`positionInArray`, `aaSeq` and `domCode` become `array_id`, `domain_type`,
`position_in_crd`, `dna_seq`, `source_directory`, `position_in_array`,
`aa_seq` and `dom_code`.

`tales()` still accepts the old spellings and renames them, so a `tell_tales()
` output directory written by an earlier version still loads.

### The distance tables

`domain_distances` and `tale_distances` previously disagreed about how to name
the pair being compared -- `RepU1`/`RepU2` in one, `TAL1`/`TAL2` in the other,
and in opposite column orders. Both now use `id1`, `id2` and `dissim`, with
`arlem_score` and `max_length` alongside for the TALE table.

`pairwise_distances()` accepts `TAL1`/`RepU1`/`Sim`/`Dissim`/`arlemScore` and
renames them, and `plot_tales_msa()` puts its `tal_sim` and `repeat_sim`
arguments through it, so passing an old-format table still works.

### On-disk output

`tell_tales()` writes `array_id` instead of `arrayID` in `arrayReport.tsv`,
`domainsReport.tsv` and `hitsReport.tsv`, and the GFFs derived from them carry
an `array_id` attribute. **Scripts that read these files by column name need
updating.**

### Fixed along the way

- `.tale_parts_from_file()` named its column `arrayIDs` when the input was
  empty and `arrayID` otherwise, so the two returns had incompatible schemas.
- `repeat_to_rvd_map_distalr()` and `tale_parts_to_rvd()` are documented as
  taking a `tales_compare()` result but read the pre-class column names, so
  both had been broken against that result.
- `.repeat_to_cluster_align()` recovered a distance by computing `100 - Sim`
  before clustering, which assumed the scale ran 0-100. It now reads the
  stored distance.
- Six `if (class(x) == "...")` comparisons became `inherits()`.

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
