# tantale: Transcription Activator-Like Effectors (TALEs) tools

Tools to find the TALE genes of *Xanthomonas* in DNA sequences, compare
and group them, align their repeat arrays and predict the plant promoter
sites they bind. The [package
website](https://scunnac.github.io/tantale/) has worked examples of each
step.

## Details

![tantale_logo](figures/tantale_logo_small.gif)

## Note

tantale runs on Linux and on Intel macOS, the platforms the pinned conda
tools exist for; it does not run on Windows. The Java tools need Java on
the PATH, which
[`tantale_setup()`](https://scunnac.github.io/tantale/reference/tantale_setup.md)
checks.

## Finding TALEs

[`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md)
finds TALE genes with HMMER profiles of the TALE domains, much as
AnnoTALE does, and is written for noisy sequences such as draft
assemblies or long reads: it can correct frameshifts against reference
TALEs before AnnoTALE reads out the RVDs.
[`run_annotale_predict()`](https://scunnac.github.io/tantale/reference/run_annotale_predict.md)
runs AnnoTALE itself, and
[`correct_tales()`](https://scunnac.github.io/tantale/reference/correct_tales.md)
runs TALEcorrection.
[`tales_from_telltales()`](https://scunnac.github.io/tantale/reference/tales_from_telltales.md)
and
[`tales_from_annotale()`](https://scunnac.github.io/tantale/reference/tales_from_annotale.md)
load either result as a
[tales](https://scunnac.github.io/tantale/reference/tales.md) object.

## Working with TALEs in R

A [tales](https://scunnac.github.io/tantale/reference/tales.md) object
holds one row per part (N-terminus, repeat or C-terminus) of each TALE,
and
[`tales_anomalies()`](https://scunnac.github.io/tantale/reference/tales_anomalies.md)
reports the TALEs whose structure is not standard. A
[tales_msa](https://scunnac.github.io/tantale/reference/tales_msa.md)
adds the alignment of the arrays.
[tale_distances](https://scunnac.github.io/tantale/reference/pairwise_distances.md)
and
[domain_distances](https://scunnac.github.io/tantale/reference/pairwise_distances.md)
objects hold distances between TALEs or between their domains.

## Comparing and grouping TALEs

[`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md)
and
[`tales_compare_functal()`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md)
reimplement the DisTAL and FuncTAL comparisons of QueTAL.
[`tales_group_kmedoids()`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md)
and
[`tales_group_hclust()`](https://scunnac.github.io/tantale/reference/tales_group_hclust.md)
group TALEs from these distances, and
[`talomes_heatmap()`](https://scunnac.github.io/tantale/reference/talomes_heatmap.md)
compares the groups across strains. AnnoTALE's own classes come from
[`run_annotale_build()`](https://scunnac.github.io/tantale/reference/run_annotale_build.md),
and its published catalogue from
[`run_annotale_load_classes()`](https://scunnac.github.io/tantale/reference/run_annotale_load_classes.md)
and
[`run_annotale_assign()`](https://scunnac.github.io/tantale/reference/run_annotale_assign.md).
The
[tale_annotations](https://scunnac.github.io/tantale/reference/tale_annotations.md)
dataset gives 128 curated TALEs of ten *X. oryzae* genomes to compare
against.

## Aligning TALEs

[`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md)
aligns the RVD or domain sequences of the arrays with MAFFT;
[`tales_consensus()`](https://scunnac.github.io/tantale/reference/tales_consensus.md)
and the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
summarise the alignment.

## Predicting targets

[`tales_predict_targets()`](https://scunnac.github.io/tantale/reference/tales_predict_targets.md)
runs Talvez
([`talvez()`](https://scunnac.github.io/tantale/reference/talvez.md)) or
PrediTALE
([`preditale()`](https://scunnac.github.io/tantale/reference/preditale.md))
on DNA sequences such as promoters, and
[`plot_target_preds()`](https://scunnac.github.io/tantale/reference/plot_target_preds.md)
draws a predicted site with the RVDs facing their bases.

## Setting up

MAFFT, HMMER, mmseqs2 and the Perl that Talvez runs on come from a conda
environment the package builds for itself, so conda (or mamba, or
micromamba) is needed for most of the package.
[`tantale_setup()`](https://scunnac.github.io/tantale/reference/tantale_setup.md)
called bare reports what is present and changes nothing;
`tantale_setup(install = TRUE)` builds or repairs the environment and
downloads the Java tools (AnnoTALE, PrediTALE, TALEcorrection) and the
example genomes. It checks tool versions as well as presence: MAFFT and
HMMER are pinned because later MAFFT versions align TALE repeat strings
differently.

Without conda,
[`reticulate::install_miniconda()`](https://rstudio.github.io/reticulate/reference/install_miniconda.html)
installs one from inside R. An existing conda, mamba or micromamba is
found through
[`reticulate::conda_binary()`](https://rstudio.github.io/reticulate/reference/conda-tools.html).

## See also

Useful links:

- <https://scunnac.github.io/tantale/>

- <https://github.com/scunnac/tantale>

- Report bugs at <https://github.com/scunnac/tantale/issues>

## Author

**Maintainer**: Sebastien Cunnac <sebastien.cunnac@ird.fr>
([ORCID](https://orcid.org/0000-0002-3695-491X))

Authors:

- Sebastien Cunnac <sebastien.cunnac@ird.fr>
  ([ORCID](https://orcid.org/0000-0002-3695-491X))

- Bao Tram Vi <vbt576@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-4319-5544))

Other contributors:

- Alvaro L. Pérez-Quintero (Author of the bundled TALVEZ 3.2 and of the
  FuncTAL table behind rvd_dna_specificity) \[copyright holder\]

- Molly Megraw (Author of the PlantTFBS Java classes bundled with TALVEZ
  (simplescancode/)) \[copyright holder\]

- Artemis G. Hatzigeorgiou (Author of the PlantTFBS Java classes bundled
  with TALVEZ (simplescancode/)) \[copyright holder\]
