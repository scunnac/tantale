# Package index

## Overview

Which object each exported function takes and returns. Boxes are object
classes, dots are functions. Hover over a node to see what it takes and
returns, click it to highlight its neighbours, and scroll to zoom. [Open
the map in a full
window](https://scunnac.github.io/tantale/function_map.md).

## TALE discovery

Find TALE genes in genomic sequence and parse them into arrays of parts.

- [`correct_tales()`](https://scunnac.github.io/tantale/reference/correct_tales.md)
  : Correct TALE ORFs in error-prone sequences

- [`tale_annotations`](https://scunnac.github.io/tantale/reference/tale_annotations.md)
  :

  Curated TALE annotations for ten published *Xanthomonas oryzae*
  genomes

- [`tales_from_annotale()`](https://scunnac.github.io/tantale/reference/tales_from_annotale.md)
  : Build a tales object from AnnoTALE's own TALE predictions

- [`tales_from_telltales()`](https://scunnac.github.io/tantale/reference/tales_from_telltales.md)
  : Build a tales object from a tell_tales run directory

- [`tell_tales()`](https://scunnac.github.io/tantale/reference/tell_tales.md)
  : Search and report on the features of TALE protein domains
  potentially encoded in subject DNA sequences

## tales objects

The central long table of TALE array parts, its constructor and
validators.

- [`as_tales()`](https://scunnac.github.io/tantale/reference/as_tales.md)
  : Coerce sequences of TALE parts to a tales object
- [`format(`*`<tales>`*`)`](https://scunnac.github.io/tantale/reference/format.tales.md)
  : Render a tales object as lines of text
- [`format(`*`<tales_msa>`*`)`](https://scunnac.github.io/tantale/reference/format.tales_msa.md)
  : Render a tales_msa object as lines of text
- [`is_tales()`](https://scunnac.github.io/tantale/reference/is_tales.md)
  : Is this a tales object?
- [`print(`*`<tales>`*`)`](https://scunnac.github.io/tantale/reference/print.tales.md)
  : Print a tales object
- [`print(`*`<tales_msa>`*`)`](https://scunnac.github.io/tantale/reference/print.tales_msa.md)
  : Print a tales_msa object
- [`` `[`( ``*`<tales>`*`)`](https://scunnac.github.io/tantale/reference/sub-.tales.md)
  : Subset a tales object
- [`summary(`*`<tales>`*`)`](https://scunnac.github.io/tantale/reference/summary.tales.md)
  [`print(`*`<summary.tales>`*`)`](https://scunnac.github.io/tantale/reference/summary.tales.md)
  : Summarise a tales object
- [`summary(`*`<tales_msa>`*`)`](https://scunnac.github.io/tantale/reference/summary.tales_msa.md)
  [`print(`*`<summary.tales_msa>`*`)`](https://scunnac.github.io/tantale/reference/summary.tales_msa.md)
  : Summarise a tales_msa object
- [`tales()`](https://scunnac.github.io/tantale/reference/tales.md) :
  Create a tales object
- [`tales_anchor_codes()`](https://scunnac.github.io/tantale/reference/tales_anchor_codes.md)
  : Codes marking a TALE array terminus
- [`tales_anomalies()`](https://scunnac.github.io/tantale/reference/tales_anomalies.md)
  : Report the biological anomalies in a tales object
- [`tales_assert_complete()`](https://scunnac.github.io/tantale/reference/tales_assert_complete.md)
  : Assert that a tales object holds complete arrays
- [`tales_bind()`](https://scunnac.github.io/tantale/reference/tales_bind.md)
  : Combine tales objects
- [`tales_names()`](https://scunnac.github.io/tantale/reference/tales_names.md)
  : The names of the TALEs in a tales object
- [`tales_namespace()`](https://scunnac.github.io/tantale/reference/tales_namespace.md)
  : The dom_code namespace of a tales object
- [`tales_requirements()`](https://scunnac.github.io/tantale/reference/tales_requirements.md)
  : Column requirements of the tales consumers
- [`validate_tales()`](https://scunnac.github.io/tantale/reference/validate_tales.md)
  : Validate a tales object

## Pairwise distances

Comparing TALEs and their domains, and the typed tables that result.

- [`as.matrix(`*`<pairwise_distances>`*`)`](https://scunnac.github.io/tantale/reference/as.matrix.pairwise_distances.md)
  : Render a distance table as a square matrix
- [`distances_assert_square()`](https://scunnac.github.io/tantale/reference/distances_assert_square.md)
  : Assert that a distance table is complete and square
- [`distances_restrict()`](https://scunnac.github.io/tantale/reference/distances_restrict.md)
  : Restrict a distance table to a set of entities
- [`is_pairwise_distances()`](https://scunnac.github.io/tantale/reference/is_pairwise_distances.md)
  : Is this a pairwise distance table?
- [`pairwise_distances()`](https://scunnac.github.io/tantale/reference/pairwise_distances.md)
  [`tale_distances()`](https://scunnac.github.io/tantale/reference/pairwise_distances.md)
  [`domain_distances()`](https://scunnac.github.io/tantale/reference/pairwise_distances.md)
  : Create a pairwise distance table
- [`rvd_dna_specificity`](https://scunnac.github.io/tantale/reference/rvd_dna_specificity.md)
  : RVD-to-DNA-binding-specificity weights
- [`` `[`( ``*`<pairwise_distances>`*`)`](https://scunnac.github.io/tantale/reference/sub-.pairwise_distances.md)
  : Subset a pairwise_distances object
- [`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_assign_domain_codes.md)
  : Assign a domain code to every distinct part sequence
- [`tales_compare_distal()`](https://scunnac.github.io/tantale/reference/tales_compare_distal.md)
  : Compute TALE and domain relatedness by domain-sequence alignment
  (DisTAL)
- [`tales_compare_functal()`](https://scunnac.github.io/tantale/reference/tales_compare_functal.md)
  : Compare TALEs by predicted DNA-binding specificity (FuncTAL)
- [`tales_domain_distances()`](https://scunnac.github.io/tantale/reference/tales_domain_distances.md)
  : Pairwise distances between distinct TALE domains
- [`tales_group_hclust()`](https://scunnac.github.io/tantale/reference/tales_group_hclust.md)
  : Group TALEs by hierarchical clustering of their pairwise distance
- [`tales_group_kmedoids()`](https://scunnac.github.io/tantale/reference/tales_group_kmedoids.md)
  : Group TALEs by k-medoids clustering of their pairwise distance
- [`tales_tale_distances()`](https://scunnac.github.io/tantale/reference/tales_tale_distances.md)
  : Pairwise distances between TALE arrays
- [`validate_pairwise_distances()`](https://scunnac.github.io/tantale/reference/validate_pairwise_distances.md)
  : Validate a pairwise distance table

## Alignment

Multiple alignment of TALE arrays and consensus over it.

- [`as.matrix(`*`<tales_msa>`*`)`](https://scunnac.github.io/tantale/reference/as.matrix.tales_msa.md)
  : Render a TALE alignment as a matrix
- [`is_tales_msa()`](https://scunnac.github.io/tantale/reference/is_tales_msa.md)
  : Is this a tales_msa object?
- [`tales_align()`](https://scunnac.github.io/tantale/reference/tales_align.md)
  : Align the repeat arrays of a tales object
- [`tales_consensus()`](https://scunnac.github.io/tantale/reference/tales_consensus.md)
  : Compute a consensus from a TALE msa
- [`tales_consensus_match()`](https://scunnac.github.io/tantale/reference/tales_consensus_match.md)
  : Do elements in a TALE msa match the consensus?
- [`tales_msa()`](https://scunnac.github.io/tantale/reference/tales_msa.md)
  : Create a tales_msa object
- [`tales_msa_width()`](https://scunnac.github.io/tantale/reference/tales_msa_width.md)
  : Width of a TALE alignment
- [`validate_tales_msa()`](https://scunnac.github.io/tantale/reference/validate_tales_msa.md)
  : Validate a tales_msa object

## Projections and conversions

Views derived from a tales object, and format conversions between them.

- [`tales_coded_strings()`](https://scunnac.github.io/tantale/reference/tales_coded_strings.md)
  : Domain-coded strings, one per TALE array
- [`tales_domain_codes()`](https://scunnac.github.io/tantale/reference/tales_domain_codes.md)
  : The domain code lookup table
- [`tales_get_dna_seq()`](https://scunnac.github.io/tantale/reference/tales_get_dna_seq.md)
  : Whole-array DNA sequence, one per TALE array
- [`tales_get_protein_seq()`](https://scunnac.github.io/tantale/reference/tales_get_protein_seq.md)
  : Whole-array protein sequence, one per TALE array
- [`tales_rvd_strings()`](https://scunnac.github.io/tantale/reference/tales_rvd_strings.md)
  : RVD strings, one per TALE array
- [`tales_to_universalmotif()`](https://scunnac.github.io/tantale/reference/tales_to_universalmotif.md)
  : Convert a tales object to a list of predicted
  DNA-binding-specificity PWMs

## Plots

Visualising alignments, array composition and talome content.

- [`plot(`*`<tales>`*`)`](https://scunnac.github.io/tantale/reference/plot.tales.md)
  : Plot the domain composition of a set of TALE arrays
- [`plot(`*`<tales_msa>`*`)`](https://scunnac.github.io/tantale/reference/plot.tales_msa.md)
  : Plot a multiple alignment of TALEs
- [`plot_target_preds()`](https://scunnac.github.io/tantale/reference/plot_target_preds.md)
  : Plot TALE RVD sequences along a potential DNA target region
- [`talomes_heatmap()`](https://scunnac.github.io/tantale/reference/talomes_heatmap.md)
  : Heatmap of RVD sequence variants across strains and TALE groups

## Target prediction

Predicting EBEs in a promoter set, and plotting the result.

- [`plot_target_preds()`](https://scunnac.github.io/tantale/reference/plot_target_preds.md)
  : Plot TALE RVD sequences along a potential DNA target region
- [`preditale()`](https://scunnac.github.io/tantale/reference/preditale.md)
  : Run TALE target predictions on DNA sequence(s) using PrediTale
- [`tales_predict_targets()`](https://scunnac.github.io/tantale/reference/tales_predict_targets.md)
  : Predict TALE target boxes
- [`talvez()`](https://scunnac.github.io/tantale/reference/talvez.md) :
  Run TALE target predictions on DNA sequence(s) using Talvez

## Setting up

Checking, building and downloading the external tools tantale drives,
and the example genomes.

- [`tantale_setup()`](https://scunnac.github.io/tantale/reference/tantale_setup.md)
  : Check, and optionally build, tantale's external dependencies
- [`tantale_genome()`](https://scunnac.github.io/tantale/reference/tantale_genome.md)
  : Path to one of the example genomes

## Package documentation

Package-level documentation.

- [`tantale`](https://scunnac.github.io/tantale/reference/tantale-package.md)
  [`tantale-package`](https://scunnac.github.io/tantale/reference/tantale-package.md)
  : tantale: Transcription Activator-Like Effectors (TALEs) tools

## External tool wrappers

Wrappers around AnnoTALE programs.

- [`run_annotale_assign()`](https://scunnac.github.io/tantale/reference/run_annotale_assign.md)
  : Assign TALEs to AnnoTALE's published classes
- [`run_annotale_build()`](https://scunnac.github.io/tantale/reference/run_annotale_build.md)
  : Run the "build" stage of AnnoTALE
- [`run_annotale_load_classes()`](https://scunnac.github.io/tantale/reference/run_annotale_load_classes.md)
  : Download AnnoTALE's catalogue of TALE classes
- [`run_annotale_predict()`](https://scunnac.github.io/tantale/reference/run_annotale_predict.md)
  : Runs the "predict" and "analyze" steps of AnnoTALE on a fasta file
