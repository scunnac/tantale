# Package index

## TALE discovery

Find TALE genes in genomic sequence and parse them into arrays of parts.

- [`correct_tales()`](https://scunnac.github.io/tantale/dev/reference/correct_tales.md)
  : Correct TALE ORFs in error-prone sequences
- [`tales_from_telltale()`](https://scunnac.github.io/tantale/dev/reference/tales_from_telltale.md)
  : Build a tales object from a tell_tales run directory
- [`tell_tales()`](https://scunnac.github.io/tantale/dev/reference/tell_tales.md)
  : Search and report on the features of TALE protein domains
  potentially encoded in subject DNA sequences

## tales objects

The central long table of TALE array parts, its constructor and
validators.

- [`as_tales()`](https://scunnac.github.io/tantale/dev/reference/as_tales.md)
  : Coerce sequences of TALE parts to a tales object
- [`format(`*`<tales>`*`)`](https://scunnac.github.io/tantale/dev/reference/format.tales.md)
  : Render a tales object as lines of text
- [`format(`*`<tales_msa>`*`)`](https://scunnac.github.io/tantale/dev/reference/format.tales_msa.md)
  : Render a tales_msa object as lines of text
- [`is_tales()`](https://scunnac.github.io/tantale/dev/reference/is_tales.md)
  : Is this a tales object?
- [`print(`*`<tales>`*`)`](https://scunnac.github.io/tantale/dev/reference/print.tales.md)
  : Print a tales object
- [`print(`*`<tales_msa>`*`)`](https://scunnac.github.io/tantale/dev/reference/print.tales_msa.md)
  : Print a tales_msa object
- [`summary(`*`<tales>`*`)`](https://scunnac.github.io/tantale/dev/reference/summary.tales.md)
  [`print(`*`<summary.tales>`*`)`](https://scunnac.github.io/tantale/dev/reference/summary.tales.md)
  : Summarise a tales object
- [`summary(`*`<tales_msa>`*`)`](https://scunnac.github.io/tantale/dev/reference/summary.tales_msa.md)
  [`print(`*`<summary.tales_msa>`*`)`](https://scunnac.github.io/tantale/dev/reference/summary.tales_msa.md)
  : Summarise a tales_msa object
- [`tales()`](https://scunnac.github.io/tantale/dev/reference/tales.md)
  : Create a tales object
- [`tales_anchor_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_anchor_codes.md)
  : Codes marking a TALE array terminus
- [`tales_anomalies()`](https://scunnac.github.io/tantale/dev/reference/tales_anomalies.md)
  : Report the biological anomalies in a tales object
- [`tales_assert_complete()`](https://scunnac.github.io/tantale/dev/reference/tales_assert_complete.md)
  : Assert that a tales object holds complete arrays
- [`tales_namespace()`](https://scunnac.github.io/tantale/dev/reference/tales_namespace.md)
  : The dom_code namespace of a tales object
- [`tales_requirements()`](https://scunnac.github.io/tantale/dev/reference/tales_requirements.md)
  : Column requirements of the tales consumers
- [`validate_tales()`](https://scunnac.github.io/tantale/dev/reference/validate_tales.md)
  : Validate a tales object

## Pairwise distances

Comparing TALEs and their domains, and the typed tables that result.

- [`as.matrix(`*`<pairwise_distances>`*`)`](https://scunnac.github.io/tantale/dev/reference/as.matrix.pairwise_distances.md)
  : Render a similarity table as a square matrix
- [`distances_assert_square()`](https://scunnac.github.io/tantale/dev/reference/distances_assert_square.md)
  : Assert that a similarity table is complete and square
- [`distances_restrict()`](https://scunnac.github.io/tantale/dev/reference/distances_restrict.md)
  : Restrict a similarity table to a set of entities
- [`is_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/is_pairwise_distances.md)
  : Is this a pairwise similarity table?
- [`pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md)
  [`tale_distances()`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md)
  [`domain_distances()`](https://scunnac.github.io/tantale/dev/reference/pairwise_distances.md)
  : Create a pairwise similarity table
- [`tales_assign_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_assign_domain_codes.md)
  : Assign a domain code to every distinct part sequence
- [`tales_compare()`](https://scunnac.github.io/tantale/dev/reference/tales_compare.md)
  : Compute TALE and repeat relatedness
- [`tales_domain_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_distances.md)
  : Pairwise distances between distinct TALE domains
- [`tales_group()`](https://scunnac.github.io/tantale/dev/reference/tales_group.md)
  : Group TALEs by similarity
- [`tales_tale_distances()`](https://scunnac.github.io/tantale/dev/reference/tales_tale_distances.md)
  : Pairwise distances between TALE arrays
- [`validate_pairwise_distances()`](https://scunnac.github.io/tantale/dev/reference/validate_pairwise_distances.md)
  : Validate a pairwise similarity table

## Alignment

Multiple alignment of TALE arrays and consensus over it.

- [`as.matrix(`*`<tales_msa>`*`)`](https://scunnac.github.io/tantale/dev/reference/as.matrix.tales_msa.md)
  : Render a TALE alignment as a matrix
- [`is_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/is_tales_msa.md)
  : Is this a tales_msa object?
- [`tales_align()`](https://scunnac.github.io/tantale/dev/reference/tales_align.md)
  : Align the repeat arrays of a tales object
- [`tales_consensus()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus.md)
  : Compute a consensus from a TALE msa
- [`tales_consensus_match()`](https://scunnac.github.io/tantale/dev/reference/tales_consensus_match.md)
  : Do elements in a TALE msa match the consensus?
- [`tales_msa()`](https://scunnac.github.io/tantale/dev/reference/tales_msa.md)
  : Create a tales_msa object
- [`tales_width()`](https://scunnac.github.io/tantale/dev/reference/tales_width.md)
  : Width of a TALE alignment
- [`validate_tales_msa()`](https://scunnac.github.io/tantale/dev/reference/validate_tales_msa.md)
  : Validate a tales_msa object

## Projections and conversions

Views derived from a tales object, and format conversions between them.

- [`repeat_to_rvd_map()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map.md)
  : Generate a mapping between Distal repeat IDs and their cognate RVD.
- [`repeat_to_rvd_map_distalr()`](https://scunnac.github.io/tantale/dev/reference/repeat_to_rvd_map_distalr.md)
  : Generate a mapping between Distal repeat IDs and their cognate RVD.
- [`tale_parts_to_rvd()`](https://scunnac.github.io/tantale/dev/reference/tale_parts_to_rvd.md)
  : Generates a RVD sequences set from a tale_parts object
- [`tales_coded_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_coded_strings.md)
  : Domain-coded strings, one per TALE array
- [`tales_domain_codes()`](https://scunnac.github.io/tantale/dev/reference/tales_domain_codes.md)
  : The domain code lookup table
- [`tales_rvd_strings()`](https://scunnac.github.io/tantale/dev/reference/tales_rvd_strings.md)
  : RVD strings, one per TALE array

## Plots

Visualising alignments, array composition and talome content.

- [`plot(`*`<tales>`*`)`](https://scunnac.github.io/tantale/dev/reference/plot.tales.md)
  : Plot the domain composition of a set of TALE arrays
- [`plot(`*`<tales_msa>`*`)`](https://scunnac.github.io/tantale/dev/reference/plot.tales_msa.md)
  : Plot a multiple alignment of TALEs
- [`plot_target_preds()`](https://scunnac.github.io/tantale/dev/reference/plot_target_preds.md)
  : Plot TALE RVD sequences along a potential DNA target region
- [`talomes_heatmap()`](https://scunnac.github.io/tantale/dev/reference/talomes_heatmap.md)
  : Heatmap plotting of rvd sequence variants

## Target prediction

Predicting EBEs in a promoter set, and plotting the result.

- [`plot_target_preds()`](https://scunnac.github.io/tantale/dev/reference/plot_target_preds.md)
  : Plot TALE RVD sequences along a potential DNA target region
- [`preditale()`](https://scunnac.github.io/tantale/dev/reference/preditale.md)
  : Run TALE target predictions on DNA sequence(s) using PrediTale.
- [`tales_predict_targets()`](https://scunnac.github.io/tantale/dev/reference/tales_predict_targets.md)
  : Predict TALE target boxes
- [`talvez()`](https://scunnac.github.io/tantale/dev/reference/talvez.md)
  : Run TALE target predictions on DNA sequence(s) using Talvez

## Setting up

Checking and building the external tools tantale drives.

- [`tantale_setup()`](https://scunnac.github.io/tantale/dev/reference/tantale_setup.md)
  : Check, and optionally build, tantale's external dependencies

## Package documentation

Package-level documentation.

- [`tantale`](https://scunnac.github.io/tantale/dev/reference/tantale-package.md)
  [`tantale-package`](https://scunnac.github.io/tantale/dev/reference/tantale-package.md)
  : tantale: Transcription Activator-Like Effectors (TALEs) tools

## External tool wrappers

Wrappers around AnnoTALE and QueTAL.

- [`functal()`](https://scunnac.github.io/tantale/dev/reference/functal.md)
  : Run functal from QueTAL to build a phylogenetic tree of TALE RVD
  sequences.
- [`run_annotale_build()`](https://scunnac.github.io/tantale/dev/reference/run_annotale_build.md)
  : Run the "build" stage of AnnoTALE.
- [`run_annotale_predict()`](https://scunnac.github.io/tantale/dev/reference/run_annotale_predict.md)
  : Runs the "predict" and "analyze" steps of AnnoTALE on a fasta file.
