#'tantale: Transcription Activator-Like Effectors (TALEs) tools
#'
#'\if{html}{\figure{tantale_logo_small.gif}{options: width=100 alt="tantale_logo"}}
#'
#'
#'@description Tools to find the TALE genes of \emph{Xanthomonas} in DNA
#'sequences, compare and group them, align their repeat arrays and predict
#'the plant promoter sites they bind. The
#'\href{https://scunnac.github.io/tantale/}{package website} has worked
#'examples of each step.
#'
#'@section Finding TALEs:
#'  [tell_tales()] finds TALE genes with HMMER profiles of the TALE domains,
#'  much as AnnoTALE does, and is written for noisy sequences such as draft
#'  assemblies or long reads: it can correct frameshifts against reference
#'  TALEs before AnnoTALE reads out the RVDs. [run_annotale_predict()] runs
#'  AnnoTALE itself, and [correct_tales()] runs TALEcorrection.
#'  [tales_from_telltales()] and [tales_from_annotale()] load either result
#'  as a [tales] object.
#'
#'@section Working with TALEs in R:
#'  A [tales] object holds one row per part (N-terminus, repeat or
#'  C-terminus) of each TALE, and [tales_anomalies()] reports the TALEs
#'  whose structure is not standard. A [tales_msa] adds the alignment of the
#'  arrays. [tale_distances] and [domain_distances] objects hold distances
#'  between TALEs or between their domains.
#'
#'@section Comparing and grouping TALEs:
#'  [tales_compare_distal()] and [tales_compare_functal()] reimplement the
#'  DisTAL and FuncTAL comparisons of QueTAL. [tales_group_kmedoids()] and
#'  [tales_group_hclust()] group TALEs from these distances, and
#'  [talomes_heatmap()] compares the groups across strains. AnnoTALE's own
#'  classes come from [run_annotale_build()], and its published catalogue
#'  from [run_annotale_load_classes()] and [run_annotale_assign()]. The
#'  [tale_annotations] dataset gives 128 curated TALEs of ten
#'  \emph{X. oryzae} genomes to compare against.
#'
#'@section Aligning TALEs:
#'  [tales_align()] aligns the RVD or domain sequences of the arrays with
#'  MAFFT; [tales_consensus()] and the \code{plot()} method summarise the
#'  alignment.
#'
#'@section Predicting targets:
#'  [tales_predict_targets()] runs Talvez ([talvez()]) or PrediTALE
#'  ([preditale()]) on DNA sequences such as promoters, and
#'  [plot_target_preds()] draws a predicted site with the RVDs facing their
#'  bases.
#'
#'@section Setting up:
#'  MAFFT, HMMER, mmseqs2 and the Perl that Talvez runs on come from a
#'  conda environment the package builds for itself, so conda (or mamba, or
#'  micromamba) is needed for most of the package. [tantale_setup()] called
#'  bare reports what is present and changes nothing;
#'  \code{tantale_setup(install = TRUE)} builds or repairs the environment
#'  and downloads the Java tools (AnnoTALE, PrediTALE, TALEcorrection) and
#'  the example genomes. It checks tool versions as well as presence: MAFFT
#'  and HMMER are pinned because later MAFFT versions align TALE repeat
#'  strings differently.
#'
#'  Without conda, [reticulate::install_miniconda()] installs one from
#'  inside R. An existing conda, mamba or micromamba is found through
#'  \code{reticulate::conda_binary()}.
#'
#'@note tantale runs on Linux and on Intel macOS, the platforms the pinned
#'  conda tools exist for; it does not run on Windows. The Java tools need
#'  Java on the PATH, which [tantale_setup()] checks.
#'
#'
#'@importFrom IRanges IRanges
#'@importFrom magrittr %>%
# Biostrings: the package calls these base-named functions unqualified, and
# they resolve to Biostrings' S4 generics -- nchar() on an XStringSet is the
# trap, since it reads as base R and only differs for S4 arguments. The list
# is every Biostrings export the package's code uses unqualified
# (codetools::findGlobals(), 2026-10-04), so each call resolves exactly as it
# did under the old import(Biostrings) (ledger 7.3, §52 F4). A new
# unqualified call to such a generic must be added here or qualified.
#'@importFrom Biostrings %in% as.data.frame as.list as.matrix duplicated
#'@importFrom Biostrings end intersect match nchar order rank setdiff
#'@importFrom Biostrings setequal sort start strsplit substr summary union
#'@importFrom dplyr mutate if_else
#'@importFrom ggplot2 ggplot aes labs geom_point geom_text facet_grid theme_light
#'@importFrom ggplot2 scale_x_continuous
#'@importFrom tidyr gather
#'@importFrom rlang .data
#'@importFrom stats hclust cutree dist as.dist as.dendrogram order.dendrogram median quantile
#'@importFrom utils read.table read.delim write.table
#'@importFrom grDevices colorRampPalette dev.off
#'@importFrom graphics abline axis box image layout mtext par strwidth text
#'@importFrom methods as new hasArg
"_PACKAGE"



