#'tantale: Transcription Activator-Like Effectors (TALEs) tools
#'
#'\if{html}{\figure{tantale_logo_small.gif}{options: width=100 alt="tantale_logo"}}
#'
#'
#'@description An integrated collection of functions for TALE mining and
#'analysis in R.
#'
#'Please take a look at the package \href{https://scunnac.github.io/tantale/}{website}
#'for further details.
#'
#'@section A TALE-oriented OOP framework:
#'
#'  \itemize{
#'    \item \code{tales}/\code{tales_msa} S3 classes, with subsetting,
#'    coercion, and plotting methods}
#'
#'
#'@section TALE mining in bacterial sequences:
#'
#'  \itemize{
#'    \item Wrapper around annotale_jar and correcTALE
#'    \item tell_tales, an R function similar to annotale_jar
#'    \item Analysis tools for RVD inventory, repeat length}
#'
#'
#'@section TALEs classification, phylogeny:
#'
#'  \itemize{
#'    \item R reimplementations of distal and functal comparisons, plus a
#'    wrapper around annotale_jar
#'    \item TALE groups inference
#'    \item Easily build Multiple alignments and generate nice plots}
#'
#'
#'@section TALE targets mining:
#'
#'  \itemize{
#'    \item Wrappers around target predictors
#'    \item General parser for results aggregation
#'    \item Connector with daTALbase (to be done)}
#'
#'@section Setting up:
#'
#'  MAFFT, HMMER, mmseqs2 and the Perl dependencies of the target predictors
#'  come from a conda environment the package builds for itself, so **conda
#'  (or mamba, or micromamba) is a prerequisite of the main workflow**. The
#'  Java tools (AnnoTALE, PrediTALE, TALE correction), which have no conda
#'  package, ship with tantale itself. The environment is built on first
#'  use, so the first call needs a network connection and takes a few
#'  minutes.
#'
#'  Start with [tantale_setup()]. Called bare it reports what is present and
#'  changes nothing; `tantale_setup(install = TRUE)` builds or repairs the
#'  environment. It checks **versions** as well as presence, which matters
#'  because MAFFT changed its `--text` mode gap handling after 7.4x and later
#'  versions align TALE repeat strings differently -- an environment left
#'  over from an older tantale gives different alignments from the same
#'  input, and nothing else would report it.
#'
#'  If you have no conda at all, `reticulate` will install one from inside R
#'  with [reticulate::install_miniconda()] (this installs miniconda). An
#'  existing conda, mamba or micromamba is found automatically through
#'  `reticulate::conda_binary()` and used instead.
#'
#'@note CAUTIONARY NOTES:
#'
#'  \itemize{
#'    \item tantale has been written with only Linux systems in mind and will very
#'     likely \strong{not work on other OS} (eg Windows)
#'    \item Some of tantale wrappers use code written in other languages :
#'     \strong{Java and Perl must be on the PATH} in your system.
#'     [tantale_setup()] checks for both.}
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



