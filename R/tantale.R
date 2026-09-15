#'tantale: Transcription Activator-Like Effectors (TALEs) tools
#'
#'\figure{tantale_logo_small.gif}{options: width=100 alt="tantale_logo"}
#'
#'
#'@description An integrated collection of functions for (IDEALLY):
#'
#'Please take a look at the package \href{https://scunnac.github.io/tantale}{website}
#'for further details.
#'
#'@section   - A TALE-oriented OOP framework:
#'
#'  \itemize{
#'    \item A TALE class and associated methods}
#'
#'
#'@section   - TALE mining in bacterial sequences:
#'
#'  \itemize{
#'    \item Wrapper around annotale_jar and correcTALE
#'    \item tell_tales, an R function similar to annotale_jar
#'    \item Analysis tools for RVD inventory, repeat lenght}
#'
#'
#'@section   - TALEs classification, phylogeny:
#'
#'  \itemize{
#'    \item Wrappers around distal, functal, annotale_jar
#'    \item TALE groups inference
#'    \item Easily build Multiple alignments and generate nice plots}
#'
#'
#'@section   - TALE targets mining:
#'
#'  \itemize{
#'    \item Wrappers around target predictors
#'    \item General parser for results aggregation
#'    \item Connector with daTALbase (to be done)}
#'
#'@note CAUTIONARY NOTES:
#'
#'  \itemize{
#'    \item tantale has been written with only Linux systems in mind and will very
#'     likely \strong{not work on other OS} (eg Windows)
#'    \item Some of tantale wrappers use code written in other languages :
#'     \strong{Java and Perl must be on the PATH} in your system.
#'    \item Furthermore, Conda and Mamba must be installed.
#'    \item For direction on how to use Conda with R, consult the
#'     \href{https://rstudio.github.io/reticulate/reference/install_miniconda.html}{install_miniconda()} help page.}
#'
#'
#'@importFrom IRanges IRanges
#'@import fs
#'@import magrittr
#'@import cli
# Biostrings is imported wholesale rather than by name. Several calls in the
# package rely on its S4 methods for base-looking generics -- nchar() on an
# XStringSet is the trap, since it reads as base R and only differs for S4
# arguments. This import used to arrive as a side effect of an S4 class
# definition in telltale.R; it is declared deliberately here instead.
# Narrowing it to importFrom() is ledger 7.3, and wants the call sites
# qualified first.
#'@import Biostrings
#'@importFrom dplyr mutate if_else
#'@importFrom ggplot2 ggplot aes labs geom_point geom_text facet_grid theme_light
#'@importFrom ggplot2 scale_x_continuous scale_fill_discrete scale_color_viridis_d
#'@importFrom tidyr gather
#'@importFrom rlang .data
#'@importFrom stats hclust cutree dist as.dist as.dendrogram order.dendrogram median quantile
#'@importFrom utils read.table read.delim write.table
#'@importFrom grDevices colorRampPalette dev.off
#'@importFrom graphics abline axis box image layout mtext par strwidth text
#'@importFrom methods as new hasArg
"_PACKAGE"



