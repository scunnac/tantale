#'tantale: Transcription Activator-Like Effectors (TALEs) tools
#'
#'\figure{tantale_logo.png}{options: width=100 alt="tantale_logo"}
#'
#'
#'@description An integrated collection of functions for (IDEALLY):
#'
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
#'    \item Wrapper around annotale
#'    \item tellTale
#'    \item Analysis tools for RVD inventory, repeat lenght}
#'
#'
#'@section   - TALEs classification, phylogeny:
#'
#'  \itemize{
#'    \item Wrappers around distal, functal, annoTALE and correcTALE
#'    \item TALE groups inference
#'    \item Easily build Multiple alignments and generate nice plots}
#'
#'
#'@section   - TALE targets mining:
#'
#'  \itemize{
#'    \item Wrappers around target predictors
#'    \item General parser for results aggregation
#'    \item Connector with daTALbase}
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
#'@docType package
#'@name tantale
#'
#'@importFrom IRanges IRanges
#'@import fs
#'@import magrittr
#'@import logger
NULL



