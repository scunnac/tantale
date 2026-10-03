#### The package's colour style ####
#
# Every plot draws its colours from here, so the figures read as one set and
# stay legible for colour-blind readers (ledger §46). The colours are Paul
# Tol's schemes (https://personal.sron.nl/~pault/), copied as hex codes:
# "muted" for categories, the light half of "YlOrBr" for sequential scales,
# the light middle of "sunset" for the diverging one, and "light" for a
# second set of categories in the same figure.
# The DNA bases keep ColorBrewer's RdYlBu, and two single-hue ramps (wine
# for ranked variants, purple for match quality) were chosen by eye
# against the alternatives.
#
# Fills behind text are kept pale wherever the text colour carries its own
# meaning (the consensus match in plot.tales_msa()), so dark text reads on
# all of them. Elsewhere the text colour follows the fill
# (.text_colour_on()).

.tol_muted <- c(rose = "#CC6677", indigo = "#332288", sand = "#DDCC77",
                green = "#117733", cyan = "#88CCEE", wine = "#882255",
                teal = "#44AA99", olive = "#999933", purple = "#AA4499")

.tol_light <- c("#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF",
                "#44BB99", "#BBCC33", "#AAAA00")

.tantale_colours <- list(
  no_value = "#DDDDDD",
  strip = "#F2F2F2",
  match = "#000000",
  mismatch = "#CC3311",
  no_consensus = "#777777",
  sequential = c("#FFFFE5", "#FFF7BC", "#FEE391", "#FEC44F", "#FB9A29"),
  diverging = c("#FDB366", "#FEDA8B", "#EAECCC", "#C2E4EF", "#98CAE1"),
  clusters = c("#8CCBBF", "#EDE5B8", "#F2F2F2"),
  variant_ranks = c("#4A1230", "#882255", "#CC6677", "#F2D7DE"),
  match_quality = c(`1` = "#6D2163", `2` = "#B565A7", `3` = "#EBD3E6"),
  # ColorBrewer RdYlBu, as Bioconductor's biovizBase uses for DNA bases
  dna_bases = c(A = "#ABD9E9", C = "#FDAE61", G = "#D7191C", T = "#2C7BB6",
                N = "#FFFFBF")
)

#' Black or white, whichever reads better on each fill
#'
#' Uses the WCAG relative luminance; the threshold 0.3 sits where black and
#' white text have about the same contrast.
#' @param fill Colours, in any form \code{col2rgb()} takes.
#' @return A character vector of \code{"black"} and \code{"white"}.
#' @noRd
.text_colour_on <- function(fill) {
  m <- grDevices::col2rgb(fill) / 255
  lin <- ifelse(m <= 0.03928, m / 12.92, ((m + 0.055) / 1.055)^2.4)
  lum <- 0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ]
  ifelse(lum > 0.3, "black", "white")
}
