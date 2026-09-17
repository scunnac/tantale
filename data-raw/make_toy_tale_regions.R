# Build the toy subject fixture used by the frameshift-correction tests.
#
# The shipped BAI3 fixture is 116 kb over two regions holding four TALEs, and
# a corrected run against it costs minutes. More to the point, nothing in it
# has a *known* right answer: the existing tests check that correction does
# not error, never that it corrects anything.
#
# This builds a smaller subject with one, from the same sequences rather than
# synthesised, so the DNA stays biologically real:
#
#   toy_intact       one complete TALE, untouched -- the negative control
#   toy_frameshift   the same TALE with one nucleotide inserted, at a
#                    position recorded in the companion TSV
#   toy_no_tale      a stretch of the same contig holding no TALE at all
#
# Two regions carrying the same TALE is deliberate: the intact copy is the
# control for the frameshifted one, so any difference in the correction is
# attributable to the inserted base and not to the two arrays differing.
#
# Run with:  Rscript data-raw/make_toy_tale_regions.R

suppressMessages({library(Biostrings); devtools::load_all(".", quiet = TRUE)})

OUT_DIR  <- "tests/testthat/data_for_tests"
FLANK    <- 500   # enough for extend_len = 300 to have somewhere to extend
SEED_ROI <- "ROI_00001"

src <- system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
                   package = "tantale", mustWork = TRUE)
genome <- readDNAStringSet(src)


## ---- locate the arrays ---------------------------------------------------

scratch <- file.path(tempdir(), "toy_build")
unlink(scratch, recursive = TRUE)
invisible(suppressWarnings(suppressMessages(
  tell_tales(subject_file = src, output_dir = scratch))))

report <- readr::read_tsv(file.path(scratch, "arrayReport.tsv"),
                          show_col_types = FALSE, progress = FALSE)
row <- as.data.frame(report[report$array_id == SEED_ROI, ])
stopifnot(nrow(row) == 1L)

contig <- genome[[row$OriginalSubjectName]]
from <- max(1L, row$Start - FLANK)
to   <- min(length(contig), row$End + FLANK)
region <- subseq(contig, from, to)
# the array's coordinates within the extracted region
arrayFrom <- row$Start - from + 1L
arrayTo   <- row$End - from + 1L

message("seed array ", SEED_ROI, " on ", row$OriginalSubjectName,
        " ", row$Start, "-", row$End, " (", row$Strand, "), ",
        row$NumberOfHits, " domain hits")
message("extracted region: ", length(region), " nt, array at ",
        arrayFrom, "-", arrayTo)


## ---- insert one nucleotide ------------------------------------------------
#
# Placed near the middle of the repeat array: far enough from either end that
# the frameshift truncates the ORF well short of the real stop, which is what
# makes it detectable. Putting it near an end would let the ORF finder simply
# pick the other side.

insertAt <- arrayFrom + (arrayTo - arrayFrom) %/% 2L
insertedBase <- "A"

frameshifted <- xscat(subseq(region, 1L, insertAt),
                      DNAString(insertedBase),
                      subseq(region, insertAt + 1L, length(region)))

stopifnot(length(frameshifted) == length(region) + 1L)
message("inserted '", insertedBase, "' after position ", insertAt,
        " of the region")


## ---- a region with no TALE ------------------------------------------------
#
# Taken from the far end of the other contig, checked below to be clear of
# every array tell_tales() found.

other <- genome[[setdiff(names(genome), row$OriginalSubjectName)[1]]]
otherName <- setdiff(names(genome), row$OriginalSubjectName)[1]
otherArrays <- as.data.frame(report[report$OriginalSubjectName == otherName, ])
# a 4 kb window after the last array on that contig
noTaleFrom <- max(otherArrays$End) + 2000L
noTaleTo <- min(length(other), noTaleFrom + 4000L)
stopifnot(noTaleFrom < noTaleTo)
noTale <- subseq(other, noTaleFrom, noTaleTo)


## ---- write ----------------------------------------------------------------

toy <- DNAStringSet(list(intact = region,
                         frameshift = frameshifted,
                         no_tale = noTale))
names(toy) <- c("toy_intact", "toy_frameshift", "toy_no_tale")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
writeXStringSet(toy, file.path(OUT_DIR, "toy_tal_regions.fasta"))

# The companion file is the point of the fixture: it records the right
# answer, so a test can assert the correction finds it rather than merely
# assert the call returned.
truth <- data.frame(
  seqname = names(toy),
  source_array = c(SEED_ROI, SEED_ROI, NA),
  source_contig = c(row$OriginalSubjectName, row$OriginalSubjectName, otherName),
  source_from = c(from, from, noTaleFrom),
  source_to = c(to, to, noTaleTo),
  array_from = c(arrayFrom, arrayFrom, NA),
  array_to = c(arrayTo, arrayTo, NA),
  insertion_at = c(NA, insertAt, NA),
  inserted_base = c(NA, insertedBase, NA),
  expected_tale = c(TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
readr::write_tsv(truth, file.path(OUT_DIR, "toy_tal_regions_truth.tsv"))

message("\nwritten to ", OUT_DIR, ":")
message("  toy_tal_regions.fasta       ", sum(width(toy)), " nt in ",
        length(toy), " regions")
message("  toy_tal_regions_truth.tsv   the known answer")
print(toy)
