# Build the frameshift-correction reference sets shipped in inst/extdata.
#
# Source: data-raw/tale_correction_ref_source.fa.gz -- the raw output of
# tell_tales() over 70 Xanthomonas oryzae genomes, assembled years ago and
# never filtered. 1057 sequences, of which 555 are exact duplicates of another
# entry and 8 are too short to be TALEs at all.
#
# Two sets come out of it. Neither is "better"; they trade coverage for size.
#
#   tale_correction_ref.fa.gz            ~494  the default
#   tale_correction_ref_representative.fa.gz ~134  a diversity-sampled subset
#
# Run with:  Rscript data-raw/make_correction_references.R

suppressMessages({library(Biostrings); library(DECIPHER)})

src <- "data-raw/tale_correction_ref_source.fa.gz"
raw <- readAAStringSet(src)
message("source: ", length(raw), " sequences, ", sum(width(raw)), " aa")


## ---- the default set: everything that carries information ----------------
#
# Deduplication is unarguable: 555 entries are byte-identical to another, the
# same TALE found in many sequenced strains. The same alignment recomputed
# adds nothing. (It does not make correction faster -- DECIPHER scores every
# reference cheaply and only aligns against the closest maxComparisons of them
# -- but it halves the file.)
#
# The length floor removes fragments that are not TALEs: the shortest entry is
# 23 aa against a median of 1198. They can never be the best reference for a
# real array.
#
# Pseudogenes are KEPT, deliberately. They come from high-quality genomes, so
# their frameshifts are real biology rather than sequencing error, and the
# point of correction is to recover the sequence as it exists in nature -- not
# to reshape every array into an intact TALE. A reference set of only intact
# TALEs risks "repairing" a genuine pseudogene into an ORF that no strain
# carries.

MIN_LENGTH <- 300

dedup <- raw[!duplicated(as.character(raw))]
full  <- dedup[width(dedup) >= MIN_LENGTH]
message("default set: ", length(full), " sequences, ", sum(width(full)), " aa",
        "  (", sum(grepl("(Pseudo)", names(full), fixed = TRUE)), " pseudogenes kept)")


## ---- the representative set: all pseudogenes + clustered intact ----------
#
# For users who want a smaller reference. Every pseudogene is carried over,
# since they are the scarce and irreplaceable part; the intact TALEs, which
# are numerous and highly similar, are clustered and one representative per
# cluster is kept.
#
# Clustering settings, and why:
#
#   includeTerminalGaps = TRUE       count length differences as differences,
#                                    so a truncated sequence stays distinct
#                                    from a full-length one.
#   penalizeGapLetterMatches = NA    the default: a run of gaps counts as one
#                                    mismatch, not one per residue. TALEs
#                                    differ chiefly by whole repeats, and
#                                    with TRUE a difference of a few repeats
#                                    (~34 aa each) swamps every other
#                                    difference. What predicts a good
#                                    correction reference is overall sequence
#                                    similarity, not matching repeat count.
#   method = "overlap"               stated for clarity but inert: the
#                                    documentation notes method only applies
#                                    when includeTerminalGaps is FALSE.
#
# The representative of each cluster is its longest member: a longer reference
# gives the aligner more to anchor on, and a truncated one can only match part
# of a candidate.

CUTOFF <- 0.04

isPseudo <- grepl("(Pseudo)", names(full), fixed = TRUE)
intact   <- full[!isPseudo]

clusters <- Clusterize(intact, cutoff = CUTOFF,
                       method = "overlap",
                       includeTerminalGaps = TRUE,
                       penalizeGapLetterMatches = NA,
                       processors = NULL, verbose = FALSE)[[1]]

pick <- vapply(split(seq_along(intact), clusters),
               function(i) i[which.max(width(intact)[i])], integer(1))
representative <- c(full[isPseudo], intact[sort(pick)])

message("representative set: ", length(representative), " sequences, ",
        sum(width(representative)), " aa",
        "  (", sum(isPseudo), " pseudogenes + ", length(pick), " intact)")


## ---- write ---------------------------------------------------------------

writeXStringSet(full, "inst/extdata/tale_correction_ref.fa.gz", compress = TRUE)
writeXStringSet(representative,
                "inst/extdata/tale_correction_ref_representative.fa.gz",
                compress = TRUE)
message("written to inst/extdata/")
