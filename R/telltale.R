
# subject_file = system.file("extdata", "bai3_sample_tal_regions.fasta", package = "tantale", mustWork = T)
# output_dir = tempdir(check = TRUE)
# hmm_dir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T)
# hmmer_path = system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = T)
# correct_array = TRUE
# correction_ref = system.file("extdata", "decipher_ref_tales_aa.fa.gz", package = "tantale", mustWork = T)
# frameshift = -11
# nterm_min_score = 300
# repeat_min_score = 20
# cterm_min_score = 200
# min_domain_hits = 4
# merge_hits = TRUE
# min_gap = 35
# taleArrayStartAnchorCode = "NTERM"
# taleArrayEndAnchorCode = "CTERM"
# extremity_codes = TRUE
# rvd_sep = "-"
# extend_len = 300
# ... = NULL

# subject_file = "/home/cunnac/TEMP/220928-8_talCor.fasta"
# output_dir = file.path("/home/cunnac/TEMP", gsub("(\\.fasta)|(\\.fa)|(\\.fna)|(\\.fsa)", "", basename(subject_file)))
# hmm_dir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T)
# hmmer_path = system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = T)
# correct_array = FALSE
# correction_ref = system.file("extdata", "decipher_ref_tales_aa.fa.gz", package = "tantale", mustWork = T)
# frameshift = -11
# nterm_min_score = 300
# repeat_min_score = 20
# cterm_min_score = 200
# min_domain_hits = 4
# merge_hits = TRUE
# min_gap = 35
# extremity_codes = TRUE
# rvd_sep = "-"
# extend_len = 300
# ... = NULL







#' Read the three TALE profile HMMs and concatenate them for nhmmer
#'
#' The search is done with one merged profile file rather than three separate
#' runs, so the three are read, checked, and written out together here.
#'
#' Each profile's \code{NAME} tag is what nhmmer reports in its
#' \code{query_name} column, which is how the hits are later told apart --
#' N-terminus from repeat from C-terminus. That makes the names part of the
#' contract between this stage and every stage that filters hits, so they are
#' returned rather than re-parsed downstream.
#'
#' @param hmm_dir Directory holding the profiles. They must carry the names
#'   below; see the TODO in \code{tell_tales()} about letting a user supply
#'   arbitrary paths.
#' @param merged_file Where to write the concatenation.
#' @return A list with the \code{nterm}, \code{repeats} and \code{cterm}
#'   profile names, in the order nhmmer will see them, plus \code{files},
#'   the three paths they were read from.
#' @noRd
.telltale_hmm_profiles <- function(hmm_dir, merged_file) {
  files <- c(file.path(hmm_dir, "Xo_TALE_Nterm_CDS_profile.hmm"),
             file.path(hmm_dir, "Xo_TALE_repeat_CDS_profile.hmm"),
             file.path(hmm_dir, "Xo_TALE_Cterm_CDS_profile.hmm"))

  lines <- plyr::llply(files, function(f) readLines(con = f))
  names <- unlist(plyr::llply(lines, function(x) {
    hmmName <- grep("NAME", x, perl = TRUE, value = TRUE)
    hmmName <- unlist(strsplit(hmmName, split = "\\s+"))
    if (length(hmmName) != 2) stop("One or several profile ",
                                   "HMM have a name with spaces. ",
                                   "Please remove them in the file at the Tag 'NAME'")
    hmmName[2]
  }))

  writeLines(text = unlist(lines), con = merged_file)
  list(nterm = names[1], repeats = names[2], cterm = names[3],
       # the paths are reported in the run log, so they are carried too
       files = stats::setNames(files, c("nterm", "repeats", "cterm")))
}


#' Run nhmmer and return the TALE domain hits that survive filtering
#'
#' The first stage of \code{tell_tales()}: search the subject sequences with
#' the merged profile, read the tabular output, and keep the hits worth
#' carrying forward.
#'
#' Filtering happens twice, for different reasons. Each domain type is scored
#' against its own threshold, because an N-terminus, a repeat and a C-terminus
#' are different lengths and score on different scales. Then whole subject
#' sequences are dropped if they carry too few hits, on the reasoning that a
#' contig with one or two stray domain matches is unlikely to hold a real
#' TALE.
#'
#' @param subject_file Sequences to search.
#' @param hmm What \code{.telltale_hmm_profiles()} returned.
#' @param paths What \code{.telltale_paths()} returned.
#' @param hmmer_path Directory holding the nhmmer binary.
#' @param nterm_min_score,repeat_min_score,cterm_min_score Per-domain score
#'   thresholds.
#' @param min_domain_hits A subject sequence is kept when it carries *more*
#'   than this many hits. Note it counts hits per subject sequence, not per
#'   TALE array; see restructuring-notes.md 8.0.
#' @return The filtered table, or \code{NULL} when nothing survives -- which
#'   the caller must treat as "stop here", since every later stage assumes at
#'   least one hit.
#' @noRd
.telltale_find_domain_hits <- function(subject_file, hmm, paths, hmmer_path,
                                       nterm_min_score, repeat_min_score,
                                       cterm_min_score, min_domain_hits) {
  .run_nhmmer_search(hmmer_path = hmmer_path,
                     subject_file = subject_file,
                     hmm_file = paths$merged_hmm,
                     search_out_file = paths$hmmer_search,
                     readable_out_file = paths$hmmer_readable)

  hits <- try(read.table(paths$hmmer_search), silent = TRUE)
  if (inherits(hits, "try-error")) {
    warning("NhmmerSearch found no TALE cds hit in ", subject_file , " Exitting...")
    return(NULL)
  }

  colnames(hits) <- c("target_name", "accession", "query_name", "accession", "hmmfrom", "hmm_to", "alifrom",
                      "ali_to", "envfrom", "env_to", "sq_len", "strand", "Evalue", "score", "bias", "description_of_target")

  ## filtering results differentially depending on the query HMM
  hits <- subset(hits,
                 query_name == hmm$nterm & score >= nterm_min_score |
                   query_name == hmm$repeats & score >= repeat_min_score |
                   query_name == hmm$cterm & score >= cterm_min_score
  )
  hits <- droplevels(hits)
  if (nrow(hits) == 0L) {
    cli::cli_warn("No record remains after filtering NhmmerSearch hits based on score. Exitting...")
    return(NULL)
  }
  ## Add a hitID column
  hits$hitID <- paste("DOM", sprintf("%05.0f", 1:nrow(hits)), sep = "_")

  ## Trick to re-order positions in an increasing order to satisfy IRanges() in preparation of creating a GRanges
  hits[, c("start", "end")] <- plyr::adply(.data = hits[, c("envfrom", "env_to")],
                                           .margins = 1, .fun = c(min, max))[, -(1:2)]
  rownames(hits) <- hits$hitID

  ## Filter out target DNA sequences that have too few repeat CDSs
  ## NB: for the sake of consistency  it would be better just to filter out
  ## from any further consideration the ARRAYS shorter than a certain value (say 5).
  ## WHAT DO WE DO ABOUT THAT?
  perSubject <- plyr::ddply(hits[, -20], ~ target_name + sq_len, nrow) # I do not know why but it fails to work if I leave the RVD column (#20)
  hits <- subset(hits, target_name %in% perSubject[perSubject$V1 > min_domain_hits, "target_name"])
  hits <- droplevels(hits)
  if (nrow(hits) == 0L) {
    # Unguarded before: the run carried on and died several stages later
    # inside Bioconductor with "Rle of type 'NULL' is not supported".
    cli::cli_warn(c("No subject sequence carries more than {min_domain_hits} TALE domain hit{?s}.",
                    "i" = "{.arg min_domain_hits} counts hits per subject sequence, not per TALE array.",
                    "x" = "Nothing left to analyse. Exitting..."))
    return(NULL)
  }

  hits
}


#' Merge hits of the same domain type that overlap each other
#'
#' nhmmer can report the same repeat twice, as two overlapping hits. Left
#' alone those become two repeats in the array, and the inferred RVD sequence
#' gains a residue that is not there.
#'
#' Merging is done per domain type, never across types: an N-terminus hit
#' overlapping a repeat hit is a real feature of where one domain ends and the
#' next begins, not a duplicate.
#'
#' The identifiers of the hits that went into each merged range are kept in
#' \code{nhmmerHitID}, separated by \code{|}, so a merged range can be traced
#' back to the raw search output.
#'
#' @param gr Hits as a \code{GRanges}, with \code{query_name} naming the
#'   domain type and \code{hitID} identifying each hit.
#' @return A \code{GRanges} of merged hits, re-identified as \code{MDOM_*}.
#' @noRd
.telltale_merge_overlapping_hits <- function(gr) {
  byDomain <- GenomicRanges::split(gr, f = gr$query_name)

  merged <- lapply(byDomain, function(g) {
    reduced <- as.data.frame(GenomicRanges::findOverlaps(
        g, g, minoverlap = 2, type = "any", ignore.strand = FALSE, select = "all")) %>%
      dplyr::group_by(queryHits) %>%
      dplyr::group_map({
        ~ GenomicRanges::reduce(g[as.numeric(.x$subjectHits)], with.revmap = FALSE)
      }) %>%
      plyranges::bind_ranges() %>%
      unique()

    # which raw hits ended up inside each merged range
    formerIDs <- as.data.frame(GenomicRanges::findOverlaps(
        g, reduced, minoverlap = 2, type = "within", ignore.strand = FALSE, select = "all")) %>%
      dplyr::group_by(subjectHits) %>%
      dplyr::group_map({
        ~ paste0(g[as.numeric(.x$queryHits)]$hitID, collapse = "|")
      }) %>%
      unlist()
    reduced$nhmmerHitID <- formerIDs
    reduced
  }) %>%
    plyranges::bind_ranges(.id = "query_name")

  merged$hitID <- paste("MDOM", sprintf("%05.0f", 1:length(merged)), sep = "_")
  names(merged) <- merged$hitID
  merged
}


#' Group neighbouring domain hits into candidate TALE arrays
#'
#' A TALE is a run of domain hits close together on the same strand: an
#' N-terminus, a series of repeats, a C-terminus. Hits separated by less than
#' \code{min_gap} are taken to belong to the same array, and each array
#' becomes a region of interest, \code{ROI_*}.
#'
#' Hits within an array should not overlap -- the merge stage exists to make
#' sure of that -- so any that still do are reported. They matter because the
#' RVD sequence is read off the repeats in order, and two overlapping hits
#' put a residue in it that is not in the protein.
#'
#' @param gr Domain hits, merged.
#' @param min_gap Largest gap, in bases, still counted as contiguous.
#' @param subject_seqs The DNA the hits were found in, for extracting each
#'   array's sequence.
#' @param hmm What \code{.telltale_hmm_profiles()} returned; used to record
#'   whether an array carries all three domain types.
#' @return A list of \code{arrays} (one range per array) and \code{by_array}
#'   (the hits, grouped, carrying the per-array metadata).
#' @noRd
.telltale_group_arrays <- function(gr, min_gap, subject_seqs, hmm) {
  ## Use reduce to obtain the regions that span "contiguous" hits
  arraysGR <- GenomicRanges::reduce(gr,
                                    drop.empty.ranges = FALSE,
                                    min.gapwidth = min_gap,
                                    with.revmap = TRUE,
                                    ignore.strand = FALSE)
  revmap <- S4Vectors::mcols(arraysGR)$revmap  # an IntegerList

  ## Use the mapping from reduced to original ranges to group the originals
  byArray <- BiocGenerics::relist(gr[unlist(revmap)], revmap)
  names(byArray) <- paste("ROI", sprintf("%05.0f", 1:length(byArray)), sep = "_")
  names(arraysGR) <- names(byArray)

  ## Make sure that hits do not overlap for some weird reason
  doHitsOverlap <- !GenomicRanges::isDisjoint(byArray)
  if (any(doHitsOverlap)) {
    warning("It appears that some hmmer hits actually overlap.\n It is thus possible that the inferred sequences of RVDs have artefactual insertions.\n")
    warning(paste0("Please check the hits in the following RegionsOfInterest:", "\n",
                   paste(names(doHitsOverlap)[doHitsOverlap], collapse = "\n"), "\n")
    )
  }

  ## Populate metadata about the elements of the list of arrays
  S4Vectors::mcols(byArray) <- S4Vectors::DataFrame(
    array_id = names(byArray),
    OriginalSubjectName = sapply(byArray,
                                 function(x) unique(as.character(GenomicRanges::seqnames(x)))),
    Start = BiocGenerics::start(arraysGR),
    End = BiocGenerics::end(arraysGR),
    Strand = BiocGenerics::strand(arraysGR),
    NumberOfHits = S4Vectors::elementNROWS(byArray),
    ArraySeq = BSgenome::getSeq(subject_seqs, arraysGR),
    AllDomains = sapply(byArray,
                        function(x) {
                          all(
                            c(hmm$nterm, hmm$repeats,
                              hmm$cterm) %in% as.character(x$query_name)
                          )
                        }
    )
  )

  list(arrays = arraysGR, by_array = byArray)
}


#' Write the run log
#'
#' Echoes the parameters the run was given and a handful of summary measures,
#' to the console and to a file. Nothing downstream reads it; it exists so
#' that a directory of results can be read months later and still say what
#' produced it.
#'
#' @param params The user-facing arguments, echoed verbatim.
#' @param log_file Where to write.
#' @param hmm What \code{.telltale_hmm_profiles()} returned.
#' @param subject_seqs The sequences that were searched.
#' @param arrays,by_array What \code{.telltale_group_arrays()} returned.
#' @param array_report The assembled array report.
#' @param gaps_below_500,gap_quartiles Gap statistics between arrays.
#' @param annotale_messages Whatever AnnoTALE complained about, if anything.
#' @return \code{NULL}, invisibly. Called for its output.
#' @noRd
.telltale_log <- function(params, log_file, hmm, subject_seqs, arrays, by_array,
                          array_report, gaps_below_500, gap_quartiles,
                          annotale_messages) {
  
  ## counts of appearance of each RVD type (excluding N- and C- terms symbols) for the log file
  #RVDtbl <- table(subset(unlist(by_array), query_name ==hmm$repeats, drop = TRUE)$RVD)
  ## Total count of repeat CDS after filtering for uniformative subject seqs for the log file
  numberOfRepeatHitsAfterFiltering <- length(subset(unlist(by_array), query_name == hmm$repeats))
  ## Distribution of the number of hits per array
  countsHitsByArrayDistri <- summary(S4Vectors::mcols(by_array)$NumberOfHits)
  ## Number of domains in arrays that display all domain types
  # completeArrayLengths <- subset(S4Vectors::mcols(by_array), AllDomains)$NumberOfHits
  
  
  ## might have been cleaner with a glue approach
  txt <- c(
    "#****************************************",
    "#**   tell_tales analysis done     **",
    
    paste("Current date:", date(), sep = "\t"),
    "#_________Provided I/O parameters __________",
    paste("File of subject DNA sequences:", params$subject_file, sep = "\t"),
    paste("TALE N-term CDS region detection HMM file:", hmm$files[["nterm"]], sep = "\t"),
    paste("TALE repeat unit CDS detection HMM file:", hmm$files[["repeats"]], sep = "\t"),
    paste("TALE C-term CDS region detection HMM file:", hmm$files[["cterm"]], sep = "\t"),
    paste("Output directory:", params$output_dir, sep = "\t"),
    
    "#____________Other parameters________________",
    paste("nterm_min_score",":", params$nterm_min_score, sep = "\t"),
    paste("repeat_min_score",":", params$repeat_min_score, sep = "\t"),
    paste("cterm_min_score",":", params$cterm_min_score, sep = "\t"),
    paste("min_domain_hits",":", params$min_domain_hits, sep = "\t"),
    paste("merge_hits",":", params$merge_hits, sep = "\t"),
    paste("min_gap",":", params$min_gap, sep = "\t"),
    paste("extend_len",":", params$extend_len, sep = "\t"),
    paste("correct_array",":", params$correct_array, sep = "\t"),
    paste("correction_ref",":", params$correction_ref, sep = "\t"),
    paste("frameshift",":", params$frameshift, sep = "\t"),
    
    "#__________Summary measures of TALE search outcome__________",
    paste("Number of analysed subject sequences :", length(subject_seqs), sep = "\t"),
    paste("Total number of TALE repeat DNA coding sequence motif hits found with the nhmmer approach:",
          numberOfRepeatHitsAfterFiltering, sep = "\t"),
    #paste("Total number of repeat HMM hits on the corresponding set of translated DNA hits:", sum(RVDtbl), sep = "\t"),
    
    paste("Total number of subject seqs with TALE motif hits after low hit number filtering:",
          length(GenomeInfoDb::seqlevels(arrays)), sep = "\t"),
    paste("Total number of distinct regions (repeat arrays) with adjacent TALE motifs :", nrow(array_report), sep = "\t"),
    paste("Total number of 'complete' arrays (with both N- and C-term flanking motifs):",
          sum(S4Vectors::mcols(by_array)$AllDomains),	sep = "\t"),
    
    #paste("Total number of distinct types of RVD:", nrow(RVDtbl), sep = "\t"),
    
    paste("Minimum array length (number of TALE domain hits):", min(array_report$NumberOfHits), sep = "\t"),
    paste("Maximum array length:", max(array_report$NumberOfHits), sep = "\t"),
    paste("Median array length:", median(array_report$NumberOfHits), sep = "\t"),
    # paste("Length of the longest 'complete' array:", max(completeArrayLengths),	sep = "\t"),
    # paste("Length of the shortest 'complete' array:", min(completeArrayLengths),	sep = "\t"),
    
    #message("Distribution of the number of TALE domain hits per array:\n")
    #message(paste(names(countsHitsByArrayDistri), countsHitsByArrayDistri, sep = "\t", collapse = "\n"))
    
    paste("Number of gaps of size below 500nt between TALE motifs arrays:", length(gaps_below_500)/2, sep = "\t"),
    
    paste("First quartile of size of gaps (below 500nt) between TALE motifs arrays:", gap_quartiles[1], sep = "\t"),
    paste("Median size of gaps (below 500nt) between TALE motifs arrays:", gap_quartiles[2], sep = "\t"),
    paste("Upper quartile of size of gaps (below 500nt) between TALE motifs arrays:", gap_quartiles[3], sep = "\t"),
    
    "#__________Noteworthy AnnoTale issues__________",
    paste("#", annotale_messages),
    
    "#*************************\n"
  )
  
  message(paste(txt, collapse = "\n"))
  logf <- file(log_file, open = "w")
  writeLines(text = txt, con = logf)
  close(logf)
  invisible(NULL)
}


#' Write before-and-after alignments of each corrected array
#'
#' For every array, three versions are aligned and written as a browsable
#' HTML page: the sequence as found, the sequence after frameshift
#' correction, and the version handed to AnnoTALE with any \code{N}
#' substituted. They exist to be looked at -- frameshift correction changes
#' the reading frame of a real sequence, and this is what lets a user see
#' what was changed and judge whether to believe it.
#'
#' Both the DNA and its translation are written, because a frameshift is
#' obvious in the protein and easy to miss in the nucleotides.
#'
#' @param raw,corrected,substituted The three versions, named by array.
#' @param dna_dir,aa_dir Where to write.
#' @return \code{NULL}, invisibly.
#' @noRd
.telltale_write_correction_alignments <- function(raw, corrected, substituted,
                                                  dna_dir, aa_dir) {
  for (n in names(corrected)) {
    rawSeq <- raw[n]
    names(rawSeq) <- paste0("raw_", n)
    correctedSeq <- corrected[n]
    names(correctedSeq) <- paste0("corrected_", n)
    substitutedSeq <- substituted[n]
    names(substitutedSeq) <- paste0("forAnnoTALE", n)

    seqToAlign <- c(rawSeq, correctedSeq, substitutedSeq)
    alignedSeqs <- DECIPHER::AlignSeqs(seqToAlign, verbose = FALSE)
    DECIPHER::BrowseSeqs(alignedSeqs,
                         htmlFile = file.path(dna_dir, glue::glue("CorrectionAlignmentDNA_{n}.html")),
                         openURL = FALSE, colWidth = 120)

    seqToAlignTranslated <- Biostrings::translate(seqToAlign, no.init.codon = TRUE,
                                                  if.fuzzy.codon = "solve")
    alignedSeqsTranslated <- DECIPHER::AlignSeqs(seqToAlignTranslated, verbose = FALSE)
    DECIPHER::BrowseSeqs(alignedSeqsTranslated,
                         htmlFile = file.path(aa_dir, glue::glue("CorrectionAlignmentAA_{n}.html")),
                         openURL = FALSE, colWidth = 120)
  }
  invisible(NULL)
}


#' Find the TALE ORF in each array, correcting frameshifts if asked
#'
#' Each extended array region is searched for its longest open reading frame,
#' which is the putative TALE coding sequence handed to AnnoTALE.
#'
#' With \code{correct_array = TRUE} the arrays are first run through
#' \code{DECIPHER::CorrectFrameshifts()} against a reference set of TALE
#' proteins. A sequencing error that shifts the reading frame truncates the
#' ORF and loses every repeat downstream of it, so a TALE that is real can
#' look like a fragment; correction restores the frame. It is expensive:
#' every array is compared against every reference, so the cost is the number
#' of arrays times the size of the reference set.
#'
#' Correction can leave \code{N} in a sequence where it inserted a base it
#' could not call. AnnoTALE will not read those, so they are substituted with
#' \code{C} before it runs -- a change to the sequence given to AnnoTALE, not
#' to the one reported, and the written alignments show exactly where it
#' happened.
#'
#' @param array_seqs The extended array sequences.
#' @param by_array The grouped hits, whose metadata gains the per-array indel
#'   counts when correction runs.
#' @param correct_array Whether to correct.
#' @param correction_ref Fasta of reference TALE proteins.
#' @param frameshift Frameshift penalty passed to DECIPHER.
#' @param paths What \code{.telltale_paths()} returned.
#' @param ... Passed to \code{DECIPHER::CorrectFrameshifts()}.
#' @return A list of \code{orf} (what AnnoTALE is given), \code{full_orf}
#'   (what is reported), and \code{by_array}, updated.
#' @noRd
.telltale_array_orfs <- function(array_seqs, by_array, correct_array,
                                 correction_ref, frameshift, paths, ...) {
  if (!correct_array) {
    orfs <- systemPipeR::predORF(x = array_seqs,
                                 n = 1, type = "gr", mode = "ORF", strand = "sense")
    fullTalOrf <- BSgenome::getSeq(array_seqs, orfs)
    names(fullTalOrf) <- as.character(GenomicRanges::seqnames(orfs))
    return(list(orf = fullTalOrf, full_orf = fullTalOrf, by_array = by_array))
  }

  # An alternative approach: https://github.com/Jstacs/Jstacs/tree/master/projects/talecorrect
  cli::cli_inform(paste0("Correcting putative TALE coding sequences. Be patient, this may take a LONG time..."))
  AAref <- Biostrings::readAAStringSet(correction_ref, seek.first.rec = TRUE, use.names = TRUE)
  ## TO SPEEDUP CORRECTION could correct only predicted ORFs that cover less than X% of the raw sequence
  correction <- DECIPHER::CorrectFrameshifts(array_seqs,
                                             AAref, type = "both",
                                             maxComparisons = length(AAref),
                                             frameShift = frameshift, ...)
  cli::cli_inform("Correction of putative TALE coding sequences is done!")
  corrected <- correction$sequences

  ####   Correction stats   ####
  indels <- function(which, name) {
    .correction_tibble(correction$indels) %>%
      dplyr::group_by(Seq) %>%
      dplyr::count(variable, name = name) %>%
      dplyr::filter(variable == which) %>%
      dplyr::select(-variable)
  }
  S4Vectors::mcols(by_array) <- merge(
    S4Vectors::mcols(by_array),
    dplyr::full_join(indels("insertions", "predicted_ins_count"),
                     indels("deletions", "predicted_dels_count"), by = "Seq"),
    by.x = "array_id", by.y = "Seq", all.x = TRUE)
  S4Vectors::mcols(by_array)[c("predicted_dels_count", "predicted_ins_count")] %<>%
    apply(., 2, function(v) ifelse(is.na(v), 0, v))

  for (n in names(corrected)[Biostrings::vcountPattern("N", corrected) > 0]) {
    cli::cli_warn(paste0("After correction, {n} sequence contains 'N's which will be substituted by 'C's in order",
                         "to run AnnoTALE analyze for RVDs prediction."))
  }
  substituted <- Biostrings::chartr("N", "C", corrected)

  .telltale_write_correction_alignments(
    raw = array_seqs, corrected = corrected, substituted = substituted,
    dna_dir = paths$correction_dna, aa_dir = paths$correction_aa)

  orfs <- systemPipeR::predORF(x = substituted,
                               n = 1, type = "gr", mode = "ORF", strand = "sense")
  talOrf <- BSgenome::getSeq(substituted, orfs)
  names(talOrf) <- as.character(GenomicRanges::seqnames(orfs))

  # The idea here was to convert back the Ns that were substituted for AnnoTALE in order to output
  # an unsubstituted orf.
  # I am not sure the code below properly does that and retrospectively,
  # it may be a problem for downstream bioinformatic analyses to have Ns in the sequence...
  #
  # insPosition <- vmatchPattern("N", corrected)
  # for (i in names(insPosition)) {
  #   if (length(width(insPosition[[i]])) == 0) next()
  #   insPositionAfCorr <- insPosition[[i]][BiocGenerics::start(insPosition[[i]]) < width(talOrf[i])]
  #   talOrf[i] <- replaceAt(talOrf[i], insPositionAfCorr, value = "N")
  # }
  list(orf = talOrf, full_orf = talOrf, by_array = by_array)
}


#' Run AnnoTALE's "analyze" stage on one putative TALE ORF
#'
#' \code{\link{run_annotale_predict}} runs AnnoTALE's predict and analyze
#' stages together, starting from a genome. \code{tell_tales()} has already
#' found the ORF by then, so it needs analyze on its own.
#'
#' @param fasta_file The ORF to analyze.
#' @param output_dir Where AnnoTALE writes its parts and RVD files.
#' @param prefix TALE name prefix; taken from the file name when absent.
#' @param annotale_jar Path to the AnnoTALE jar.
#' @return AnnoTALE's exit status, invisibly.
#' @noRd
.run_annotale_analyze <- function(fasta_file,
                                  output_dir = getwd(),
                                  prefix = NULL,
                                  annotale_jar = system.file("tools", "AnnoTALEcli-1.5.jar",
                                                             package = "tantale", mustWork = TRUE)) {
  stopifnot(dir.exists(output_dir) || dir.create(path = output_dir, showWarnings = TRUE,
                                                 recursive = TRUE, mode = "775"))
  # Define a prefix for TALEs (assembly ID) derived from the genome file name.
  if (is.null(prefix)) {
    prefix <- gsub(pattern = "^(.*)\\.(fasta|fa|fas)$",
                   replacement = "\\1", basename(fasta_file),
                   perl = TRUE)
  }
  comAnalyze <- paste0(
    "java -jar ", annotale_jar,
    " analyze ",
    " t=", fasta_file,
    " outdir=", output_dir
  )
  invisible(system(comAnalyze, ignore.stdout = TRUE, ignore.stderr = TRUE))
}


#' Read each array's RVD sequence and domain composition from AnnoTALE
#'
#' AnnoTALE is run once per array, in its own directory, and its output read
#' back. It is the step that turns a coding sequence into the biology: which
#' parts are the N-terminus, the repeats and the C-terminus, and what RVD
#' each repeat carries.
#'
#' It fails on some ORFs, and does so in several ways -- the jar errors, or
#' it writes no protein parts, or it writes an empty RVD file. All of them
#' mean the same thing here, that this array yielded nothing, so each is
#' caught, the empty file removed so nothing downstream reads it, and the
#' array reported and skipped rather than aborting the run.
#'
#' @param orfs The putative TALE ORFs, named by array.
#' @param by_array The grouped hits, read for each array's source sequence
#'   name.
#' @param annotale_dir Directory to create the per-array subdirectories in.
#' @return A list of \code{rvds} (one sequence per array that worked),
#'   \code{domains} (their domain composition) and \code{messages} (what
#'   failed, for the run log).
#' @noRd
.telltale_run_annotale <- function(orfs, by_array, annotale_dir) {
  messages <- character()

  out <- lapply(names(orfs), function(talOrfID) {
    AnnotaleDir <- file.path(annotale_dir, talOrfID)
    dir.create(AnnotaleDir)
    TalOrf <- orfs[talOrfID]
    correctedTalOrfFile <- file.path(AnnotaleDir, "putativeTalOrf.fasta")
    Biostrings::writeXStringSet(TalOrf, correctedTalOrfFile)

    cli::cli_inform("Now running AnnoTALE analyze for {talOrfID}")
    checkAnnoTale <- try(.run_annotale_analyze(correctedTalOrfFile, AnnotaleDir), silent = TRUE)

    prot_parts_files <- file.path(AnnotaleDir, "TALE_Protein_parts.fasta")
    annoTaleRVD_file <- file.path(AnnotaleDir, "TALE_RVDs.fasta")
    seqOfRVDs <- try(Biostrings::readAAStringSet(annoTaleRVD_file,
                                                 seek.first.rec = TRUE,
                                                 use.names = TRUE),
                     silent = TRUE)
    prot_parts <- try(Biostrings::readAAStringSet(prot_parts_files), silent = TRUE)

    if (any(
      inherits(checkAnnoTale, "try-error"), # in case annotale does not work
      if (inherits(prot_parts, "try-error")) { # in case annotale does not return a prot_parts file or if it is empty.
        file.exists(prot_parts_files) && file.remove(prot_parts_files)
        TRUE
      } else {
        if (file.exists(prot_parts_files) && length(prot_parts) == 0L) {
          file.remove(prot_parts_files) # should also return TRUE
        }
      },
      if (inherits(seqOfRVDs, "try-error")) { # in case annotale works but cannot find rvds or rvd seq file is empty.
        file.exists(annoTaleRVD_file) && file.remove(annoTaleRVD_file)
        TRUE
      } else {
        if (Biostrings::width(seqOfRVDs) == 0) file.remove(annoTaleRVD_file) # should also return TRUE
      }
    )) {
      messages <<- c(messages,
                     (m <- glue::glue("Annotale failed to parse TALE domains for {talOrfID}.")))
      cli::cli_warn(m)
      return(list(rvds = Biostrings::AAStringSet(), domains = data.frame()))
    }
    names(seqOfRVDs) <- talOrfID

    ## domains report
    stops <- Biostrings::vcountPattern("*", prot_parts)
    domainsReport <- tibble::tibble(
      "array_id" = talOrfID,
      "seqnames" = S4Vectors::mcols(by_array)$OriginalSubjectName[S4Vectors::mcols(by_array)$array_id == talOrfID],
      "query_name" = gsub("(.+\\: )|( \\d+)", "", names(prot_parts)),
      "codon_count" = Biostrings::width(prot_parts) - stops
    )

    list(rvds = seqOfRVDs, domains = domainsReport)
  })

  # Each element carries an AAStringSet and a data frame. That pairing used to
  # be an exported S4 class whose only purpose was to let sapply() return both
  # at once; a list does it without putting an implementation detail in the
  # package's API.
  list(rvds = unlist(Biostrings::AAStringSetList(lapply(out, `[[`, "rvds"))),
       domains = do.call(rbind, lapply(out, `[[`, "domains")),
       messages = messages)
}


#' Align the N- and C-termini of every array found
#'
#' The termini are the parts of a TALE that do not vary with its target: the
#' repeats differ from one TALE to the next by design, while the flanking
#' regions are near-identical across a strain's TALEs. Aligning them across
#' arrays is therefore a way to see whether a predicted array is a plausible
#' TALE at all -- a terminus that does not align with the others is a sign
#' the prediction is wrong, or that the sequence is.
#'
#' Nothing downstream reads the alignments; they are written for a person to
#' look at. Fewer than two arrays makes an alignment meaningless, and that
#' case is reported rather than attempted.
#'
#' @param annotale_dir Directory holding one AnnoTALE output per array.
#' @param output_dir Where the HTML goes.
#' @param type \code{"DNA"} or \code{"AA"}.
#' @return The collected parts, one element per terminus.
#' @noRd
.telltale_align_termini <- function(annotale_dir, output_dir, type = c("DNA", "AA")) {
  type <- match.arg(type)
  spec <- switch(type,
    DNA = list(parts = "TALE_DNA_parts.fasta", label = "DNA",
               read = Biostrings::readDNAStringSet, setlist = Biostrings::DNAStringSetList,
               suffix = "DNAAlignment.html"),
    AA  = list(parts = "TALE_Protein_parts.fasta", label = "protein",
               read = Biostrings::readAAStringSet, setlist = Biostrings::AAStringSetList,
               suffix = "AAAlignment.html"))

  partFiles <- list.files(annotale_dir, spec$parts, recursive = TRUE, full.names = TRUE)

  sapply(c("N-terminus", "C-terminus"), function(part) {
    allpart <- sapply(partFiles, function(p) {
      allpart <- spec$read(p, seek.first.rec = TRUE)
      onepart <- allpart[grepl(part, names(allpart))]
      names(onepart) <- basename(dirname(p))   # the ROI this part came from
      onepart
    }, simplify = "array", USE.NAMES = FALSE) %>%
      spec$setlist() %>%
      unlist()

    if (length(allpart) > 1) {
      alignment <- DECIPHER::AlignSeqs(allpart, verbose = FALSE)
      DECIPHER::BrowseSeqs(alignment,
                           htmlFile = file.path(output_dir, glue::glue("{part}{spec$suffix}")),
                           openURL = FALSE, colWidth = 120)
    } else {
      cli::cli_warn("Skipping {part} TALE {spec$label} regions alignment because the input sequence has less than 2 putative TALEs.")
    }
    allpart
  }, USE.NAMES = TRUE)
}


#' Finish the RVD strings and attach them to the arrays
#'
#' AnnoTALE reports the RVDs of an array as a dash-separated string. Three
#' things are done to it here.
#'
#' A lowercase letter in an RVD is AnnoTALE's way of flagging a repeat whose
#' length departs from the canonical ~34 aa. Such an array is marked
#' \code{aberrantRepeat}, because an aberrant repeat changes how the array
#' should be read and is worth knowing about before the RVDs are used to
#' predict targets.
#'
#' The separator becomes \code{rvd_sep}, whatever the caller asked for.
#'
#' Finally, with \code{extremity_codes}, each string is bracketed by codes
#' standing for the termini, so that a string of RVDs and a string of repeat
#' codes describe the same number of parts. Where a terminus was detected the
#' code names it; where the array simply ends without one, the code is
#' \code{XXXXX} -- a terminus is presumed present but was not identified,
#' which is a different statement from its absence.
#'
#' @param rvds The RVD strings as AnnoTALE reported them.
#' @param by_array The grouped hits; gains \code{SeqOfRVD} and
#'   \code{aberrantRepeat}.
#' @param hmm What \code{.telltale_hmm_profiles()} returned, for recognising
#'   which termini a given array actually has.
#' @param rvd_sep Separator between RVDs.
#' @param extremity_codes Whether to bracket with terminus codes.
#' @return A list of the finished \code{rvds} and the updated
#'   \code{by_array}.
#' @noRd
.telltale_finish_rvd_strings <- function(rvds, by_array, hmm, rvd_sep,
                                         extremity_codes) {
  # a lowercase letter marks a repeat of non-canonical length
  aberrantRepeat <- sapply(rvds, function(s) {
    ifelse(length(s) > 0, grepl("[a-z]", s), NA)
  })

  rvds <- gsub("\\-", rvd_sep, rvds) %>% Biostrings::AAStringSet()

  if (extremity_codes) {
    # This is necessary for other tantale utilities that can operate on 'full'
    # domains sequences, ie downstream of distal, for TALE domains sequences
    # alignments.
    anchors <- tales_anchor_codes()   # NTERM, CTERM, XXXXX
    for (s in names(rvds)) {
      present <- as.character(by_array[[s]]$query_name)
      rvds[s] <- paste(ifelse(hmm$nterm %in% present, anchors[1], anchors[3]),
                       rvds[s], sep = rvd_sep)
      rvds[s] <- paste(rvds[s],
                       ifelse(hmm$cterm %in% present, anchors[2], anchors[3]),
                       sep = rvd_sep)
    }
  }

  S4Vectors::mcols(by_array) <- merge(
    S4Vectors::mcols(by_array),
    data.frame(SeqOfRVD = rvds, aberrantRepeat = aberrantRepeat,
               array_id = names(rvds)),
    by = "array_id", all.x = TRUE)
  S4Vectors::mcols(by_array)$SeqOfRVD[is.na(S4Vectors::mcols(by_array)$SeqOfRVD)] <- ""

  list(rvds = rvds, by_array = by_array)
}


#' Write the run's reports
#'
#' Three tables and two GFFs, each answering a different question about the
#' same run: \code{hitsReport} one row per domain hit, \code{domainsReport}
#' one row per part AnnoTALE named, \code{arrayReport} one row per array with
#' its RVD sequence. The GFFs put the same ranges where a genome browser can
#' show them against the original sequence.
#'
#' The RVD fasta is the file the rest of the package reads: it is what
#' \code{\link{tales_from_telltale}} and the target predictors start from.
#' Arrays that yielded no RVDs are left out of it rather than written empty.
#'
#' @param by_array The grouped hits, carrying the per-array metadata.
#' @param domains_report What AnnoTALE reported, per part.
#' @param unmerged The hits before merging, or \code{NULL}.
#' @param paths What \code{.telltale_paths()} returned.
#' @return The array report, which the run log also summarises.
#' @noRd
.telltale_write_reports <- function(by_array, domains_report, unmerged, paths) {
  ## Report on the hmmer hits, both as tab-delimited and as gff
  hitsReport <- lapply(as.list(by_array, use.names = TRUE),
                       function(gr) tibble::as_tibble(as.data.frame(gr))) %>%
    dplyr::bind_rows(.id = "array_id")
  readr::write_tsv(x = hitsReport, file = paths$hits_report)
  .hits_report_to_gff(paths$hits_report) # saving to gff format
  readr::write_tsv(x = domains_report, file = paths$domains_report)

  ##   Report with info on arrays, including the seq of RVD
  arrayReport <- as.data.frame(
    S4Vectors::mcols(by_array)[
      order(S4Vectors::mcols(by_array)$OriginalSubjectName,
            S4Vectors::mcols(by_array)$NumberOfHits), ]
  )
  readr::write_tsv(x = arrayReport, file = paths$array_report)

  ## Write a gff with everything collated
  allGR <- c(unlist(by_array),
             GenomicRanges::makeGRangesFromDataFrame(
               S4Vectors::mcols(by_array),
               seqnames.field = "OriginalSubjectName",
               keep.extra.columns = TRUE),
             # the unmerged hits are included only when merging happened, so
             # that a merged range can be compared against what went into it
             unmerged)
  rtracklayer::export.gff3(allGR, paths$all_ranges_gff)

  ## Write a fasta file of the seqs of RVDs
  rvds <- Biostrings::BStringSet(S4Vectors::mcols(by_array)$SeqOfRVD)
  names(rvds) <- S4Vectors::mcols(by_array)$array_id
  rvds <- rvds[!Biostrings::width(rvds) == 0]
  Biostrings::writeXStringSet(x = rvds, paths$rvd_sequences)

  arrayReport
}


#' Add the per-array measures that come from the ORF and the termini
#'
#' Four numbers per array, all of them ways of asking "does this look like a
#' whole TALE?": the length of each terminus, the length of the longest ORF,
#' and what fraction of the array that ORF covers. An array whose ORF covers
#' only part of its length is one where something -- a frameshift, a
#' premature stop, a mis-called region -- interrupts the coding sequence.
#'
#' @param by_array The grouped hits.
#' @param ends_aa The terminus sequences, from
#'   \code{.telltale_align_termini()}.
#' @param full_orf The longest ORF per array.
#' @param array_seqs The extended array sequences the ORFs were found in.
#' @return \code{by_array}, with the measures in its metadata.
#' @noRd
.telltale_add_array_measures <- function(by_array, ends_aa, full_orf, array_seqs) {
  endsAAlength <- lapply(names(ends_aa), function(e) {
    stringset <- ends_aa[e] %>% Biostrings::AAStringSetList(., use.names = FALSE) %>% unlist()
    df <- data.frame(names(stringset), BiocGenerics::width(stringset))
    colnames(df) <- c("array_id", paste0(e, "AAlength"))
    df
  })
  S4Vectors::mcols(by_array) <- merge(S4Vectors::mcols(by_array),
                                      do.call(merge, endsAAlength),
                                      by = "array_id", all.x = TRUE)

  moreInfo <- merge(
    S4Vectors::mcols(by_array),
    data.frame(array_id = names(full_orf),
               LongestOrfLength = Biostrings::nchar(full_orf),
               OrfCovOverArrayLength = round(100 * Biostrings::nchar(full_orf) /
                                               GenomicRanges::width(array_seqs[names(full_orf)])),
               LongestORFSeq = full_orf),
    by = "array_id", all.x = TRUE, sort = FALSE)
  rownames(moreInfo) <- moreInfo$array_id
  S4Vectors::mcols(by_array) <- moreInfo[rownames(S4Vectors::mcols(by_array)), ]
  by_array
}


#' Distances between neighbouring arrays on the same sequence
#'
#' Summarised in the run log only. A short gap between two arrays can mean
#' they are really one TALE whose middle was missed, so the distribution is
#' worth a glance when a strain's TALE count looks wrong.
#'
#' @param arrays The array ranges.
#' @return A list of the gaps under 500 nt and their quartiles.
#' @noRd
.telltale_array_gaps <- function(arrays) {
  bySeqlevel <- split(arrays, GenomicRanges::seqnames(arrays))
  gaps <- sapply(bySeqlevel, function(x) {
    t(as.data.frame(GenomicRanges::distanceToNearest(x)))[3, ]
  })
  gaps <- unlist(gaps)
  gaps <- gaps[!is.na(gaps)]
  below500 <- gaps[gaps <= 500]
  list(below_500 = below500,
       quartiles = quantile(below500, probs = c(0.25, 0.50, 0.75)))
}


#' Put the domain hits on the genome, with their sequences
#'
#' Turns the tabular nhmmer output into ranges on the original sequences, and
#' records each hit's DNA alongside two counts derived from it.
#'
#' \code{frameshift_count} is the hit's length modulo three. A domain that
#' codes for protein should be a whole number of codons, so a non-zero value
#' says the hit is not in frame -- which is the signal frameshift correction
#' exists to act on.
#'
#' The sequence names are restored to the originals here: they were
#' simplified on the way in because spaces in fasta headers break the
#' downstream parsing.
#'
#' @param hits The filtered tabular hits.
#' @param subject_seqs The sequences that were searched.
#' @param seqlevels,seqinfo The original names and sequence information.
#' @return A \code{GRanges} of hits carrying their sequence.
#' @noRd
.telltale_hits_to_ranges <- function(hits, subject_seqs, seqlevels, seqinfo) {
  gr <- GenomicRanges::makeGRangesFromDataFrame(
    df = hits, keep.extra.columns = TRUE, seqnames.field = "target_name")
  names(gr) <- gr$hitID

  ## Updating seqinfo with original seqinfo from the sequences before renaming
  gr <- GenomeInfoDb::renameSeqlevels(gr, value = seqlevels)
  GenomeInfoDb::seqinfo(gr, pruning.mode = "coarse") <- seqinfo[GenomeInfoDb::seqlevels(gr)]
  gr
}


#' Attach each hit's DNA sequence and its codon arithmetic
#' @inheritParams .telltale_hits_to_ranges
#' @param gr Hits as ranges.
#' @return \code{gr}, with \code{seq}, \code{codon_count} and
#'   \code{frameshift_count} in its metadata.
#' @noRd
.telltale_add_hit_seqs <- function(gr, subject_seqs) {
  hitSeqs <- BSgenome::getSeq(subject_seqs, gr)
  S4Vectors::mcols(gr) %<>% cbind(
    data.frame(
      "seq" = as.character(hitSeqs),
      "codon_count" = Biostrings::nchar(hitSeqs) %/% 3,
      # not a whole number of codons: the hit is out of frame
      "frameshift_count" = Biostrings::nchar(hitSeqs) %% 3,
      row.names = NULL,
      check.rows = TRUE)
  )
  gr
}


#' Every file and directory a tell_tales() run writes
#'
#' Computed once, up front, so that the rest of the function reads as a
#' pipeline over data rather than over paths. These are the only values in
#' \code{tell_tales()} that are written once and then carried, unchanged,
#' through every stage to the end.
#'
#' The directories are created here too: a path this returns can be written
#' to without checking.
#'
#' @param output_dir Directory the run writes into.
#' @param correct_array Whether frameshift correction is on. The two
#'   correction alignment directories exist only then, so the corresponding
#'   entries are absent when it is \code{FALSE} rather than naming a
#'   directory that was never created.
#' @return A named list of absolute paths.
#' @noRd
.telltale_paths <- function(output_dir, correct_array) {
  dir.create(output_dir, recursive = TRUE, mode = "755", showWarnings = FALSE)
  p <- list(
    output         = output_dir,
    # one directory per region of interest is created under this one, later,
    # by the AnnoTALE stage
    annotale       = file.path(output_dir, "annotale"),
    ## Tabular file reporting on individual TALE domain hits
    hits_report    = file.path(output_dir, "hitsReport.tsv"),
    domains_report = file.path(output_dir, "domainsReport.tsv"),
    ## Tabular file reporting on putative TALEs (contiguous arrays of domain hits)
    array_report   = file.path(output_dir, "arrayReport.tsv"),
    ## Gff file with all the identified domains and arrays and their associated data
    all_ranges_gff = file.path(output_dir, "allRanges.gff"),
    ## fasta of tal orfs that have rvds, and of those predicted not to
    tale_orf_fasta = file.path(output_dir, "putativeTalOrf.fasta"),
    pseudo_tal     = file.path(output_dir, "pseudoTalCds.fasta"),
    ## the selected seqs of RVDs, without the separator
    rvd_sequences  = file.path(output_dir, "rvdSequences.fas"),
    ## the three HMM profiles concatenated, which is what nhmmer is given
    merged_hmm     = file.path(output_dir, "TALE_CDS_all_diagnostic_regions_hmmfile.out"),
    hmmer_search   = file.path(output_dir, "hmmerSearchOut.txt"),
    hmmer_readable = file.path(output_dir, "nhmmerHumanReadableOutputOfLastRun.txt"),
    ## logging info and some general analysis measures
    log            = file.path(output_dir, "tell_tales.log")
  )
  dir.create(p$annotale)
  if (correct_array) {
    p$correction_dna <- file.path(output_dir, "CorrectionAlignmentDNA")
    p$correction_aa  <- file.path(output_dir, "CorrectionAlignmentAA")
    dir.create(p$correction_dna, showWarnings = FALSE)
    dir.create(p$correction_aa, showWarnings = FALSE)
  }
  p
}


#' Search and report on the features of TALE protein domains potentially encoded
#' in subject DNA sequences
#'
#' \code{tell_tales} has been primarily written to report on 'corrected' TALE RVD
#' sequences in indels prone, noisy DNA sequences (suboptimally polished genomes
#' assembly, raw reads of long read sequencing technologies [eg PacBio, ONT])
#' that would otherwise be missed by conventional tools (eg AnnoTALE).
#'
#' The approach is first to use \href{http://hmmer.org/}{HMMER} to find and
#' categorize regions in the input DNA sequence that are related to the coding
#' sequence of canonical TALE protein domains (N-Term, repeats, C-term). Hits
#' that are (nearly [see the min_gap parameter]) adjacent are grouped in
#' "taleArrays" which are considered as potential tal genes.
#'
#'
#' If the \code{correct_array} parameter is turned off, the longest
#' predicted open reading frame (+extend_len) for each talArray is fed to
#' \href{http://www.jstacs.de/index.php/AnnoTALE}{AnnoTALE} to detect TALE
#' domains in the predicted translation product. The Results should hence be
#' very similar to what would be obtained with AnnoTALE, plus many additional
#' informative output files such as tabular reports.
#'
#' If \code{correct_array} is turned on, these talearrays are passed to the
#' \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}} function that
#' attemps to 'correct' potential frameshifts in the taleArray sequences. This
#' conveniently removes many artefactual indels but bear in mind that this may
#' also \strong{erroneously} 'correct' genuine frame shifts which can be highly
#' relevant especially for truncTALEs or iTALES. The resulting 'corrected'
#' taleArray open reading frames are then passed to AnnoTALE.
#'
#' Note that occasionally, when a putative open reading frame does not encode a
#' canonical TALE protein (early frame shift, incomplete ORF, etc...), the
#' "analyze" module of AnnoTALE outputs DNA parts but no protein parts and/or
#' RVD sequence. This should be detected and reported in the tell_tales log.
#'
#'
#' @param subject_file Fasta file with DNA sequence(s) to be searched for the
#'   presence of TALE coding sequences (CDS).
#' @param output_dir Path of the output directory. If not specified, results will
#'   be written to current working folder.
#' @param hmm_dir Specify the path to a folder holding the hmmfiles if you
#'   do not want to use the ones provided with tantale.
#' @param nterm_min_score Minimal nhmmer score cut_off value to
#'   consider the hit as genuine
#' @param repeat_min_score Minimal nhmmer score cut_off value to consider
#'   the hit as genuine
#' @param cterm_min_score Minimal nhmmer score cut_off value to
#'   consider the hit as genuine
#' @param min_domain_hits Minimum number of nhmmer hits for a subject
#'   sequence to be reported as having TALE diagnostic regions. This is a way to
#'   simplify output a little by getting ride of uninformative sequences
#' @param merge_hits Perform overlapping hits merging per domain type. Should not
#'   be modified.
#' @param min_gap Minimum gap in base pairs between two tale domain hits for
#'   them to be considered distinct. If the length of the gap is below this
#'   value, domains are considered "contiguous" and grouped in the same array.
#' @param extremity_codes Set this to \code{FALSE} if you do not want the
#'   N- and C-TREM anchor codes in the output sequences of RVD
#' @param rvd_sep Symbol acting as a separator in RVD sequences
#' @param hmmer_path Specify the path to a directory holding the HMMER executable
#'   if you do not want to use the ones provided with tantale.
#' @param extend_len number of nucleotides to extend in 3'-end at the tal
#'   ORF prediction stage.
#' @param correct_array True or False
#' @param correction_ref Reference AA sequences for tal array
#'   predicted ORF correction if you do not want to use the ones provided with
#'   tantale.
#' @param frameshift This is an internal parameter of the
#'   \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}} function. The
#'   default is 11 and fiddle with this at your own risk...
#' @param ... Additional parameters for the
#'   \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}} function.
#' @return This functions has only side effects (writing files, mostly).
#'   However, if everything ran smoothly, it will invisibly return the path of
#'   the directory where output files were written.
#'
#'
#'   List of output files:
#'   \itemize{
#'   \item allRanges.gff: gff file of all Tal arrays detected by HMMer
#'   \item arrayReport.tsv: report of all Tal arrays. In the arrayReport.tsv,
#'   column \emph{predicted_dels_count}/\emph{predicted_ins_count} shows the
#'   number of putative deletions/insertions in the raw sequences that have been
#'   corrected in the corrected sequences with the
#'   function \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}}.
#'   \item hitsReport.tsv: report of all hits detected by HMMer
#'   \item hitsReport.gff: gff file of all hits detected by HMMer
#'   \item domainsReport.tsv: report of all Tal amino acid domains detected by AnnoTALE analyze
#'   \item putativeTalOrf.fasta: Tal putative ORFs
#'   \item pseudoTalCds.fasta: pseudo Tal CDS, putative Tal array ORFs detected by HMMer for whch
#'    AnnoTALE analyze failed to find RVD(s).
#'   \item rvdSequences.fas: Sequence of RVDs (separated by rvd_sep) predicted to be encoded in the Tal array
#'    ORFs by AnnoTALE. Note that if extremity_codes is \code{TRUE} (by default),
#'    the N- and C-TREM anchor codes will be appended at the beginning and end of the sequences
#'    if the corresponding domain coding sequence was wound by HMMer at the DNA level.
#'    If no such HMMer hits were found, the "XXXXX" string will be appended
#'    to denote that AA sequences outside of the RVD array are likely to be atypical.
#'   \item C-terminusAAAlignment.html: protein alignment of all C-termini
#'   \item C-terminusDNAAlignment.html: DNA alignment of all C-termini
#'   \item N-terminusAAAlignment.html: protein alignment of all N-termini
#'   \item N-terminusDNAAlignment.html: DNA alignment of all N-termini
#'   \item TALE_CDS_all_diagnostic_regions_hmmfile.out: HMMER profile used for tale cds search.
#'   \item hmmerSearchOut.txt: ignore
#'   \item nhmmerHumanReadableOutputOfLastRun.txt: primary HMMER output file.
#'   \item tell_tales.log: a log file
#'   \item annotale folder: folder containing result of AnnoTALE analyze for all Tal arrays
#'   \item CorrectionAlignmentAA folder: folder containing protein alignment of Tal array detected by HMMer and corrected Tal array if \code{correct_array} = TRUE
#'   \item CorrectionAlignmentDNA folder: folder containing DNA alignment of Tal array detected by HMMer and corrected Tal array if \code{correct_array = TRUE}
#'   }
#' @export
#' @family TALE discovery
tell_tales <- function(
  subject_file,
  output_dir = getwd(),
  hmm_dir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T),
  nterm_min_score = 300,
  repeat_min_score = 20,
  cterm_min_score = 200,
  min_domain_hits = 4,
  merge_hits = TRUE,
  min_gap = 35,
  extremity_codes = TRUE,
  rvd_sep = "-",
  hmmer_path = system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = T),
  extend_len = 300,
  correct_array = FALSE,
  correction_ref = system.file("extdata", "decipher_ref_tales_aa.fa.gz", package = "tantale", mustWork = T),
  frameshift = -11,
  ...
) {

  # @param taleArrayStartAnchorCode This scalar character vector will symbolize a
  #   TALE N-TERM CDS hit in the RVD sequence
  # @param taleArrayEndAnchorCode This scalar character vector will symbolize a
  #   TALE C-TERM CDS hit in the RVD sequence

  #### TODO ####
  # Add an ooptional argument that olds the circularity status of molecules in genome
  # update the seqinfo objects accrodingly. This may solve some issues if a tale is located
  # at the junction of extremities in a circular molecule.
  # Could also implement an autotmated mecanisms "findCircular" that would
  # find a temr (eg 'circular') in the sequence title and act accrodingly.
  
  
  ####   Paths of output files   ####
  paths <- .telltale_paths(output_dir, correct_array)
  
  
  ####   Checks for parameters and other things   ####

  ## Deal with spaces in sequence names because this messes up parsing of HMMER output
  cli::cli_inform("HMMER is very picky about forbiden characters in sequence name. Renaming sequences in {subject_file}.")
  originalSeqs <- Biostrings::readDNAStringSet(filepath = subject_file)
  Rsamtools::indexFa(subject_file)
  originalSeqInfo <- Rsamtools::seqinfo(Rsamtools::FaFile(subject_file))
  originalSeqlevels <- names(originalSeqs)
  foolproofSeqlevels <- paste0("seq", 1:length(originalSeqlevels))
  names(originalSeqlevels) <- foolproofSeqlevels
  names(originalSeqs) <- foolproofSeqlevels
  cli::cli_inform(paste0("Original seq names : {glue::glue_collapse(originalSeqlevels, sep = ' ; ')}"))
  cli::cli_inform(paste0("Dummy seq names : {glue::glue_collapse(names(originalSeqlevels), sep = ' ; ')}"))
  subject_file <- tempfile()
  Biostrings::writeXStringSet(originalSeqs, filepath = subject_file)
  
  ####   Read the TALE profile HMMs and concatenate them for nhmmer   ####
  ## TODO come up with a mechanism for the user to be able to provide the FULL PATH
  ## to custom hmm !!! hmm_dir parameter is useless unless custom hmm are named
  ## as specified in .telltale_hmm_profiles()
  hmm <- .telltale_hmm_profiles(hmm_dir, paths$merged_hmm)
  
  ####   Find the TALE domain CDS hits   #####
  nhmmerTabularOutput <- .telltale_find_domain_hits(
    subject_file = subject_file, hmm = hmm, paths = paths,
    hmmer_path = hmmer_path,
    nterm_min_score = nterm_min_score, repeat_min_score = repeat_min_score,
    cterm_min_score = cterm_min_score, min_domain_hits = min_domain_hits)
  # Every stage below assumes at least one hit.
  if (is.null(nhmmerTabularOutput)) return(invisible(output_dir))

  #####   Put the hits on the genome   ####
  ## Load in R the DNA sequences that are queried for TALE CDS
  subjectDNASequences <- Biostrings::readDNAStringSet(filepath = subject_file)
  names(subjectDNASequences) <- originalSeqlevels[match(names(subjectDNASequences), names(originalSeqlevels))]

  nhmmerOutputGRBeforeMerge <- .telltale_hits_to_ranges(
    nhmmerTabularOutput, subjectDNASequences, originalSeqlevels, originalSeqInfo)

  #####   Domain-wise merge of overlapping hits  #####
  nhmmerOutputGR <- if (merge_hits) {
    .telltale_merge_overlapping_hits(nhmmerOutputGRBeforeMerge)
  } else {
    nhmmerOutputGRBeforeMerge
  }

  #####   Record each hit's DNA sequence   #####
  nhmmerOutputGR <- .telltale_add_hit_seqs(nhmmerOutputGR, subjectDNASequences)

  #####   Group (nearly) adjacent hits in "TALE array" regions   #####
  grouped <- .telltale_group_arrays(nhmmerOutputGR, min_gap, subjectDNASequences, hmm)
  arraysGR <- grouped$arrays
  hitsByArraysLst <- grouped$by_array

  #####   Extend DNA Tal arrays   #####
  ## Extract the genomic sequence of arrays +-bp on the borders
  completeArraysGR <- arraysGR
  extdCompleteArraysGR <- GenomicRanges::resize(completeArraysGR,
                                                width = GenomicRanges::width(completeArraysGR) + extend_len,
                                                fix = "start", ignore.strand = FALSE) %>%
    GenomicRanges::trim(use.names = TRUE)
  extdCompleteArraysSeqs <- BSgenome::getSeq(subjectDNASequences, extdCompleteArraysGR)
  
  

  orfResult <- .telltale_array_orfs(
    array_seqs = extdCompleteArraysSeqs, by_array = hitsByArraysLst,
    correct_array = correct_array, correction_ref = correction_ref,
    frameshift = frameshift, paths = paths, ...)
  TalOrfForAnnoTALE <- orfResult$orf
  fullTalOrf <- orfResult$full_orf
  hitsByArraysLst <- orfResult$by_array

  #### AnnoTALE analyze on tal ORFs  ####
  # Shall we also run the predict stage of annotale? May be it will do a better job at
  # finding orf and/or filtering out "pseudo tales" because some times analyse output a RVD from
  # a domain that does not look like a repeat....
  annotale <- .telltale_run_annotale(TalOrfForAnnoTALE, hitsByArraysLst, paths$annotale)
  seqsOfRVDs <- annotale$rvds
  domainsReport <- annotale$domains
  annoTaleMessages <- annotale$messages

  #### TODO  #####
  # The exact content of the files below needs to be reassesed and 
  # we need to determine if this is really what we want.
  # save tals orfs that have rvds
  Biostrings::writeXStringSet(fullTalOrf[names(fullTalOrf) %in% names(seqsOfRVDs)], paths$tale_orf_fasta)
  
  # save tals that DO NOT have rvds
  Biostrings::writeXStringSet(extdCompleteArraysSeqs[!names(extdCompleteArraysSeqs) %in% names(seqsOfRVDs)], paths$pseudo_tal)
  
  rvdResult <- .telltale_finish_rvd_strings(seqsOfRVDs, hitsByArraysLst, hmm,
                                            rvd_sep, extremity_codes)
  seqsOfRVDs <- rvdResult$rvds
  hitsByArraysLst <- rvdResult$by_array

  #### Align N-term and C-term ####
  .telltale_align_termini(paths$annotale, output_dir, type = "DNA")
  endsAA <- .telltale_align_termini(paths$annotale, output_dir, type = "AA")

  #### Per-array measures from the ORF and the termini ####
  hitsByArraysLst <- .telltale_add_array_measures(
    hitsByArraysLst, endsAA, fullTalOrf, extdCompleteArraysSeqs)

  ####   Gaps between neighbouring arrays, for the log   ####
  arrayGaps <- .telltale_array_gaps(arraysGR)
  gaplengthBetweenHitDomainsbelow500 <- arrayGaps$below_500
  quartilesGapLength <- arrayGaps$quartiles

  ####   Write tabulated reports and sequence files   #####
  arrayReport <- .telltale_write_reports(
    by_array = hitsByArraysLst, domains_report = domainsReport,
    unmerged = if (merge_hits) nhmmerOutputGRBeforeMerge else NULL,
    paths = paths)

  ####   Generate info messages and log file about the analysis   #####
  .telltale_log(
    params = list(subject_file = subject_file, output_dir = output_dir,
                  nterm_min_score = nterm_min_score,
                  repeat_min_score = repeat_min_score,
                  cterm_min_score = cterm_min_score,
                  min_domain_hits = min_domain_hits, merge_hits = merge_hits,
                  min_gap = min_gap, extend_len = extend_len,
                  correct_array = correct_array, correction_ref = correction_ref,
                  frameshift = frameshift),
    log_file = paths$log, hmm = hmm, subject_seqs = subjectDNASequences,
    arrays = arraysGR, by_array = hitsByArraysLst, array_report = arrayReport,
    gaps_below_500 = gaplengthBetweenHitDomainsbelow500,
    gap_quartiles = quartilesGapLength,
    annotale_messages = annoTaleMessages)

  return(invisible(output_dir))
}


#### Helpers for tell_tales() ####
#
# Formerly R/tellTale_utilities.R. Folded in here because tell_tales() is the
# only caller of every one of them.

.get_hmmer <- function() {
  pathOfHmmerBinsDir <- system.file("tools", "hmmer-3.3", "bin", package = "tantale", mustWork = TRUE)
  return(pathOfHmmerBinsDir)
}


.check_hmmer <- function(hmmer_path) {
  cmd <- file.path(hmmer_path, "hmmsearch -h | grep \"^#\"")
  if (system(command = cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)) {
    stop("HMMER is not in PATH. Follow instructions at http://hmmer.org/documentation.html to install it.")
  } else {
    out <- system(command = cmd,intern = TRUE)
    cli::cli_inform(gsub("^#[ ]?", "", out[2:3]))
  }
}






.run_nhmmer_search <- function(hmmer_path = NULL, subject_file, hmm_file, search_out_file, readable_out_file) {
  if (is.null(hmmer_path)) hmmer_path <- .get_hmmer()
  .check_hmmer(hmmer_path)
  searchCmd <- paste(file.path(hmmer_path, "nhmmer"),
                     "--tblout",
                     search_out_file,
                     hmm_file,
                     subject_file,
                     ">",
                     readable_out_file,
                     sep = " "
  )
  system(command = searchCmd, ignore.stderr = FALSE, intern = TRUE)
}




.correction_tibble <- function(indels) {
  indelsTble <- lapply(indels, function(lst) {
    info <- tibble::tibble()
    colnames(info) <- c("variable","value")
    for (talOrfID in 1:length(lst)) {
      element <- lst[talOrfID]
      if (length(unlist(element)) == 0L) {
        next()
      } else {
        info %<>% dplyr::bind_rows(tibble::tibble(variable = names(element), value = as.numeric(unlist(element))))
      }
    }
    return(info)
  }
  ) %>% dplyr::bind_rows(.id = "Seq")
  return(indelsTble)
} 




.hits_report_to_gff <- function(f = "hitsReport.csv") {
  # Convert the info contained in a HitReport file into a GFF file for display by
  # a genome viewer.
  # The f parameter corresponds to the path to a hitsReport file.
  # Read the file as a data.frame
  hitsReport <- read.delim(f)
  # Create a GenomicRange that will be converted.
  hitsGR <- GenomicRanges::makeGRangesFromDataFrame(hitsReport, keep.extra.columns=TRUE)
  # Write a gff3 file to disk with this info.
  rtracklayer::export.gff3(hitsGR,
                           con = file.path(dirname(f), paste0(sub("\\..*$", "", basename(f)), ".gff"))
  )
  
}
