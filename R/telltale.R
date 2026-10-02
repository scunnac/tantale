
# subject_file = system.file("extdata", "bai3_sample_tal_regions.fasta", package = "tantale", mustWork = T)
# output_dir = tempdir(check = TRUE)
# hmm_dir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T)
# hmmer_path = NULL   # NULL -> the tantale conda env
# correct_array = TRUE
# correction_ref = system.file("extdata", "tale_correction_ref.fa.gz", package = "tantale", mustWork = T)
# frameshift = -11
# nterm_min_score = 300
# repeat_min_score = 20
# cterm_min_score = 200
# min_dna_hits = 4
# merge_hits = TRUE
# min_gap = 35
# taleArrayStartAnchorCode = "NTERM"
# taleArrayEndAnchorCode = "CTERM"
# extremity_codes = TRUE
# rvd_sep = "-"
# extend_len = 300
# ... = NULL

# subject_file = system.file("extdata", "bai3_sample_tal_genomic_regions.fasta", package = "tantale", mustWork = T)
# output_dir = file.path(tempdir(), gsub("(\\.fasta)|(\\.fa)|(\\.fna)|(\\.fsa)", "", basename(subject_file)))
# hmm_dir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T)
# hmmer_path = NULL   # NULL -> the tantale conda env
# correct_array = FALSE
# correction_ref = system.file("extdata", "tale_correction_ref.fa.gz", package = "tantale", mustWork = T)
# frameshift = -11
# nterm_min_score = 300
# repeat_min_score = 20
# cterm_min_score = 200
# min_dna_hits = 4
# merge_hits = TRUE
# min_gap = 35
# extremity_codes = TRUE
# rvd_sep = "-"
# extend_len = 300
# ... = NULL






#### Helpers for tell_tales() ####
#


.get_hmmer <- function() {
  # HMMER comes from the tantale conda environment rather than being shipped
  # inside the package. 3.3.2 there was checked against the 3.3 that used to
  # be bundled: identical hits, only the version banner in the raw output
  # differs (restructuring-notes.md 7.4).
  file.path(.tantale_env_prefix(), "bin")
}


.check_hmmer <- function(hmmer_path) {
  cmd <- file.path(hmmer_path, "hmmsearch -h | grep \"^#\"")
  if (system(command = cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)) {
    cli::cli_abort(
      c("HMMER is not on the {.envvar PATH}.",
        "i" = "{.run tantale_setup(install = TRUE)} installs it in the {.val tantale} environment.",
        "i" = "Or install it yourself: {.url http://hmmer.org/documentation.html}"),
      class = c("tantale_error_hmmer_missing", "tantale_error"))
  } else {
    out <- system(command = cmd,intern = TRUE)
    cli::cli_inform(gsub("^#[ ]?", "", out[2:3]))
  }
}

.run_nhmmer_search <- function(hmmer_path = NULL,
                               subject_file, hmm_file,
                               search_out_file, readable_out_file) {
  if (is.null(hmmer_path)) hmmer_path <- .get_hmmer()
  .check_hmmer(hmmer_path)
  searchCmd <- paste(shQuote(file.path(hmmer_path, "nhmmer")),
                     "--tblout",
                     shQuote(search_out_file),
                     shQuote(hmm_file),
                     shQuote(subject_file),
                     ">",
                     shQuote(readable_out_file),
                     sep = " "
  )
  .tantale_exec(searchCmd, what = "nHMMER search")
}



.hits_report_to_gff <- function(f = "hits_report.tsv") {
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

  lines <- lapply(files, function(f) readLines(con = f))
  names <- unlist(lapply(lines, function(x) {
    hmmName <- grep("NAME", x, perl = TRUE, value = TRUE)
    hmmName <- unlist(strsplit(hmmName, split = "\\s+"))
    if (length(hmmName) != 2) {
      cli::cli_abort(
        c("A profile HMM has spaces in its name.",
          "i" = "Remove them from the {.field NAME} tag in the HMM file."),
        class = c("tantale_error_hmm_name", "tantale_error"))
    }
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
#' @param hmmer_path Directory holding the \code{nhmmer} binary.
#'   \code{NULL}, the default, uses the \code{tantale} conda environment,
#'   creating it on first use.
#' @param nterm_min_score,repeat_min_score,cterm_min_score Per-domain score
#'   thresholds.
#' @param min_dna_hits A subject sequence is kept when it carries at least
#'   this many hits. Counts hits per subject sequence, not per TALE array --
#'   it is a cheap pre-filter that discards whole contigs carrying nothing but
#'   stray matches. Short *arrays* are filtered separately, after grouping.
#' @return The filtered table, or \code{NULL} when nothing survives -- which
#'   the caller must treat as "stop here", since every later stage assumes at
#'   least one hit.
#' @noRd
.telltale_find_domain_hits <- function(subject_file, hmm, paths, hmmer_path,
                                       nterm_min_score, repeat_min_score,
                                       cterm_min_score, min_dna_hits) {
  .run_nhmmer_search(hmmer_path = hmmer_path,
                     subject_file = subject_file,
                     hmm_file = paths$merged_hmm,
                     search_out_file = paths$hmmer_search,
                     readable_out_file = paths$hmmer_readable)

  hits <- try(read.table(paths$hmmer_search), silent = TRUE)
  if (inherits(hits, "try-error")) {
    cli::cli_warn("No TALE CDS hits found in {.file {subject_file}}.",
                  class = "tantale_warning_no_hits")
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
  ## Add a hit_id column
  hits$hit_id <- paste("DOM", sprintf("%05.0f", 1:nrow(hits)), sep = "_")

  ## nhmmer does not guarantee envfrom <= env_to (a reverse-strand hit can
  ## report them the other way round), and IRanges() requires start <= end.
  hits$start <- pmin(hits$envfrom, hits$env_to)
  hits$end   <- pmax(hits$envfrom, hits$env_to)
  rownames(hits) <- hits$hit_id

  ## Filter out target DNA sequences that have too few repeat CDSs
  ## NB: for the sake of consistency  it would be better just to filter out
  ## from any further consideration the ARRAYS shorter than a certain value (say 5).
  ## WHAT DO WE DO ABOUT THAT?
  perSubject <- dplyr::count(hits, target_name, sq_len, name = "V1")
  # >=, not >: the argument is documented as a minimum, and a sequence
  # carrying exactly that many hits used to be dropped.
  hits <- subset(hits, target_name %in% perSubject[perSubject$V1 >= min_dna_hits, "target_name"])
  hits <- droplevels(hits)
  if (nrow(hits) == 0L) {
    # Unguarded before: the run carried on and died several stages later
    # inside Bioconductor with "Rle of type 'NULL' is not supported".
    cli::cli_warn(c("No subject sequence carries at least {min_dna_hits} nhmmer DNA hit{?s}.",
                    "i" = "{.arg min_dna_hits} counts hits per subject sequence, not per TALE array.",
                    "x" = "Nothing left to analyse. Exitting..."))
    return(NULL)
  }

  hits
}


#' Merge hits of the same domain type that overlap each other
#'
#' nhmmer can report the same repeat twice, as two overlapping hits. Left
#' alone, that repeat is counted twice in \code{n_dna_hits} and by the
#' \code{min_array_length} filter, and listed twice in \code{hits_report.tsv}.
#' The RVDs are unaffected: AnnoTALE reads them from the array's ORF.
#'
#' Merging is done per domain type, never across types: an N-terminus hit
#' overlapping a repeat hit is a real feature of where one domain ends and the
#' next begins, not a duplicate.
#'
#' The identifiers of the hits that went into each merged range are kept in
#' \code{nhmmer_hit_id}, separated by \code{|}, so a merged range can be traced
#' back to the raw search output.
#'
#' @param gr Hits as a \code{GRanges}, with \code{query_name} naming the
#'   domain type and \code{hit_id} identifying each hit.
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
        ~ paste0(g[as.numeric(.x$queryHits)]$hit_id, collapse = "|")
      }) %>%
      unlist()
    reduced$nhmmer_hit_id <- formerIDs
    reduced
  }) %>%
    plyranges::bind_ranges(.id = "query_name")

  merged$hit_id <- paste("MDOM", sprintf("%05.0f", 1:length(merged)), sep = "_")
  names(merged) <- merged$hit_id
  merged
}


#' Group neighbouring domain hits into candidate TALE arrays
#'
#' A TALE is a run of domain hits close together on the same strand: an
#' N-terminus, a series of repeats, a C-terminus. Hits separated by less than
#' \code{min_gap} are taken to belong to the same array, and each array
#' becomes a region of interest, \code{ROI_*}.
#'
#' Hits of the same domain type should not overlap within an array (the
#' merge stage joins them), so any that still do are reported: they inflate
#' \code{n_dna_hits}. A terminus hit overlapping the adjacent repeat hit by a
#' few nucleotides is where one domain ends and the next begins, and is not
#' reported.
#'
#' @param gr Domain hits, merged.
#' @param min_gap Largest gap, in bases, still counted as contiguous.
#' @param subject_seqs The DNA the hits were found in, for extracting each
#'   array's sequence.
#' @param hmm What \code{.telltale_hmm_profiles()} returned; used to record
#'   whether an array carries all three domain types, and to recognise which
#'   hits are repeats.
#' @param min_array_length Arrays with fewer repeat units than this are
#'   dropped. \code{0} keeps everything.
#' @return A list of \code{arrays} (one range per array) and \code{by_array}
#'   (the hits, grouped, carrying the per-array metadata), or \code{NULL} if
#'   the length filter removed every array.
#' @noRd
.telltale_group_arrays <- function(gr, min_gap, subject_seqs, hmm,
                                   min_array_length = 0) {
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

  ## Drop arrays with too few repeats, before anything is computed about them.
  ## Repeats, not all hits: an array should not be penalised for having had a
  ## terminus missed, and the repeat count is what "array length" means for a
  ## TALE.
  if (min_array_length > 0) {
    repeatCount <- vapply(byArray, function(x) sum(as.character(x$query_name) == hmm$repeats),
                          integer(1))
    tooShort <- repeatCount < min_array_length
    if (any(tooShort)) {
      cli::cli_inform(c("Dropping {sum(tooShort)} array{?s} with fewer than {min_array_length} repeat{?s}.",
                        "i" = "Array{?s}: {.val {names(byArray)[tooShort]}}"))
      byArray <- byArray[!tooShort]
      arraysGR <- arraysGR[!tooShort]
    }
    if (length(byArray) == 0L) {
      cli::cli_warn(c("No TALE array has at least {min_array_length} repeat{?s}.",
                      "x" = "Nothing left to analyse. Exitting..."))
      return(NULL)
    }
    # renumber so the ROI ids stay contiguous
    names(byArray) <- paste("ROI", sprintf("%05.0f", seq_along(byArray)), sep = "_")
    names(arraysGR) <- names(byArray)
  }

  ## Overlaps between hits of the same domain type only (§36)
  doHitsOverlap <- vapply(byArray, function(x) {
    !all(GenomicRanges::isDisjoint(GenomicRanges::split(x, as.character(x$query_name))))
  }, logical(1))
  if (any(doHitsOverlap)) {
    cli::cli_warn(
      c("Some nhmmer hits of the same domain type overlap, so {.field n_dna_hits} counts these domains twice.",
        "i" = "Region{?s}: {.val {names(doHitsOverlap)[doHitsOverlap]}}",
        "i" = "{.fn tell_tales} merges such hits unless {.code merge_hits = FALSE}."),
      class = "tantale_warning_overlapping_hits")
  }

  ## Populate metadata about the elements of the list of arrays
  S4Vectors::mcols(byArray) <- S4Vectors::DataFrame(
    array_id = names(byArray),
    seqnames = sapply(byArray,
                      function(x) unique(as.character(GenomicRanges::seqnames(x)))),
    start = BiocGenerics::start(arraysGR),
    end = BiocGenerics::end(arraysGR),
    strand = BiocGenerics::strand(arraysGR),
    n_dna_hits = S4Vectors::elementNROWS(byArray),
    array_seq = BSgenome::getSeq(subject_seqs, arraysGR),
    nterm_dna_hit = sapply(byArray, function(x) hmm$nterm %in% as.character(x$query_name)),
    cterm_dna_hit = sapply(byArray, function(x) hmm$cterm %in% as.character(x$query_name))
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
  countsHitsByArrayDistri <- summary(S4Vectors::mcols(by_array)$n_dna_hits)
  ## Number of domains in arrays that display all domain types
  # completeArrayLengths <- subset(S4Vectors::mcols(by_array), nterm_dna_hit & cterm_dna_hit)$n_dna_hits
  
  
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
    paste("nterm_min_score:", params$nterm_min_score, sep = "\t"),
    paste("repeat_min_score:", params$repeat_min_score, sep = "\t"),
    paste("cterm_min_score:", params$cterm_min_score, sep = "\t"),
    paste("terminus_max_evalue:", params$terminus_max_evalue, sep = "\t"),
    paste("min_dna_hits:", params$min_dna_hits, sep = "\t"),
    paste("min_array_length:", params$min_array_length, sep = "\t"),
    paste("merge_hits:", params$merge_hits, sep = "\t"),
    paste("min_gap:", params$min_gap, sep = "\t"),
    paste("extend_len:", params$extend_len, sep = "\t"),
    paste("correct_array:", params$correct_array, sep = "\t"),
    paste("correction_ref:", params$correction_ref, sep = "\t"),
    paste("max_comparisons:",
          if (is.null(params$max_comparisons)) "all" else params$max_comparisons,
          sep = "\t"),
    paste("frameshift:", params$frameshift, sep = "\t"),
    
    "#__________Summary measures of TALE search outcome__________",
    paste("Number of analysed subject sequences :", length(subject_seqs), sep = "\t"),
    paste("Total number of TALE repeat DNA coding sequence motif hits found with the nhmmer approach:",
          numberOfRepeatHitsAfterFiltering, sep = "\t"),
    #paste("Total number of repeat HMM hits on the corresponding set of translated DNA hits:", sum(RVDtbl), sep = "\t"),
    
    paste("Total number of subject seqs with TALE motif hits after low hit number filtering:",
          length(GenomeInfoDb::seqlevels(arrays)), sep = "\t"),
    paste("Total number of distinct regions (repeat arrays) with adjacent TALE motifs :", nrow(array_report), sep = "\t"),
    paste("Number of arrays with nhmmer DNA hits for both termini:",
          sum(S4Vectors::mcols(by_array)$nterm_dna_hit & S4Vectors::mcols(by_array)$cterm_dna_hit),
          sep = "\t"),
    paste("Number of arrays whose AnnoTALE N-terminus matches the TALE N-terminal protein profile:",
          sum(S4Vectors::mcols(by_array)$nterm_aa_hit, na.rm = TRUE), sep = "\t"),
    paste("Number of arrays whose AnnoTALE C-terminus matches the TALE C-terminal protein profile:",
          sum(S4Vectors::mcols(by_array)$cterm_aa_hit, na.rm = TRUE), sep = "\t"),
    
    #paste("Total number of distinct types of RVD:", nrow(RVDtbl), sep = "\t"),
    
    paste("Minimum array length (number of nhmmer DNA hits):", min(array_report$n_dna_hits), sep = "\t"),
    paste("Maximum array length:", max(array_report$n_dna_hits), sep = "\t"),
    paste("Median array length:", median(array_report$n_dna_hits), sep = "\t"),
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
  
  # verbatim: cli_inform() would rewrap the lines into one paragraph and
  # read any {...} in a path as code to evaluate
  cli::cli_verbatim(txt)
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
                         htmlFile = file.path(dna_dir, glue::glue("correction_alignment_dna_{n}.html")),
                         openURL = FALSE, colWidth = 120)

    seqToAlignTranslated <- Biostrings::translate(seqToAlign, no.init.codon = TRUE,
                                                  if.fuzzy.codon = "solve")
    alignedSeqsTranslated <- DECIPHER::AlignSeqs(seqToAlignTranslated, verbose = FALSE)
    DECIPHER::BrowseSeqs(alignedSeqsTranslated,
                         htmlFile = file.path(aa_dir, glue::glue("correction_alignment_aa_{n}.html")),
                         openURL = FALSE, colWidth = 120)
  }
  invisible(NULL)
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
#' @param max_comparisons How many references each array may be aligned
#'   against. \code{NULL} means all of them.
#' @param paths What \code{.telltale_paths()} returned.
#' @param ... Passed to \code{DECIPHER::CorrectFrameshifts()}.
#' @return A list of \code{orf} (what AnnoTALE is given), \code{full_orf}
#'   (what is reported), and \code{by_array}, updated.
#' @noRd
.telltale_array_orfs <- function(array_seqs, by_array, correct_array,
                                 correction_ref, frameshift, paths,
                                 max_comparisons = NULL, ...) {
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
  if (is.null(max_comparisons)) max_comparisons <- length(AAref)
  correction <- DECIPHER::CorrectFrameshifts(array_seqs,
                                             AAref, type = "both",
                                             maxComparisons = max_comparisons,
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
#' @return \code{0}, invisibly. Stops with
#'   \code{tantale_error_annotale_failed} if AnnoTALE exits with a non-zero
#'   status; \code{.telltale_run_annotale()} catches that and skips the
#'   array.
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
    "java -jar ", shQuote(annotale_jar),
    " analyze ",
    " t=", shQuote(fasta_file),
    " outdir=", shQuote(output_dir)
  )
  .annotale_exec(comAnalyze, "analyze", quiet = TRUE)
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
    correctedTalOrfFile <- file.path(AnnotaleDir, "putative_tal_orf.fasta")
    Biostrings::writeXStringSet(TalOrf, correctedTalOrfFile)

    cli::cli_inform("Now running AnnoTALE analyze for {talOrfID}")
    checkAnnoTale <- try(.run_annotale_analyze(correctedTalOrfFile, AnnotaleDir), silent = TRUE)

    prot_parts_files <- file.path(AnnotaleDir, "TALE_Protein_parts.fasta")
    dna_parts_file <- file.path(AnnotaleDir, "TALE_DNA_parts.fasta")
    annoTaleRVD_file <- file.path(AnnotaleDir, "TALE_RVDs.fasta")
    seqOfRVDs <- try(Biostrings::readAAStringSet(annoTaleRVD_file,
                                                 seek.first.rec = TRUE,
                                                 use.names = TRUE),
                     silent = TRUE)
    prot_parts <- try(Biostrings::readAAStringSet(prot_parts_files), silent = TRUE)

    if (any(
      inherits(checkAnnoTale, "try-error"), # in case annotale does not work
      if (inherits(prot_parts, "try-error") || length(prot_parts) == 0L) {
        # No protein parts: AnnoTALE splits the DNA before translating, and
        # still writes the DNA parts when the protein cannot be split. Those
        # DNA parts are not trustworthy (a 9-nt "repeat" was seen), so both
        # files go, and the array stays out of the tales object.
        file.exists(prot_parts_files) && file.remove(prot_parts_files)
        file.exists(dna_parts_file) && file.remove(dna_parts_file)
        TRUE
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
      cli::cli_warn(m, parent = if (inherits(checkAnnoTale, "try-error")) {
        attr(checkAnnoTale, "condition")
      })
      return(list(rvds = Biostrings::AAStringSet(), domains = data.frame()))
    }
    names(seqOfRVDs) <- talOrfID

    ## domains report
    stops <- Biostrings::vcountPattern("*", prot_parts)
    domainsReport <- tibble::tibble(
      "array_id" = talOrfID,
      "seqnames" = S4Vectors::mcols(by_array)$seqnames[S4Vectors::mcols(by_array)$array_id == talOrfID],
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
#' @param paths The `.telltale_paths()` list, for where the HTML goes --
#'   `n_terminus_dna`/`c_terminus_dna`/`n_terminus_aa`/`c_terminus_aa`.
#' @param type \code{"DNA"} or \code{"AA"}.
#' @return The collected parts, one element per terminus.
#' @noRd
.telltale_align_termini <- function(annotale_dir, paths, type = c("DNA", "AA")) {
  type <- match.arg(type)
  spec <- switch(type,
    DNA = list(parts = "TALE_DNA_parts.fasta", label = "DNA",
               read = Biostrings::readDNAStringSet, setlist = Biostrings::DNAStringSetList,
               html = c(`N-terminus` = paths$n_terminus_dna, `C-terminus` = paths$c_terminus_dna)),
    AA  = list(parts = "TALE_Protein_parts.fasta", label = "protein",
               read = Biostrings::readAAStringSet, setlist = Biostrings::AAStringSetList,
               html = c(`N-terminus` = paths$n_terminus_aa, `C-terminus` = paths$c_terminus_aa)))

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
                           htmlFile = spec$html[[part]],
                           openURL = FALSE, colWidth = 120)
    } else {
      cli::cli_warn("Skipping {part} TALE {spec$label} regions alignment because the input sequence has less than 2 putative TALEs.")
    }
    allpart
  }, USE.NAMES = TRUE)
}


# How many profile positions a terminus match may stop short of the end that
# adjoins the repeats. Genuine termini in the article genomes reach within 2;
# the frameshifted N-termini of BAI3-1-1's raw assembly stop 138 short
# (ledger §42).
.terminus_max_profile_gap <- 10L


#' Does each terminus look like a canonical TALE terminal domain?
#'
#' AnnoTALE calls "N-terminus" whatever the ORF encodes upstream of the first
#' repeat, and "C-terminus" whatever it encodes downstream of the last one.
#' The ORF is taken to be translated, so these segments exist, but nothing
#' says they resemble the terminal domains of a TALE: an ORF that starts or
#' ends inside a frameshifted region yields segments of unrelated sequence.
#' Each segment is therefore searched with \code{hmmsearch} against the TALE
#' N- or C-terminal protein profile.
#'
#' A match must also reach the end of the profile that adjoins the repeats
#' (the last positions of the N-terminal profile, the first ones of the
#' C-terminal profile), within \code{max_profile_gap} positions. A frameshift
#' inside a terminus puts its repeat-side part in another reading frame, so
#' such a segment matches only up to the frameshift; a genuine truncated
#' terminus is shortened at its far end and still reaches the repeats (ledger
#' §42).
#'
#' Kept as a function because it is meant to serve every reader of AnnoTALE
#' output, \code{tell_tales()} first.
#'
#' @param termini A list of two \code{AAStringSet}s named \code{"N-terminus"}
#'   and \code{"C-terminus"}, each named by array, as
#'   \code{.telltale_align_termini(type = "AA")} returns them. A stop codon
#'   is removed before the search.
#' @param max_evalue A terminus matches when \code{hmmsearch}'s per-sequence
#'   E-value is at most this. Only domains whose own (independent) E-value is
#'   at most this count towards the profile coverage.
#' @param max_profile_gap How many profile positions a match may stop short
#'   of the end adjoining the repeats.
#' @param hmm_dir Directory holding \code{Xo_TALE_Nterm_AA_profile.hmm} and
#'   \code{Xo_TALE_Cterm_AA_profile.hmm}.
#' @param hmmer_path Directory holding the \code{hmmsearch} binary;
#'   \code{NULL} uses the tantale conda environment.
#' @return A tibble with one row per array holding at least one terminus:
#'   \code{array_id}, \code{nterm_aa_evalue}, \code{cterm_aa_evalue} (\code{NA}
#'   when there is no segment or \code{hmmsearch} reports no match),
#'   \code{nterm_aa_profile_gap}, \code{cterm_aa_profile_gap} (profile
#'   positions between the match and the end adjoining the repeats; \code{NA}
#'   without a significant domain), \code{nterm_aa_hit}, \code{cterm_aa_hit}
#'   (\code{TRUE} for a match, \code{FALSE} for a segment without one,
#'   \code{NA} without a segment).
#' @noRd
.tale_termini_hmmsearch <- function(termini, max_evalue, hmm_dir, hmmer_path = NULL,
                                    max_profile_gap = .terminus_max_profile_gap) {
  if (is.null(hmmer_path)) hmmer_path <- .get_hmmer()
  profiles <- c(`N-terminus` = file.path(hmm_dir, "Xo_TALE_Nterm_AA_profile.hmm"),
                `C-terminus` = file.path(hmm_dir, "Xo_TALE_Cterm_AA_profile.hmm"))
  if (!all(file.exists(profiles))) {
    cli::cli_abort(
      c("Cannot find the TALE terminus protein profile{?s} {.file {profiles[!file.exists(profiles)]}}.",
        "i" = "{.arg hmm_dir} must hold them next to the DNA profiles."),
      class = c("tantale_error_hmm_missing", "tantale_error"))
  }

  evalues <- lapply(c("N-terminus", "C-terminus"), function(part) {
    isNterm <- part == "N-terminus"
    seqs <- Biostrings::AAStringSet(gsub("*", "", as.character(termini[[part]]), fixed = TRUE))
    found <- tibble::tibble(array_id = as.character(names(seqs)), evalue = NA_real_,
                            profile_gap = NA_integer_)
    # hmmsearch cannot take an empty sequence; such a segment simply matches nothing
    seqs <- seqs[Biostrings::width(seqs) > 0L]
    if (length(seqs) > 0L) {
      seqFile <- tempfile(fileext = ".fasta")
      tblFile <- tempfile(fileext = ".tbl")
      domFile <- tempfile(fileext = ".domtbl")
      on.exit(unlink(c(seqFile, tblFile, domFile)), add = TRUE)
      Biostrings::writeXStringSet(seqs, seqFile)
      searchCmd <- paste(shQuote(file.path(hmmer_path, "hmmsearch")),
                         "--noali --tblout", shQuote(tblFile),
                         "--domtblout", shQuote(domFile),
                         shQuote(profiles[[part]]), shQuote(seqFile),
                         "> /dev/null")
      .tantale_exec(searchCmd, what = glue::glue("hmmsearch of the {part} profile"))
      readFields <- function(f) {
        strsplit(grep("^#", readLines(f), value = TRUE, invert = TRUE), "\\s+")
      }
      # --tblout: target name in field 1, full-sequence E-value in field 5
      fields <- readFields(tblFile)
      hits <- tibble::tibble(array_id = vapply(fields, `[`, character(1), 1),
                             evalue = as.numeric(vapply(fields, `[`, character(1), 5)))
      found$evalue <- hits$evalue[match(found$array_id, hits$array_id)]
      # --domtblout: target name in field 1, profile length in field 6, the
      # domain's independent E-value in field 13, its profile coordinates in
      # fields 16 and 17
      fields <- readFields(domFile)
      domains <- tibble::tibble(array_id = vapply(fields, `[`, character(1), 1),
                                qlen = as.integer(vapply(fields, `[`, character(1), 6)),
                                i_evalue = as.numeric(vapply(fields, `[`, character(1), 13)),
                                hmm_from = as.integer(vapply(fields, `[`, character(1), 16)),
                                hmm_to = as.integer(vapply(fields, `[`, character(1), 17))) %>%
        dplyr::filter(i_evalue <= max_evalue) %>%
        dplyr::mutate(gap = if (isNterm) qlen - hmm_to else hmm_from - 1L) %>%
        dplyr::group_by(array_id) %>%
        dplyr::summarise(gap = min(gap), .groups = "drop")
      found$profile_gap <- domains$gap[match(found$array_id, domains$array_id)]
    }
    found %>%
      dplyr::mutate(hit = !is.na(evalue) & evalue <= max_evalue &
                      !is.na(profile_gap) & profile_gap <= max_profile_gap) %>%
      dplyr::rename_with(~ paste0(if (isNterm) "nterm" else "cterm", "_aa_", .x),
                         c(evalue, profile_gap, hit))
  })

  dplyr::full_join(evalues[[1]], evalues[[2]], by = "array_id") %>%
    dplyr::arrange(array_id)
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
#' \code{\link{tales_from_telltales}} and the target predictors start from.
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
      order(S4Vectors::mcols(by_array)$seqnames,
            S4Vectors::mcols(by_array)$n_dna_hits), ]
  )
  readr::write_tsv(x = arrayReport, file = paths$array_report)

  ## Write a gff with everything collated
  allGR <- c(unlist(by_array),
             GenomicRanges::makeGRangesFromDataFrame(
               S4Vectors::mcols(by_array),
               seqnames.field = "seqnames",
               keep.extra.columns = TRUE),
             # the unmerged hits are included only when merging happened, so
             # that a merged range can be compared against what went into it
             unmerged)
  rtracklayer::export.gff3(allGR, paths$all_ranges_gff)

  ## Write a fasta file of the seqs of RVDs
  rvds <- Biostrings::BStringSet(S4Vectors::mcols(by_array)$rvd_string)
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
    df <- data.frame(names(stringset), .aa_residue_count(stringset))
    # "N-terminus"/"C-terminus" (from .telltale_align_termini()) to the
    # snake_case column names array_report.tsv actually carries.
    lengthCol <- if (e == "N-terminus") "nterm_aa_length" else "cterm_aa_length"
    colnames(df) <- c("array_id", lengthCol)
    df
  })
  S4Vectors::mcols(by_array) <- merge(S4Vectors::mcols(by_array),
                                      do.call(merge, endsAAlength),
                                      by = "array_id", all.x = TRUE)

  moreInfo <- merge(
    S4Vectors::mcols(by_array),
    data.frame(array_id = names(full_orf),
               longest_orf_length = Biostrings::nchar(full_orf),
               orf_coverage = round(100 * Biostrings::nchar(full_orf) /
                                       GenomicRanges::width(array_seqs[names(full_orf)])),
               longest_orf_seq = full_orf),
    by = "array_id", all.x = TRUE, sort = FALSE)
  rownames(moreInfo) <- moreInfo$array_id
  S4Vectors::mcols(by_array) <- moreInfo[rownames(S4Vectors::mcols(by_array)), ]
  by_array
}


#' Amino acid residues in each sequence, not counting stop codons
#'
#' AnnoTALE ends a terminus record with \code{*} when the stop codon falls
#' inside the part (a truncated C-terminus). The \code{tales} object strips
#' it (\code{tales_ingest.R}), so the report's lengths must not count it
#' either, or the two disagree by one on exactly the truncated arrays
#' (ledger §32.4).
#'
#' @param x An \code{AAStringSet}.
#' @return An integer vector, one count per sequence.
#' @noRd
.aa_residue_count <- function(x) {
  as.integer(Biostrings::width(x) -
               Biostrings::letterFrequency(x, "*", as.prob = FALSE)[, 1])
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
  names(gr) <- gr$hit_id

  ## Updating seqinfo with original seqinfo from the sequences before renaming.
  ## Only the sequences carrying hits are in gr; renameSeqlevels() warns about
  ## any other name it is given (§36).
  gr <- GenomeInfoDb::renameSeqlevels(
    gr, value = seqlevels[names(seqlevels) %in% GenomeInfoDb::seqlevels(gr)])
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


#' Rewrite the subject sequences under names HMMER will accept
#'
#' nhmmer truncates a sequence name at the first space, so two contigs whose
#' headers differ only after a space become indistinguishable in its output
#' and their hits are silently pooled. Every sequence is therefore rewritten
#' to a temporary file under a placeholder name, and the originals restored
#' once the hits are back.
#'
#' @param subject_file The user's fasta.
#' @return A list of the rewritten \code{file}, the \code{seqlevels} map
#'   from placeholder back to original, and the original \code{seqinfo}.
#' @noRd
.telltale_prepare_subject <- function(subject_file) {
  cli::cli_inform("HMMER is very picky about forbiden characters in sequence name. Renaming sequences in {subject_file}.")
  originalSeqs <- Biostrings::readDNAStringSet(filepath = subject_file)
  seqNames <- names(originalSeqs)
  dupNames <- unique(seqNames[duplicated(seqNames)])
  if (length(dupNames) > 0L || !all(nzchar(seqNames))) {
    cli::cli_abort(
      c("Every sequence in {.file {subject_file}} needs a unique, non-empty name.",
        "x" = if (length(dupNames) > 0L) "Duplicated: {.val {dupNames}}",
        "x" = if (!all(nzchar(seqNames))) "{sum(!nzchar(seqNames))} sequence{?s} without a name."),
      class = c("tantale_error_seqnames", "tantale_error"))
  }
  # From the sequences in memory, with their full headers (the names the hits
  # are restored to). A samtools index would key on the first word of each
  # header only, and would be written next to the user's file (§36).
  originalSeqInfo <- GenomeInfoDb::Seqinfo(seqnames = names(originalSeqs),
                                           seqlengths = Biostrings::width(originalSeqs))

  originalSeqlevels <- names(originalSeqs)
  foolproofSeqlevels <- paste0("seq", 1:length(originalSeqlevels))
  names(originalSeqlevels) <- foolproofSeqlevels
  names(originalSeqs) <- foolproofSeqlevels

  cli::cli_inform(paste0("Original seq names : {glue::glue_collapse(originalSeqlevels, sep = ' ; ')}"))
  cli::cli_inform(paste0("Dummy seq names : {glue::glue_collapse(names(originalSeqlevels), sep = ' ; ')}"))

  renamed <- tempfile()
  Biostrings::writeXStringSet(originalSeqs, filepath = renamed)
  list(file = renamed, seqlevels = originalSeqlevels, seqinfo = originalSeqInfo)
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
    hits_report    = file.path(output_dir, "hits_report.tsv"),
    domains_report = file.path(output_dir, "domains_report.tsv"),
    ## Tabular file reporting on putative TALEs (contiguous arrays of domain hits)
    array_report   = file.path(output_dir, "array_report.tsv"),
    ## Gff file with all the identified domains and arrays and their associated data
    all_ranges_gff = file.path(output_dir, "all_ranges.gff"),
    ## fasta of tal orfs that have rvds, and of those predicted not to
    tale_orf_fasta = file.path(output_dir, "putative_tal_orf.fasta"),
    pseudo_tal     = file.path(output_dir, "pseudo_tal_cds.fasta"),
    ## the selected seqs of RVDs, without the separator
    rvd_sequences  = file.path(output_dir, "rvd_sequences.fas"),
    ## the three HMM profiles concatenated, which is what nhmmer is given
    merged_hmm     = file.path(output_dir, "tale_cds_all_diagnostic_regions_hmmfile.out"),
    hmmer_search   = file.path(output_dir, "hmmer_search_out.txt"),
    hmmer_readable = file.path(output_dir, "nhmmer_human_readable_output_of_last_run.txt"),
    ## logging info and some general analysis measures
    log            = file.path(output_dir, "tell_tales.log"),
    ## termini alignments, written by .telltale_align_termini() -- only when
    ## at least 2 putative TALEs are found for that terminus
    n_terminus_dna = file.path(output_dir, "n_terminus_dna_alignment.html"),
    c_terminus_dna = file.path(output_dir, "c_terminus_dna_alignment.html"),
    n_terminus_aa  = file.path(output_dir, "n_terminus_aa_alignment.html"),
    c_terminus_aa  = file.path(output_dir, "c_terminus_aa_alignment.html")
  )
  dir.create(p$annotale)
  if (correct_array) {
    p$correction_dna <- file.path(output_dir, "correction_alignment_dna")
    p$correction_aa  <- file.path(output_dir, "correction_alignment_aa")
    dir.create(p$correction_dna, showWarnings = FALSE)
    dir.create(p$correction_aa, showWarnings = FALSE)
  }
  p
}


#### Actual tell_tales() function definition ####


#' Search and report on the features of TALE protein domains potentially encoded
#' in subject DNA sequences
#'
#' \code{tell_tales} has been primarily written to report on 'corrected' TALE RVD
#' sequences in indels prone, noisy DNA sequences (suboptimally polished genomes
#' assembly, raw reads of long read sequencing technologies such as PacBio or ONT)
#' that would otherwise be missed by conventional tools (eg AnnoTALE).
#'
#' The approach is first to use \href{http://hmmer.org/}{HMMER} to find and
#' categorize regions in the input DNA sequence that are related to the coding
#' sequence of canonical TALE protein domains (N-Term, repeats, C-term). Hits
#' that are (nearly -- see the `min_gap` parameter) adjacent are grouped in
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
#' attempts to 'correct' potential frameshifts in the taleArray sequences. This
#' conveniently removes many artefactual indels but bear in mind that this may
#' also \strong{erroneously} 'correct' genuine frame shifts which can be highly
#' relevant especially for truncTALEs or iTALEs. The resulting 'corrected'
#' taleArray open reading frames are then passed to AnnoTALE.
#'
#' Note that occasionally, when a putative open reading frame does not encode a
#' canonical TALE protein (early frame shift, incomplete ORF, etc...), the
#' "analyze" module of AnnoTALE outputs DNA parts but no protein parts and/or
#' RVD sequence. This should be detected and reported in the tell_tales log.
#'
#' Each sequence of \code{subject_file} is treated as linear. A \emph{tal}
#' gene that spans the junction of a circular molecule (the two ends of an
#' assembled chromosome or plasmid) is cut in two, and is reported, if at
#' all, as two partial arrays at the ends of the sequence. Rotating the
#' sequence so that it starts elsewhere avoids this.
#'
#'
#' @param subject_file Fasta file with DNA sequence(s) to be searched for the
#'   presence of TALE coding sequences (CDS).
#' @param output_dir Path of the output directory. If not specified, results will
#'   be written to current working folder.
#' @param hmm_dir Folder holding the profile HMMs, if you do not want the
#'   ones provided with tantale. It must hold files with the same names: the
#'   three DNA profiles of the nhmmer search
#'   (\code{Xo_TALE_Nterm_CDS_profile.hmm},
#'   \code{Xo_TALE_repeat_CDS_profile.hmm},
#'   \code{Xo_TALE_Cterm_CDS_profile.hmm}) and the two protein profiles of
#'   the terminus check (\code{Xo_TALE_Nterm_AA_profile.hmm},
#'   \code{Xo_TALE_Cterm_AA_profile.hmm}).
#' @param nterm_min_score Minimal nhmmer score cut_off value to
#'   consider the hit as genuine
#' @param repeat_min_score Minimal nhmmer score cut_off value to consider
#'   the hit as genuine
#' @param cterm_min_score Minimal nhmmer score cut_off value to
#'   consider the hit as genuine
#' @param terminus_max_evalue Maximum \code{hmmsearch} E-value for the
#'   segment AnnoTALE reports on either side of the repeats to count as a TALE
#'   N- or C-terminus. The segments are searched with the TALE terminal-domain
#'   protein profiles of \code{hmm_dir}; this decides the \code{NTERM},
#'   \code{CTERM} and \code{XXXXX} codes (see \code{\link{tales_anchor_codes}}).
#'   Genuine termini truncated to about 40 residues still match with E-values
#'   below 1e-18. The match must also reach, within 10 positions, the end of
#'   the profile that adjoins the repeats: a terminus whose repeat-side part
#'   is in another reading frame after a frameshift matches only up to the
#'   frameshift, and is coded \code{XXXXX}.
#' @param min_dna_hits Minimum number of nhmmer hits for a subject
#'   sequence (a contig, a chromosome) to be considered further. A cheap way
#'   to discard whole sequences that carry nothing but stray matches, before
#'   any expensive work is done on them. It says nothing about the length of
#'   the TALE arrays found within a sequence that passes -- see
#'   \code{min_array_length} for that.
#' @param min_array_length Minimum number of \strong{repeat units} for a TALE
#'   array to be kept. Defaults to \code{0}, which keeps everything.
#'
#'   Counting repeats rather than all hits means an array is not penalised for
#'   having had its termini missed, and matches what "array length" usually
#'   means for a TALE: the number of repeats is what determines how long a
#'   target box it recognises.
#'
#'   Whether a short array is noise or a genuinely truncated TALE is a
#'   judgement about the biology, which is why nothing is discarded unless you
#'   ask. A pseudogene with three surviving repeats is real, and may be what
#'   you are looking for.
#' @param merge_hits Merge overlapping nhmmer hits of the same domain type,
#'   since nhmmer can report one repeat as two overlapping hits. With
#'   \code{FALSE}, such a repeat is counted twice in \code{n_dna_hits} and
#'   by \code{min_array_length}, and a warning names the arrays concerned.
#' @param min_gap Minimum gap in base pairs between two tale domain hits for
#'   them to be considered distinct. If the length of the gap is below this
#'   value, domains are considered "contiguous" and grouped in the same array.
#' @param extremity_codes Set this to \code{FALSE} if you do not want the
#'   terminus codes in the RVD strings of \code{rvd_sequences.fas} and
#'   \code{array_report.tsv}.
#' @param rvd_sep Symbol acting as a separator in RVD sequences
#' @param hmmer_path Specify the path to a directory holding the HMMER executable
#'   if you do not want to use the ones provided with tantale.
#' @param extend_len number of nucleotides to extend in 3'-end at the tal
#'   ORF prediction stage.
#' @param correct_array Whether to pass each array through
#'   \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}} before
#'   AnnoTALE sees it. \code{FALSE} by default; see Details for what turning
#'   it on buys (removing artefactual indels) and risks (erroneously
#'   "correcting" a genuine frameshift).
#' @param correction_ref Fasta of reference TALE proteins to correct
#'   against. Two are shipped, both built by
#'   \code{data-raw/make_correction_references.R} from the same source:
#'
#'   \itemize{
#'     \item \code{tale_correction_ref.fa.gz} (default, 494 sequences) --
#'       every distinct TALE protein of at least 300 aa found across 70
#'       \emph{Xanthomonas oryzae} genomes.
#'     \item \code{tale_correction_ref_representative.fa.gz} (136) -- a
#'       diversity-sampled subset, for a smaller footprint.
#'   }
#'
#'   Both keep the pseudogenes. Their frameshifts came from high-quality
#'   genomes and so are real biology, and correction is meant to recover a
#'   sequence as it exists in nature rather than reshape every array into an
#'   intact TALE. A reference of only intact TALEs risks "repairing" a
#'   genuine pseudogene into an ORF no strain carries.
#' @param max_comparisons How many reference proteins each array may be
#'   aligned against during frameshift correction, and \strong{the main
#'   control on how long correction takes}. \code{NULL}, the default, allows
#'   all of them.
#'
#'   \code{DECIPHER::CorrectFrameshifts()} scores every reference with a
#'   cheap distance first, sorts them, and only then aligns against the
#'   closest \code{max_comparisons} of them -- stopping sooner if one is
#'   close enough. So the references never reached cost almost nothing, and
#'   lowering this does not change \emph{which} references are preferred,
#'   only how deep the search goes before settling for the best seen.
#'
#'   Measured on four arrays against the 1057-sequence source set, all giving
#'   byte-identical corrected sequences:
#'
#'   \tabular{lr}{
#'     \strong{max_comparisons} \tab \strong{seconds} \cr
#'     all (1057) \tab 252 \cr
#'     400 \tab 179 \cr
#'     100 \tab 47 \cr
#'     50 \tab 23 \cr
#'     20 \tab 10 \cr
#'   }
#'
#'   \strong{The trade-off.} A cap risks a divergent array whose only good
#'   reference lies outside the closest \code{max_comparisons} by the cheap
#'   pre-screen. That pre-screen is an approximation, so a low cap trusts it
#'   to rank the truly best reference near the top.
#'
#'   When it fails, it corrects the array against a poor reference. That is
#'   worse than leaving the array uncorrected, because the result still looks
#'   like a corrected ORF. Against a deliberately small 20-sequence
#'   reference, the same four arrays give:
#'
#'   \tabular{ll}{
#'     \strong{max_comparisons} \tab \strong{indels called per array} \cr
#'     all (20), 20, 10 \tab 2, 2, 0, 1 \cr
#'     5 \tab 2, 2, 0, 2 \cr
#'     2 \tab 9, 11, 0, 15 \cr
#'   }
#'
#'   At 2 the aligner cannot reach a decent reference and invents indels
#'   wholesale. What matters is whether the closest \code{max_comparisons}
#'   references are genuinely close, whatever the size of the reference set:
#'   20 of 1057 is ample, 5 of 20 is not. With a large reference set a cap in
#'   the tens is safe and very much faster; with a small or a poorly matched
#'   one, prefer the default and pay for the full search.
#' @param frameshift Frameshift penalty passed to
#'   \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}}'s
#'   \code{frameShift}. tantale's default is \code{-11}, overriding
#'   DECIPHER's own \code{-15}; change it only with a reason.
#' @param ... Additional parameters for the
#'   \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}} function.
#' @return Called for its side effects (writing files). If everything runs
#'   smoothly, it invisibly returns the path of the directory where output
#'   files were written.
#'
#'
#'   List of output files:
#'   \itemize{
#'   \item all_ranges.gff: gff file of all Tal arrays detected by HMMer
#'   \item array_report.tsv: one row per candidate TALE array. Columns
#'   named \code{*_dna_*} describe the nhmmer search of the subject DNA;
#'   columns named \code{*_aa_*} describe the protein segments AnnoTALE
#'   extracted from the array's longest ORF.
#'   \itemize{
#'     \item \emph{array_id}, \emph{seqnames}, \emph{start}, \emph{end},
#'     \emph{strand}: the array's identifier and the span of its nhmmer hits
#'     on the subject sequence.
#'     \item \emph{n_dna_hits}: number of nhmmer hits (N-terminus, repeats and
#'     C-terminus profiles together) grouped in the array. A terminus hit
#'     usually overlaps the adjacent repeat hit by a few nucleotides, at the
#'     boundary between the two domains.
#'     \item \emph{array_seq}: DNA sequence of that span.
#'     \item \emph{nterm_dna_hit}, \emph{cterm_dna_hit}: whether an nhmmer hit
#'     of the N- (C-) terminus DNA profile is part of the array, anywhere in
#'     it.
#'     \item \emph{rvd_string}: the RVDs AnnoTALE read, separated by
#'     \code{rvd_sep}, with the terminus codes described under
#'     \emph{rvd_sequences.fas}. Empty when AnnoTALE found no RVD.
#'     \item \emph{has_aberrant_repeat}: whether AnnoTALE flagged a repeat of
#'     non-canonical length (a lowercase letter in its RVD).
#'     \item \emph{nterm_aa_evalue}, \emph{cterm_aa_evalue}: E-value of the
#'     \code{hmmsearch} match between the segment AnnoTALE reported upstream
#'     (downstream) of the repeats and the TALE N- (C-) terminal protein
#'     profile. \code{NA} when there is no segment, or no match with an
#'     E-value up to 10.
#'     \item \emph{nterm_aa_profile_gap}, \emph{cterm_aa_profile_gap}:
#'     number of profile positions between the end of that match and the end
#'     of the profile that adjoins the repeats (the last position of the
#'     N-terminal profile, the first of the C-terminal one). \code{0} for a
#'     match that reaches the repeats, \code{NA} when there is no match.
#'     \item \emph{nterm_aa_hit}, \emph{cterm_aa_hit}: \code{TRUE} when that
#'     E-value is at most \code{terminus_max_evalue} and the profile gap at
#'     most 10, \code{FALSE} for a segment that does not match,
#'     \code{NA} when AnnoTALE reported no segment on that side.
#'     \item \emph{nterm_aa_length}, \emph{cterm_aa_length}: length of those
#'     segments in amino acid residues, excluding a stop codon, as in the
#'     \code{tales} object's \code{aa_seq}.
#'     \item \emph{longest_orf_length}, \emph{longest_orf_seq}: the longest
#'     ORF found in the array region extended by \code{extend_len}
#'     nucleotides at its 3' end.
#'     \item \emph{orf_coverage}: that ORF's length as a percentage of the
#'     extended region's length.
#'     \item \emph{predicted_dels_count}, \emph{predicted_ins_count} (with
#'     \code{correct_array = TRUE}): number of putative deletions/insertions
#'     in the raw sequence that
#'     \code{\link[DECIPHER:CorrectFrameshifts]{CorrectFrameshifts}}
#'     corrected.
#'   }
#'   \item hits_report.tsv: report of all hits detected by HMMer
#'   \item hits_report.gff: gff file of all hits detected by HMMer
#'   \item domains_report.tsv: report of all Tal amino acid domains detected by AnnoTALE analyze
#'   \item putative_tal_orf.fasta: for each array in which AnnoTALE found
#'    RVDs, the DNA of the longest ORF of the array region extended by
#'    \code{extend_len} nucleotides at its 3' end (after frameshift
#'    correction with \code{correct_array = TRUE}). This is the putative
#'    TALE coding sequence AnnoTALE analysed.
#'   \item pseudo_tal_cds.fasta: for each array in which AnnoTALE found no
#'    RVD, the DNA of the array region extended by \code{extend_len}
#'    nucleotides at its 3' end, as found in \code{subject_file}: candidate
#'    pseudogenes, assembly errors or false detections, kept for inspection.
#'   \item rvd_sequences.fas: the RVDs (separated by \code{rvd_sep}) of each
#'    array for which AnnoTALE found at least one. With \code{extremity_codes
#'    = TRUE} (the default), each string is bracketed by terminus codes (see
#'    \code{\link{tales_anchor_codes}}): \code{NTERM} (\code{CTERM}) when the
#'    segment AnnoTALE reported upstream (downstream) of the repeats matches
#'    the TALE N- (C-) terminal protein profile, \code{XXXXX} when it does
#'    not, and no code when AnnoTALE reported no segment on that side.
#'   \item c_terminus_aa_alignment.html: protein alignment of all C-termini
#'   (only written when at least 2 were found; skipped with a warning otherwise)
#'   \item c_terminus_dna_alignment.html: DNA alignment of all C-termini (same condition)
#'   \item n_terminus_aa_alignment.html: protein alignment of all N-termini (same condition)
#'   \item n_terminus_dna_alignment.html: DNA alignment of all N-termini (same condition)
#'   \item tale_cds_all_diagnostic_regions_hmmfile.out: HMMER profile used for tale cds search.
#'   \item hmmer_search_out.txt: ignore
#'   \item nhmmer_human_readable_output_of_last_run.txt: primary HMMER output file.
#'   \item tell_tales.log: a log file
#'   \item annotale folder: folder containing result of AnnoTALE analyze for all Tal arrays
#'   \item correction_alignment_aa folder: folder containing protein alignment of Tal array detected by HMMer and corrected Tal array if \code{correct_array} = TRUE
#'   \item correction_alignment_dna folder: folder containing DNA alignment of Tal array detected by HMMer and corrected Tal array if \code{correct_array = TRUE}
#'   }
#' @export
#' @family TALE discovery
#' @examples
#' \donttest{
#' # Needs nhmmer and AnnoTALE, resolved from the tantale conda environment
#' # (and a Java runtime) on first use.
#' subj <- system.file("extdata", "bai3_sample_tal_genomic_regions.fasta",
#'                     package = "tantale")
#' out <- tempfile("tell_tales_example")
#' tell_tales(subject_file = subj, output_dir = out)
#' tales_from_telltales(out)
#' }
tell_tales <- function(
  subject_file,
  output_dir = getwd(),
  hmm_dir = system.file("extdata", "hmmProfile", package = "tantale", mustWork = T),
  nterm_min_score = 300,
  repeat_min_score = 20,
  cterm_min_score = 200,
  terminus_max_evalue = 1e-5,
  min_dna_hits = 4,
  min_array_length = 0,
  merge_hits = TRUE,
  min_gap = 35,
  extremity_codes = TRUE,
  rvd_sep = "-",
  hmmer_path = NULL,
  extend_len = 300,
  correct_array = FALSE,
  correction_ref = system.file("extdata", "tale_correction_ref.fa.gz", package = "tantale", mustWork = T),
  max_comparisons = NULL,
  frameshift = -11,
  ...
) {

  # @param taleArrayStartAnchorCode This scalar character vector will symbolize a
  #   TALE N-TERM CDS hit in the RVD sequence
  # @param taleArrayEndAnchorCode This scalar character vector will symbolize a
  #   TALE C-TERM CDS hit in the RVD sequence

  # Sequences are treated as linear (see Details). A `circular` argument
  # could set the seqinfo so that a tal gene spanning the junction of a
  # circular molecule is found whole (ledger §43).

  ####   Paths of output files   ####
  paths <- .telltale_paths(output_dir, correct_array)
  
  
  ####   Checks for parameters and other things   ####

  ## Deal with spaces in sequence names because this messes up parsing of hmmer output
  ## original_subject_file is kept for the log (.telltale_log()) -- subject_file
  ## itself is about to be reassigned to a renamed temp copy, and the log is
  ## meant to say what the caller actually ran on.
  original_subject_file <- subject_file
  subject <- .telltale_prepare_subject(subject_file)
  subject_file <- subject$file
  originalSeqlevels <- subject$seqlevels
  originalSeqInfo <- subject$seqinfo

  hmm <- .telltale_hmm_profiles(hmm_dir, paths$merged_hmm)
  
  ####   Find the TALE domain CDS hits   #####
  nhmmerTabularOutput <- .telltale_find_domain_hits(
    subject_file = subject_file, hmm = hmm, paths = paths,
    hmmer_path = hmmer_path,
    nterm_min_score = nterm_min_score, repeat_min_score = repeat_min_score,
    cterm_min_score = cterm_min_score, min_dna_hits = min_dna_hits)
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
  grouped <- .telltale_group_arrays(nhmmerOutputGR, min_gap, subjectDNASequences, hmm,
                                    min_array_length = min_array_length)
  if (is.null(grouped)) return(invisible(output_dir))
  arraysGR <- grouped$arrays
  hitsByArraysLst <- grouped$by_array

  #####   Extend DNA Tal arrays   #####
  ## Extract the genomic sequence of arrays +-bp on the borders
  ## An array near a sequence end is extended past it; resize() warns about
  ## that, and trim() clips it right after, so only that warning is muffled
  ## (§36)
  completeArraysGR <- arraysGR
  extdCompleteArraysGR <- withCallingHandlers(
    GenomicRanges::resize(completeArraysGR,
                          width = GenomicRanges::width(completeArraysGR) + extend_len,
                          fix = "start", ignore.strand = FALSE),
    warning = function(w) {
      if (grepl("out-of-bound range", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }) %>%
    GenomicRanges::trim(use.names = TRUE)
  extdCompleteArraysSeqs <- BSgenome::getSeq(subjectDNASequences, extdCompleteArraysGR)
  
  

  orfResult <- .telltale_array_orfs(
    array_seqs = extdCompleteArraysSeqs, by_array = hitsByArraysLst,
    correct_array = correct_array, correction_ref = correction_ref,
    frameshift = frameshift, max_comparisons = max_comparisons,
    paths = paths, ...)
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

  # save tals orfs that have rvds
  Biostrings::writeXStringSet(fullTalOrf[names(fullTalOrf) %in% names(seqsOfRVDs)], paths$tale_orf_fasta)
  
  # save tals that DO NOT have rvds
  Biostrings::writeXStringSet(extdCompleteArraysSeqs[!names(extdCompleteArraysSeqs) %in% names(seqsOfRVDs)], paths$pseudo_tal)
  
  #### Align N-term and C-term ####
  .telltale_align_termini(paths$annotale, paths, type = "DNA")
  endsAA <- .telltale_align_termini(paths$annotale, paths, type = "AA")

  #### Do the termini look like TALE terminal domains? ####
  terminiHits <- .tale_termini_hmmsearch(endsAA, max_evalue = terminus_max_evalue,
                                         hmm_dir = hmm_dir, hmmer_path = hmmer_path)

  #### RVD strings ####
  # A lowercase letter in an RVD is AnnoTALE's way of flagging a repeat whose
  # length departs from the canonical ~34 aa. An aberrant repeat changes how
  # the array should be read, which is worth knowing before the RVDs are used
  # to predict targets.
  rvdStrings <- as.character(seqsOfRVDs)
  hasAberrantRepeat <- ifelse(nzchar(rvdStrings), grepl("[a-z]", rvdStrings), NA) %>%
    stats::setNames(names(rvdStrings))
  rvdStrings <- gsub("-", rvd_sep, rvdStrings, fixed = TRUE)

  # Terminus codes, so that a string of RVDs and a string of repeat codes
  # describe the same number of parts: NTERM/CTERM for a segment matching the
  # TALE terminus profile, XXXXX for a segment that does not, nothing when
  # AnnoTALE reported no segment on that side of the repeats.
  if (extremity_codes) {
    anchors <- unname(tales_anchor_codes())   # NTERM, CTERM, XXXXX
    endHits <- terminiHits[match(names(rvdStrings), terminiHits$array_id), ]
    ntermCode <- ifelse(endHits$nterm_aa_hit, anchors[1], anchors[3])
    ctermCode <- ifelse(endHits$cterm_aa_hit, anchors[2], anchors[3])
    rvdStrings <- paste0(ifelse(is.na(ntermCode), "", paste0(ntermCode, rvd_sep)),
                         rvdStrings,
                         ifelse(is.na(ctermCode), "", paste0(rvd_sep, ctermCode))) %>%
      stats::setNames(names(rvdStrings))
  }

  arrayMeta <- S4Vectors::mcols(hitsByArraysLst)
  arrayMeta$rvd_string <- unname(rvdStrings[arrayMeta$array_id])
  arrayMeta$rvd_string[is.na(arrayMeta$rvd_string)] <- ""
  arrayMeta$has_aberrant_repeat <- unname(hasAberrantRepeat[arrayMeta$array_id])
  endHits <- terminiHits[match(arrayMeta$array_id, terminiHits$array_id), ]
  arrayMeta$nterm_aa_evalue <- endHits$nterm_aa_evalue
  arrayMeta$cterm_aa_evalue <- endHits$cterm_aa_evalue
  arrayMeta$nterm_aa_profile_gap <- endHits$nterm_aa_profile_gap
  arrayMeta$cterm_aa_profile_gap <- endHits$cterm_aa_profile_gap
  arrayMeta$nterm_aa_hit <- endHits$nterm_aa_hit
  arrayMeta$cterm_aa_hit <- endHits$cterm_aa_hit
  S4Vectors::mcols(hitsByArraysLst) <- arrayMeta

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
    params = list(subject_file = original_subject_file, output_dir = output_dir,
                  nterm_min_score = nterm_min_score,
                  repeat_min_score = repeat_min_score,
                  cterm_min_score = cterm_min_score,
                  terminus_max_evalue = terminus_max_evalue,
                  min_dna_hits = min_dna_hits,
                  min_array_length = min_array_length, merge_hits = merge_hits,
                  min_gap = min_gap, extend_len = extend_len,
                  correct_array = correct_array, correction_ref = correction_ref,
                  max_comparisons = max_comparisons,
                  frameshift = frameshift),
    log_file = paths$log, hmm = hmm, subject_seqs = subjectDNASequences,
    arrays = arraysGR, by_array = hitsByArraysLst, array_report = arrayReport,
    gaps_below_500 = gaplengthBetweenHitDomainsbelow500,
    gap_quartiles = quartilesGapLength,
    annotale_messages = annoTaleMessages)

  return(invisible(output_dir))
}
