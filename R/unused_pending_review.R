#### Functions with no caller in the package ####
#
# Nothing in R/ calls any of these. They are parked here rather than deleted
# so that the question "is this actually useless?" can be answered one at a
# time, with the code in front of us, instead of guessed at now.
#
# They are still part of the package: sourced, checked, and callable with
# `tantale:::`. Parking is a bookkeeping move, not a retirement -- retired
# code goes to inst/legacy/ instead, where it is no longer compiled.
#
# Per function, what is known:
#
# .write_hmm_file()        wraps hmmbuild.  Never called. tell_tales() uses the
#                          pre-built profile in inst/extdata/hmmProfile, so
#                          nothing in the package ever builds one.
# .run_hmmer_search()      wraps hmmsearch (protein). Never called; tell_tales()
#                          uses .run_nhmmer_search(), the nucleotide sibling,
#                          which stays in tellTale_utilities.R.
# .run_hmmalign()          wraps hmmalign. Never called.
# .extract_seqs_from_hits() Never called. Carried the comment "THIS SHOULD BE
#                          MADE OBSOLETE AND CODE USING IT SHOULD BE MODIFIED";
#                          that apparently happened, without the function being
#                          revisited.
# .rvd_to_repeat_align()   Not called anywhere in R/, but tests/testthat
#                          exercises it, so unlike the others it still has
#                          coverage. Its sibling .rvd_to_match_align() is live.
#
# The first three are a matched set: together they are a complete HMMER
# profile-building workflow that the package does not currently perform.

.write_hmm_file <- function(hmmer_path = NULL, alignment_file, hmm_out_file) {
  if (is.null(hmmer_path)) hmmer_path <- .get_hmmer()
  .check_hmmer(hmmer_path)
  buildCmd <- paste(file.path(hmmer_path,"hmmbuild"),
                    hmm_out_file,
                    alignment_file,
                    sep = " ")
  commandOut <- system(command = buildCmd, ignore.stderr = FALSE, intern = TRUE)
  return(commandOut)
}

.run_hmmer_search <- function(hmmer_path = NULL, subject_file, hmm_file, search_out_file, readable_out_file) {
  if (is.null(hmmer_path)) hmmer_path <- .get_hmmer()
  .check_hmmer(hmmer_path)
  searchCmd <- paste(file.path(hmmer_path, "hmmsearch"),
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

.run_hmmalign <- function(hmmer_path = NULL, hmm_file, seqs_file, align_out_file) {
  if (is.null(hmmer_path)) hmmer_path <- .get_hmmer()
  .check_hmmer(hmmer_path)
  alignCmd <- paste(file.path(hmmer_path, "hmmalign"),
                    "--outformat Phylip", #Stockholm, SELEX, Clustal, Phylip, Pfam, A2M, PSIBLAST.
                    "--trim",
                    hmm_file,
                    seqs_file,
                    ">", align_out_file,
                    sep = " "
  )
  system(command = alignCmd, ignore.stderr = FALSE, intern = TRUE)
}

## !! THIS SHOULD BE MADE OBSOLETE AND CODE USING IT SHOULD BE MODIFIED
.extract_seqs_from_hits <- function(nhmmer_hits, dna_seqs){
  repeatSeqsSetList <- mapply(
    function(hitID, start, end, strand, subjectID, sequences) {
      seq <- XVector::subseq(sequences[subjectID], start, end)
      if (strand == "-") {seq <- Biostrings::reverseComplement(seq)}
      names(seq) <- hitID
      return(seq)
    },
    hitID = nhmmer_hits$hitID,
    start = nhmmer_hits$start,
    end = nhmmer_hits$end,
    strand = nhmmer_hits$strand,
    subjectID = nhmmer_hits$target_name,
    MoreArgs = list(sequences = dna_seqs),
    USE.NAMES = FALSE)
  do.call(c, repeatSeqsSetList)
}

#' View an RVD alignment in terms of repeat codes
#'
#' The inverse direction of \code{repeat_to_rvd_align()}: given an alignment
#' computed on RVDs, substitute each non-gap cell with the corresponding repeat
#' code. The two alphabets differ greatly -- a few dozen RVDs against hundreds
#' of mostly-singleton repeat codes -- so aligning on one and viewing as the
#' other is a genuinely different result from aligning on the other directly.
#'
#' The mapping is **positional**: the k-th non-gap cell of a row is taken to be
#' the k-th element of that row's repeat vector. That is only well defined when
#' the two agree in length, which is now checked.
#'
#' @param rvd_msa_by_group A character matrix of aligned RVDs, rows named by array.
#' @param repeat_vecs A named list of repeat-code vectors, one per row of
#'   \code{rvd_msa_by_group}.
#' @return A character matrix with the dimensions and dimnames of \code{rvd_msa_by_group}.
#' @keywords internal
.rvd_to_repeat_align <- function(rvd_msa_by_group, repeat_vecs) {
  missingRows <- setdiff(rownames(rvd_msa_by_group), names(repeat_vecs))
  if (length(missingRows) > 0L) {
    cli::cli_abort(
      c("Every aligned row needs a matching entry in {.arg repeat_vecs}.",
        "x" = "Missing: {.val {missingRows}}"),
      class = c("tantale_error_rvd_repeat_missing", "tantale_error"))
  }
  nonGap <- rowSums(!is.na(rvd_msa_by_group))
  lens <- lengths(repeat_vecs[rownames(rvd_msa_by_group)])
  bad <- which(nonGap != lens)
  if (length(bad) > 0L) {
    cli::cli_abort(
      c("The back-mapping is positional, so each row must have as many non-gap \\
         cells as it has repeat codes.",
        "x" = "Mismatched row{?s}: {.val {rownames(rvd_msa_by_group)[bad]}}",
        "i" = "non-gap cells {nonGap[bad]} vs {lens[bad]} repeat codes"),
      class = c("tantale_error_rvd_repeat_length", "tantale_error"))
  }
  repSeqs <- lapply(rownames(rvd_msa_by_group), function(r) {
    rvdSeq <- rvd_msa_by_group[r,]
    repSeq <- repeat_vecs[[r]]

    n = 1
    for (i in 1:length(rvdSeq)) {
      if (is.na(rvdSeq[i])) {
        next()
      } else {
        rvdSeq[i] <- repSeq[n]
        n <- n + 1
      }
    }
    rvdSeq <- matrix(rvdSeq, nrow = 1)
    return(rvdSeq)
  })
  repeatMsaByGroup <- do.call(rbind, repSeqs)
  repeatMsaByGroup <- matrix(repeatMsaByGroup, nrow = nrow(rvd_msa_by_group))
  rownames(repeatMsaByGroup) <- rownames(rvd_msa_by_group)
  colnames(repeatMsaByGroup) <- colnames(rvd_msa_by_group)
  return(repeatMsaByGroup)
}



# ---------------------------------------------------------------------------
# repeat_to_rvd_align() -- unexported 2026-09. Nothing in R/ calls it. It
# turned an alignment matrix of repeat codes into one of RVDs, which was
# needed when plotting took two separate matrices. A tales_msa carries both
# layers at once, so plot() names them (fill = "dom_code", label = "rvd")
# and the conversion has no remaining purpose. Tests still exercise it.
# ---------------------------------------------------------------------------

#' Substitute Distal repeat IDs for RVDs in a TALE alignment matrix.
#'
#'
#' @param repeat_align A multiple TALE repeat sequences alignment in the form of
#'   a matrix as returned by
#'   \code{\link{tales_align}}.
#' @param rvd_map The return value of the
#'   \code{\link[tantale:repeat_to_rvd_map]{repeat_to_rvd_map}} function or the 
#'   \code{\link[tantale:repeat_to_rvd_map_distalr]{repeat_to_rvd_map_distalr}} function
#'   if you used the \code{\link{tales_compare}} function.
#'   
#'
#' @return A TALE alignment matrix made up of RVD sequences.
#' @noRd
repeat_to_rvd_align <- function(repeat_align , rvd_map) {
  states <- unique(as.vector(repeat_align))
  ##### TODO: check that all values in states are present in the rvd_map df ####
  # If not, error
  rvd_align <- t(
    apply(repeat_align, 1,
          function(repeatSeq){
            rvdSeq <- rvd_map$RVD[match(repeatSeq, rvd_map$repeatID)]
          }
    )
  )
  rvd_align <- matrix(rvd_align, nrow = nrow(repeat_align)) # in case of 1-row matrix
  rownames(rvd_align) <- rownames(repeat_align)
  colnames(rvd_align) <- colnames(repeat_align)
  return(rvd_align)
}
