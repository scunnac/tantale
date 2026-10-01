# Retired: repeat_to_rvd_map(), repeat_to_rvd_map_distalr(),
# tale_parts_to_rvd(), .rvd_to_repeat_align()
#
# Superseded by the tales object (ledger §37, maintainer's road map
# dev/notes_for_claude.md):
#   repeat_to_rvd_map()          the dom_code -> rvd mapping is two columns of a
#   repeat_to_rvd_map_distalr()  tales; its one-RVD-per-code assertion is the
#                                "aa_seq_rvd_inconsistent" anomaly of tales()
#                                together with the dom_code/aa_seq bijection
#                                check (ledger §2, §6)
#   tale_parts_to_rvd()          tales_rvd_strings(x, repeats_only = FALSE)
#   .rvd_to_repeat_align()       never called; a tales_msa carries both layers
#
# Their tests (test_untested_exports.R, test_error_conditions.R) were removed
# with them; see git history.
#
# Moved out of R/ on 2026-10-01.

#' Generate a mapping between Distal repeat IDs and their cognate RVD
#'
#' Uses Distal repeat sequences and RVD sequences from a set of TALEs to return
#' the association between repeat ID and RVD.
#'
#' Care must be taken that TALEs in the two sets of sequences have the same name.
#' The function checks that the two sets of sequences have the same structure
#' (names and lengths), so they must also agree on whether they include the
#' N- and C-terminal codes.
#'
#' @param repeat_vecs Expects a list of Distal repeat IDs character
#'   vectors. Each \strong{named} element corresponding to a TALE.
#' @param rvd_vecs A named list of RVD vectors, one per array, parallel to
#'   \code{repeat_vecs}.
#' @return A two columns repeatID - RVD data frame.
#' @export
#' @family tales projections
#' @examples
#' repeat_vecs <- list(A1 = c("12", "45", "12"), A2 = c("45", "78"))
#' rvd_vecs <- list(A1 = c("HD", "NI", "HD"), A2 = c("NI", "NG"))
#' repeat_to_rvd_map(repeat_vecs, rvd_vecs)
repeat_to_rvd_map <- function(repeat_vecs, rvd_vecs) {
  # Making sure, these objects are indentical in every ways but the actual values of the vectors
  stopifnot(setequal(names(repeat_vecs), names(rvd_vecs)))
  lrep <- sapply(repeat_vecs, length)
  lrep <- lrep[order(names(lrep))]
  lrvd <- sapply(rvd_vecs, length)
  lrvd <- lrvd[order(names(lrvd))]
  stopifnot(names(lrvd) == names(lrep))
  stopifnot(apply(cbind(lrep, lrvd), 1, function(x) x[1] == x[2]))
  
  # Function to melt the lists
  .l2df <- function(l) {
    dplyr::bind_rows(lapply(l, function(v) data.frame(idx = 1:length(v),
                                                      repeats = v,
                                                      stringsAsFactors = FALSE)
    ),
    .id = "Name")
  }
  talesRepeatDf <- .l2df(repeat_vecs)
  talesRvdDf <- .l2df(rvd_vecs)
  stopifnot(nrow(talesRepeatDf) == nrow(talesRvdDf))
  # Merging to have the repeat ID vs RVDs
  repeat2rvd <- dplyr::full_join(talesRepeatDf, talesRvdDf, by = c("idx" = "idx", "Name" = "Name"))
  colnames(repeat2rvd) <- c("names", "idx", "repeatID", "RVD")
  # Check that for each repeat there is only one corresponding RVD (the converse is NOT true)
  repeat2rvd <- repeat2rvd %>% dplyr::group_by(repeatID, RVD) %>%
    dplyr::summarise(count = dplyr::n())
  check <- repeat2rvd %>%
    dplyr::arrange(repeatID) %>%
    dplyr::group_by(repeatID) %>%
    dplyr::summarise(count = dplyr::n())
  stopifnot(sum(check$count) == length(unique(repeat2rvd$repeatID)))
  # Simplifiy the df
  repeat2rvd <- repeat2rvd %>% dplyr::select(repeatID, RVD) %>% dplyr::ungroup()
  
  return(repeat2rvd)
}


#' Generate a mapping between Distal repeat IDs and their cognate RVD
#'
#' Returns the association between each domain code (\code{dom_code}) and
#' its RVD, one row per distinct pair.
#'
#' @param tale_parts A \code{\link{tales}} object or data frame with
#'   \code{dom_code} and \code{rvd} columns, such as the \code{tales} element
#'   of \code{\link{tales_compare_distal}}'s result.
#' @return A two columns repeatID - RVD data frame.
#' @export
#' @family tales projections
#' @examples
#' parts <- data.frame(dom_code = c(1, 2, 1, 3), rvd = c("HD", "NI", "HD", "NG"))
#' repeat_to_rvd_map_distalr(parts)
repeat_to_rvd_map_distalr <- function(tale_parts) {
  .tales_require(tale_parts, "repeat_to_rvd_map_distalr")
  tale_parts %>% 
    dplyr::select(dom_code, rvd) %>%
    dplyr::distinct() %>%
    dplyr::rename(repeatID = dom_code,  RVD = rvd) %>%
    dplyr::arrange(repeatID)
}
























#' Generates a RVD sequences set from a tale_parts object
#'
#' Returns a \code{\link[Biostrings]{BStringSet}} of RVD sequences, one per
#' array, ordered by \code{position_in_array} and joined with \code{sep}.
#'
#' @param tale_parts A \code{\link{tales}} object or data frame with
#'   \code{array_id}, \code{position_in_array} and \code{rvd} columns.
#' @param sep Separator joining consecutive RVDs within an array's string.
#' @param rvd_only Return only RVDs and omit N- and C-terminal domains.
#' @return A \code{\link[Biostrings]{BStringSet}}, one element per array,
#'   named by \code{array_id}.
#' @export
#' @family tales projections
#' @examples
#' parts <- data.frame(
#'   array_id = c("A1", "A1", "A1", "A2", "A2"),
#'   position_in_array = c(1, 2, 3, 1, 2),
#'   rvd = c("NTERM", "HD", "CTERM", "NTERM", "NI")
#' )
#' tale_parts_to_rvd(parts)
tale_parts_to_rvd <- function(tale_parts, sep = "-", rvd_only = FALSE) {

  if(rvd_only) {
    # tales_anchor_codes() rather than a retyped list: the set has THREE
    # members, and the hardcoded pair here silently kept "XXXXX" -- a
    # terminus detected but not identified -- in a repeats-only string.
    tale_parts %<>% dplyr::filter(!rvd %in% tales_anchor_codes())
  }

  rvdStrings <- tale_parts %>%
    dplyr::group_by(array_id) %>%
    dplyr::arrange(position_in_array) %>%
    dplyr::summarise(rvdString = paste(rvd, collapse = sep))
  rvdStringsSet <- Biostrings::BStringSet(rvdStrings$rvdString)
  names(rvdStringsSet) <- rvdStrings$array_id
  return(rvdStringsSet)
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


