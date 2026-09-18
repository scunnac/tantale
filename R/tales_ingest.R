#### Building a tales object from a tell_tales run ####
#
# tales_from_telltale() is the public entry point; everything else here is a
# private step of it. Kept together, and out of distalr.R, so that file is
# about comparing TALEs rather than reading them off disk.
#
# The layout being read is what tell_tales() writes: one directory per
# region of interest under annotale/, holding AnnoTALE's split of the CDS
# into N-terminus, repeats and C-terminus, plus a separate file of RVD
# sequences. The two are cross-checked against each other here, because
# AnnoTALE and the RVD caller can disagree about how many repeats an array
# has.

.tale_parts_from_file <- function(fasta) {
  if (grepl("TALE_Protein_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readAAStringSet(fasta)
  if (grepl("TALE_DNA_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readDNAStringSet(fasta)
  if (length(taleStrings) == 0L) {
    cli::cli_warn("No part sequence found in: {fasta}. Returning an empty tibble.")
    tbl <- tibble::tibble(array_id = character(),
                   domain_type = character(),
                   position_in_array = character(),
                   position_in_crd = character(),
                   string = character(),
                   source_directory = character())
    return(tbl)
  } 
  tbl <- tibble::tibble(array_id = gsub("(.*): .*", "\\1", names(taleStrings)),
                        domain_type = gsub(".*: (.*?)[ ]?[0-9]{0,}$", "\\1", names(taleStrings)),
                        position_in_crd = gsub(".*: repeat[ ]([0-9]{0,})$", "\\1", names(taleStrings)) %>%
                          as.integer() %>%
                          suppressWarnings(),
                        string = as.character(taleStrings) %>% as.vector(),
                        source_directory = dirname(fasta)
  )
  missingTerm <- setdiff(c("N-terminus", "C-terminus"), unique(tbl$domain_type))
  if (length(missingTerm) != 0L) {
    cli::cli_warn("Array {unique(tbl$array_id)} is missing a {missingTerm} domain in {fasta}")
    missingTerm <- tibble::tibble(array_id = unique(tbl$array_id),
                                  domain_type = missingTerm,
                                  position_in_crd = NA,
                                  string = NA,
                                  source_directory = dirname(fasta))
    tbl <- dplyr::bind_rows(tbl, missingTerm)
  }
  tbl %<>% 
    dplyr::rowwise() %>%
    dplyr::mutate(position_in_array = switch(domain_type,
                                                 `N-terminus` = 1,
                                                 `repeat` = position_in_crd + 1,
                                                 `C-terminus` = nrow(tbl)
                                                 ))
  return(tbl)
}

.rvds_from_annotale_file <- function(fasta) {
  if (!grepl("TALE_RVDs.fasta", basename(fasta))) {
    cli::cli_abort("The provided file does not seem to be an AnnoTALE RVDs file: {fasta}", class = c("tantale_error"))
  } else {
    rvdTble <- .split_list(fasta) %>%
      lapply(function(x) tibble::tibble(string = x,
                                        position_in_crd = 1:length(x))
             ) %>%
      dplyr::bind_rows(.id = "array_id")
    rvdTble <- rvdTble %>% dplyr::mutate(source_directory = dirname(fasta),
                                         domain_type = "repeat")
  }
  return(rvdTble)
}

#' Read TALE parts from a tell_tales output directory
#'
#' Implementation behind \code{\link{tales_from_telltale}} and the deprecated
#' \code{\link{tales_from_telltale}}. Internal so that package code can call it without
#' tripping the deprecation warning.
#' @noRd
.tale_parts <- function(telltale_dir) {
  # Get info from telltale output dir
  # !!!! array_id are assumed to be unique !!!!
  protPartsFiles <- list.files(telltale_dir, "TALE_Protein_parts.fasta", recursive = T, full.names = T)
  dnaPartsFiles <- list.files(telltale_dir, "TALE_DNA_parts.fasta", recursive = T, full.names = T)
  if (telltale_dir %>% dirname() %>% unique() %>% length() != 1L) {
    cli::cli_warn("The provided path most likely does not correspond to a SINGLE tell_tales output directory.")
  }
  # Fetch info from annotale/telltale files with .tale_parts_from_file
  taleProtString <- lapply(protPartsFiles, .tale_parts_from_file) %>% dplyr::bind_rows()
  taleDnaString <- lapply(dnaPartsFiles, .tale_parts_from_file) %>% dplyr::bind_rows()
  #stopifnot(nrow(taleProtString) == nrow(taleDnaString))
  # Join info in a table with one domain per row
  tale_parts <- dplyr::full_join(taleDnaString %>% dplyr::rename(dna_seq = string),
                                taleProtString %>% dplyr::rename(aa_seq = string),
                                by = c("array_id", "domain_type", "position_in_array", "position_in_crd", "source_directory"),
                                relationship = "one-to-one") %>%
    dplyr::mutate(aa_seq = gsub("[*]", "", aa_seq))

  # Get RVDs
  # NOTE: could be easier to get the RVDs directly from AnnoTALE output with
  # .rvds_from_annotale_file() but I currently feel that it is good to
  # be aware of disagreements between AnnoTALE diagnostic on terminal domains presence in AA seqs
  # and nhmmer diagnostic on terminal domains CDS presence on DNA.
  rvds <- .split_list(list.files(path = telltale_dir,
                                        pattern = "rvd_sequences.fas",
                                        recursive = F,
                                        full.names = T)
                             ) %>%
    lapply(function(x) tibble::tibble(rvd = x, position_in_array = 1:length(x) )) %>%
    dplyr::bind_rows(.id = "array_id")  
  anchorCodes <- tales_anchor_codes()
  
  
  # Some checks on the consistency between parts and rvd sequences
  # if nhmmer did not report on a C-Term CDS, the corresponding domain
  # "CTERM" tag will not be written in the rvd slot of the table.
  arraysConsistency <- dplyr::full_join(tale_parts %>% dplyr::count(array_id, name = "AnnoTALELength"),
                                        rvds %>% dplyr::count(array_id, name = "rvdFileLength"),
                                        by = dplyr::join_by(array_id)) %>%
    dplyr::mutate(sameLength = AnnoTALELength == rvdFileLength)
  
  if (any(is.na(arraysConsistency$sameLength))) {
    cli::cli_warn("There are mismatches in array IDs between rvd seq file and AnnoTALE parts files:")
    cli::cli_warn("There are mismatches in array IDs between rvd seq file and AnnoTALE parts files.")
  } else if (!all(arraysConsistency$sameLength, na.rm = TRUE)) {
    cli::cli_abort(
      c("Array lengths are inconsistent between the rvd seq file and the AnnoTALE parts files.",
        "i" = "Affected arrays: {.val {arraysConsistency$array_id[!arraysConsistency$sameLength]}}"),
      class = c("tantale_error_parts_inconsistent", "tantale_error"))
  }
  
  # Include RVDs in the talParts tibble
  tale_parts <- dplyr::left_join(tale_parts, 
                                rvds,
                                by = c("array_id", "position_in_array"),
                                unmatched = "drop", relationship = "one-to-one")
  # Include seqnames in the talParts tibble
  tale_parts %<>% dplyr::left_join(
    readr::read_tsv(list.files(telltale_dir, "hits_report.tsv", recursive = T, full.names = T),
                    show_col_types = FALSE) %>%
      dplyr::select(array_id, seqnames) %>%
      dplyr::distinct(),
    by = "array_id", relationship = "many-to-one"
  )
  # Check talparts
  partsWithMissingAaSeq <- tale_parts %>% dplyr::filter(is.na(aa_seq)) %>% dplyr::pull(array_id) %>% unique()
  partsWithMissingDnaSeq <- tale_parts %>% dplyr::filter(is.na(dna_seq)) %>% dplyr::pull(array_id) %>% unique()
  partsWithMissingRvdSeq <- tale_parts %>% dplyr::filter(is.na(rvd)) %>% dplyr::pull(array_id) %>% unique()
  if (any(sapply(list(partsWithMissingAaSeq, partsWithMissingDnaSeq, partsWithMissingRvdSeq), length) != 0L)) {
    cli::cli_warn(c("The returned tale_parts has records with missing sequences.",
                    "i" = "Affected array{?s}: {.val {unique(c(partsWithMissingAaSeq, partsWithMissingDnaSeq, partsWithMissingRvdSeq))}}"))
  }
  return(tale_parts)
}


#' Build a tales object from a tell_tales run directory
#'
#' Reads the AnnoTALE/telltale part files of a single
#' \code{\link[tantale:tell_tales]{tell_tales}} output directory and returns a
#' validated \code{\link{tales}} object.
#'
#' The result carries no \code{dom_code}: that surrogate key is minted later,
#' by the relatedness computation, over the whole set of parts being analysed.
#'
#' @param sanitize If \code{TRUE}, arrays carrying biological anomalies are
#'   removed with a warning naming them and why; if \code{FALSE} (default) they
#'   are kept and merely warned about. See \code{\link{tales_anomalies}}.
#' @param telltale_dir Path to a single \code{\link[tantale:tell_tales]{tell_tales}}
#'   output directory.
#' @return A validated \code{tales} object.
#' @export
#' @family TALE discovery
#' @examples
#' tales_from_telltale(system.file("extdata", "tellTaleExampleOutput",
#'                                 package = "tantale"))
tales_from_telltale <- function(telltale_dir, sanitize = FALSE) {
  tales(.tale_parts(telltale_dir), sanitize = sanitize)
}
