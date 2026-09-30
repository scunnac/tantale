#### Building a tales object from a tell_tales run directory ####
#
# tales_from_telltale() is the public entry point; everything else here is a
# private step of it. Kept together, and out of distalr.R, so that file is
# about comparing TALEs rather than reading them off disk.
#
# The layout being read is what tell_tales() writes: one directory per
# region of interest under annotale/, holding AnnoTALE's split of the ORF
# into N-terminus, repeats and C-terminus (as protein and as DNA) and its
# RVDs, plus array_report.tsv, which says whether each terminus matches the
# TALE terminal-domain protein profile.

.tale_parts_from_file <- function(fasta) {
  if (grepl("TALE_Protein_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readAAStringSet(fasta)
  if (grepl("TALE_DNA_parts.fasta", basename(fasta))) taleStrings <- Biostrings::readDNAStringSet(fasta)
  if (length(taleStrings) == 0L) {
    cli::cli_warn("No part sequence found in: {fasta}. Returning an empty tibble.")
    tbl <- tibble::tibble(array_id = character(),
                   domain_type = character(),
                   position_in_crd = integer(),
                   string = character(),
                   source_directory = character())
    return(tbl)
  } 
  tibble::tibble(array_id = gsub("(.*): .*", "\\1", names(taleStrings)),
                 domain_type = gsub(".*: (.*?)[ ]?[0-9]{0,}$", "\\1", names(taleStrings)),
                 position_in_crd = gsub(".*: repeat[ ]([0-9]{0,})$", "\\1", names(taleStrings)) %>%
                   as.integer() %>%
                   suppressWarnings(),
                 string = as.character(taleStrings) %>% as.vector(),
                 source_directory = dirname(fasta)
  )
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
#' Implementation behind \code{\link{tales_from_telltale}}.
#' @noRd
.tale_parts <- function(telltale_dir) {
  # Get info from telltale output dir
  # !!!! array_id are assumed to be unique !!!!
  protPartsFiles <- list.files(telltale_dir, "TALE_Protein_parts.fasta", recursive = T, full.names = T)
  dnaPartsFiles <- list.files(telltale_dir, "TALE_DNA_parts.fasta", recursive = T, full.names = T)
  rvdFiles <- list.files(telltale_dir, "TALE_RVDs.fasta", recursive = T, full.names = T)
  if (telltale_dir %>% dirname() %>% unique() %>% length() != 1L) {
    cli::cli_warn("The provided path most likely does not correspond to a SINGLE tell_tales output directory.")
  }
  arrayReport <- readr::read_tsv(file.path(telltale_dir, "array_report.tsv"), show_col_types = FALSE)
  if (!all(c("nterm_aa_hit", "cterm_aa_hit") %in% names(arrayReport))) {
    cli::cli_abort(
      c("{.file {telltale_dir}} was written by an older version of {.fn tell_tales}.",
        "x" = "Its {.file array_report.tsv} has no {.field nterm_aa_hit}/{.field cterm_aa_hit} column, which the terminus codes are read from.",
        "i" = "Run {.fn tell_tales} again on the same sequences."),
      class = c("tantale_error_telltale_outdated", "tantale_error"))
  }

  # Fetch info from annotale/telltale files with .tale_parts_from_file
  taleProtString <- lapply(protPartsFiles, .tale_parts_from_file) %>% dplyr::bind_rows()
  taleDnaString <- lapply(dnaPartsFiles, .tale_parts_from_file) %>% dplyr::bind_rows()
  # Join info in a table with one domain per row
  tale_parts <- dplyr::full_join(taleDnaString %>% dplyr::rename(dna_seq = string),
                                 taleProtString %>% dplyr::rename(aa_seq = string),
                                 by = c("array_id", "domain_type", "position_in_crd", "source_directory"),
                                 relationship = "one-to-one") %>%
    dplyr::mutate(aa_seq = gsub("[*]", "", aa_seq))

  # A part in only one of AnnoTALE's two files means the files disagree about
  # that array. The whole array goes: dropping only the part would leave a
  # hole in the array and shift every position after it.
  incomplete <- tale_parts %>%
    dplyr::filter(is.na(dna_seq) | is.na(aa_seq)) %>%
    dplyr::pull(array_id) %>%
    unique()
  if (length(incomplete) != 0L) {
    cli::cli_warn(
      c("AnnoTALE's protein and DNA parts disagree for {length(incomplete)} array{?s}, left out of the result.",
        "i" = "Affected array{?s}: {.val {incomplete}}"),
      class = c("tantale_warning_parts_inconsistent", "tantale_warning"))
    tale_parts <- tale_parts %>% dplyr::filter(!array_id %in% incomplete)
  }

  # When AnnoTALE reports no terminus on one side, there is probably nothing
  # there, and the array has no part on that side.
  sides <- tale_parts %>%
    dplyr::group_by(array_id) %>%
    dplyr::summarise(no_nterm = !any(domain_type == "N-terminus"),
                     no_cterm = !any(domain_type == "C-terminus"))
  noNterm <- sides$array_id[sides$no_nterm]
  noCterm <- sides$array_id[sides$no_cterm]
  if (length(noNterm) + length(noCterm) != 0L) {
    cli::cli_warn(
      c("AnnoTALE reported no terminus on one side of the repeats for some arrays.",
        if (length(noNterm)) c("i" = "No N-terminus: {.val {noNterm}}"),
        if (length(noCterm)) c("i" = "No C-terminus: {.val {noCterm}}")),
      class = c("tantale_warning_terminus_absent", "tantale_warning"))
  }

  # Positions are counted, since a terminus may be absent
  tale_parts <- tale_parts %>%
    dplyr::arrange(array_id,
                   match(domain_type, c("N-terminus", "repeat", "C-terminus")),
                   position_in_crd) %>%
    dplyr::group_by(array_id) %>%
    dplyr::mutate(position_in_array = dplyr::row_number()) %>%
    dplyr::ungroup()

  # RVDs of the repeats, as AnnoTALE read them. AnnoTALE writes one RVD per
  # repeat part, so any repeat without its RVD, or the reverse, is a bug.
  rvds <- lapply(rvdFiles, .rvds_from_annotale_file) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(!array_id %in% incomplete) %>%
    dplyr::select(array_id, domain_type, position_in_crd, rvd = string)
  repeats <- tale_parts %>% dplyr::filter(domain_type == "repeat")
  rvdKeys <- c("array_id", "domain_type", "position_in_crd")
  unmatched <- unique(c(dplyr::anti_join(repeats, rvds, by = rvdKeys)$array_id,
                        dplyr::anti_join(rvds, repeats, by = rvdKeys)$array_id))
  if (length(unmatched) != 0L) {
    cli::cli_abort(
      c("AnnoTALE's RVDs and repeat parts do not correspond.",
        "i" = "Affected array{?s}: {.val {unmatched}}"),
      class = c("tantale_error_parts_inconsistent", "tantale_error"))
  }
  repeats <- dplyr::left_join(repeats, rvds, by = rvdKeys, relationship = "one-to-one")

  # Terminus codes, decided by tell_tales() from the protein profile search
  anchors <- unname(tales_anchor_codes())   # NTERM, CTERM, XXXXX
  termini <- tale_parts %>%
    dplyr::filter(domain_type != "repeat") %>%
    dplyr::left_join(dplyr::select(arrayReport, array_id, nterm_aa_hit, cterm_aa_hit),
                     by = "array_id", relationship = "many-to-one") %>%
    dplyr::mutate(
      is_hit = dplyr::if_else(domain_type == "N-terminus", nterm_aa_hit, cterm_aa_hit),
      rvd = dplyr::case_when(is.na(is_hit) ~ NA_character_,
                             !is_hit ~ anchors[3],
                             domain_type == "N-terminus" ~ anchors[1],
                             TRUE ~ anchors[2])) %>%
    dplyr::select(-nterm_aa_hit, -cterm_aa_hit, -is_hit)

  tale_parts <- dplyr::bind_rows(repeats, termini) %>%
    dplyr::arrange(array_id, position_in_array) %>%
    dplyr::relocate(position_in_array, .before = aa_seq)

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
    cli::cli_warn(c("The returned {.fn tantale::tales} object has records with missing sequences.",
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
#' One row per part: the N-terminus, each repeat and the C-terminus, as
#' AnnoTALE split the array's longest ORF. The \code{rvd} column holds the RVD
#' of a repeat, or a terminus code (see \code{\link{tales_anchor_codes}}):
#' \code{NTERM}/\code{CTERM} when the terminus matches the TALE
#' terminal-domain protein profile, \code{XXXXX} when it does not. When
#' AnnoTALE reported no terminus on one side, the array has no part there, with
#' a warning. An array whose protein and DNA parts disagree is left out, with
#' a warning.
#'
#' The directory must have been written by the current version of
#' \code{tell_tales()}, whose \code{array_report.tsv} holds the terminus
#' check; an older one is an error.
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
