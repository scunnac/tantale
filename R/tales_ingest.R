#### Building a tales object from AnnoTALE output ####
#
# tales_from_telltales() and tales_from_annotale() are the public entry points;
# everything else here is a private step of them. Kept together, and out of
# distalr.R, so that file is about comparing TALEs rather than reading them
# off disk.
#
# Both read AnnoTALE's split of each ORF into N-terminus, repeats and
# C-terminus (as protein and as DNA) and its RVDs. tell_tales() writes one
# directory per region of interest under annotale/, plus array_report.tsv,
# which says whether each terminus matches the TALE terminal-domain protein
# profile. run_annotale_predict() writes one set of files for all TALEs in
# Analyze/, and the termini are searched with the profiles on reading.

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
  # AnnoTALE's own names carry the TALE's location after a space
  # ("MAI1-tempTALE1 [624136-627961:1]"); the id is the first word
  tibble::tibble(array_id = sub(" .*$", "", gsub("(.*): .*", "\\1", names(taleStrings))),
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
    cli::cli_abort("{.file {fasta}} is not an AnnoTALE RVD file ({.file TALE_RVDs.fasta}).",
                   class = c("tantale_error_annotale_file", "tantale_error"))
  } else {
    # AnnoTALE writes an empty record for a TALE in which it found no repeat
    # (seen on pseudogenes); that TALE simply has no RVD
    rvdTble <- withCallingHandlers(
      .split_list(fasta),
      tantale_warning_empty_element = function(w) invokeRestart("muffleWarning"))
    names(rvdTble) <- sub(" .*$", "", names(rvdTble))
    rvdTble <- rvdTble %>%
      lapply(function(x) tibble::tibble(string = x,
                                        position_in_crd = seq_along(x))
             ) %>%
      dplyr::bind_rows(.id = "array_id")
    rvdTble <- rvdTble %>% dplyr::mutate(source_directory = dirname(fasta),
                                         domain_type = "repeat")
  }
  return(rvdTble)
}


#' Read TALE parts from a tell_tales output directory
#'
#' Implementation behind \code{\link{tales_from_telltales}}.
#' @noRd
.tale_parts <- function(telltale_dir) {
  # !!!! array_id are assumed to be unique !!!!
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
  tale_parts <- .tale_parts_assemble(telltale_dir)
  seqnames <- readr::read_tsv(list.files(telltale_dir, "hits_report.tsv", recursive = T, full.names = T),
                              show_col_types = FALSE) %>%
    dplyr::select(array_id, seqnames) %>%
    dplyr::distinct()
  .tale_parts_finish(tale_parts,
                     aa_hits = dplyr::select(arrayReport, array_id, nterm_aa_hit, cterm_aa_hit),
                     seqnames = seqnames)
}


#' Read TALE parts from AnnoTALE's analyze output
#'
#' Implementation behind \code{\link{tales_from_annotale}}.
#' @noRd
.tale_parts_annotale <- function(annotale_dir, terminus_max_evalue, hmm_dir) {
  tale_parts <- .tale_parts_assemble(annotale_dir)
  termini <- lapply(c(`N-terminus` = "N-terminus", `C-terminus` = "C-terminus"), function(part) {
    x <- tale_parts[tale_parts$domain_type == part, ]
    stats::setNames(Biostrings::AAStringSet(x$aa_seq), x$array_id)
  })
  aaHits <- .tale_termini_hmmsearch(termini, max_evalue = terminus_max_evalue, hmm_dir = hmm_dir)
  # the contig of each TALE, from predict's GFF3 when it is there
  gffFiles <- list.files(annotale_dir, "^GFF__.*\\.gff3$", recursive = TRUE, full.names = TRUE)
  seqnames <- if (length(gffFiles) > 0L) {
    lapply(gffFiles, function(f) {
      gff <- utils::read.delim(f, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
      gff <- gff[gff$V3 == "mRNA", ]
      tibble::tibble(array_id = sub("^.*Id=([^;]+).*$", "\\1", gff$V9), seqnames = gff$V1)
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::distinct()
  }
  .tale_parts_finish(tale_parts, aa_hits = aaHits, seqnames = seqnames)
}


#' Assemble AnnoTALE's parts and RVDs into one table
#'
#' Reads every \code{TALE_Protein_parts.fasta}, \code{TALE_DNA_parts.fasta}
#' and \code{TALE_RVDs.fasta} under \code{dir}, joins them part by part and
#' numbers the parts along each array. The termini get no code yet.
#'
#' @param dir A tell_tales output directory, or a directory holding
#'   AnnoTALE's analyze output.
#' @return A tibble, one row per part.
#' @noRd
.tale_parts_assemble <- function(dir) {
  protPartsFiles <- list.files(dir, "TALE_Protein_parts.fasta", recursive = T, full.names = T)
  dnaPartsFiles <- list.files(dir, "TALE_DNA_parts.fasta", recursive = T, full.names = T)
  rvdFiles <- list.files(dir, "TALE_RVDs.fasta", recursive = T, full.names = T)
  if (length(protPartsFiles) == 0L) {
    cli::cli_abort(
      c("No {.file TALE_Protein_parts.fasta} under {.file {dir}}.",
        "i" = "Expected the output of {.fn tell_tales} or of AnnoTALE's analyze stage."),
      class = c("tantale_error_annotale_missing", "tantale_error"))
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

  dplyr::bind_rows(repeats, tale_parts %>% dplyr::filter(domain_type != "repeat"))
}


#' Code the termini, attach the contigs and check the result
#'
#' @param tale_parts What \code{.tale_parts_assemble()} returned.
#' @param aa_hits One row per array: \code{array_id}, \code{nterm_aa_hit},
#'   \code{cterm_aa_hit}, as \code{.tale_termini_hmmsearch()} returns them.
#' @param seqnames \code{array_id} and \code{seqnames}, or \code{NULL}.
#' @return The parts, as \code{tales()} takes them.
#' @noRd
.tale_parts_finish <- function(tale_parts, aa_hits, seqnames = NULL) {
  # NTERM/CTERM for a terminus matching its protein profile, XXXXX otherwise
  anchors <- unname(tales_anchor_codes())   # NTERM, CTERM, XXXXX
  repeats <- tale_parts %>% dplyr::filter(domain_type == "repeat")
  termini <- tale_parts %>%
    dplyr::filter(domain_type != "repeat") %>%
    dplyr::select(-dplyr::any_of("rvd")) %>%
    dplyr::left_join(dplyr::select(aa_hits, array_id, nterm_aa_hit, cterm_aa_hit),
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

  if (!is.null(seqnames)) {
    tale_parts <- dplyr::left_join(tale_parts, seqnames, by = "array_id",
                                   relationship = "many-to-one")
  }
  # Check talparts
  partsWithMissingAaSeq <- tale_parts %>% dplyr::filter(is.na(aa_seq)) %>% dplyr::pull(array_id) %>% unique()
  partsWithMissingDnaSeq <- tale_parts %>% dplyr::filter(is.na(dna_seq)) %>% dplyr::pull(array_id) %>% unique()
  partsWithMissingRvdSeq <- tale_parts %>% dplyr::filter(is.na(rvd)) %>% dplyr::pull(array_id) %>% unique()
  if (any(lengths(list(partsWithMissingAaSeq, partsWithMissingDnaSeq, partsWithMissingRvdSeq)) != 0L)) {
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
#' tales_from_telltales(system.file("extdata", "tellTaleExampleOutput",
#'                                  package = "tantale"))
tales_from_telltales <- function(telltale_dir, sanitize = FALSE) {
  tales(.tale_parts(telltale_dir), sanitize = sanitize)
}


#' Build a tales object from AnnoTALE's own TALE predictions
#'
#' Reads the output of AnnoTALE's analyze stage, as written by
#' \code{\link{run_annotale_predict}}, and returns a validated
#' \code{\link{tales}} object with the same columns as
#' \code{\link{tales_from_telltales}}.
#'
#' AnnoTALE predict finds TALE genes in a genome; analyze splits each one into
#' its N-terminal region, repeats and C-terminal region, and reads the RVD of
#' every repeat. \code{\link{tell_tales}} runs analyze too, on ORFs it
#' delimits itself from nhmmer hits of the TALE DNA profiles, so the two need
#' not report the same set of TALEs.
#'
#' As in \code{tell_tales()}, each terminal segment is searched with the TALE
#' N- or C-terminal protein profile (\code{hmmsearch}, from the tantale
#' environment; see \code{\link{tantale_setup}}). The \code{rvd} column holds
#' \code{NTERM}/\code{CTERM} for a segment that matches its profile and
#' \code{XXXXX} for one that does not (see \code{\link{tales_anchor_codes}}).
#'
#' \code{array_id} is AnnoTALE's name for the TALE (\code{MAI1-tempTALE1}),
#' without the location AnnoTALE appends to it. \code{seqnames} is filled
#' when predict's GFF3 file is found under \code{annotale_dir}.
#'
#' @param annotale_dir A directory holding AnnoTALE analyze's
#'   \code{TALE_Protein_parts.fasta}, \code{TALE_DNA_parts.fasta} and
#'   \code{TALE_RVDs.fasta}, in it or in a subdirectory: the
#'   \code{output_dir} of \code{run_annotale_predict()} will do.
#' @param terminus_max_evalue Maximum \code{hmmsearch} E-value for a
#'   terminal segment to be coded \code{NTERM}/\code{CTERM}. As in
#'   \code{\link{tell_tales}}, the match must also reach the end of the
#'   profile that adjoins the repeats.
#' @param hmm_dir Directory holding the TALE terminus protein profiles.
#' @inheritParams tales_from_telltales
#' @return A validated \code{tales} object.
#' @export
#' @family TALE discovery
#' @examples
#' \donttest{
#' # Needs the tantale environment for hmmsearch (see tantale_setup()).
#' tales_from_annotale(system.file("extdata", "annotaleExampleOutput",
#'                                 package = "tantale"))
#' }
tales_from_annotale <- function(annotale_dir, terminus_max_evalue = 1e-5, sanitize = FALSE,
                                hmm_dir = system.file("extdata", "hmmProfile", package = "tantale",
                                                      mustWork = TRUE)) {
  tales(.tale_parts_annotale(annotale_dir, terminus_max_evalue = terminus_max_evalue,
                             hmm_dir = hmm_dir),
        sanitize = sanitize)
}
