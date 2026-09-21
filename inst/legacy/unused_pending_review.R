# Retired: functions with no caller in the package
#
# Moved here from R/unused_pending_review.R, where they had sat parked
# (sourced, checked, callable with `tantale:::`, but not called by anything)
# pending a one-at-a-time review of whether each was actually useless. That
# review is now done: none of the six below has any caller in R/ or any test
# dependency, so they move from "parked" to "retired" -- no longer sourced or
# compiled as part of the package (dev/restructuring-notes.md ledger §17).
#
# Two siblings that started out in the same parked file, repeat_to_rvd_align()
# and .rvd_to_repeat_align(), did NOT move here: both are still real test
# dependencies (test_plot_tales_msa.R's fixture setup, and
# test_error_conditions.R's own tests of the second one), so they were kept
# as real, active code in R/conversion.R instead.
#
# Per function, what was known when this was still R/unused_pending_review.R:
#
# .write_hmm_file()        wraps hmmbuild. Never called. tell_tales() uses the
#                          pre-built profile in inst/extdata/hmmProfile, so
#                          nothing in the package ever builds one.
# .run_hmmer_search()      wraps hmmsearch (protein). Never called; tell_tales()
#                          uses .run_nhmmer_search(), the nucleotide sibling,
#                          which stays in R/tellTale_utilities.R.
# .run_hmmalign()          wraps hmmalign. Never called.
# .extract_seqs_from_hits() Never called. Its one call, `XVector::subseq()`,
#                          is a gap: `XVector` was dropped from DESCRIPTION's
#                          Imports 2026-09-20 (ledger §14) since this was its
#                          only reference anywhere in the package. Re-add it
#                          before ever reviving this function.
# .tales_compare_core()    Superseded by tales_compare_distal(), which now
#                          composes tales_assign_domain_codes(),
#                          tales_domain_distances() and tales_tale_distances()
#                          instead of calling this (ledger §8.5b). Kept as a
#                          reference for what the composed version had to
#                          reproduce; the golden baseline is what confirmed it
#                          does.
# .run_in_conda()          Superseded by .tantale_bin()/.tantale_exec()
#                          (ledger §12): every caller now resolves the
#                          executable to an absolute path inside the
#                          environment instead of going through `conda run`,
#                          which had two real faults -- `-n` by name (ledger
#                          §7.4a) and only the first command of a compound
#                          string running inside the environment.
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

.tales_compare_core <- function(tale_parts, ncores = 1,
                                aln_method = "DECIPHER", conda_bin = "auto") {

  #### Reality checks ####

  ## Make sure we are dealing only with parts that have defined protein sequences.
  if (any(is.na(tale_parts$aa_seq) | tale_parts$aa_seq == "")) {
    # must match the guard above, or an empty-string part lists nothing
    badArrays <- unique(tale_parts$array_id[is.na(tale_parts$aa_seq) | tale_parts$aa_seq == ""])
    cli::cli_abort(
      c("Some of the provided TALE parts have no amino acid sequence.",
        "i" = "Affected array{?s}: {.val {badArrays}}"),
      class = c("tantale_error_parts_no_aa", "tantale_error"))
  }
  if (any(is.na(tale_parts$dna_seq) | tale_parts$dna_seq == "")) {
    cli::cli_warn("It seems that some of the provided TALE parts miss the DNA sequence!")
  }

  # Assign domain codes
  tale_parts %<>% dplyr::group_by(aa_seq) %>%
    dplyr::mutate(dom_code = dplyr::cur_group_id() %>% unlist() %>% as.character()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dom_code = dplyr::if_else(is.na(aa_seq), as.character(NA), dom_code))


  #### Assemble repeat code strings and write in a file for arlem ####
  cli::cli_inform("Assemble repeat code TALE strings and write in a file for ARLEM")

  codesSeqsfile <- tempfile(fileext = ".fasta")
  repeatStrings <- tale_parts %>%
    dplyr::group_by(array_id) %>%
    dplyr::arrange(position_in_array) %>%
    dplyr::summarise(repeatString = paste(dom_code, collapse = " "),
                     posString = paste(position_in_array, collapse = " "))
  codesSeqSet <- Biostrings::BStringSet(repeatStrings$repeatString)
  names(codesSeqSet) <- repeatStrings$array_id

  # # Must use seqinr because Biostrings wraps sequences in fasta file which messes up Arlem...
  # codesSeqLst <- as.list(repeatStrings$repeatString)
  # names(codesSeqLst) <- repeatStrings$array_id
  # seqinr::write.fasta(codesSeqLst, names = names(codesSeqLst),
  #                     file.out = codesSeqsfile, as.string = TRUE, nbchar = 10000)

  Biostrings::writeXStringSet(x = codesSeqSet,
                              filepath = codesSeqsfile,
                              format = "fasta", width = 20000L)



  #### Compute systematic pairwise dissimilarities (distances) between 'repeat' units. ####
  # Get unique domains sequences
  taleAaParts <- Biostrings::AAStringSet(tale_parts$aa_seq)
  names(taleAaParts) <- tale_parts$dom_code
  uniqueTaleAaParts <- unique(taleAaParts)
  stopifnot(!anyDuplicated(names(uniqueTaleAaParts)))
  stopifnot(!anyDuplicated(names(unique(taleAaParts))))

  # Get pairwise repeat aa sequence dissimilarity scores in a long tibble
  cli::cli_inform(paste0("Computing a distance matrix between TALE parts amino acid sequences ",
                   "using: {aln_method}"))
  if (aln_method == "mmseq2") {
    dissimLong <- .pairwise_align_mmseq2(part_aa_set = uniqueTaleAaParts, ncores = ncores,
                                        conda_bin = conda_bin)
    #saveRDS(dissimLong, file = "/home/cunnac/TEMP/dissimLong")
  } else if (aln_method == "Biostrings") {
    dissimLong <- .pairwise_align_biostrings(part_aa_set = uniqueTaleAaParts, ncores = ncores)
  } else if (aln_method == "DECIPHER") {
    dissimLong <- .pairwise_align_decipher(part_aa_set = uniqueTaleAaParts, ncores = ncores)
  } else {
    cli::cli_abort("{.arg aln_method} must be one of {.val Biostrings}, {.val mmseq2} or {.val DECIPHER}, not {.val {aln_method}}.",
                   class = c("tantale_error_aln_method", "tantale_error"))
  }
  dissimLong %<>% dplyr::mutate(sim = 100 - dissim)
  # Convert Distance (dissimilarity) measures to Similarity with a four-parameter logistic function
  # pair_align_scores %<>% dplyr::mutate(Sim = 100/(1+exp(-1*-0.9*(Dissim-3))))


  # Convert to square matrix
  dissimMat <- .pairwise_long_to_matrix(dissimLong, "dissim")
  stopifnot(nrow(dissimMat) == ncol(dissimMat))
  # reorder row and colnames because I suspect arlem expect them in increasing order
  dissimMat <- dissimMat[rownames(dissimMat) %>% as.numeric() %>% order(),
                         colnames(dissimMat) %>% as.numeric() %>% order()]

  #### Generate an ARLEM cost matrix ####
  if (TRUE) {
    method <- "minkowski"
    cli::cli_inform(paste0("Generate an ARLEM cost matrix which meets triangle inequality criteria by computing ",
                     "the {method} distance between pairwise distance vectors."))
    dissimMat <- as.matrix(stats::dist(dissimMat, method = method, p = 3.5, diag = TRUE, upper = TRUE))
    dissimMat <- dissimMat/max(dissimMat) * 100
  }
  # if (!fossil::tri.ineq(dissimMat)) {
  #   logger::log_error("TALE domains dissimilarity (distance) matrix does not respect the triangle inequality",
  #                     "Arlem will fail. Aborting...")
  #   stop()
  # }

  #### Prepare Arlem cfile with systematic pairwise distances between 'repeat' units. ####
  # Get parameters for arlem
  TypeNo <- glue::glue("# Type no. ", nrow(dissimMat))
  Types <- glue::glue("# Types ", paste(1:nrow(dissimMat), collapse = " "))
  # Convert mat to character, beware of the ceiling in conversion...
  dissimMat <- matrix(ceiling(dissimMat) %>% format(),
                      ncol = ncol(dissimMat),
                      dimnames = list(rownames(dissimMat), colnames(dissimMat))
                      )
  # 'erase' lower triangle and diag
  dissimMat[lower.tri(dissimMat)] <- ""
  diag(dissimMat) <- ""

  # Convert mat rows to strings of space separated values
  dissimMatLines <- apply(dissimMat, 1, function(row) {
    string <- paste(row, collapse = " ")
    gsub("^[ ]+", "", string)
  })
  # remove empty last line
  dissimMatLines <- dissimMatLines[1:(length(dissimMatLines) - 1)]
  # Add the arlem stuff
  dissimMatLines <- c(TypeNo, Types,
                      "# Indel align 10", "# Indel hist 10", "# Dup 10",
                      "# matrix",
                      dissimMatLines)

  # write lines in a temp cfile
  cfile <- tempfile()
  writeLines(cfile, text = dissimMatLines)

  #### run arlem with a system call and parse std output ####
  arlemPath <- system.file("tools", "arlem", "arlem", package = "tantale", mustWork = T)
  arlemCmd <- glue::glue("{shQuote(arlemPath)} -f {shQuote(codesSeqsfile)} -cfile {shQuote(cfile)} -align -insert -showalign")
  cli::cli_inform("Running ARLEM version 1.0 : ")
  cli::cli_inform("Copyright by Mohamed I. Abouelhoda")
  cli::cli_inform(paste0("Plz. cite Abouelhoda, Giegerich, Behzadi, and Steyaert"))
  arlemRawRes <- system(arlemCmd, intern = TRUE)
  arlemSelfRes <- grep("Processed Seq[.]:", arlemRawRes, value = TRUE)
  arlemSelfScores <- gsub("Processed Seq[.]: ([0-9]{1,}) Score: ([0-9]{1,}),.*", "\\1|\\2",
                          substring(arlemSelfRes, 1, 35)) %>%
    strsplit(split = "\\|") %>%
    lapply(function(s) t(as.matrix(as.numeric(s)))) %>%
    do.call(rbind, .) %>% tibble::as_tibble(.name_repair = "minimal")
  arlemRes <- grep("Score of aligning Seq:",
                   arlemRawRes, value = TRUE)
  arlemScores <- gsub("Score of aligning Seq:([0-9]+), Seq:([0-9]+) =([0-9]+)", "\\1|\\2|\\3",
                      arlemRes)
  arlemScores <- strsplit(arlemScores, split = "\\|")
  arlemScores <- lapply(arlemScores, function(s) {t(as.matrix(as.numeric(s)))}) %>%
    do.call(rbind, .) %>%
    tibble::as_tibble(.name_repair = "minimal")
  colnames(arlemScores) <- c("id1", "id2", "arlem_score")
  # Shaping into matrix to have scores in both directions (fill diag and triangle)
  arlemScoresMat <- .pairwise_long_to_matrix(arlemScores, "arlem_score")
  arlemScoresMat <- cbind("0" = NA, arlemScoresMat)
  arlemScoresMat <- rbind(arlemScoresMat, NA)
  rownames(arlemScoresMat)[length(codesSeqSet)] <- length(codesSeqSet) - 1
  arlemScores <- stats::as.dist(t(arlemScoresMat), diag = TRUE, upper = TRUE) %>% as.matrix() %>%
    .matrix_to_long(value_name = "arlem_score") %>%
    tibble::as_tibble()
  colnames(arlemScores) <- c("id1", "id2", "arlem_score")

  #### Compute normalized arlem scores and include array_ids rather than arlem index ####
  arrayLengths <- tale_parts %>% dplyr::group_by(array_id) %>% dplyr::count()

  normArlemScoresTble <- arlemScores %>%
    dplyr::mutate(
      id1 = names(codesSeqSet)[id1 + 1],
      id2 = names(codesSeqSet)[id2 + 1]
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      max_length = max(arrayLengths$n[arrayLengths$array_id == id1],
                       arrayLengths$n[arrayLengths$array_id == id2]),
      norm_arlem_score = arlem_score/max_length
    ) %>%
    dplyr::ungroup()

  #### Check features of the Arlem results table
  arraysCount <- codesSeqSet %>% length()
  if (nrow(normArlemScoresTble) != arraysCount^2) {
    allCombs <- expand.grid(names(codesSeqSet), names(codesSeqSet), stringsAsFactors = FALSE) %>% tibble::as_tibble()
    colnames(allCombs) <- c("id1", "id2")
    absentCombs <- dplyr::left_join(allCombs, normArlemScoresTble) %>%
      dplyr::filter(is.na(norm_arlem_score))
      cli::cli_abort(
        c("The TALE similarity table does not have the expected number of comparisons.",
          "x" = "Expected {arraysCount^2}, got {nrow(normArlemScoresTble)}; {nrow(absentCombs)} missing.",
          "i" = "First missing pair{?s}: {.val {paste(utils::head(absentCombs$id1, 3), utils::head(absentCombs$id2, 3), sep = \"/\")}}"),
        class = c("tantale_error_arlem_incomplete", "tantale_error"))
  }


  #### return the raw pieces; callers class and assemble them ####
  cli::cli_inform("Finished computing TALE and repeat relatedness.")
  list(
    tale_parts = tale_parts,
    dissim_long = dissimLong,
    tal_sim = normArlemScoresTble,
    coded_seq_set = codesSeqSet
  )
}

.run_in_conda <- function(env_name, command,
                             conda_bin = "auto",
                             cwd = getwd(),
                             ...) {
  conda_bin <- reticulate::conda_binary(conda_bin)
  # activateEnvCmd <- glue::glue("eval \"$({conda_bin} shell hook -s posix)\"",
  #                              "; micromamba activate {env_name}")
  # fullCommand <- glue::glue_collapse(c(activateEnvCmd, command), sep = "; ")
  fullCommand <- glue::glue("eval \"$({conda_bin} shell hook -s posix)\"",
                            "{conda_bin} run --cwd {cwd} -n {env_name} {command}",
                            .sep = "; ")
  system(command = fullCommand, ...)
}
