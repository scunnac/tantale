##### The 'tales_msa' class #####
# A 'tales' that additionally carries an alignment coordinate. Gaps are
# *implicit*: a gap at (array, alignment_position) is the absence of a row, so
# every row remains a real part and every 'tales' invariant holds unchanged.
# See dev/class-design.md §4.


#### Constructors ####

#' Low-level constructor for a tales_msa object
#'
#' Attaches the class without validating. Use \code{\link{tales_msa}} unless the
#' invariants are already established.
#'
#' @param x A \code{tales} object carrying an \code{alignment_position} column.
#' @param alignment_width Integer width of the alignment. Stored as an
#'   attribute rather than derived, because subsetting arrays can empty the last
#'   column and would silently shrink a derived value.
#' @return A \code{tales_msa} object.
#' @keywords internal
new_tales_msa <- function(x, alignment_width = NULL) {
  stopifnot(is.data.frame(x))
  if (!is_tales(x)) x <- new_tales(x)
  if (!is.null(alignment_width)) {
    attr(x, "alignment_width") <- as.integer(alignment_width)
  }
  class(x) <- unique(c("tales_msa", class(x)))
  x
}

#' Is this a tales_msa object?
#' @param x An object.
#' @return A logical scalar.
#' @export
#' @family TALE alignment
is_tales_msa <- function(x) inherits(x, "tales_msa")

#' Width of a TALE alignment
#'
#' The number of columns in the alignment, including those that are all gaps in
#' the object at hand. Stored rather than derived: subsetting arrays can empty
#' the last column, which would silently shrink \code{max(alignment_position)}.
#'
#' @param x A \code{tales_msa} object.
#' @return An integer scalar, or \code{NULL} if unset.
#' @export
#' @family TALE alignment
tales_width <- function(x) {
  attr(x, "alignment_width", exact = TRUE)
}

#' Create a tales_msa object
#'
#' A \code{\link{tales}} object plus an \code{alignment_position} column. Gaps
#' are implicit: a gap is simply the absence of a row at that
#' (\code{array_id}, \code{alignment_position}).
#'
#' Normally produced by \code{\link{tales_align}} rather than called directly.
#'
#' @param x A data frame with the \code{\link{tales}} columns plus
#'   \code{alignment_position}.
#' @param alignment_width Integer alignment width; defaults to
#'   \code{max(alignment_position)}.
#' @param dom_code_namespace Optional scalar string, see
#'   \code{\link{tales_namespace}}.
#' @return A validated \code{tales_msa} object.
#' @export
#' @family TALE alignment
tales_msa <- function(x, alignment_width = NULL, dom_code_namespace = NULL) {
  x <- tales(x, dom_code_namespace = dom_code_namespace)
  if ("alignment_position" %in% names(x) && is.numeric(x$alignment_position)) {
    x$alignment_position <- as.integer(x$alignment_position)
  }
  if (is.null(alignment_width) && "alignment_position" %in% names(x) &&
      nrow(x) > 0L) {
    alignment_width <- max(x$alignment_position, na.rm = TRUE)
  }
  validate_tales_msa(new_tales_msa(x, alignment_width = alignment_width))
}


#### Validator ####

#' Validate a tales_msa object
#'
#' Checks every \code{\link{validate_tales}} invariant, then those specific to
#' an alignment. As for \code{tales}, only properties closed under row
#' subsetting are checked here; grid completeness is a precondition of the
#' functions that need it.
#'
#' @param x A \code{tales_msa} object.
#' @return \code{x}, invisibly, if valid; otherwise an error.
#' @export
#' @family TALE alignment
validate_tales_msa <- function(x) {
  validate_tales(x)

  if (!"alignment_position" %in% names(x)) {
    cli::cli_abort(
      "A {.cls tales_msa} object requires an {.field alignment_position} column.",
      class = c("tantale_error_tales_missing_column", "tantale_error")
    )
  }
  .tales_check_type(x, "alignment_position", is.integer, "an integer vector")

  if (nrow(x) == 0L) return(invisible(x))

  ## 11. positive, never NA
  if (anyNA(x$alignment_position) || any(x$alignment_position < 1L)) {
    cli::cli_abort(
      "{.field alignment_position} must be a positive integer without {.val NA}.",
      class = c("tantale_error_msa_position", "tantale_error")
    )
  }

  ## 12. unique within an array
  if (any(duplicated(data.frame(a = x$array_id, p = x$alignment_position)))) {
    cli::cli_abort(
      "{.field alignment_position} must be unique within an array.",
      class = c("tantale_error_msa_duplicate", "tantale_error")
    )
  }

  ## 13. an alignment may insert gaps but never reorder parts
  .tales_msa_check_order(x)

  ## 14. within the declared width
  w <- tales_width(x)
  if (!is.null(w) && any(x$alignment_position > w)) {
    cli::cli_abort(
      c("{.field alignment_position} must not exceed the alignment width ({w}).",
        "x" = "Found up to {max(x$alignment_position)}."),
      class = c("tantale_error_msa_width", "tantale_error")
    )
  }

  invisible(x)
}

#' Check that the alignment preserves part order
#'
#' Within an array, ranking by \code{alignment_position} must agree with ranking
#' by \code{position_in_array}: an alignment inserts gaps, it never reorders.
#' @noRd
.tales_msa_check_order <- function(x) {
  o <- x[order(x$array_id, x$position_in_array), c("array_id", "alignment_position")]
  bad <- unique(o$array_id[
    stats::ave(o$alignment_position, o$array_id, FUN = function(z) c(0L, diff(z))) < 0L
  ])
  if (length(bad) > 0L) {
    cli::cli_abort(
      c("An alignment may insert gaps but must not reorder parts.",
        "x" = "{.field alignment_position} disagrees with {.field position_in_array} in {.val {utils::head(bad, 5)}}."),
      class = c("tantale_error_msa_order", "tantale_error")
    )
  }
  invisible(NULL)
}

#' @noRd
.tales_msa_contract_holds <- function(x) {
  "alignment_position" %in% names(x) && is.integer(x$alignment_position)
}


#### Views ####

#' Render a TALE alignment as a matrix
#'
#' Materialises the rectangular form: one row per array, one column per
#' alignment position, cells holding the requested layer. This is the only place
#' the gapped matrix is built — the long object is the canonical storage.
#'
#' @param x A \code{tales_msa} object.
#' @param value Name of the column to fill cells with. Defaults to the first
#'   available of \code{rvd}, \code{dom_code}.
#' @param gap Value to use for gaps. Defaults to \code{NA}.
#' @param ... Ignored.
#' @return A character matrix with arrays as rows and alignment positions as
#'   columns.
#' @method as.matrix tales_msa
#' @export
#' @family TALE alignment
as.matrix.tales_msa <- function(x, value = NULL, gap = NA, ...) {
  if (is.null(value)) {
    value <- intersect(TALES_RESIDUE_COLS, names(x))[1]
  }
  if (!value %in% names(x)) {
    cli::cli_abort(
      c("This alignment has no {.field {value}} layer.",
        "i" = "Available: {.field {setdiff(names(x), c('array_id', 'alignment_position'))}}"),
      class = c("tantale_error_msa_layer", "tantale_error")
    )
  }
  arrays <- unique(x$array_id)
  width <- tales_width(x) %||% max(x$alignment_position, na.rm = TRUE)

  m <- matrix(gap, nrow = length(arrays), ncol = width,
              dimnames = list(arrays, as.character(seq_len(width))))
  idx <- cbind(match(x$array_id, arrays), x$alignment_position)
  m[idx] <- as.character(x[[value]])
  m
}

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#### Alignment ####

#' Align the repeat arrays of a tales object
#'
#' Aligns TALE arrays on one of their residue layers with MAFFT's text mode,
#' returning the alignment as a \code{\link{tales_msa}} — the input object with
#' an \code{alignment_position} column added, so every other layer
#' (\code{rvd}, \code{dom_code}, \code{aa_seq}, ...) remains available.
#'
#' The mapping back from MAFFT's output is **positional**: the k-th non-gap cell
#' of an aligned row is the k-th part fed in. That is well defined only because
#' this function builds MAFFT's input from \code{x} itself, which is why
#' \code{\link{tales_assert_complete}} is enforced first.
#'
#' @param x A \code{\link{tales}} object holding complete arrays.
#' @param residue_col Which layer to align on: \code{"rvd"} (default) or
#'   \code{"dom_code"}. Given explicitly rather than guessed from the values.
#' @param repeat_sims Scoring matrix for the residues being aligned.
#'   \code{NULL} (default) means none. Pass \code{"rvd"} to opt in to the
#'   built-in RVD similarity matrix when aligning RVDs, or a
#'   \code{\link{domain_distances}} object when aligning repeat codes. Optional similarity table passed to MAFFT as a scoring
#'   matrix, as accepted by \code{\link{tales_align}}.
#' @param ... Passed to \code{\link{tales_align}} (e.g. \code{mafft_opts}).
#' @return A \code{tales_msa} object.
#' @export
#' @family TALE alignment
tales_align <- function(x, residue_col = c("rvd", "dom_code"),
                        repeat_sims = NULL, ...) {
  residue_col <- match.arg(residue_col)
  if (!is_tales(x)) {
    cli::cli_abort("{.arg x} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  if (!residue_col %in% names(x)) {
    cli::cli_abort(
      c("{.arg x} has no {.field {residue_col}} column to align on.",
        "i" = "Available residue column{?s}: {.field {intersect(TALES_RESIDUE_COLS, names(x))}}"),
      class = c("tantale_error_msa_layer", "tantale_error")
    )
  }
  if (nrow(x) == 0L) {
    cli::cli_abort("{.arg x} has no parts to align.",
                   class = c("tantale_error_msa_empty", "tantale_error"))
  }
  tales_assert_complete(x, arg = "x")

  ## One space-separated string per array, in position_in_array order.
  ## A space is safe for both layers: neither RVDs nor repeat codes contain one.
  ord <- x[order(x$array_id, x$position_in_array), ]
  seqs <- split(as.character(ord[[residue_col]]), ord$array_id)
  seqs <- lapply(seqs, paste, collapse = " ")

  m <- .build_repeat_msa(input_seqs = seqs, sep = " ", repeat_sims = repeat_sims,
                        gap_symbol = NA, ...)

  if (!setequal(rownames(m), unique(x$array_id))) {
    cli::cli_abort(
      "MAFFT returned a different set of arrays than was submitted.",
      class = c("tantale_error_msa_backmap", "tantale_error")
    )
  }

  ## Positional back-mapping. MAFFT is run with --reorder by default, so rows
  ## must be matched by name, never by position.
  mapping <- lapply(rownames(m), function(a) {
    nz <- which(!is.na(m[a, ]))
    parts <- sort(x$position_in_array[x$array_id == a])
    if (length(nz) != length(parts)) {
      cli::cli_abort(
        c("Alignment of {.val {a}} returned {length(nz)} residue{?s} for {length(parts)} part{?s}.",
          "i" = "The back-mapping is positional and requires them to agree."),
        class = c("tantale_error_msa_backmap", "tantale_error")
      )
    }
    tibble::tibble(array_id = a, position_in_array = parts,
                   alignment_position = as.integer(nz))
  })
  mapping <- dplyr::bind_rows(mapping)

  out <- dplyr::left_join(x, mapping, by = c("array_id", "position_in_array"))
  tales_msa(out, alignment_width = ncol(m),
            dom_code_namespace = tales_namespace(x))
}


#### Running MAFFT to produce the alignment ####
#
# tales_align() is the only caller of .build_repeat_msa(), and
# .build_repeat_msa() the only caller of the two score-table helpers, so the
# three live here beside it rather than in msa.R, which is about drawing an
# alignment rather than building one.

#' Align TALE sequences with MAFFT text mode
#'
#' Implementation behind \code{\link{tales_align}} and the deprecated
#' \code{\link{tales_align}}. Internal so that package code can call it
#' without tripping the deprecation warning.
#' Normalise a pairwise table to what MAFFT's --textmatrix wants
#'
#' MAFFT scores matches, so it wants a *similarity*: higher means more alike.
#' Accepts either the canonical \code{\link{pairwise_distances}} vocabulary
#' (\code{id1}, \code{id2}, \code{dissim}) or the legacy one
#' (\code{RepU1}, \code{RepU2}, \code{Sim}), inverting the distance on the
#' way in.
#'
#' @param x A data frame of pairwise scores.
#' @return A three-column data frame named \code{id1}, \code{id2}, \code{sim}.
#' @keywords internal
.as_mafft_score_table <- function(x) {
  nms <- names(x)
  if (all(c("id1", "id2") %in% nms)) {
    if ("dissim" %in% nms) {
      out <- data.frame(id1 = x$id1, id2 = x$id2, sim = 100 - x$dissim,
                        stringsAsFactors = FALSE)
    } else if ("sim" %in% nms) {
      out <- data.frame(id1 = x$id1, id2 = x$id2, sim = x$sim,
                        stringsAsFactors = FALSE)
    } else {
      cli::cli_abort(
        c("A pairwise table must carry {.field dissim} or {.field sim}.",
          "i" = "Got: {.field {nms}}"),
        class = c("tantale_error_msa_sim_table", "tantale_error")
      )
    }
    return(out)
  }
  if (all(c("RepU1", "RepU2", "Sim") %in% nms)) {
    out <- as.data.frame(x[, c("RepU1", "RepU2", "Sim")])
    names(out) <- c("id1", "id2", "sim")
    return(out)
  }
  cli::cli_abort(
    c("Cannot read {.arg repeat_sims}.",
      "i" = "Expected {.field id1}/{.field id2}/{.field dissim}, or the legacy {.field RepU1}/{.field RepU2}/{.field Sim}.",
      "x" = "Got: {.field {nms}}"),
    class = c("tantale_error_msa_sim_table", "tantale_error")
  )
}

#' RVD similarity matrix for aligning RVD sequences
#'
#' @description
#' The repeat-level similarity matrix is keyed by \code{dom_code}, which is
#' meaningless for an RVD sequence, so RVD alignments have had no scoring matrix
#' at all -- MAFFT treated \code{NI} and \code{NN} as no more alike than
#' \code{NI} and \code{HD}. This supplies the missing one, from the Spearman
#' correlation of each RVD's A/C/G/T preference profile
#' (\code{rvdSimDf}, derived from TALVEZ's \code{mat1}).
#'
#' @details
#' \code{XX} -- a terminus detected but not identified -- has a uniform base
#' profile, so its correlation with everything is undefined. Those cells are
#' filled as **neutral** (0): we know nothing about it, so it should neither
#' attract nor repel.
#'
#' The one exception is \code{XX} against itself, which is set to the maximum.
#' That is forced, not chosen: MAFFT produces unusable output when the diagonal
#' is not high, and a low diagonal would anyway assert that an \code{XX} must
#' *not* align with an \code{XX}, which is a stronger claim than ignorance.
#'
#' Scale is irrelevant -- MAFFT normalises the matrix, so a linear rescale or
#' offset leaves the alignment unchanged. Only relative structure
#' matters, so the correlations are used as they are.
#'
#' @param residues Character vector of the RVDs present in the alignment.
#' @return A data frame of \code{id1}, \code{id2}, \code{sim} covering
#'   every ordered pair of \code{residues}.
#' @keywords internal
.rvd_score_table <- function(residues) {
  residues <- unique(as.character(residues))
  grid <- expand.grid(id1 = residues, id2 = residues,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  key <- paste(rvdSimDf$rvd1, rvdSimDf$rvd2)
  grid$sim <- rvdSimDf$Cor[match(paste(grid$id1, grid$id2), key)]
  # unknown pairings, and anything involving XX, are neutral
  grid$sim[is.na(grid$sim)] <- 0
  # ... except a symbol against itself, which must stay high for MAFFT
  grid$sim[grid$id1 == grid$id2] <- 1
  grid
}

#' @noRd
.build_repeat_msa <- function(input_seqs, sep = " ", repeat_sims = NULL,
                           mafft_opts = "--localpair --maxiterate 1000 --reorder --op 0 --ep 5 --thread 1",
                           mafft_path = system.file("tools", "mafft-linux64",package = "tantale", mustWork = TRUE),
                           gap_symbol = NA) {
  # A bunch of tempfiles
  simMatHexFile <- tempfile(pattern = "simMatHexFile")
  simMatAsciiFile <- tempfile(pattern = "simMatAsciiFile")
  hexFile <- tempfile(pattern = "hexFile")
  asciFile <- tempfile(pattern = "asciFile")
  mafftAsciiOutFile <- tempfile(pattern = "mafftAsciiOutFile")
  mafftHexOutFile <- tempfile(pattern = "mafftHexOutFile")

  # This is the data frame that will enable conversion of RVDs to Hexadecimal codes
  asciitable = data.frame(hex = as.raw(1:255),
                          printable =rawToChar(as.raw(1:255),multiple=TRUE),
                          stringsAsFactors = FALSE)
  mafftExcludedHex <- as.raw(c(0x0, 0x3E, 0x3D, 0x3C, 0x2D, 0x20, 0x0d, 0x0a))
  asciitableForMafft <- asciitable[! asciitable$hex %in% mafftExcludedHex, ]
  #Encoding(asciitableForMafft$printable) <- "bytes"
  
  # Load repeat/RVD sequences
  seqsAsVectors <- suppressWarnings(.split_list(input_seqs, sep = sep))
  residues <- unique(unlist(seqsAsVectors))
  
  # Deals with cases where the nomber of sequences is < 2
  if (length(seqsAsVectors) == 0L) {
    cli::cli_warn("There are no sequences to align. Returning an empty matrix.")
    return(matrix())
  }
  if (length(seqsAsVectors) == 1L) {
    cli::cli_inform("Only one sequence to align. Returning it as a one-row matrix.")
    msaOfResiduesAsMatrix <- as.matrix(as.data.frame(seqsAsVectors))
    msaOfResiduesAsMatrix <- matrix(msaOfResiduesAsMatrix, nrow = 1)
    rownames(msaOfResiduesAsMatrix) <- colnames(as.data.frame(seqsAsVectors))
    colnames(msaOfResiduesAsMatrix) <- 1:length(msaOfResiduesAsMatrix)
    return(msaOfResiduesAsMatrix)
  }
  
  # Determine the type of 'elements' (rvd or repeat) contained in the sequences
  frequentRvds <- c("NN", "NG", "HD", "NI", "N*", "NS")
  if(! any(residues %in% frequentRvds)) {
    cli::cli_inform(paste0("Will be assuming sequences contain repeat unit codes because ",
                     "none of the RVDs obtained from input sequences matches ",
                     "a list of 'frequent RVDs': {paste(frequentRvds, collapse = ' ')}"))
    repeatType <- "repeatUnit"
  } else {
    cli::cli_inform("Input sequences are detected as RVD sequences.")
    repeatType <- "rvds"
  }
  if( length(residues) > nrow(asciitableForMafft) ) {
    cli::cli_warn("The number of distinct residues (RVDs or repeat units) must be 248 or fewer.")
    cli::cli_abort(
      c("These sequences contain {length(residues)} distinct residues.",
        "i" = "MAFFT's text mode encodes each residue as one byte, which caps the alphabet at 248."),
      class = c("tantale_error_msa_alphabet", "tantale_error"))
  }
  

  
  # Coding residues in hexadecimal representations and concatenating them for mafft --text
  seqsOfHex <- sapply(seqsAsVectors, function(x) {
    idxs <- match(x, residues)
    paste(asciitableForMafft$hex[idxs], collapse = " ")
    }
  )
  seqsOfHex <- Biostrings::BStringSet(seqsOfHex)
  names(seqsOfHex) <- names(seqsAsVectors)
  # Write to a temp file
  Biostrings::writeXStringSet(seqsOfHex, filepath = hexFile)


  # If provided recode also the distance matrix
  # "rvd" opts in to the built-in RVD matrix; NULL and FALSE mean no matrix.
  # Opt-in rather than on-by-default: the matrix demonstrably makes alignments
  # more compact (fewer gaps, narrower) but there is no evidence it makes them
  # biologically better, and defaulting it on would silently change every
  # existing RVD alignment.
  if (is.null(repeat_sims) || isFALSE(repeat_sims)) {
    maffMatOpt <- ""
  } else if (identical(repeat_sims, "rvd")) {
    if (repeatType != "rvds") {
      cli::cli_abort(
        c('{.code repeat_sims = "rvd"} only applies to an RVD alignment.',
          "i" = "These sequences look like repeat unit codes."),
        class = c("tantale_error_msa_sim_table", "tantale_error")
      )
    }
    cli::cli_inform("Scoring the RVD alignment with the built-in RVD similarity matrix.")
    repeatSims <- .rvd_score_table(residues)
    repeatSims$id1 <- asciitableForMafft$hex[match(repeatSims$id1, residues)]
    repeatSims$id2 <- asciitableForMafft$hex[match(repeatSims$id2, residues)]
    colnames(repeatSims) <- NULL
    write.table(repeatSims, file = simMatHexFile, row.names = FALSE, fileEncoding = "ASCII")
    maffMatOpt <- glue::glue("--textmatrix {shQuote(simMatHexFile)}")
  } else if (repeatType == "rvds") {
    cli::cli_abort(
      c("A repeat similarity table cannot score an RVD alignment.",
        "i" = "It is keyed by {.field dom_code}, which has no meaning for RVDs.",
        "i" = 'Use {.code repeat_sims = "rvd"} for the built-in RVD matrix, or leave it {.code NULL} for none.'),
      class = c("tantale_error_msa_sim_table", "tantale_error")
    )
  } else {
    cli::cli_inform("The provided similarity matrix file will be used to compute msa.")
    if (length(repeat_sims) > 1 &&
        (is.data.frame(repeat_sims) | tibble::is_tibble(repeat_sims))
    ) {
      repeatSims <- .as_mafft_score_table(repeat_sims)
    } else if (length(repeat_sims) == 1 && is.character(repeat_sims)) {
      repeatSims <- .format_repeat_dist_mat(repeat_sims)
    } else {
      cli::cli_abort(
        c("Cannot read {.arg repeat_sims}.",
          "i" = "Expected a table of pairwise scores, or the path to a {.file *_Repeatmatrix.mat} file written by Distal.",
          "i" = 'Use {.code repeat_sims = "rvd"} for the built-in RVD matrix, or leave it {.code NULL} for none.',
          "x" = "Got {.obj_type_friendly {repeat_sims}}."),
        class = c("tantale_error_msa_sim_table", "tantale_error"))
    }

    stopifnot(all(residues %in% unique(repeatSims$id1)))
    repeatSims <- subset(repeatSims, id1 %in% residues & id2 %in% residues)
    stopifnot(all.equal(nrow(repeatSims), length(residues)^2))
    repeatSims$id1 <- asciitableForMafft$hex[match(repeatSims$id1, residues)]
    repeatSims$id2 <- asciitableForMafft$hex[match(repeatSims$id2, residues)]
    colnames(repeatSims) <-  NULL
    write.table(repeatSims, file = simMatHexFile, row.names = FALSE, fileEncoding = "ASCII")
    maffMatOpt <- glue::glue("--textmatrix {shQuote(simMatHexFile)}")
  }

  # Running mafft msa
  cli::cli_inform("Now running MAFFT (Copyright 2002-2007 Kazutaka Katoh) on TALE array sequences.")
  hex2text <- shQuote(file.path(mafft_path, "mafftdir", "libexec", "hex2maffttext"))
  text2hex <- shQuote(file.path(mafft_path, "mafftdir", "libexec", "maffttext2hex"))
  mafftBin <- shQuote(file.path(mafft_path, "mafft.bat"))
  asciiConverstionCmd <- glue::glue("{hex2text} {shQuote(hexFile)} > {shQuote(asciFile)}")
  mafftCmd <-  glue::glue("{mafftBin} {maffMatOpt} --text {mafft_opts} {shQuote(asciFile)} > {shQuote(mafftAsciiOutFile)}")
  MsaConversionToHexCmd <- glue::glue("{text2hex} {shQuote(mafftAsciiOutFile)} > {shQuote(mafftHexOutFile)}")
  res <- system(command = paste(asciiConverstionCmd, mafftCmd, MsaConversionToHexCmd, sep = "; "),
         ignore.stdout = FALSE, ignore.stderr = FALSE, intern = FALSE)

  # Getting msa output and converting back to alignment of residues
  msaOfHex <- Biostrings::readBStringSet(mafftHexOutFile)
  if (length(msaOfHex) == 0L) {
    cli::cli_warn("MAFFT did not complete successfully.")
    cli::cli_abort("MAFFT exit status: {res}", class = c("tantale_error"))
  } else {
  }
  #cat(as.character(msaOfHex), sep = "\n")
  msaAsHexVectors <- stringr::str_split(as.character(msaOfHex), pattern = " ")
  msaAsHexVectors <- lapply(msaAsHexVectors, function(x) x[-length(x)]) # Remove last "" element
  #cat(knitr::kable(t(matrix(msaAsHexVectors))), sep = " ")
  msaOfResiduesAsMatrix <- t(
    sapply(msaAsHexVectors, function(x) {
      idxs <- match(x, asciitableForMafft$hex, nomatch = NA)
      #cat("idx in asciiTable: ", idxs, "\n")
      seqOfResidues <- residues[idxs]
      seqOfResidues[is.na(seqOfResidues)] <- gap_symbol
      #cat("Seq of residues: ", seqOfResidues, "\n")
      seqOfResidues
    }
    )
  )
  rownames(msaOfResiduesAsMatrix) <- names(msaOfHex)
  colnames(msaOfResiduesAsMatrix) <- 1:ncol(msaOfResiduesAsMatrix)
  return(msaOfResiduesAsMatrix)
}
