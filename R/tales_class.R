##### The 'tales' class #####
# A long table of TALE array parts: one row per (array, slot).
# See dev/class-design.md §2 for the identity, column contract and invariants
# this file implements. Conditions are signalled with cli and carry classes so
# they can be asserted on (restructuring-notes.md §9.5).


#### Column schema ####

# Required: the primary key.
TALES_KEY_COLS <- c("array_id", "position_in_array")

# At least one of these must be present (the residue vocabulary).
TALES_RESIDUE_COLS <- c("rvd", "dom_code")

# Validated only if present.
TALES_OPTIONAL_COLS <- c(
  "domain_type", "position_in_crd", "aa_seq", "dna_seq",
  "seqnames", "source_directory"
)

TALES_DOMAIN_TYPES <- c("N-terminus", "repeat", "C-terminus")

# Legacy camelCase -> target snake_case. Kept so the class can be used against
# output of the not-yet-renamed pipeline; see class-design.md §1.1.
# `seqnames` is deliberately absent: it keeps its Bioconductor spelling.
TALES_LEGACY_NAMES <- c(
  arrayID         = "array_id",
  positionInArray = "position_in_array",
  positionInCrd   = "position_in_crd",
  domainType      = "domain_type",
  aaSeq           = "aa_seq",
  dnaSeq          = "dna_seq",
  domCode         = "dom_code",
  sourceDirectory = "source_directory"
)


#' Codes marking a TALE array terminus
#'
#' The values a \code{rvd} column takes on non-repeat parts. \code{"NTERM"} and
#' \code{"CTERM"} mark identified termini; \code{"XXXXX"} marks a terminus whose
#' CDS was detected but for which no HMMer hit was found, so its identity is
#' unknown (see \code{\link[tantale:tell_tales]{tell_tales}}).
#'
#' These share the \code{rvd} column with real RVDs, so code that distinguishes
#' repeats from termini by value should use this function rather than spelling
#' the codes out.
#'
#' @return A character vector of anchor codes.
#' @export
#' @examples
#' tales_anchor_codes()
#' @family tales objects
tales_anchor_codes <- function() {
  c("NTERM", "CTERM", "XXXXX")
}


#### Constructors ####

#' Low-level constructor for a tales object
#'
#' Attaches the class without validating. Use \code{\link{tales}} unless you
#' have already established the invariants.
#'
#' @param x A tibble.
#' @param dom_code_namespace Optional scalar string identifying the
#'   \code{dom_code} namespace this object belongs to (see
#'   \code{dev/class-design.md} §3.5). Carried, never recomputed.
#' @return A \code{tales} object.
#' @keywords internal
#' @family tales objects
new_tales <- function(x, dom_code_namespace = NULL) {
  stopifnot(is.data.frame(x))
  x <- tibble::as_tibble(x)
  if (!is.null(dom_code_namespace)) {
    attr(x, "dom_code_namespace") <- dom_code_namespace
  }
  class(x) <- unique(c("tales", class(x)))
  x
}

#' Is this a tales object?
#' @param x An object.
#' @return A logical scalar.
#' @export
#' @family tales objects
is_tales <- function(x) inherits(x, "tales")

#' The dom_code namespace of a tales object
#'
#' Identifies the run whose \code{dom_code} values this object carries. Objects
#' from different runs must not be joined: \code{dom_code} is minted per run, so
#' the codes collide across runs and a join would silently succeed against the
#' wrong repeats.
#'
#' @param x A \code{tales} object.
#' @return A scalar string, or \code{NULL} if the object is not stamped.
#' @export
#' @family tales objects
#' @examples
#' rvd_fasta <- system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",
#'                          package = "tantale")
#' x <- as_tales(rvd_fasta, sep = "-")
#' tales_namespace(x) # NULL -- dom_code has not been minted yet
tales_namespace <- function(x) {
  attr(x, "dom_code_namespace", exact = TRUE)
}

#' Compute the dom_code namespace identifier for a set of parts
#'
#' A content hash of the sorted unique amino acid sequences of a run.
#'
#' Hashing exactly that set is not arbitrary: \code{dom_code} is
#' \code{cur_group_id()} over \code{aa_seq}, so the code assignment is
#' a deterministic function of the sorted unique sequences and nothing else. Two
#' runs over identical parts therefore hash alike and are correctly treated as
#' compatible, while any change to the part set yields a different namespace.
#'
#' \code{xxhash64} is used rather than a cryptographic digest: this is an
#' identity tag guarding against accidental cross-run joins, not a security
#' boundary, so a short fast hash is the right trade.
#'
#' Called once, where a run is defined; never recomputed on a subset — that is
#' the whole point of the tag.
#' @noRd
.tales_dom_code_namespace <- function(aa_seq) {
  digest::digest(sort(unique(as.character(aa_seq))), algo = "xxhash64")
}


#' Create a tales object
#'
#' Builds a \code{tales} object from a data frame of TALE array parts: one row
#' per (array, slot). Column names still using the legacy camelCase spelling
#' (\code{arrayID}, \code{positionInArray}, ...) are renamed on the way in.
#'
#' @param x A data frame with at least \code{array_id} and
#'   \code{position_in_array} columns, plus at least one of \code{rvd} or
#'   \code{dom_code}. Other recognised columns (\code{domain_type},
#'   \code{position_in_crd}, \code{aa_seq}, \code{dna_seq}, \code{seqnames},
#'   \code{source_directory}) are validated if present. Any further column is
#'   preserved untouched.
#' @param dom_code_namespace Optional scalar string, see
#'   \code{\link{tales_namespace}}.
#' @param sanitize If \code{TRUE}, arrays carrying biological anomalies --
#'   missing sequences, impossible terminus arrangements, coordinate
#'   disagreements -- are removed, with a warning naming them and why. If
#'   \code{FALSE} (default) they are kept and merely warned about, so odd
#'   predictions can still be loaded and inspected. Structural corruption is
#'   an error either way. See \code{\link{tales_anomalies}}.
#' @return A validated \code{tales} object.
#' @export
#' @family tales objects
#' @examples
#' parts <- data.frame(
#'   array_id = c("A1", "A1", "A1", "A2", "A2", "A2"),
#'   position_in_array = c(1L, 2L, 3L, 1L, 2L, 3L),
#'   domain_type = c("N-terminus", "repeat", "C-terminus",
#'                   "N-terminus", "repeat", "C-terminus"),
#'   rvd = c("NTERM", "HD", "CTERM", "NTERM", "NI", "CTERM")
#' )
#' tales(parts)
tales <- function(x, dom_code_namespace = NULL, sanitize = FALSE) {
  if (!is.data.frame(x)) {
    cli::cli_abort(
      "{.arg x} must be a data frame, not {.obj_type_friendly {x}}.",
      class = c("tantale_error_tales_type", "tantale_error")
    )
  }
  x <- .tales_rename_legacy(x)
  if ("position_in_array" %in% names(x) && is.numeric(x$position_in_array)) {
    x$position_in_array <- as.integer(x$position_in_array)
  }
  out <- validate_tales(new_tales(x, dom_code_namespace = dom_code_namespace))
  # Structural problems have already aborted above. Biological anomalies are
  # reported here, and removed if asked for.
  .tales_report_anomalies(out, sanitize = sanitize)
}

#' Rename legacy camelCase columns to the target schema
#' @param x A data frame.
#' @return The same data frame with recognised legacy names replaced.
#' @noRd
.tales_rename_legacy <- function(x) {
  hit <- intersect(names(x), names(TALES_LEGACY_NAMES))
  if (length(hit) == 0L) return(x)
  clash <- intersect(unname(TALES_LEGACY_NAMES[hit]), names(x))
  if (length(clash) > 0L) {
    cli::cli_abort(
      c("Cannot rename legacy columns: the target name{?s} {.field {clash}} {?is/are} already present.",
        "i" = "Drop or rename {cli::qty(clash)}the duplicate{?s} before calling {.fn tales}."),
      class = c("tantale_error_tales_name_clash", "tantale_error")
    )
  }
  names(x)[match(hit, names(x))] <- unname(TALES_LEGACY_NAMES[hit])
  x
}


#### Construction from external sources ####

#' Coerce sequences of TALE parts to a tales object
#'
#' Turns \code{sep}-separated TALE sequences — RVD strings, or repeat-code
#' strings — into a \code{\link{tales}} object, taking sequence order as
#' \code{position_in_array} and sequence names as \code{array_id}.
#'
#' The result is deliberately column-poor: a bare sequence file carries no
#' domain types, amino acid sequences or source contigs, so only
#' \code{array_id}, \code{position_in_array} and the chosen residue column are
#' produced. That is a valid \code{tales} — see \code{dev/class-design.md}
#' §2.3 for why \code{seqnames} and the rest are optional.
#'
#' @param x A path to a fasta file, a \code{BStringSet}/\code{AAStringSet}, a
#'   list of strings, or a data frame (which is passed to \code{\link{tales}}).
#' @param sep Separator between elements of a sequence. Use \code{"-"} for RVD
#'   sequences and \code{" "} for repeat-code strings.
#' @param residue_col Which residue column the parsed elements become:
#'   \code{"rvd"} (default) or \code{"dom_code"}. Given explicitly rather than
#'   guessed from the values.
#' @param ... Passed to methods.
#' @return A validated \code{tales} object.
#' @export
#' @family tales objects
#' @examples
#' rvd_fasta <- system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",
#'                          package = "tantale")
#' as_tales(rvd_fasta, sep = "-")
as_tales <- function(x, ...) {
  UseMethod("as_tales")
}

#' @rdname as_tales
#' @export
as_tales.data.frame <- function(x, ...) {
  tales(x, ...)
}

#' @rdname as_tales
#' @export
as_tales.default <- function(x, sep = "-", residue_col = c("rvd", "dom_code"), ...) {
  residue_col <- match.arg(residue_col)
  seqs <- .split_list(x, sep = sep)

  if (is.null(names(seqs)) || anyNA(names(seqs)) || !all(nzchar(names(seqs)))) {
    cli::cli_abort(
      c("Every sequence must be named; the names become {.field array_id}.",
        "i" = "A fasta file supplies these from its headers."),
      class = c("tantale_error_tales_unnamed", "tantale_error")
    )
  }

  out <- tibble::tibble(
    array_id = rep(names(seqs), lengths(seqs)),
    position_in_array = unlist(lapply(lengths(seqs), seq_len), use.names = FALSE),
    residue = unlist(seqs, use.names = FALSE)
  )
  names(out)[names(out) == "residue"] <- residue_col
  tales(out)
}


#### Combining ####

#' Combine tales objects
#'
#' Row-binds two or more \code{\link{tales}} objects into one, reconciling
#' the invariants that a plain \code{\link[dplyr]{bind_rows}} would not check:
#' \code{array_id} uniqueness across inputs, the \code{dom_code} namespace,
#' and the \code{group} column.
#'
#' @details
#' Not named \code{c.tales()} for two reasons: base \code{c()} dispatch is
#' leaky (mixing a \code{tales} with an unrelated object can silently drop
#' attributes rather than error, the wrong failure mode for a class whose
#' point is invariants that must not go silent), and
#' \code{on_namespace_mismatch} has no room in \code{c()}'s signature.
#'
#' \code{tales_msa} inputs are refused rather than silently demoted --
#' \code{alignment_width}/\code{alignment_position} are a coordinate system
#' specific to one alignment run, and binding two runs' matrices would
#' produce an object that misdescribes what a given column means in each
#' row. Demote explicitly with \code{\link{as_tales}} first if that is really
#' what you want.
#'
#' \code{group} is dropped from the result whenever present on any input, not
#' reconciled: it is a clustering result over one specific distance matrix
#' and one specific set of arrays, so two "group 1"s from separate
#' \code{\link{tales_group}} calls are not comparable, not merely at risk of
#' colliding. Recompute it with \code{\link{tales_group}} after
#' \code{\link{tales_compare}} on the bound result.
#'
#' \code{tale_distances}/\code{domain_distances} are untouched: this function
#' only binds \code{tales} data. A companion distance table from either input
#' does not describe the bound object -- re-run \code{\link{tales_compare}}
#' if you need one.
#'
#' @param ... Two or more \code{\link{tales}} objects (not \code{tales_msa}).
#' @param on_namespace_mismatch What to do when the inputs' \code{dom_code}
#'   namespaces (see \code{\link{tales_namespace}}) disagree. \code{"recode"}
#'   (default) drops the stale \code{dom_code} column and calls
#'   \code{\link{tales_assign_domain_codes}} on the bound result, over the
#'   union of \code{aa_seq}; \code{"error"} aborts instead, for callers who
#'   want that stricter behaviour.
#' @param sanitize Passed to \code{\link{tales}} for the final construction.
#' @return A validated \code{\link{tales}} object.
#' @export
#' @family tales objects
#' @examples
#' a <- tales(data.frame(
#'   array_id = "A1", position_in_array = 1:2,
#'   rvd = c("NTERM", "HD")
#' ))
#' b <- tales(data.frame(
#'   array_id = "A2", position_in_array = 1:2,
#'   rvd = c("NTERM", "NI")
#' ))
#' tales_bind(a, b)
tales_bind <- function(..., on_namespace_mismatch = c("recode", "error"), sanitize = FALSE) {
  on_namespace_mismatch <- match.arg(on_namespace_mismatch)
  inputs <- rlang::list2(...)

  if (length(inputs) == 0L) {
    cli::cli_abort(
      "{.fn tales_bind} needs at least one {.cls tales} object.",
      class = c("tantale_error_bind_type", "tantale_error")
    )
  }

  for (x in inputs) {
    if (!is_tales(x)) {
      cli::cli_abort(
        "Every argument to {.fn tales_bind} must be a {.cls tales} object.",
        class = c("tantale_error_bind_type", "tantale_error")
      )
    }
    if (is_tales_msa(x)) {
      cli::cli_abort(
        c("{.fn tales_bind} does not accept {.cls tales_msa} input.",
          "i" = "Demote with {.fn as_tales} first if binding across alignment runs is really what you want."),
        class = c("tantale_error_bind_tales_msa", "tantale_error")
      )
    }
  }

  .tales_bind_check_array_ids(inputs)

  has_group <- any(vapply(inputs, function(x) "group" %in% names(x), logical(1)))
  if (has_group) {
    cli::cli_inform(
      c("Dropping {.field group} from the bound result.",
        "i" = "It is a clustering result over each input's own distances, not a comparable label. Recompute with {.fn tales_group} after {.fn tales_compare}."),
      class = c("tantale_message_bind_group_dropped", "tantale_message")
    )
    inputs <- lapply(inputs, function(x) {
      x$group <- NULL
      x
    })
  }

  bound <- dplyr::bind_rows(lapply(inputs, tibble::as_tibble))

  namespace <- .tales_bind_reconcile_namespace(inputs, bound, on_namespace_mismatch)
  bound <- namespace$x

  tales(bound, dom_code_namespace = namespace$namespace, sanitize = sanitize)
}

#' @keywords internal
.tales_bind_check_array_ids <- function(inputs) {
  ids <- lapply(inputs, function(x) unique(x$array_id))
  all_ids <- unlist(ids, use.names = FALSE)
  dup <- unique(all_ids[duplicated(all_ids)])
  if (length(dup) > 0L) {
    cli::cli_abort(
      c("{.field array_id} must be disjoint across the inputs to {.fn tales_bind}.",
        "x" = "Present in more than one input: {.val {utils::head(dup, 8)}}"),
      class = c("tantale_error_bind_array_id_collision", "tantale_error")
    )
  }
  invisible(NULL)
}

#' @keywords internal
.tales_bind_reconcile_namespace <- function(inputs, bound, on_namespace_mismatch) {
  namespaces <- lapply(inputs, tales_namespace)
  agree <- length(unique(namespaces)) <= 1L
  if (agree) {
    return(list(x = bound, namespace = namespaces[[1L]]))
  }
  if (on_namespace_mismatch == "error") {
    cli::cli_abort(
      c("The inputs to {.fn tales_bind} carry different {.field dom_code} namespaces.",
        "i" = "Pass {.code on_namespace_mismatch = \"recode\"} to reassign {.field dom_code} over their union instead."),
      class = c("tantale_error_bind_namespace_mismatch", "tantale_error")
    )
  }
  cli::cli_inform(
    c("Reassigning {.field dom_code} over the union of inputs.",
      "i" = "The inputs' {.field dom_code} namespaces disagreed, so any distance table from either is now stale."),
    class = c("tantale_message_bind_namespace_recoded", "tantale_message")
  )
  bound$dom_code <- NULL
  # new_tales(), not tales(): this is an intermediate shape on its way to the
  # single official construction in tales_bind() itself (step 6 of the
  # ledger's plan, §5.2) -- validating here too would just double the
  # anomaly warnings the final tales() call already reports.
  recoded <- tales_assign_domain_codes(new_tales(bound))
  list(x = tibble::as_tibble(recoded), namespace = tales_namespace(recoded))
}


#### Validator ####

#' Validate a tales object
#'
#' Checks the column contract and every invariant that is closed under row
#' subsetting. Properties that hold only of a *complete* object — that an array
#' carries all its parts, numbered contiguously from 1 — are deliberately not
#' checked here; they are preconditions of the functions that need them, such
#' as alignment.
#'
#' @param x A \code{tales} object.
#' @return \code{x}, invisibly, if valid; otherwise an error.
#' @export
#' @family tales objects
#' @examples
#' parts <- data.frame(
#'   array_id = c("A1", "A1"), position_in_array = c(1L, 2L),
#'   rvd = c("NTERM", "HD")
#' )
#' validate_tales(tales(parts))
validate_tales <- function(x) {
  cols <- names(x)

  ## Column contract -------------------------------------------------------
  missing <- setdiff(TALES_KEY_COLS, cols)
  if (length(missing) > 0L) {
    cli::cli_abort(
      "A {.cls tales} object requires the column{?s} {.field {missing}}.",
      class = c("tantale_error_tales_missing_column", "tantale_error")
    )
  }
  if (!any(TALES_RESIDUE_COLS %in% cols)) {
    cli::cli_abort(
      c("A {.cls tales} object requires at least one residue column.",
        "i" = "Expected one of {.field {TALES_RESIDUE_COLS}}."),
      class = c("tantale_error_tales_missing_column", "tantale_error")
    )
  }
  .tales_check_type(x, "array_id", is.character, "a character vector")
  .tales_check_type(x, "position_in_array", is.integer, "an integer vector")
  for (nm in intersect(c(TALES_RESIDUE_COLS, "domain_type", "aa_seq",
                         "dna_seq", "seqnames", "source_directory"), cols)) {
    .tales_check_type(x, nm, is.character, "a character vector")
  }
  if ("position_in_crd" %in% cols) {
    .tales_check_type(x, "position_in_crd", is.integer, "an integer vector")
  }

  ## An empty tales is valid; the row invariants are vacuous ----------------
  if (nrow(x) == 0L) return(invisible(x))

  ## Hard invariants -------------------------------------------------------
  if (anyNA(x$array_id)) {
    cli::cli_abort("{.field array_id} must not contain {.val NA}.",
                   class = c("tantale_error_tales_na", "tantale_error"))
  }
  if (anyNA(x$position_in_array) || any(x$position_in_array < 1L)) {
    cli::cli_abort(
      "{.field position_in_array} must be a positive integer without {.val NA}.",
      class = c("tantale_error_tales_position", "tantale_error")
    )
  }
  .tales_check_key(x)

  ## Structural, conditional -------------------------------------------------
  # dom_code is a surrogate key into both similarity tables, so a broken
  # bijection makes those joins wrong rather than merely odd.
  if (all(c("aa_seq", "dom_code") %in% cols)) .tales_check_bijection(x)
  if ("position_in_crd" %in% cols) .tales_check_crd_unique(x)

  # Everything else -- missing sequences, impossible terminus arrangements,
  # coordinate disagreements, attributes varying within an array -- is a
  # biological anomaly rather than a structural one. See .tales_anomalies().

  invisible(x)
}

#' @keywords internal
.tales_check_type <- function(x, nm, predicate, expected) {
  if (!predicate(x[[nm]])) {
    cli::cli_abort(
      "{.field {nm}} must be {expected}, not {.obj_type_friendly {x[[nm]]}}.",
      class = c("tantale_error_tales_type", "tantale_error")
    )
  }
  invisible(NULL)
}


#' Report the biological anomalies in a tales object
#'
#' @description
#' Lists the arrays that are *odd* rather than *unreadable*: missing sequence
#' data, impossible domain-type arrangements, coordinate disagreements, an
#' amino acid sequence paired with more than one RVD, or an attribute that
#' varies within an array when it should not.
#'
#' Such arrays are accepted by \code{\link{tales}} -- real TALE predictions are
#' messy, and refusing to load them would force cleaning outside the package and
#' destroy the diagnostic signal. Construction warns about them; this function
#' tells you which and why; \code{tales(x, sanitize = TRUE)} removes them.
#'
#' @param x A \code{\link{tales}} object.
#' @return A tibble of \code{array_id}, \code{check} and \code{detail}, one
#'   row per anomaly. Zero rows if the object is clean.
#' @export
#' @family tales objects
#' @examples
#' # A2 carries two N-termini -- an impossible arrangement.
#' odd <- data.frame(
#'   array_id = c("A1", "A1", "A2", "A2", "A2"),
#'   position_in_array = c(1L, 2L, 1L, 2L, 3L),
#'   domain_type = c("N-terminus", "repeat",
#'                   "N-terminus", "N-terminus", "repeat"),
#'   rvd = c("NTERM", "HD", "NTERM", "NTERM", "NI")
#' )
#' x <- suppressWarnings(tales(odd))
#' tales_anomalies(x)
#' tales(odd, sanitize = TRUE) # drops A2 instead of merely warning
tales_anomalies <- function(x) {
  .tales_anomalies(x)
}


#' Biological anomalies in a tales object
#'
#' Collects the array-level anomalies that make an object *odd* rather than
#' *unreadable*: missing sequence data, impossible domain-type arrangements,
#' coordinate disagreements, an amino acid sequence paired with more than one
#' RVD, attributes that should be constant within an array but are not.
#'
#' These are deliberately **not** errors. Real TALE predictions are messy, and a
#' class that refuses to load them forces cleaning outside the package and
#' destroys exactly the diagnostic signal a user wants. They are reported as a
#' warning on construction and can be removed with \code{sanitize = TRUE}.
#'
#' Structural violations -- a duplicated key, an \code{NA} \code{array_id}, a
#' broken \code{aa_seq}/\code{dom_code} bijection -- are a different matter and
#' remain hard errors: the table cannot be interpreted at all, and downstream
#' code would silently compute wrong answers rather than merely odd ones.
#'
#' @param x A data frame with at least the tales key columns.
#' @return A tibble with one row per (array, anomaly): \code{array_id},
#'   \code{check} and \code{detail}. Zero rows if the object is clean.
#' @keywords internal
.tales_anomalies <- function(x) {
  cols <- names(x)
  out <- list()
  add <- function(ids, check, detail) {
    ids <- unique(ids[!is.na(ids)])
    if (length(ids)) out[[length(out) + 1L]] <<-
      tibble::tibble(array_id = ids, check = check, detail = detail)
  }
  if (nrow(x) == 0L) {
    return(tibble::tibble(array_id = character(), check = character(),
                          detail = character()))
  }

  ## missing sequence data ---------------------------------------------------
  for (nm in intersect(c(TALES_RESIDUE_COLS, "aa_seq", "dna_seq"), cols)) {
    bad <- is.na(x[[nm]]) | !nzchar(x[[nm]])
    add(x$array_id[bad], paste0("missing_", nm),
        paste0("part(s) with no ", nm))
  }

  ## aa_seq -> rvd must be a function -----------------------------------------
  # Not a full bijection: many distinct aa_seq legitimately share one rvd
  # (different proteins, same DNA-binding specificity -- the whole reason
  # dom_code and rvd are separate layers, dev/restructuring-notes.md's
  # tales-class article). What must not happen is the same aa_seq pairing
  # with two different rvd values: dom_code is minted from aa_seq alone
  # (cur_group_id()), so that would silently give tales_domain_codes() two
  # rows for one dom_code, breaking its "one row per distinct code" contract.
  if (all(c("aa_seq", "rvd") %in% cols)) {
    known <- !is.na(x$aa_seq) & nzchar(x$aa_seq) & !is.na(x$rvd) & nzchar(x$rvd)
    if (any(known)) {
      n_rvd <- tapply(x$rvd[known], x$aa_seq[known], function(z) length(unique(z)))
      bad_aa_seq <- names(n_rvd)[n_rvd > 1L]
      if (length(bad_aa_seq)) {
        add(x$array_id[known & x$aa_seq %in% bad_aa_seq], "aa_seq_rvd_inconsistent",
            "the same amino acid sequence is paired with more than one rvd")
      }
    }
  }

  ## domain_type arrangement -------------------------------------------------
  if ("domain_type" %in% cols) {
    unknown <- !x$domain_type %in% TALES_DOMAIN_TYPES
    add(x$array_id[unknown], "domain_type_unknown",
        "domain_type outside the expected vocabulary")
    for (type in c("N-terminus", "C-terminus")) {
      n <- tapply(x$domain_type == type, x$array_id, sum)
      add(names(n)[!is.na(n) & n > 1L], "terminus_duplicated",
          paste0("more than one ", type))
    }
    nterm <- x$domain_type == "N-terminus" & x$position_in_array != 1L
    add(x$array_id[nterm], "terminus_misplaced", "N-terminus not at position 1")
    cterm <- x$domain_type == "C-terminus"
    if (any(cterm)) {
      maxpos <- tapply(x$position_in_array, x$array_id, max)
      cpos <- tapply(x$position_in_array[cterm], x$array_id[cterm], max)
      shared <- intersect(names(cpos), names(maxpos))
      add(shared[cpos[shared] != maxpos[shared]], "terminus_misplaced",
          "C-terminus not at the end of its array")
    }
  }

  ## coordinate agreement ----------------------------------------------------
  if ("position_in_crd" %in% cols) {
    if ("domain_type" %in% cols) {
      wrong <- is.na(x$position_in_crd) != (x$domain_type != "repeat")
      add(x$array_id[wrong], "crd_placement",
          "position_in_crd is not NA on exactly the non-repeat parts")
    }
    keep <- !is.na(x$position_in_crd)
    if (any(keep)) {
      add(x$array_id[keep][x$position_in_crd[keep] < 1L], "crd_placement",
          "position_in_crd is not positive")
      o <- x[keep, c("array_id", "position_in_array", "position_in_crd")]
      o <- o[order(o$array_id, o$position_in_array), ]
      add(o$array_id[stats::ave(o$position_in_crd, o$array_id,
                                FUN = function(z) c(0L, diff(z))) < 0L],
          "crd_order", "position_in_crd does not increase with position_in_array")
    }
  }

  ## attributes that should be constant within an array ----------------------
  for (nm in intersect(c("seqnames", "group"), cols)) {
    n <- tapply(x[[nm]], x$array_id, function(z) length(unique(z)))
    add(names(n)[!is.na(n) & n > 1L], paste0(nm, "_inconsistent"),
        paste0(nm, " varies within the array"))
  }

  if (!length(out)) {
    return(tibble::tibble(array_id = character(), check = character(),
                          detail = character()))
  }
  unique(do.call(rbind, out))
}


#' Warn about, or drop, the arrays flagged by .tales_anomalies()
#' @keywords internal
.tales_report_anomalies <- function(x, sanitize = FALSE, arg = "x") {
  an <- .tales_anomalies(x)
  if (nrow(an) == 0L) return(x)
  ids <- unique(an$array_id)
  reasons <- unique(an$check)
  if (isTRUE(sanitize)) {
    cli::cli_warn(
      c("Dropped {length(ids)} array{?s} with biological anomalies.",
        "x" = "Array{?s}: {.val {utils::head(ids, 8)}}",
        "i" = "Reason{?s}: {.field {reasons}}"),
      class = c("tantale_warning_tales_sanitized", "tantale_warning")
    )
    return(x[!x$array_id %in% ids, , drop = FALSE])
  }
  cli::cli_warn(
    c("{length(ids)} array{?s} {?has/have} biological anomalies.",
      "x" = "Array{?s}: {.val {utils::head(ids, 8)}}",
      "i" = "Reason{?s}: {.field {reasons}}",
      "i" = "Inspect with {.fn tales_anomalies}, or drop with {.code sanitize = TRUE}."),
    class = c("tantale_warning_tales_anomalous", "tantale_warning")
  )
  x
}



#' @keywords internal
.tales_check_key <- function(x) {
  dup <- duplicated(x[TALES_KEY_COLS])
  if (any(dup)) {
    offenders <- unique(x$array_id[dup])
    cli::cli_abort(
      c("{.field array_id} and {.field position_in_array} must together be unique.",
        "x" = "{sum(dup)} duplicated row{?s} in {length(offenders)} array{?s}: {.val {utils::head(offenders, 5)}}"),
      class = c("tantale_error_tales_duplicate_key", "tantale_error")
    )
  }
  invisible(NULL)
}


#' Structural part of the CRD contract: the coordinate must be a key
#'
#' Placement (NA on exactly the non-repeat parts) and order agreement are
#' *biological* checks and live in \code{.tales_anomalies()}; only uniqueness
#' is structural, since a repeated coordinate makes the array unindexable.
#' @keywords internal
.tales_check_crd_unique <- function(x) {
  keep <- !is.na(x$position_in_crd)
  if (!any(keep)) return(invisible(NULL))
  if (any(duplicated(data.frame(a = x$array_id[keep], p = x$position_in_crd[keep])))) {
    cli::cli_abort("{.field position_in_crd} must be unique within an array.",
                   class = c("tantale_error_tales_crd", "tantale_error"))
  }
  invisible(NULL)
}



#### Preconditions ####
# Properties that hold only of a *complete* object. Not invariants: ordinary
# subsetting breaks them legitimately. Checked by the functions that need them.

#' Assert that a tales object holds complete arrays
#'
#' Checks that every array carries all of its parts: \code{position_in_array}
#' running contiguously from 1, and — when \code{domain_type} and
#' \code{position_in_crd} are present — \code{position_in_crd} equal to
#' \code{position_in_array} minus the number of non-repeat parts before it.
#'
#' This is a **precondition**, not an invariant: \code{filter(x, domain_type ==
#' "repeat")} legitimately produces a valid \code{tales} that is no longer
#' complete. Alignment requires completeness, because the mapping back from a
#' MAFFT alignment is positional — the k-th aligned residue is the k-th part
#' fed in.
#'
#' @param x A \code{tales} object.
#' @param arg Name of the argument being checked, for the error message.
#' @return \code{x}, invisibly.
#' @export
#' @family tales objects
#' @examples
#' rvd_fasta <- system.file("extdata", "TalA_RVDSeqs_AnnoTALE.fasta",
#'                          package = "tantale")
#' x <- as_tales(rvd_fasta, sep = "-")
#' tales_assert_complete(x) # every array is 1..n already -- no error
#'
#' \dontrun{
#' # Keeping only the repeats breaks completeness, which tales_align() needs.
#' repeats_only <- x[x$position_in_array > 1, ]
#' tales_assert_complete(repeats_only)
#' }
tales_assert_complete <- function(x, arg = "x") {
  if (!is_tales(x)) {
    cli::cli_abort("{.arg {arg}} must be a {.cls tales} object.",
                   class = c("tantale_error_tales_type", "tantale_error"))
  }
  if (nrow(x) == 0L) return(invisible(x))

  n <- tapply(x$position_in_array, x$array_id, function(p) identical(sort(p), seq_along(p)))
  if (any(!n)) {
    cli::cli_abort(
      c("{.arg {arg}} must hold complete arrays.",
        "x" = "{.field position_in_array} is not 1..n in {.val {utils::head(names(n)[!n], 5)}}.",
        "i" = "Subsetting parts (e.g. keeping only repeats) breaks completeness."),
      class = c("tantale_error_tales_incomplete", "tantale_error")
    )
  }

  if (all(c("domain_type", "position_in_crd") %in% names(x))) {
    o <- x[order(x$array_id, x$position_in_array), ]
    offset <- stats::ave(as.integer(o$domain_type != "repeat"), o$array_id, FUN = cumsum)
    keep <- o$domain_type == "repeat"
    if (!identical(o$position_in_crd[keep], as.integer((o$position_in_array - offset)[keep]))) {
      bad <- unique(o$array_id[keep][o$position_in_crd[keep] != (o$position_in_array - offset)[keep]])
      cli::cli_abort(
        c("{.field position_in_crd} is inconsistent with {.field position_in_array} in {.arg {arg}}.",
          "i" = "Expected {.code position_in_array - (non-repeat parts before it)}.",
          "x" = "Affected array{?s}: {.val {utils::head(bad, 5)}}"),
        class = c("tantale_error_tales_incomplete", "tantale_error")
      )
    }
  }
  invisible(x)
}

#' @keywords internal
.tales_check_bijection <- function(x) {
  per_seq <- tapply(x$dom_code, x$aa_seq, function(z) length(unique(z)))
  per_code <- tapply(x$aa_seq, x$dom_code, function(z) length(unique(z)))
  if (any(per_seq > 1L) || any(per_code > 1L)) {
    cli::cli_abort(
      c("{.field aa_seq} and {.field dom_code} must be in one-to-one correspondence.",
        "i" = "{.field dom_code} is a surrogate key for the part's amino acid sequence."),
      class = c("tantale_error_tales_dom_code", "tantale_error")
    )
  }
  invisible(NULL)
}



#### dplyr integration ####
# The column contract is re-checked on every verb; the row invariants are not,
# because all of them are closed under row subsetting. `mutate()` can still
# collide the key by overwriting it in place, so that one check lives in
# dplyr_col_modify(). See class-design.md §2.6.

#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.tales <- function(data, template) {
  if (!.tales_contract_holds(data)) {
    return(.tales_declass(data))
  }
  out <- NextMethod()
  .tales_regrade(out, template)
}

#' @exportS3Method dplyr::dplyr_col_modify
dplyr_col_modify.tales <- function(data, cols) {
  out <- NextMethod()
  if (is_tales(out) && nrow(out) > 0L) .tales_check_key(out)
  out
}

# dplyr's dplyr_col_select() only calls dplyr_reconstruct() for plain
# data.frame/data.table (verified, dplyr 1.2.1); for a tibble subclass it relies
# on the class's own `[`. So column-dropping is caught here, not in
# dplyr_reconstruct() -- this is what makes select(x, -array_id) degrade.
#' @export
`[.tales` <- function(x, ...) {
  template <- x
  out <- NextMethod()
  if (!is.data.frame(out)) return(out)
  if (!.tales_contract_holds(out)) return(.tales_declass(out))
  .tales_regrade(out, template)
}

#' Drop the tales classes, leaving a plain tibble
#' @param x A data frame.
#' @return The same data without the \code{tales}/\code{tales_msa} classes.
#' @noRd
.tales_declass <- function(x) {
  class(x) <- setdiff(class(x), c("tales_msa", "tales"))
  attr(x, "dom_code_namespace") <- NULL
  attr(x, "alignment_width") <- NULL
  tibble::as_tibble(x)
}

#' Restore class and attributes after an operation, demoting if needed
#'
#' Degradation is graded: an object that loses \code{alignment_position} but
#' keeps the \code{tales} contract becomes a plain \code{tales} rather than
#' dropping all the way to a tibble.
#' @noRd
.tales_regrade <- function(out, template) {
  attr(out, "dom_code_namespace") <- tales_namespace(template)
  if (inherits(out, "tales_msa")) {
    if (.tales_msa_contract_holds(out)) {
      attr(out, "alignment_width") <- attr(template, "alignment_width", exact = TRUE)
    } else {
      class(out) <- setdiff(class(out), "tales_msa")
      attr(out, "alignment_width") <- NULL
    }
  }
  out
}

#' @keywords internal
.tales_contract_holds <- function(x) {
  all(TALES_KEY_COLS %in% names(x)) &&
    any(TALES_RESIDUE_COLS %in% names(x)) &&
    is.character(x$array_id) &&
    is.integer(x$position_in_array)
}
