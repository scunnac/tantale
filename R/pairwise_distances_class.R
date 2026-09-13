##### The 'pairwise_distances' class #####
# The long form of a square pairwise similarity over one entity set: one row per
# ordered pair. Subclasses 'tale_distances' and 'domain_distances' add no structure, only
# entity semantics. See dev/class-design.md §3.
#
# Squareness, a diagonal of 100 and symmetry are *preconditions*, not
# invariants: the package legitimately filters these tables asymmetrically
# (conversion.R, inside a loop over alignment columns) and in two steps
# (msa.R), passing through non-square intermediates on purpose.


#### Column schema ####

PAIRWISE_DISTANCES_ID_COLS <- c("id1", "id2")
PAIRWISE_DISTANCES_VALUE_COL <- "dissim"

PAIRWISE_DISTANCES_OPTIONAL_COLS <- c(
  "arlem_score", "max_length"
)

# Legacy spellings -> canonical. Both entity vocabularies collapse onto the
# same id columns; renaming by name also fixes repeat.similarity listing its
# ids in the order RepU2, RepU1.
PAIRWISE_DISTANCES_LEGACY_NAMES <- c(
  TAL1           = "id1",
  TAL2           = "id2",
  RepU1          = "id1",
  RepU2          = "id2",
  Sim            = "sim",
  Dissim         = "dissim",
  arlemScore     = "arlem_score",
  maxLength      = "max_length",
  normArlemScore = "norm_arlem_score"
)


#### Constructors ####

#' Low-level constructor for a pairwise_distances object
#'
#' @param x A tibble.
#' @param subclass Optional entity subclass, \code{"tale_distances"} or
#'   \code{"domain_distances"}.
#' @param dom_code_namespace Optional scalar string, see
#'   \code{\link{tales_namespace}}.
#' @return A \code{pairwise_distances} object.
#' @keywords internal
new_pairwise_distances <- function(x, subclass = NULL, dom_code_namespace = NULL) {
  stopifnot(is.data.frame(x))
  x <- tibble::as_tibble(x)
  if (!is.null(dom_code_namespace)) {
    attr(x, "dom_code_namespace") <- dom_code_namespace
  }
  class(x) <- unique(c(subclass, "pairwise_distances", class(x)))
  x
}

#' Is this a pairwise similarity table?
#' @param x An object.
#' @return A logical scalar.
#' @export
#' @family pairwise distances
is_pairwise_distances <- function(x) inherits(x, "pairwise_distances")

#' Create a pairwise similarity table
#'
#' The long form of a square pairwise similarity over one entity set: one row
#' per ordered pair of entities, with a \code{sim} score.
#'
#' \code{tale_distances()} and \code{domain_distances()} are the entity-specific flavours:
#' similarity between whole TALE arrays, and between individual repeat units.
#' They add no structure, only semantics — every method is written once on the
#' parent.
#'
#' Legacy column spellings (\code{TAL1}/\code{TAL2}, \code{RepU1}/\code{RepU2},
#' \code{Sim}, \code{Dissim}, \code{arlemScore}, ...) are renamed on the way in.
#'
#' @param x A data frame with \code{id1}, \code{id2} and \code{dissim} columns.
#'   Legacy spellings and the similarity vocabulary are accepted and converted:
#'   a table carrying \code{Sim} or \code{normArlemScore} instead of a distance
#'   is folded into \code{dissim}, and those restatements are then dropped so
#'   only one copy of the quantity is stored. Further columns
#'   (\code{arlem_score}, \code{max_length}, or anything else) are preserved.
#' @param dom_code_namespace Optional scalar string identifying the run whose
#'   \code{dom_code} values this table is keyed by, see
#'   \code{\link{tales_namespace}}. Relevant for \code{domain_distances()}, whose ids
#'   *are* \code{dom_code}s.
#' @return A validated \code{pairwise_distances} object.
#' @export
#' @family pairwise distances
pairwise_distances <- function(x, dom_code_namespace = NULL) {
  .new_validated_distances(x, subclass = NULL, dom_code_namespace = dom_code_namespace)
}

#' @rdname pairwise_distances
#' @export
tale_distances <- function(x, dom_code_namespace = NULL) {
  .new_validated_distances(x, subclass = "tale_distances", dom_code_namespace = dom_code_namespace)
}

#' @rdname pairwise_distances
#' @export
domain_distances <- function(x, dom_code_namespace = NULL) {
  .new_validated_distances(x, subclass = "domain_distances", dom_code_namespace = dom_code_namespace)
}

#' @noRd
.new_validated_distances <- function(x, subclass, dom_code_namespace) {
  if (!is.data.frame(x)) {
    cli::cli_abort(
      "{.arg x} must be a data frame, not {.obj_type_friendly {x}}.",
      class = c("tantale_error_distances_type", "tantale_error")
    )
  }
  x <- .pairwise_distances_rename_legacy(x)
  # Only the distance is stored. Keeping several copies of one quantity invites
  # them drifting out of step (ledger 9.6), so sim and norm_arlem_score -- both
  # exact restatements of the distance -- are folded in and dropped.
  # norm_arlem_score is preferred as the source: it *is* the distance, whereas
  # deriving from sim round-trips through 100 - x for no gain.
  if (!"dissim" %in% names(x)) {
    if ("norm_arlem_score" %in% names(x)) {
      x[["dissim"]] <- as.numeric(x[["norm_arlem_score"]])
    } else if ("sim" %in% names(x)) {
      x[["dissim"]] <- 100 - as.numeric(x[["sim"]])
    }
  }
  x <- x[setdiff(names(x), c("sim", "norm_arlem_score"))]
  for (nm in intersect(PAIRWISE_DISTANCES_ID_COLS, names(x))) {
    x[[nm]] <- as.character(x[[nm]])
  }
  # Canonical column order, not just canonical names: repeat.similarity lists
  # its ids as RepU2, RepU1, so renaming alone would leave the two tables
  # ordered differently -- the rule-5 wart this class exists partly to fix.
  lead <- intersect(c(PAIRWISE_DISTANCES_ID_COLS, PAIRWISE_DISTANCES_VALUE_COL), names(x))
  x <- x[c(lead, setdiff(names(x), lead))]
  validate_pairwise_distances(
    new_pairwise_distances(x, subclass = subclass, dom_code_namespace = dom_code_namespace)
  )
}

#' @noRd
.pairwise_distances_rename_legacy <- function(x) {
  hit <- intersect(names(x), names(PAIRWISE_DISTANCES_LEGACY_NAMES))
  if (length(hit) == 0L) return(x)
  target <- unname(PAIRWISE_DISTANCES_LEGACY_NAMES[hit])
  clash <- intersect(target, names(x))
  if (length(clash) > 0L) {
    cli::cli_abort(
      c("Cannot rename legacy columns: the target name{?s} {.field {clash}} {?is/are} already present.",
        "i" = "Drop or rename the duplicate{?s} first."),
      class = c("tantale_error_distances_name_clash", "tantale_error")
    )
  }
  names(x)[match(hit, names(x))] <- target
  x
}


#### Validator ####

#' Validate a pairwise similarity table
#'
#' Checks the column contract, and nothing else. Squareness, the diagonal and
#' symmetry are deliberately *not* checked here: the package filters these
#' tables asymmetrically on purpose, so those properties are preconditions of
#' the methods that need them — see \code{\link{distances_assert_square}}.
#'
#' @param x A \code{pairwise_distances} object.
#' @return \code{x}, invisibly, if valid; otherwise an error.
#' @export
#' @family pairwise distances
validate_pairwise_distances <- function(x) {
  needed <- c(PAIRWISE_DISTANCES_ID_COLS, PAIRWISE_DISTANCES_VALUE_COL)
  missing <- setdiff(needed, names(x))
  if (length(missing) > 0L) {
    cli::cli_abort(
      "A {.cls pairwise_distances} object requires the column{?s} {.field {missing}}.",
      class = c("tantale_error_distances_missing_column", "tantale_error")
    )
  }
  for (nm in PAIRWISE_DISTANCES_ID_COLS) {
    if (!is.character(x[[nm]])) {
      cli::cli_abort(
        "{.field {nm}} must be a character vector, not {.obj_type_friendly {x[[nm]]}}.",
        class = c("tantale_error_distances_type", "tantale_error")
      )
    }
  }
  if (!is.numeric(x[[PAIRWISE_DISTANCES_VALUE_COL]])) {
    cli::cli_abort(
      "{.field {PAIRWISE_DISTANCES_VALUE_COL}} must be numeric, not {.obj_type_friendly {x[[PAIRWISE_DISTANCES_VALUE_COL]]}}.",
      class = c("tantale_error_distances_type", "tantale_error")
    )
  }
  if (nrow(x) == 0L) return(invisible(x))
  if (anyNA(x$id1) || anyNA(x$id2)) {
    cli::cli_abort("{.field id1} and {.field id2} must not contain {.val NA}.",
                   class = c("tantale_error_distances_na", "tantale_error"))
  }
  invisible(x)
}


#### Preconditions ####

#' Assert that a similarity table is complete and square
#'
#' Checks that the table holds every ordered pair of the ids it contains —
#' \code{n^2} rows for \code{n} ids. Phrased over the ids *present*, so that a
#' symmetric subset stays square: filtering both id columns to the same set of
#' entities preserves this, while filtering one of them does not.
#'
#' A **precondition**, not an invariant. \code{conversion.R} filters on
#' \code{id1} alone inside a loop over alignment columns, and \code{msa.R}
#' filters the two id columns in succession, passing through a non-square
#' intermediate. Both are correct; enforcing squareness everywhere would
#' outlaw them.
#'
#' @param x A \code{pairwise_distances} object.
#' @param arg Name of the argument being checked, for the error message.
#' @return \code{x}, invisibly.
#' @export
#' @family pairwise distances
distances_assert_square <- function(x, arg = "x") {
  if (!is_pairwise_distances(x)) {
    cli::cli_abort("{.arg {arg}} must be a {.cls pairwise_distances} object.",
                   class = c("tantale_error_distances_type", "tantale_error"))
  }
  ids <- union(x$id1, x$id2)
  n <- length(ids)
  if (nrow(x) != n^2) {
    missing1 <- setdiff(ids, x$id1)
    missing2 <- setdiff(ids, x$id2)
    cli::cli_abort(
      c("{.arg {arg}} must hold every pair of the {n} id{?s} it contains ({n^2} row{?s}), but has {nrow(x)}.",
        if (length(missing1)) c("x" = "Absent from {.field id1}: {.val {utils::head(missing1, 5)}}"),
        if (length(missing2)) c("x" = "Absent from {.field id2}: {.val {utils::head(missing2, 5)}}"),
        "i" = "Filtering both id columns to the same entities keeps a table square; filtering one does not."),
      class = c("tantale_error_distances_not_square", "tantale_error")
    )
  }
  invisible(x)
}


#### Views ####

#' Render a similarity table as a square matrix
#'
#' Materialises the wide form, with ids as both row and column names, sorted.
#' This replaces the hand-written \code{acast(x, id1 ~ id2, value.var = "sim")}
#' that appears at four call sites in the package.
#'
#' @param x A \code{pairwise_distances} object.
#' @param value Name of the column to fill cells with. Defaults to \code{sim}.
#' @param ... Ignored.
#' @return A numeric matrix, square, with sorted ids as dimnames.
#' @method as.matrix pairwise_distances
#' @export
as.matrix.pairwise_distances <- function(x, value = PAIRWISE_DISTANCES_VALUE_COL, ...) {
  if (!value %in% names(x)) {
    cli::cli_abort(
      c("This table has no {.field {value}} column.",
        "i" = "Available: {.field {setdiff(names(x), PAIRWISE_DISTANCES_ID_COLS)}}"),
      class = c("tantale_error_distances_layer", "tantale_error")
    )
  }
  distances_assert_square(x)
  ids <- sort(union(x$id1, x$id2))
  m <- matrix(NA_real_, nrow = length(ids), ncol = length(ids),
              dimnames = list(ids, ids))
  m[cbind(match(x$id1, ids), match(x$id2, ids))] <- as.numeric(x[[value]])
  m
}

#' Restrict a similarity table to a set of entities
#'
#' Filters **both** id columns, which is what keeps the result square. The
#' package currently does this by hand at three call sites, in two steps.
#'
#' @param x A \code{pairwise_distances} object.
#' @param ids A character vector of entity ids to keep.
#' @return A \code{pairwise_distances} over \code{ids} only.
#' @export
#' @family pairwise distances
distances_restrict <- function(x, ids) {
  if (!is_pairwise_distances(x)) {
    cli::cli_abort("{.arg x} must be a {.cls pairwise_distances} object.",
                   class = c("tantale_error_distances_type", "tantale_error"))
  }
  ids <- as.character(ids)
  absent <- setdiff(ids, union(x$id1, x$id2))
  if (length(absent) > 0L) {
    cli::cli_abort(
      c("{length(absent)} requested id{?s} {?is/are} absent from the table.",
        "x" = "{.val {utils::head(absent, 5)}}"),
      class = c("tantale_error_distances_unknown_id", "tantale_error")
    )
  }
  x[x$id1 %in% ids & x$id2 %in% ids, , drop = FALSE]
}


#### dplyr integration ####
# Same policy as 'tales': the column contract is re-checked cheaply on every
# verb, nothing else is. Squareness is not re-checked, by design.

#' @exportS3Method dplyr::dplyr_reconstruct
dplyr_reconstruct.pairwise_distances <- function(data, template) {
  if (!.pairwise_sim_contract_holds(data)) {
    return(.pairwise_sim_declass(data))
  }
  out <- NextMethod()
  attr(out, "dom_code_namespace") <- tales_namespace(template)
  out
}

# As for 'tales', dplyr's dplyr_col_select() does not call dplyr_reconstruct()
# for a tibble subclass, so column-dropping is caught here.
#' @export
`[.pairwise_distances` <- function(x, ...) {
  ns <- tales_namespace(x)
  out <- NextMethod()
  if (!is.data.frame(out)) return(out)
  if (!.pairwise_sim_contract_holds(out)) return(.pairwise_sim_declass(out))
  attr(out, "dom_code_namespace") <- ns
  out
}

#' @noRd
.pairwise_sim_contract_holds <- function(x) {
  all(c(PAIRWISE_DISTANCES_ID_COLS, PAIRWISE_DISTANCES_VALUE_COL) %in% names(x)) &&
    is.character(x$id1) && is.character(x$id2) &&
    is.numeric(x[[PAIRWISE_DISTANCES_VALUE_COL]])
}

#' @noRd
.pairwise_sim_declass <- function(x) {
  class(x) <- setdiff(class(x), c("tale_distances", "domain_distances", "pairwise_distances"))
  attr(x, "dom_code_namespace") <- NULL
  tibble::as_tibble(x)
}
