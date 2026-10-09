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

#' Is this a pairwise distance table?
#' @param x An object.
#' @return A logical scalar.
#' @examples
#' d <- data.frame(
#'   id1 = c("A1", "A1", "A2", "A2"),
#'   id2 = c("A1", "A2", "A1", "A2"),
#'   dissim = c(0, 35, 35, 0)
#' )
#' is_pairwise_distances(d)
#' is_pairwise_distances(pairwise_distances(d))
#' @export
#' @family pairwise distances
is_pairwise_distances <- function(x) inherits(x, "pairwise_distances")

#' Create a pairwise distance table
#'
#' The long form of a square pairwise distance over one entity set: one row
#' per ordered pair of entities, with a \code{dissim} value (0 for identical
#' entities, larger for more different ones).
#'
#' \code{tale_distances()} and \code{domain_distances()} are the
#' entity-specific flavours: distances between whole TALE arrays, and between
#' distinct domains, repeats and the two termini alike. For
#' \code{domain_distances}, \code{dissim} is the percentage of amino acids
#' that differ between two domains (see \code{\link{tales_domain_distances}});
#' for \code{tale_distances}, it is the array alignment cost between two TALEs
#' (see \code{\link{tales_tale_distances}}). The two flavours add no
#' structure, only meaning: every method is written once on the parent.
#'
#' @param x A data frame with \code{id1}, \code{id2} and \code{dissim} columns.
#'   The function also accepts a similarity: it converts a table carrying
#'   \code{sim} or \code{norm_arlem_score} instead of a distance into
#'   \code{dissim} and drops the similarity columns, so that only one copy
#'   of the quantity remains. It keeps any further column
#'   (\code{arlem_score}, \code{max_length}, or anything else).
#' @param dom_code_namespace Optional scalar string identifying the run whose
#'   \code{dom_code} values this table is keyed by, see
#'   \code{\link{tales_namespace}}. Relevant for \code{domain_distances()}, whose ids
#'   *are* \code{dom_code}s.
#' @return A validated \code{pairwise_distances} object.
#' @export
#' @family pairwise distances
#' @examples
#' d <- data.frame(
#'   id1 = c("A1", "A1", "A2", "A2"),
#'   id2 = c("A1", "A2", "A1", "A2"),
#'   dissim = c(0, 35, 35, 0)
#' )
#' pairwise_distances(d)
#' domain_distances(d, dom_code_namespace = "example")
#'
#' # A similarity is folded into the distance.
#' similar <- data.frame(id1 = "A1", id2 = "A2", sim = 65)
#' tale_distances(similar)
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
  # Canonical column order: ids first, then the distance, whatever order the
  # input had them in.
  lead <- intersect(c(PAIRWISE_DISTANCES_ID_COLS, PAIRWISE_DISTANCES_VALUE_COL), names(x))
  x <- x[c(lead, setdiff(names(x), lead))]
  validate_pairwise_distances(
    new_pairwise_distances(x, subclass = subclass, dom_code_namespace = dom_code_namespace)
  )
}



#### Validator ####

#' Validate a pairwise distance table
#'
#' Checks the column contract, and nothing else. Squareness, the diagonal and
#' symmetry are preconditions of the methods that need them, checked by
#' \code{\link{distances_assert_square}}: a table filtered on one id column is
#' a legitimate intermediate step.
#'
#' @param x A \code{pairwise_distances} object.
#' @return \code{x}, invisibly, if valid; otherwise an error.
#' @examples
#' d <- data.frame(
#'   id1 = c("A1", "A1", "A2", "A2"),
#'   id2 = c("A1", "A2", "A1", "A2"),
#'   dissim = c(0, 35, 35, 0)
#' )
#' validate_pairwise_distances(pairwise_distances(d))
#'
#' # Without its dissim column, a table is not a pairwise_distances
#' try(validate_pairwise_distances(d[, c("id1", "id2")]))
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

#' Assert that a distance table is complete and square
#'
#' Checks that the table holds every ordered pair of the ids it contains —
#' \code{n^2} rows for \code{n} ids. Phrased over the ids *present*, so that a
#' symmetric subset stays square: filtering both id columns to the same set of
#' entities preserves this, while filtering one of them does not.
#'
#' Squareness is a **precondition** of the functions that need it. The class
#' itself does not require it, because filtering the two id columns one
#' after the other passes through a legitimate non-square intermediate.
#'
#' @param x A \code{pairwise_distances} object.
#' @param arg Name of the argument being checked, for the error message.
#' @return \code{x}, invisibly.
#' @export
#' @family pairwise distances
#' @examples
#' d <- pairwise_distances(data.frame(
#'   id1 = c("A1", "A1", "A2", "A2"),
#'   id2 = c("A1", "A2", "A1", "A2"),
#'   dissim = c(0, 35, 35, 0)
#' ))
#' distances_assert_square(d)
#'
#' try(distances_assert_square(d[1:3, ])) # missing the A2-A2 pair -- errors
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

#' Render a distance table as a square matrix
#'
#' Materialises the wide form, with ids as both row and column names, sorted.
#'
#' \code{x} must be square (see \code{\link{distances_assert_square}}, called
#' here for you): every id present in both \code{id1} and \code{id2}. A table
#' filtered on one id column only is not; use
#' \code{\link{distances_restrict}} to filter both.
#'
#' @param x A \code{pairwise_distances} object.
#' @param value Name of the column to fill cells with. Defaults to
#'   \code{dissim}.
#' @param ... Ignored.
#' @return A numeric matrix, square, with sorted ids as dimnames.
#' @method as.matrix pairwise_distances
#' @export
#' @family pairwise distances
#' @examples
#' d <- pairwise_distances(data.frame(
#'   id1 = c("A1", "A1", "A2", "A2"),
#'   id2 = c("A1", "A2", "A1", "A2"),
#'   dissim = c(0, 35, 35, 0)
#' ))
#' as.matrix(d)
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

#' Restrict a distance table to a set of entities
#'
#' Filters **both** id columns, which is what keeps the result square.
#'
#' @param x A \code{pairwise_distances} object.
#' @param ids A character vector of entity ids to keep.
#' @return A \code{pairwise_distances} over \code{ids} only.
#' @export
#' @family pairwise distances
#' @examples
#' d <- pairwise_distances(data.frame(
#'   id1 = c("A1", "A1", "A1", "A2", "A2", "A2", "A3", "A3", "A3"),
#'   id2 = c("A1", "A2", "A3", "A1", "A2", "A3", "A1", "A2", "A3"),
#'   dissim = c(0, 35, 60, 35, 0, 40, 60, 40, 0)
#' ))
#' distances_restrict(d, c("A1", "A2"))
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

#' Subset a pairwise_distances object
#'
#' Subsets rows and columns like an ordinary tibble, but the class travels
#' only while the result still satisfies the contract: \code{id1}/\code{id2}
#' (both character) and \code{dissim} (numeric) all present. Dropping any of
#' them leaves something that can no longer be described as a pairwise
#' distance table, and the result is a plain tibble, without an error.
#'
#' Squareness is not checked here (see \code{\link{distances_assert_square}}):
#' subsetting rows is an everyday way to end up with a non-square table.
#'
#' @param x A \code{pairwise_distances} object.
#' @param ... Passed on to the tibble/data frame method.
#' @return The subset: still a \code{pairwise_distances} (or its
#'   \code{tale_distances}/\code{domain_distances} subclass) if the contract
#'   holds, otherwise a plain tibble.
#' @export
#' @family pairwise distances
#' @examples
#' d <- pairwise_distances(data.frame(
#'   id1 = c("A1", "A1", "A2", "A2"),
#'   id2 = c("A1", "A2", "A1", "A2"),
#'   dissim = c(0, 35, 35, 0)
#' ))
#' is_pairwise_distances(d[1:2, ]) # row subsetting never breaks the contract
#' class(d[, "id1"])               # dropping dissim/id2: the class steps aside
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
