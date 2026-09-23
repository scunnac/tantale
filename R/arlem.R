

#### ARLEM's alignment score, computed in R ####
#
# A re-implementation of the ARLEM 1.0 program that tales_tale_distances()
# used to run (as `arlem -align -insert`), so that the TALE-level distances
# do not depend on a Linux-x86-64-only executable whose licence does not
# cover redistribution (ledger 28). The executable is gone; the code that
# drove it is in inst/legacy/arlem_binary.R.
#
# The model is ARLEM's own, from Abouelhoda, Giegerich, Behzadi & Steyaert,
# "Alignment of minisatellite maps: a minimum spanning tree-based approach",
# APBC 2008, pp. 261-272 (the JBCB 2009 version is
# https://doi.org/10.1142/S0219720009004060). An array is a "map" of units
# (here, domain codes). Two maps are aligned with three kinds of event:
#
#   - a match of two units, costing M(a, b), the substitution matrix;
#   - a duplication, where a unit gives rise to a tandem copy of itself that
#     may then diverge: d(a, b) = Dup + M(a, b);
#   - an insertion, a unit that arose from neither (cost `Indel hist`).
#
# Runs of units that one array has and the other lacks must be explained as
# a duplication history grown from the unit at their left or right edge
# (section 3 of the paper: optimal left/right "ordered directed spanning
# trees", Cl and Cr below), and the alignment recurrence (section 4.2)
# chooses between matching units and absorbing a run into such a history.
# A sentinel unit "$" is prepended to every map so that a leading run can be
# explained by insertions; mutating "$" into anything costs 99999, as in the
# binary.
#
# Checked against the binary: identical scores on ~6000 random map pairs
# (alphabets of 2-60 types, maps of 1-40 units, a range of Dup/indel costs,
# with and without insertions) and on real TALE arrays (ledger 33). The
# binary's own answers for a subset are kept as a test fixture
# (tests/testthat/data_for_tests/arlem_reference_scores.rds). Two things
# observed in the binary that the paper does not say: of its cost file's
# two indel costs only `Indel hist` was used, and a score is the plain
# integer sum of these costs.

#' Column minima of a matrix
#'
#' Most of the run time is spent here; `matrixStats::colMins()` is about
#' seven times faster than a base-R equivalent (ledger 33).
#' @noRd
.col_min <- function(x) matrixStats::colMins(x)

#' Optimal duplication histories over every interval of one map
#'
#' @param u Integer type indices of the map's units, "$" first.
#' @param M Substitution cost matrix, indexed by type.
#' @param dup,ins Duplication and insertion costs.
#' @param insert Allow insertions (ARLEM's `-insert`).
#' @return `list(l, r, r_upper)`: `l[i, j]` is the cheapest history
#'   producing units `i..j` from unit `i` alone (a left history), `r[i, j]`
#'   the same grown from unit `j` (a right history); `r_upper` is `r` with
#'   `Inf` below the diagonal.
#' @noRd
.arlem_histories <- function(u, M, dup, ins, insert = TRUE) {
  n <- length(u)
  # edge[a, b]: the cost of unit a giving rise to unit b
  edge <- dup + M[u, u, drop = FALSE]
  if (insert) edge <- pmin(edge, ins)
  cl <- cr <- matrix(0, n, n)
  if (n < 2L) return(list(l = cl, r = cr, r_upper = cr))
  for (len in seq_len(n - 1L)) {
    for (i in seq_len(n - len)) {
      j <- i + len
      if (len == 1L) {
        cl[i, j] <- edge[i, j]
        cr[i, j] <- edge[j, i]
        next
      }
      # j is the root's last child, rooting a right history of its own
      k <- i:(j - 1L)
      split <- min(cl[i, k] + cr[k + 1L, j])
      # or the tree splits at an inner node k into two histories
      k <- (i + 1L):(j - 1L)
      cl[i, j] <- min(cl[i, k] + cl[k, j], split + edge[i, j])
      cr[i, j] <- min(cr[i, k] + cr[k, j], split + edge[j, i])
    }
  }
  # r with the cells below the diagonal (no such interval) made unusable,
  # so .arlem_align() can take a whole column's minimum in one call
  r_upper <- cr
  r_upper[lower.tri(r_upper)] <- Inf
  list(l = cl, r = cr, r_upper = r_upper)
}

#' ARLEM's alignment score of two maps
#'
#' `A[i, j]` is the cost of aligning the first `i` units of `s` with the
#' first `j` of `r` ("$" at index 1). `B` is the paper's A' table: the best
#' way to end in `r[tr..j]` grown rightward-from `r[j]`, stored so that
#' simultaneous right duplications in both maps cost O(n^3), not O(n^4).
#' @param s,r Integer type indices, "$" first.
#' @param hs,hr Their histories, from `.arlem_histories()`.
#' @noRd
.arlem_align <- function(s, r, hs, hr, M) {
  n <- length(s)
  m <- length(r)
  Msr <- M[s, r, drop = FALSE]
  Msym <- pmin(Msr, t(M[r, s, drop = FALSE]))
  A <- B <- matrix(Inf, n, m)
  jj <- seq_len(m)[-1L]
  for (i in seq_len(n)) {
    row <- rep(Inf, m)
    if (i == 1L) {
      row[1L] <- 0
    } else {
      # s[i] matched with r[j]
      row[jj] <- Msr[i, jj] + A[i - 1L, jj - 1L]
      # s[(l + 1):i] grown from s[l]
      l <- seq_len(i - 1L)
      row <- pmin(row, .col_min(A[l, , drop = FALSE] + hs$l[l, i]))
      # s[ts:i] and r[tr:j] grown from s[i] and r[j], which match
      if (m > 1L) {
        ts <- 2:i
        row[jj] <- pmin(row[jj],
                        .col_min(B[ts - 1L, jj, drop = FALSE] + hs$r[ts, i]) +
                          Msym[i, jj])
      }
    }
    # r[(k + 1):j] grown from r[k]: needs this row's own earlier cells
    for (j in jj) {
      k <- seq_len(j - 1L)
      row[j] <- min(row[j], row[k] + hr$l[k, j])
    }
    A[i, ] <- row
    if (m > 1L) B[i, jj] <- .col_min(row[-m] + hr$r_upper[jj, jj, drop = FALSE])
  }
  A[n, m]
}

#' Every pairwise ARLEM score between arrays, computed in R
#'
#' The scores the ARLEM program printed, in the shape its parser returned:
#' `id1`/`id2` are 0-based record indices of `coded`, `id1 < id2`.
#'
#' @param coded The arrays as space-separated domain codes,
#'   from `.coded_seq_set()`.
#' @param cost The substitution matrix, from `.arlem_cost_matrix()`, with
#'   domain codes as dimnames.
#' @return A tibble `id1`, `id2`, `arlem_score`.
#' @noRd
.arlem_scores_r <- function(coded, cost, dup = .arlem_dup_cost,
                            ins = .arlem_indel_cost, insert = TRUE) {
  types <- rownames(cost)
  k <- length(types)
  dollar <- k + 1L
  M <- matrix(99999, dollar, dollar)
  M[seq_len(k), seq_len(k)] <- cost
  M[dollar, dollar] <- 0

  maps <- strsplit(as.character(coded), " ", fixed = TRUE)
  u <- lapply(maps, function(x) c(dollar, match(x, types)))
  unknown <- setdiff(unlist(maps), types)
  if (length(unknown)) {
    cli::cli_abort(
      c("Some domain codes in the arrays have no row in the cost matrix.",
        "x" = "Missing: {.val {utils::head(unknown, 5)}}."),
      class = c("tantale_error_arlem_codes", "tantale_error"))
  }

  if (length(u) < 2L) {
    return(tibble::tibble(id1 = integer(), id2 = integer(),
                          arlem_score = numeric()))
  }
  cli::cli_inform("Aligning {length(u)} TALE arrays pairwise ({choose(length(u), 2)} pair{?s}).")
  hist <- lapply(u, .arlem_histories, M = M, dup = dup, ins = ins,
                 insert = insert)
  pairs <- utils::combn(length(u), 2L)
  score <- vapply(seq_len(ncol(pairs)), function(p) {
    a <- pairs[1L, p]
    b <- pairs[2L, p]
    .arlem_align(u[[a]], u[[b]], hist[[a]], hist[[b]], M)
  }, numeric(1))
  tibble::tibble(id1 = pairs[1L, ] - 1L, id2 = pairs[2L, ] - 1L,
                 arlem_score = score)
}
