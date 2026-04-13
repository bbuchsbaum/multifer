#' Form latent units from a vector of roots
#'
#' Converts a sorted vector of singular values (or eigenvalues) into a
#' `multifer_units` table. The default behaviour treats every root as an
#' independent singleton component. Near-tie grouping into subspace units
#' is available as an explicit opt-in (see `group_near_ties`).
#'
#' @param roots Numeric vector of roots (singular values / eigenvalues),
#'   sorted in descending order. All values must be non-negative.
#' @param selected Logical vector of the same length as `roots`, marking
#'   which roots were declared significant by the test ladder. If `NULL`
#'   (the default), all roots are treated as selected.
#' @param group_near_ties Logical scalar (default `FALSE`). When `TRUE`,
#'   consecutive roots whose relative gap falls below `tie_threshold` are
#'   merged into a single `"subspace"` unit. When `FALSE` (the v1
#'   default), every root becomes its own `"component"` unit regardless
#'   of proximity.
#' @param tie_threshold Single non-negative numeric. The relative-gap
#'   threshold used when `group_near_ties = TRUE`. A pair `(r_i, r_{i+1})`
#'   is considered tied when
#'   `(r_i - r_{i+1}) / max(r_i, .Machine$double.eps) < tie_threshold`.
#'   Default `0.01`.
#'
#' @return A `multifer_units` table built via [infer_units()].
#'
#' @details
#' Near-tie grouping algorithm (opt-in only): walk roots in order.
#' Compare each consecutive pair via the relative gap. While the gap is
#' below `tie_threshold`, extend the current group. When the gap meets or
#' exceeds the threshold, close the current group and open a new one.
#' A group of size 1 becomes a `"component"` unit with `identifiable =
#' TRUE`; a group of size >= 2 becomes a `"subspace"` unit with
#' `identifiable = FALSE`. A group is `selected = TRUE` only when every
#' member root is individually selected in the caller-supplied `selected`
#' vector (a partially-selected subspace carries `selected = FALSE`).
#'
#' @seealso [infer_units()], [infer_result()]
#'
#' @examples
#' # Default: every root becomes its own component unit.
#' form_units(c(5, 3, 1))
#'
#' # Opt-in grouping: roots 1 and 2 are near-tied, root 3 is isolated.
#' form_units(c(5.0001, 5, 1), group_near_ties = TRUE, tie_threshold = 0.01)
#'
#' @export
form_units <- function(roots,
                       selected        = NULL,
                       group_near_ties = FALSE,
                       tie_threshold   = 0.01) {

  # ---- input validation ------------------------------------------------------

  if (!is.numeric(roots)) {
    stop("`roots` must be a numeric vector.", call. = FALSE)
  }
  if (any(roots < 0)) {
    stop("`roots` must be non-negative.", call. = FALSE)
  }
  if (length(roots) > 1L && any(diff(roots) > 0)) {
    stop("`roots` must be sorted in descending order.", call. = FALSE)
  }
  if (is.null(selected)) {
    selected <- rep(TRUE, length(roots))
  } else {
    if (!is.logical(selected)) {
      stop("`selected` must be a logical vector.", call. = FALSE)
    }
    if (length(selected) != length(roots)) {
      stop("`selected` must have the same length as `roots`.", call. = FALSE)
    }
  }
  if (!(is.numeric(tie_threshold) && length(tie_threshold) == 1L &&
        !is.na(tie_threshold) && tie_threshold >= 0)) {
    stop("`tie_threshold` must be a single non-negative numeric.", call. = FALSE)
  }

  n <- length(roots)

  # ---- empty case ------------------------------------------------------------

  if (n == 0L) {
    return(infer_units())
  }

  # ---- build groups ----------------------------------------------------------

  if (!group_near_ties) {
    # Default: every root is its own singleton group.
    groups <- as.list(seq_len(n))
  } else {
    # Opt-in near-tie grouping.
    groups        <- list()
    current_group <- integer(0)

    for (i in seq_len(n)) {
      current_group <- c(current_group, i)

      # Decide whether to close the current group:
      #   - always close on the last root
      #   - close when the gap to the next root exceeds the threshold
      close_group <- (i == n)
      if (!close_group) {
        gap_rel <- (roots[i] - roots[i + 1L]) /
                   max(roots[i], .Machine$double.eps)
        close_group <- (gap_rel >= tie_threshold)
      }

      if (close_group) {
        groups        <- c(groups, list(current_group))
        current_group <- integer(0)
      }
    }
  }

  # ---- build unit table columns from groups ----------------------------------

  n_units    <- length(groups)
  unit_id    <- sprintf("u%d", seq_len(n_units))
  unit_type  <- character(n_units)
  members    <- vector("list", n_units)
  identif    <- logical(n_units)
  sel        <- logical(n_units)

  for (k in seq_len(n_units)) {
    g          <- groups[[k]]          # integer vector of root indices
    members[[k]] <- as.integer(g)

    if (length(g) == 1L) {
      unit_type[k] <- "component"
      identif[k]   <- TRUE
    } else {
      unit_type[k] <- "subspace"
      identif[k]   <- FALSE
    }

    # A group is selected only when ALL member roots are selected.
    sel[k] <- all(selected[g])
  }

  # ---- construct and return --------------------------------------------------

  infer_units(
    unit_id      = unit_id,
    unit_type    = unit_type,
    members      = members,
    identifiable = identif,
    selected     = sel
  )
}
