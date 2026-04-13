#' Variable stability from a bootstrap artifact
#'
#' Consumes a \code{multifer_bootstrap_artifact} produced by
#' \code{\link{bootstrap_fits}()} and produces an
#' \code{\link{infer_variable_stability}()} table keyed by
#' \code{unit_id}, \code{domain}, and \code{variable}.
#'
#' For each unit, each domain, and each variable in that domain, we
#' summarize the column of the unit's loadings across bootstrap
#' replicates as:
#' \itemize{
#'   \item \code{estimate}: the mean across reps of the (signed)
#'     loading value at this variable.
#'   \item \code{stability}: a normalized stability measure in
#'     \code{[0, 1]}, defined as \code{1 / (1 + sd / (|mean| + eps))}.
#'     Equals 1 when sd = 0 (perfect agreement) and trends to 0 as
#'     variation across reps overwhelms the mean.
#' }
#'
#' Multi-member subspace units pool members: stability is computed on
#' the per-variable Frobenius norm across the unit's columns rather
#' than on a single column, because the individual loading vectors
#' inside a tied subspace are not identifiable.
#'
#' This is a \strong{stability measure, not a p-value}. Variable
#' significance is deferred to Phase 3 per Part 5 section 38.
#'
#' @param artifact A \code{multifer_bootstrap_artifact}.
#' @param units A \code{multifer_units} table.
#'
#' @return A \code{multifer_variable_stability} table.
#' @export
variable_stability_from_bootstrap <- function(artifact, units) {

  if (!inherits(artifact, "multifer_bootstrap_artifact")) {
    stop("`artifact` must be a multifer_bootstrap_artifact.", call. = FALSE)
  }
  if (!inherits(units, "multifer_units")) {
    stop("`units` must be a multifer_units table.", call. = FALSE)
  }

  domains <- artifact$domains
  R       <- artifact$R
  reps    <- artifact$reps
  eps     <- 1e-12

  if (R < 1L || length(reps) == 0L) {
    return(infer_variable_stability())
  }

  unit_ids <- units$unit_id
  members  <- attr(units, "members")
  if (is.null(members)) members <- vector("list", nrow(units))

  template_dims <- vapply(
    domains,
    function(d) nrow(reps[[1L]]$aligned_loadings[[d]]),
    integer(1L)
  )
  n_rows <- length(unit_ids) * sum(template_dims)

  rows_unit_id  <- character(n_rows)
  rows_domain   <- character(n_rows)
  rows_variable <- character(n_rows)
  rows_estimate <- numeric(n_rows)
  rows_stab     <- numeric(n_rows)
  row_pos <- 1L

  for (u in seq_along(unit_ids)) {
    uid     <- unit_ids[u]
    member_cols <- as.integer(members[[u]])
    if (length(member_cols) == 0L) next

    for (d in domains) {

      # First rep tells us how many variables in this domain.
      first_L <- reps[[1L]]$aligned_loadings[[d]]
      n_var   <- nrow(first_L)
      var_names <- if (!is.null(rownames(first_L))) rownames(first_L)
                   else as.character(seq_len(n_var))

      rep_norms <- matrix(NA_real_, nrow = n_var, ncol = R)
      first_member <- matrix(NA_real_, nrow = n_var, ncol = R)
      for (b in seq_len(R)) {
        Lb <- reps[[b]]$aligned_loadings[[d]]
        if (any(member_cols > ncol(Lb))) next
        unit_loadings <- Lb[, member_cols, drop = FALSE]
        rep_norms[, b] <- sqrt(rowSums(unit_loadings * unit_loadings))
        first_member[, b] <- unit_loadings[, 1L]
      }

      stats_norm <- .row_mean_sd(rep_norms)
      mean_norm <- stats_norm$mean
      sd_norm   <- stats_norm$sd
      stability <- 1 / (1 + sd_norm / (abs(mean_norm) + eps))

      estimate <- rowMeans(first_member, na.rm = TRUE)

      idx <- row_pos:(row_pos + n_var - 1L)
      rows_unit_id[idx]  <- uid
      rows_domain[idx]   <- d
      rows_variable[idx] <- var_names
      rows_estimate[idx] <- estimate
      rows_stab[idx]     <- stability
      row_pos <- row_pos + n_var
    }
  }

  if (row_pos == 1L) {
    return(infer_variable_stability())
  }

  keep <- seq_len(row_pos - 1L)
  n <- length(keep)
  infer_variable_stability(
    unit_id            = rows_unit_id[keep],
    domain             = rows_domain[keep],
    variable           = rows_variable[keep],
    estimate           = rows_estimate[keep],
    stability          = rows_stab[keep],
    selection_freq     = rep(NA_real_, n),
    weight_sensitivity = rep(NA_real_, n)
  )
}
