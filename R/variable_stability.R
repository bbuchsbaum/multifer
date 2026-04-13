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

  rows_unit_id  <- character(0)
  rows_domain   <- character(0)
  rows_variable <- character(0)
  rows_estimate <- numeric(0)
  rows_stab     <- numeric(0)

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

      # Collect the unit's loading SLAB for every replicate:
      # array of shape (n_var, length(member_cols), R).
      slab <- array(NA_real_,
                    dim = c(n_var, length(member_cols), R))
      for (b in seq_len(R)) {
        Lb <- reps[[b]]$aligned_loadings[[d]]
        if (any(member_cols > ncol(Lb))) next
        slab[, , b] <- Lb[, member_cols, drop = FALSE]
      }

      # Per-variable stability: pool columns by row Frobenius norm so
      # tied subspaces are not penalized for in-subspace rotation.
      # rep_norms[i, b] = sqrt(sum slab[i, , b]^2)
      rep_norms <- apply(slab, c(1L, 3L), function(v) sqrt(sum(v^2)))
      # rep_norms is (n_var x R)

      mean_norm <- rowMeans(rep_norms, na.rm = TRUE)
      sd_norm   <- apply(rep_norms, 1L, function(v) {
        v <- v[is.finite(v)]
        if (length(v) < 2L) return(0)
        stats::sd(v)
      })
      stability <- 1 / (1 + sd_norm / (abs(mean_norm) + eps))

      # estimate: signed mean of the FIRST member column per variable
      first_col_slab <- slab[, 1L, , drop = TRUE]
      if (is.null(dim(first_col_slab))) {
        first_col_slab <- matrix(first_col_slab, nrow = n_var)
      }
      estimate <- rowMeans(first_col_slab, na.rm = TRUE)

      rows_unit_id  <- c(rows_unit_id,  rep(uid, n_var))
      rows_domain   <- c(rows_domain,   rep(d,   n_var))
      rows_variable <- c(rows_variable, var_names)
      rows_estimate <- c(rows_estimate, estimate)
      rows_stab     <- c(rows_stab,     stability)
    }
  }

  n <- length(rows_unit_id)
  infer_variable_stability(
    unit_id            = rows_unit_id,
    domain             = rows_domain,
    variable           = rows_variable,
    estimate           = rows_estimate,
    stability          = rows_stab,
    selection_freq     = rep(NA_real_, n),
    weight_sensitivity = rep(NA_real_, n)
  )
}
