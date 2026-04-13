#' Score stability from a bootstrap artifact
#'
#' \strong{Design lock (Part 5 review).} Score stability is computed by
#' projecting the \strong{original} observations through each bootstrap
#' fit, NOT by summarizing scores on the bootstrap sample itself.
#' Bootstrap samples have duplicates and omissions, which makes
#' per-observation summaries on the resampled data interpretable only
#' for the random subset of observations that happen to be present.
#' Projecting the original data through each rep's loadings yields a
#' clean per-observation distribution of scores under the rep's
#' fitted basis.
#'
#' Concretely, for replicate \code{b} the score of observation \code{i}
#' on the unit's first member column is
#' \deqn{s_{i,b} = (X_i - \bar X)^t v_b}
#' where \code{v_b} is the rep's aligned loading column and the
#' centering uses the column means of the ORIGINAL data.
#'
#' For multi-member subspace units we report stability of the
#' per-observation Frobenius norm across the unit's columns, which is
#' rotation-invariant within the subspace.
#'
#' @param artifact A \code{multifer_bootstrap_artifact}.
#' @param original_data For oneblock, a numeric matrix. For cross, a
#'   list with \code{X} and \code{Y}. Must be the same data that was
#'   passed to \code{\link{bootstrap_fits}()}.
#' @param units A \code{multifer_units} table.
#' @param quantiles Length-2 numeric in \code{[0, 1]} for the lower and
#'   upper interval bounds. Default \code{c(0.025, 0.975)}.
#'
#' @return A \code{multifer_score_stability} table keyed by
#'   \code{(unit_id, domain, observation)}.
#' @export
score_stability_from_bootstrap <- function(artifact,
                                           original_data,
                                           units,
                                           quantiles = c(0.025, 0.975)) {

  if (!inherits(artifact, "multifer_bootstrap_artifact")) {
    stop("`artifact` must be a multifer_bootstrap_artifact.", call. = FALSE)
  }
  if (!inherits(units, "multifer_units")) {
    stop("`units` must be a multifer_units table.", call. = FALSE)
  }
  if (!is.numeric(quantiles) || length(quantiles) != 2L) {
    stop("`quantiles` must be a length-2 numeric vector.", call. = FALSE)
  }

  domains <- artifact$domains
  R       <- artifact$R
  reps    <- artifact$reps

  unit_ids <- units$unit_id
  members  <- attr(units, "members")
  if (is.null(members)) members <- vector("list", nrow(units))

  ## Resolve per-domain data matrices and centers.
  data_by_domain <- list()
  if (identical(domains, "X")) {
    if (!is.matrix(original_data)) {
      stop("For oneblock, `original_data` must be a matrix.", call. = FALSE)
    }
    data_by_domain$X <- sweep(original_data, 2L, colMeans(original_data), "-")
  } else {
    if (!is.list(original_data) ||
        is.null(original_data$X) || is.null(original_data$Y)) {
      stop("For cross, `original_data` must be a list with X and Y.", call. = FALSE)
    }
    data_by_domain$X <- sweep(original_data$X, 2L, colMeans(original_data$X), "-")
    data_by_domain$Y <- sweep(original_data$Y, 2L, colMeans(original_data$Y), "-")
  }

  rows_unit_id     <- character(0)
  rows_domain      <- character(0)
  rows_observation <- integer(0)
  rows_estimate    <- numeric(0)
  rows_lower       <- numeric(0)
  rows_upper       <- numeric(0)
  rows_leverage    <- numeric(0)

  for (u in seq_along(unit_ids)) {
    uid         <- unit_ids[u]
    member_cols <- as.integer(members[[u]])
    if (length(member_cols) == 0L) next

    for (d in domains) {
      Xd  <- data_by_domain[[d]]
      n   <- nrow(Xd)

      # Per-rep score matrix for this unit, shape n x length(member_cols).
      # We collapse over members via Frobenius norm so subspace units
      # are rotation-invariant.
      rep_norms <- matrix(NA_real_, nrow = n, ncol = R)
      for (b in seq_len(R)) {
        Lb <- reps[[b]]$aligned_loadings[[d]]
        if (any(member_cols > ncol(Lb))) next
        S_b <- Xd %*% Lb[, member_cols, drop = FALSE]
        rep_norms[, b] <- sqrt(rowSums(S_b * S_b))
      }

      estimate <- rowMeans(rep_norms, na.rm = TRUE)
      lower    <- apply(rep_norms, 1L, function(v) {
        v <- v[is.finite(v)]
        if (length(v) == 0L) return(NA_real_)
        as.numeric(stats::quantile(v, quantiles[1L]))
      })
      upper    <- apply(rep_norms, 1L, function(v) {
        v <- v[is.finite(v)]
        if (length(v) == 0L) return(NA_real_)
        as.numeric(stats::quantile(v, quantiles[2L]))
      })

      rows_unit_id     <- c(rows_unit_id,     rep(uid, n))
      rows_domain      <- c(rows_domain,      rep(d,   n))
      rows_observation <- c(rows_observation, seq_len(n))
      rows_estimate    <- c(rows_estimate,    estimate)
      rows_lower       <- c(rows_lower,       lower)
      rows_upper       <- c(rows_upper,       upper)
      rows_leverage    <- c(rows_leverage,    rep(NA_real_, n))
    }
  }

  infer_score_stability(
    unit_id     = rows_unit_id,
    domain      = rows_domain,
    observation = rows_observation,
    estimate    = rows_estimate,
    lower       = rows_lower,
    upper       = rows_upper,
    leverage    = rows_leverage
  )
}
