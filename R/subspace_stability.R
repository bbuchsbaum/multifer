#' Subspace stability from a bootstrap artifact
#'
#' For each unit, summarizes the principal angles between the unit's
#' loading subspace in the original fit and in each bootstrap rep,
#' across all domains. Singleton component units have a 1D subspace
#' and yield one angle per rep; subspace units (grouped tied roots
#' from \code{\link{form_units}()} opt-in mode) have a multi-d
#' subspace and yield multiple angles per rep.
#'
#' Per-unit summary across reps:
#' \itemize{
#'   \item \code{principal_angle_mean}: average across reps of the
#'     per-rep mean angle.
#'   \item \code{principal_angle_max}: maximum across reps of the
#'     per-rep maximum angle.
#' }
#'
#' Stability label heuristic (max angle in radians):
#' \itemize{
#'   \item \code{"stable"} if max < 0.1
#'   \item \code{"tied"} if the unit is multi-dimensional and max
#'     in \code{[0.1, pi/4]}
#'   \item \code{"unstable"} if max > pi/4
#' }
#'
#' Reuses \code{\link{principal_angles}()} (which prefers
#' \code{multivarious::prinang} when installed).
#'
#' @param artifact A \code{multifer_bootstrap_artifact}.
#' @param original_fit The output of \code{adapter$refit} on the
#'   original data, used as the alignment reference.
#' @param adapter The same \code{multifer_adapter} that produced
#'   \code{original_fit} -- needed to extract the original loadings
#'   per domain.
#' @param units A \code{multifer_units} table.
#'
#' @return A \code{multifer_subspace_stability} table.
#' @export
subspace_stability_from_bootstrap <- function(artifact,
                                               original_fit,
                                               adapter,
                                               units) {

  if (!inherits(artifact, "multifer_bootstrap_artifact")) {
    stop("`artifact` must be a multifer_bootstrap_artifact.", call. = FALSE)
  }
  if (!inherits(units, "multifer_units")) {
    stop("`units` must be a multifer_units table.", call. = FALSE)
  }
  if (!inherits(adapter, "multifer_adapter")) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }

  domains <- artifact$domains
  R       <- artifact$R
  reps    <- artifact$reps
  am      <- artifact$method_align

  unit_ids <- units$unit_id
  unit_type <- units$unit_type
  members  <- attr(units, "members")
  if (is.null(members)) members <- vector("list", nrow(units))

  orig_loadings <- stats::setNames(
    lapply(domains, function(d) adapter$loadings(original_fit, d)),
    domains
  )

  rows_unit_id  <- character(length(unit_ids))
  rows_pa_mean  <- numeric(length(unit_ids))
  rows_pa_max   <- numeric(length(unit_ids))
  rows_align    <- character(length(unit_ids))
  rows_label    <- character(length(unit_ids))
  row_pos <- 1L

  for (u in seq_along(unit_ids)) {
    uid         <- unit_ids[u]
    member_cols <- as.integer(members[[u]])
    if (length(member_cols) == 0L) next

    per_rep_mean_all <- rep(NA_real_, R)
    per_rep_max_all  <- rep(NA_real_, R)

    for (b in seq_len(R)) {
      per_domain_angles <- numeric(0)
      for (d in domains) {
        Vrep_all <- reps[[b]]$aligned_loadings[[d]]
        Vref_all <- orig_loadings[[d]]
        if (any(member_cols > ncol(Vref_all)) || any(member_cols > ncol(Vrep_all))) {
          next
        }
        Vref <- Vref_all[, member_cols, drop = FALSE]
        Vrep <- Vrep_all[, member_cols, drop = FALSE]
        ang  <- principal_angles(Vref, Vrep)
        per_domain_angles <- c(per_domain_angles, ang)
      }
      if (length(per_domain_angles) == 0L) {
        next
      }
      per_rep_mean_all[b] <- mean(per_domain_angles)
      per_rep_max_all[b]  <- max(per_domain_angles)
    }

    finite_mean <- per_rep_mean_all[is.finite(per_rep_mean_all)]
    finite_max  <- per_rep_max_all[is.finite(per_rep_max_all)]
    pa_mean <- if (length(finite_mean) > 0L) mean(finite_mean) else NA_real_
    pa_max  <- if (length(finite_max)  > 0L) max(finite_max)   else NA_real_

    label <- if (is.na(pa_max)) {
      "unknown"
    } else if (pa_max < 0.1) {
      "stable"
    } else if (unit_type[u] == "subspace" && pa_max <= pi / 4) {
      "tied"
    } else if (pa_max > pi / 4) {
      "unstable"
    } else {
      "stable"
    }

    rows_unit_id[row_pos] <- uid
    rows_pa_mean[row_pos] <- pa_mean
    rows_pa_max[row_pos]  <- pa_max
    # alignment_method maps to the schema's allowed set
    align_label <- if (length(member_cols) > 1L) "subspace" else am
    if (align_label == "sign") align_label <- "sign"
    if (align_label == "procrustes") align_label <- "procrustes"
    rows_align[row_pos] <- align_label
    rows_label[row_pos] <- label
    row_pos <- row_pos + 1L
  }

  if (row_pos == 1L) {
    return(infer_subspace_stability())
  }

  keep <- seq_len(row_pos - 1L)
  infer_subspace_stability(
    unit_id              = rows_unit_id[keep],
    principal_angle_mean = rows_pa_mean[keep],
    principal_angle_max  = rows_pa_max[keep],
    alignment_method     = rows_align[keep],
    stability_label      = rows_label[keep]
  )
}
