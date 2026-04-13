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

  rows_unit_id  <- character(0)
  rows_pa_mean  <- numeric(0)
  rows_pa_max   <- numeric(0)
  rows_align    <- character(0)
  rows_label    <- character(0)

  for (u in seq_along(unit_ids)) {
    uid         <- unit_ids[u]
    member_cols <- as.integer(members[[u]])
    if (length(member_cols) == 0L) next

    per_rep_mean_all <- numeric(0)
    per_rep_max_all  <- numeric(0)

    for (b in seq_len(R)) {
      per_domain_angles <- numeric(0)
      for (d in domains) {
        Vref <- orig_loadings[[d]][, member_cols, drop = FALSE]
        Vrep <- reps[[b]]$aligned_loadings[[d]][, member_cols, drop = FALSE]
        ang  <- principal_angles(Vref, Vrep)
        per_domain_angles <- c(per_domain_angles, ang)
      }
      if (length(per_domain_angles) == 0L) next
      per_rep_mean_all <- c(per_rep_mean_all, mean(per_domain_angles))
      per_rep_max_all  <- c(per_rep_max_all,  max(per_domain_angles))
    }

    pa_mean <- if (length(per_rep_mean_all) > 0L) mean(per_rep_mean_all) else NA_real_
    pa_max  <- if (length(per_rep_max_all)  > 0L) max(per_rep_max_all)   else NA_real_

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

    rows_unit_id <- c(rows_unit_id, uid)
    rows_pa_mean <- c(rows_pa_mean, pa_mean)
    rows_pa_max  <- c(rows_pa_max,  pa_max)
    # alignment_method maps to the schema's allowed set
    align_label <- if (length(member_cols) > 1L) "subspace" else am
    if (align_label == "sign") align_label <- "sign"
    if (align_label == "procrustes") align_label <- "procrustes"
    rows_align <- c(rows_align, align_label)
    rows_label <- c(rows_label, label)
  }

  infer_subspace_stability(
    unit_id              = rows_unit_id,
    principal_angle_mean = rows_pa_mean,
    principal_angle_max  = rows_pa_max,
    alignment_method     = rows_align,
    stability_label      = rows_label
  )
}
