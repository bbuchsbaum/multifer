#' Row bootstrap perturbation loop
#'
#' Runs R bootstrap replicates of an adapter's fit on resampled data,
#' aligning loadings to the original fit on each replicate so that
#' downstream stability consumers (variable / score / subspace) see
#' per-replicate aligned fits without having to re-do alignment.
#'
#' For (oneblock, variance), uses ordinary row bootstrap with
#' replacement. For (cross, covariance | correlation), uses PAIRED row
#' bootstrap -- a single index vector is sampled and applied to BOTH
#' X and Y to preserve the observation pairing that the design
#' assumes.
#'
#' @section Score alignment for procrustes:
#' When \code{method_align = "procrustes"}, after the column permutation is
#' applied to both loadings and scores, the loadings are further rotated by
#' the orthogonal matrix Q from the Procrustes solution
#' \code{M = t(Vref) \%*\% Vb_perm}, \code{svd(M) = U D V^t},
#' \code{Q = V \%*\% t(U)}, and the SAME orthogonal matrix Q is applied to
#' the scores (right-multiplied). This keeps loadings and scores consistent:
#' \code{scores \%*\% Q} represents the same projections as
#' \code{loadings \%*\% Q}.
#'
#' @param recipe A compiled \code{multifer_infer_recipe}.
#' @param adapter The \code{multifer_adapter} object (from
#'   \code{recipe$adapter} or looked up from the registry).
#' @param data For oneblock, a numeric matrix. For cross, a list with
#'   elements \code{X} and \code{Y}.
#' @param original_fit The output of \code{adapter$refit} applied to the
#'   original (unperturbed) data. Used as the alignment reference.
#' @param units A \code{multifer_units} table produced by
#'   \code{\link{form_units}()}.
#' @param R Integer >= 1, number of bootstrap replicates.
#' @param method_align Character, one of \code{"sign"} or
#'   \code{"procrustes"} (default \code{"sign"}). Passed to
#'   \code{\link{align_loadings}()}. See section on score alignment for
#'   procrustes details.
#' @param seed Integer or NULL for reproducibility. When non-NULL, the
#'   current RNG state is saved and restored on exit so the caller's RNG
#'   stream is unchanged.
#'
#' @return An object of class \code{multifer_bootstrap_artifact}, a list
#'   with:
#'   \describe{
#'     \item{\code{reps}}{List of length R. Each element is a list with:
#'       \code{fit} (the adapter fit for this replicate),
#'       \code{aligned_loadings} (named list by domain, each a p x r
#'       matrix aligned to original_fit),
#'       \code{aligned_scores} (named list by domain, each an n x r matrix
#'       of scores from the replicate fit with the same column permutation
#'       and sign/rotation applied as was applied to loadings),
#'       \code{resample_indices} (integer vector of length n).}
#'     \item{\code{R}}{Integer.}
#'     \item{\code{method_align}}{Character.}
#'     \item{\code{domains}}{Character vector of domain labels.}
#'     \item{\code{seed}}{Integer or NULL.}
#'   }
#'
#' @export
bootstrap_fits <- function(recipe,
                           adapter,
                           data,
                           original_fit,
                           units,
                           R = 500L,
                           method_align = c("sign", "procrustes"),
                           seed = NULL) {

  ## --- input validation -------------------------------------------------------

  if (!inherits(recipe, "multifer_infer_recipe")) {
    stop("`recipe` must be a multifer_infer_recipe.", call. = FALSE)
  }
  if (!inherits(adapter, "multifer_adapter")) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }
  if (!is.numeric(R) || length(R) != 1L || is.na(R) ||
      R != as.integer(R) || R < 1L) {
    stop("`R` must be a positive integer.", call. = FALSE)
  }
  R <- as.integer(R)

  method_align <- match.arg(method_align)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be a single integer value or NULL.", call. = FALSE)
    }
    seed <- as.integer(seed)
  }

  ## --- determine geometry and domains -----------------------------------------

  geom_kind <- recipe$shape$geometry$kind
  rel_kind  <- recipe$shape$relation$kind

  if (geom_kind == "oneblock") {
    if (!is.matrix(data)) {
      stop(
        "For oneblock geometry, `data` must be a numeric matrix.",
        call. = FALSE
      )
    }
    domains <- "X"
    n       <- nrow(data)
  } else if (geom_kind == "cross") {
    if (!is.list(data) || is.null(data$X) || is.null(data$Y)) {
      stop(
        "For cross geometry, `data` must be a list with elements `X` and `Y`.",
        call. = FALSE
      )
    }
    domains <- c("X", "Y")
    n       <- nrow(data$X)
    if (nrow(data$Y) != n) {
      stop(
        "For cross geometry, `data$X` and `data$Y` must have the same number of rows.",
        call. = FALSE
      )
    }
  } else {
    stop(
      sprintf(
        "bootstrap_fits only supports 'oneblock' and 'cross' geometries. Got: '%s'.",
        geom_kind
      ),
      call. = FALSE
    )
  }

  ## --- extract original loadings per domain for alignment reference -----------

  orig_loadings <- stats::setNames(
    lapply(domains, function(d) adapter$loadings(original_fit, d)),
    domains
  )

  ## --- Phase 1.5 fast-path detection -----------------------------------------
  # If the adapter exposes core() and update_core(), precompute the core once
  # and use update_core per replicate (O(n*k^2 + k_x*k_y) per rep) instead of
  # refit (O(n*p^2) or O(p_x*p_y*min(n,p,q))). The cross correlation branch
  # does not have a clean fast path (per-replicate QR rescaling) and falls
  # back to refit by gating on rel_kind.
  has_fast_path <- !is.null(adapter$core) &&
                   !is.null(adapter$update_core) &&
                   (geom_kind == "oneblock" ||
                    (geom_kind == "cross" && rel_kind == "covariance"))

  core_obj <- if (has_fast_path) {
    adapter$core(original_fit, data)
  } else {
    NULL
  }

  ## --- RNG state: save and restore if seed supplied ---------------------------

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      get(".Random.seed", envir = .GlobalEnv)
    else
      NULL
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  ## --- bootstrap loop ---------------------------------------------------------

  reps <- vector("list", R)

  for (b in seq_len(R)) {

    ## -- resample ------------------------------------------------------------

    indices <- base::sample.int(n, n, replace = TRUE)

    if (geom_kind == "oneblock") {
      rep_data <- data[indices, , drop = FALSE]
    } else {
      # Paired row bootstrap: SAME indices applied to both X and Y.
      rep_data <- list(
        X        = data$X[indices, , drop = FALSE],
        Y        = data$Y[indices, , drop = FALSE],
        relation = rel_kind
      )
    }

    ## -- refit (slow path) or core-update (fast path) ------------------------

    rep_fit <- if (has_fast_path) {
      adapter$update_core(core_obj, indices = indices)
    } else {
      adapter$refit(original_fit, rep_data)
    }

    ## -- align loadings and scores per domain --------------------------------

    rep_aligned_loadings <- vector("list", length(domains))
    names(rep_aligned_loadings) <- domains

    rep_aligned_scores <- vector("list", length(domains))
    names(rep_aligned_scores) <- domains

    for (d in domains) {
      Vref <- orig_loadings[[d]]
      Vb   <- adapter$loadings(rep_fit, d)
      Sb   <- adapter$scores(rep_fit, d)

      ## Step 1: column permutation
      perm    <- match_components(Vref, Vb)
      Vb_perm <- Vb[, perm, drop = FALSE]
      Sb_perm <- Sb[, perm, drop = FALSE]

      if (method_align == "sign") {
        ## Step 2a: sign flips
        ## signs[j] = sign of diag(t(Vref) %*% Vb_perm)[j]
        ## Use colSums(Vref * Vb_perm) as the cheap diagonal equivalent.
        signs       <- sign(base::colSums(Vref * Vb_perm))
        signs[signs == 0L] <- 1L

        aligned_L   <- base::sweep(Vb_perm, 2L, signs, `*`)
        aligned_S   <- base::sweep(Sb_perm, 2L, signs, `*`)

      } else {
        ## Step 2b: Procrustes orthogonal rotation
        ## M = t(Vref) %*% Vb_perm; svd(M) = U D V^t; Q = V %*% t(U)
        M         <- base::crossprod(Vref, Vb_perm)
        sv        <- svd(M)
        Q         <- sv$v %*% t(sv$u)

        aligned_L <- Vb_perm %*% Q
        aligned_S <- Sb_perm %*% Q
      }

      rep_aligned_loadings[[d]] <- aligned_L
      rep_aligned_scores[[d]]   <- aligned_S
    }

    reps[[b]] <- list(
      fit               = rep_fit,
      aligned_loadings  = rep_aligned_loadings,
      aligned_scores    = rep_aligned_scores,
      resample_indices  = indices
    )
  }

  ## --- assemble artifact ------------------------------------------------------

  structure(
    list(
      reps           = reps,
      R              = R,
      method_align   = method_align,
      domains        = domains,
      seed           = seed,
      used_fast_path = has_fast_path
    ),
    class = "multifer_bootstrap_artifact"
  )
}


#' @export
print.multifer_bootstrap_artifact <- function(x, ...) {
  cat("<multifer_bootstrap_artifact>\n")
  cat("  R:            ", x$R,                           "\n", sep = "")
  cat("  domains:      ", paste(x$domains, collapse = ", "), "\n", sep = "")
  cat("  method_align: ", x$method_align,                "\n", sep = "")
  if (!is.null(x$seed)) {
    cat("  seed:         ", x$seed, "\n", sep = "")
  }
  invisible(x)
}
