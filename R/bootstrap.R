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
#' @param parallel One of `"sequential"`, `"mirai"`, `"auto"`.
#'   `"sequential"` runs the replicate loop in-process using a single
#'   RNG stream (Phase 1 behavior). `"mirai"` / `"auto"` fan out per-rep
#'   tasks to a mirai daemon pool with deterministic per-rep seeds.
#' @param fast_path One of `"auto"`, `"off"`. `"auto"` uses the
#'   adapter's core()/update_core() hooks when available; `"off"`
#'   forces the refit path.
#' @param core_rank Positive integer cap on the core rank used by the
#'   fast path, or NULL to leave the full thin SVD untruncated.
#'   Truncation is the knob that converts the §30 thin-SVD cache into
#'   actual wall-clock wins; defaults to `50L`.
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
                           seed = NULL,
                           parallel = c("sequential", "mirai", "auto"),
                           fast_path = c("auto", "off"),
                           core_rank = 50L) {

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
  parallel     <- match.arg(parallel)
  fast_path    <- match.arg(fast_path)
  if (!is.null(core_rank)) {
    if (!is.numeric(core_rank) || length(core_rank) != 1L ||
        is.na(core_rank) || core_rank < 1L ||
        core_rank != as.integer(core_rank)) {
      stop("`core_rank` must be a positive integer scalar or NULL.",
           call. = FALSE)
    }
    core_rank <- as.integer(core_rank)
  }

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
  has_fast_path <- fast_path == "auto" &&
                   !is.null(adapter$core) &&
                   !is.null(adapter$update_core) &&
                   (geom_kind == "oneblock" ||
                    (geom_kind == "cross" && rel_kind == "covariance"))

  core_obj <- if (has_fast_path) {
    co <- adapter$core(original_fit, data)
    # Truncate to core_rank so the fast-path inner SVD is n x k instead
    # of n x full_rank. Without this cap the core-update work equals the
    # refit work and the §30 speedup does not materialize.
    .truncate_core(co, core_rank)
  } else {
    NULL
  }

  ## --- trim original loadings to match the core rank on the fast path -------
  # The alignment reference must have the same column count as the updated
  # rep fits, or match_components() will error on shape mismatch.
  if (has_fast_path) {
    k_trunc <- .core_effective_rank(core_obj)
    orig_loadings <- lapply(orig_loadings, function(V) {
      k_use <- min(k_trunc, ncol(V))
      V[, seq_len(k_use), drop = FALSE]
    })
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

  ## --- per-replicate closure --------------------------------------------------
  # Pure function: takes one replicate index and returns the full per-rep
  # artifact. Everything it reads is captured by value in the enclosing
  # bootstrap_fits() frame; everything it writes lives in the returned list.
  # This is the mirai-serializable fan-out unit.
  #
  # Bind all package helpers that the closure calls as locals so they live
  # in the function frame (= closure's enclosing env) and get serialized
  # with the task. Without this, mirai workers that do not have multifer
  # attached to their search path cannot resolve bare-name references.
  .fn_match_components   <- match_components
  .fn_truncate_fit_to_core <- .truncate_fit_to_core
  .k_trunc_local <- if (has_fast_path) .core_effective_rank(core_obj) else NA_integer_

  rep_fn <- function(b) {
    indices <- base::sample.int(n, n, replace = TRUE)

    rep_data <- if (geom_kind == "oneblock") {
      data[indices, , drop = FALSE]
    } else {
      list(
        X        = data$X[indices, , drop = FALSE],
        Y        = data$Y[indices, , drop = FALSE],
        relation = rel_kind
      )
    }

    rep_fit <- if (has_fast_path) {
      adapter$update_core(core_obj, indices = indices)
    } else {
      adapter$refit(original_fit, rep_data)
    }

    # Truncate the rep fit's loadings/scores to core_rank on the fast path
    # so the alignment loop sees the same column count as the reference.
    # This is a no-op on the slow path.
    if (has_fast_path) {
      rep_fit <- .fn_truncate_fit_to_core(rep_fit, .k_trunc_local)
    }

    rep_aligned_loadings <- stats::setNames(vector("list", length(domains)), domains)
    rep_aligned_scores   <- stats::setNames(vector("list", length(domains)), domains)

    for (d in domains) {
      Vref <- orig_loadings[[d]]
      Vb   <- adapter$loadings(rep_fit, d)
      Sb   <- adapter$scores(rep_fit, d)

      perm    <- .fn_match_components(Vref, Vb)
      Vb_perm <- Vb[, perm, drop = FALSE]
      Sb_perm <- Sb[, perm, drop = FALSE]

      if (method_align == "sign") {
        signs       <- sign(base::colSums(Vref * Vb_perm))
        signs[signs == 0L] <- 1L
        aligned_L   <- base::sweep(Vb_perm, 2L, signs, `*`)
        aligned_S   <- base::sweep(Sb_perm, 2L, signs, `*`)
      } else {
        M  <- base::crossprod(Vref, Vb_perm)
        sv <- base::svd(M)
        Q  <- sv$v %*% t(sv$u)
        aligned_L <- Vb_perm %*% Q
        aligned_S <- Sb_perm %*% Q
      }

      rep_aligned_loadings[[d]] <- aligned_L
      rep_aligned_scores[[d]]   <- aligned_S
    }

    list(
      fit               = rep_fit,
      aligned_loadings  = rep_aligned_loadings,
      aligned_scores    = rep_aligned_scores,
      resample_indices  = indices
    )
  }

  ## --- dispatch: sequential (legacy single-stream) or parallel (per-rep seed) ---

  reps <- if (parallel == "sequential") {
    # Single-stream RNG loop: exactly the Phase 1 behavior. Every reproducibility
    # test that hard-codes replicate indices relies on this stream ordering.
    out <- vector("list", R)
    for (b in seq_len(R)) out[[b]] <- rep_fn(b)
    out
  } else {
    # Parallel path: each task gets a deterministic per-rep seed so results
    # are reproducible within the same backend for a given master seed.
    task_seeds <- multifer_task_seeds(R, master_seed = seed)
    multifer_parallel_lapply(
      X       = seq_len(R),
      FUN     = rep_fn,
      seeds   = task_seeds,
      backend = parallel
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
      used_fast_path = has_fast_path,
      parallel       = parallel
    ),
    class = "multifer_bootstrap_artifact"
  )
}

# ---------------------------------------------------------------------------
# Core-object helpers. Handle both oneblock (U, d, V) and cross (Ux/dx/Vx,
# Uy/dy/Vy) key schemes. Truncation is the knob that turns §30's
# thin-SVD-cached fast path into actual wall-clock wins: without it, the
# inner SVD cost equals the refit cost.
# ---------------------------------------------------------------------------

.truncate_core <- function(core_obj, rank_cap) {
  if (is.null(core_obj) || is.null(rank_cap)) return(core_obj)
  if (!is.null(core_obj$U) && !is.null(core_obj$d) && !is.null(core_obj$V)) {
    # oneblock schema.
    k <- min(rank_cap, length(core_obj$d))
    if (k < length(core_obj$d)) {
      idx <- seq_len(k)
      core_obj$U <- core_obj$U[, idx, drop = FALSE]
      core_obj$d <- core_obj$d[idx]
      core_obj$V <- core_obj$V[, idx, drop = FALSE]
      core_obj$k <- k
    }
  } else if (!is.null(core_obj$Ux) && !is.null(core_obj$Uy)) {
    # cross covariance schema.
    kx <- min(rank_cap, length(core_obj$dx))
    ky <- min(rank_cap, length(core_obj$dy))
    if (kx < length(core_obj$dx)) {
      ix <- seq_len(kx)
      core_obj$Ux <- core_obj$Ux[, ix, drop = FALSE]
      core_obj$dx <- core_obj$dx[ix]
      core_obj$Vx <- core_obj$Vx[, ix, drop = FALSE]
    }
    if (ky < length(core_obj$dy)) {
      iy <- seq_len(ky)
      core_obj$Uy <- core_obj$Uy[, iy, drop = FALSE]
      core_obj$dy <- core_obj$dy[iy]
      core_obj$Vy <- core_obj$Vy[, iy, drop = FALSE]
    }
  }
  core_obj
}

.core_effective_rank <- function(core_obj) {
  if (is.null(core_obj)) return(NA_integer_)
  if (!is.null(core_obj$d)) return(length(core_obj$d))
  if (!is.null(core_obj$dx) && !is.null(core_obj$dy)) {
    return(min(length(core_obj$dx), length(core_obj$dy)))
  }
  NA_integer_
}

.truncate_fit_to_core <- function(fit, k) {
  if (is.na(k) || is.null(fit)) return(fit)
  # prcomp-shaped (rotation / x / sdev)
  if (!is.null(fit$rotation) && !is.null(fit$sdev)) {
    k_use <- min(k, length(fit$sdev))
    fit$sdev     <- fit$sdev[seq_len(k_use)]
    fit$rotation <- fit$rotation[, seq_len(k_use), drop = FALSE]
    if (!is.null(fit$x)) fit$x <- fit$x[, seq_len(k_use), drop = FALSE]
    return(fit)
  }
  # svd-shaped (u / d / v)
  if (!is.null(fit$u) && !is.null(fit$d) && !is.null(fit$v)) {
    k_use <- min(k, length(fit$d))
    fit$u <- fit$u[, seq_len(k_use), drop = FALSE]
    fit$d <- fit$d[seq_len(k_use)]
    fit$v <- fit$v[, seq_len(k_use), drop = FALSE]
    return(fit)
  }
  # cross-shaped (Wx, Wy, d, Tx, Ty)
  if (!is.null(fit$Wx) && !is.null(fit$Wy)) {
    k_use <- min(k, length(fit$d))
    fit$d  <- fit$d[seq_len(k_use)]
    fit$Wx <- fit$Wx[, seq_len(k_use), drop = FALSE]
    fit$Wy <- fit$Wy[, seq_len(k_use), drop = FALSE]
    if (!is.null(fit$Tx)) fit$Tx <- fit$Tx[, seq_len(k_use), drop = FALSE]
    if (!is.null(fit$Ty)) fit$Ty <- fit$Ty[, seq_len(k_use), drop = FALSE]
    return(fit)
  }
  fit
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
