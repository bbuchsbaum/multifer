#' multifer adapter for pls::plsr
#'
#' A Tier-1 refit-only adapter wrapping [pls::plsr()] for the
#' `(cross, predictive)` family. This adapter is intentionally
#' predictive rather than covariance-based: component significance is
#' understood as held-out incremental predictive gain, not as a
#' singular-value tail ratio.
#'
#' @param adapter_id Character, default `"plsr_refit"`.
#' @param adapter_version Character, default `"0.0.1"`.
#' @param ncomp Optional integer cap on fitted components. When `NULL`,
#'   the refit hook chooses `min(nrow(X) - 1, ncol(X), ncol(Y), 50L)`.
#' @param method Character PLS algorithm passed to [pls::plsr()].
#'   Default `"simpls"`.
#'
#' @return A `multifer_adapter`.
#' @export
adapter_plsr_refit <- function(adapter_id = "plsr_refit",
                               adapter_version = "0.0.1",
                               ncomp = NULL,
                               method = "simpls") {

  require_pls <- function() {
    if (!requireNamespace("pls", quietly = TRUE)) {
      stop("pls is required by adapter_plsr_refit but is not installed.",
           call. = FALSE)
    }
  }

  fit_plsr <- function(new_data) {
    require_pls()
    .validate_plsr_data(new_data)
    k <- if (is.null(ncomp)) {
      max(1L, min(nrow(new_data$X) - 1L,
                  ncol(new_data$X),
                  ncol(new_data$Y),
                  50L))
    } else {
      as.integer(ncomp)
    }
    pls::plsr(
      Y ~ X,
      data = list(X = new_data$X, Y = new_data$Y),
      ncomp = k,
      validation = "none",
      method = method
    )
  }

  infer_adapter(
    adapter_id      = adapter_id,
    adapter_version = adapter_version,
    shape_kinds     = "cross",
    capabilities    = capability_matrix(
      list(
        geometry = "cross",
        relation = "predictive",
        targets  = c("component_significance", "variable_stability",
                     "score_stability", "subspace_stability")
      )
    ),

    roots = function(x, ...) {
      require_pls()
      as.numeric(x$Xvar)
    },

    scores = function(x, domain = c("X", "Y"), ...) {
      require_pls()
      domain <- match.arg(domain)
      if (domain == "X") {
        unclass(x$scores)
      } else {
        unclass(x$Yscores)
      }
    },

    loadings = function(x, domain = c("X", "Y"), ...) {
      require_pls()
      domain <- match.arg(domain)
      if (domain == "X") {
        unclass(x$loadings)
      } else {
        unclass(x$Yloadings)
      }
    },

    truncate = function(x, k, ...) {
      require_pls()
      .truncate_plsr_fit(x, k)
    },

    residualize = function(x, k, data, ...) {
      require_pls()
      .validate_plsr_data(data)

      kk <- min(max(as.integer(k), 0L), x$ncomp)
      if (kk <= 0L) {
        return(list(
          X = sweep(data$X, 2L, x$Xmeans, "-"),
          Y = sweep(data$Y, 2L, x$Ymeans, "-")
        ))
      }

      T_scores <- unclass(x$scores)[, seq_len(kk), drop = FALSE]
      P_load   <- unclass(x$loadings)[, seq_len(kk), drop = FALSE]
      Xc <- sweep(data$X, 2L, x$Xmeans, "-")
      Yc <- sweep(data$Y, 2L, x$Ymeans, "-")

      Yhat <- adapter_plsr_predict(x, data, k = kk) -
        matrix(x$Ymeans, nrow = nrow(data$Y), ncol = length(x$Ymeans), byrow = TRUE)

      list(
        X = Xc - T_scores %*% t(P_load),
        Y = Yc - Yhat
      )
    },

    refit = function(x, new_data, ...) {
      fit_plsr(new_data)
    },

    predict_response = function(x, new_data, k = NULL, ...) {
      require_pls()
      adapter_plsr_predict(x, new_data, k = k)
    },

    null_action = function(x, data, ...) {
      perm <- sample.int(nrow(data$Y))
      list(X = data$X, Y = data$Y[perm, , drop = FALSE])
    },

    component_stat = function(x, data, k, split = NULL, ...) {
      require_pls()
      if (is.null(split) || is.null(split$train) || is.null(split$test)) {
        stop(
          "`component_stat` for plsr_refit requires a split with train/test indices.",
          call. = FALSE
        )
      }
      .plsr_split_incremental_gain(
        data = data,
        split = split,
        k = k,
        fit_fun = fit_plsr
      )
    },

    validity_level       = "conditional",
    declared_assumptions = c("paired_rows", "continuous_response"),
    checked_assumptions  = .plsr_checks()
  )
}

adapter_plsr_predict <- function(x, new_data, k = NULL) {
  .validate_plsr_data(new_data)

  if (is.null(k)) {
    k <- x$ncomp
  }
  kk <- min(max(as.integer(k), 0L), x$ncomp)
  if (kk <= 0L) {
    return(matrix(
      x$Ymeans,
      nrow = nrow(new_data$X),
      ncol = length(x$Ymeans),
      byrow = TRUE
    ))
  }

  pred <- stats::predict(
    x,
    newdata = list(X = new_data$X),
    ncomp = kk
  )
  pred <- drop(pred[, , 1L, drop = FALSE])
  if (is.vector(pred)) {
    pred <- matrix(pred, ncol = length(x$Ymeans))
  }
  pred
}

.truncate_plsr_fit <- function(x, k) {
  kk <- min(max(as.integer(k), 0L), x$ncomp)
  out <- x
  out$ncomp <- kk

  if (kk <= 0L) {
    out$coefficients <- array(
      0,
      dim = c(nrow(x$coefficients), ncol(x$coefficients), 0L),
      dimnames = list(dimnames(x$coefficients)[[1L]],
                      dimnames(x$coefficients)[[2L]],
                      character(0))
    )
    out$scores <- x$scores[, 0L, drop = FALSE]
    out$loadings <- x$loadings[, 0L, drop = FALSE]
    out$Yscores <- x$Yscores[, 0L, drop = FALSE]
    out$Yloadings <- x$Yloadings[, 0L, drop = FALSE]
    out$projection <- x$projection[, 0L, drop = FALSE]
    out$fitted.values <- x$fitted.values[, , 0L, drop = FALSE]
    out$residuals <- x$residuals[, , 0L, drop = FALSE]
    out$Xvar <- numeric(0)
    return(out)
  }

  out$coefficients <- x$coefficients[, , seq_len(kk), drop = FALSE]
  out$scores <- x$scores[, seq_len(kk), drop = FALSE]
  out$loadings <- x$loadings[, seq_len(kk), drop = FALSE]
  out$Yscores <- x$Yscores[, seq_len(kk), drop = FALSE]
  out$Yloadings <- x$Yloadings[, seq_len(kk), drop = FALSE]
  out$projection <- x$projection[, seq_len(kk), drop = FALSE]
  out$fitted.values <- x$fitted.values[, , seq_len(kk), drop = FALSE]
  out$residuals <- x$residuals[, , seq_len(kk), drop = FALSE]
  out$Xvar <- x$Xvar[seq_len(kk)]
  out
}

.validate_plsr_data <- function(data) {
  if (!is.list(data) || is.null(data$X) || is.null(data$Y)) {
    stop("plsr_refit expects data as a list with X and Y.", call. = FALSE)
  }
  if (!is.matrix(data$X) || !is.numeric(data$X)) {
    stop("plsr_refit expects X to be a numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(data$Y) || !is.numeric(data$Y)) {
    stop("plsr_refit expects Y to be a numeric matrix.", call. = FALSE)
  }
  if (nrow(data$X) != nrow(data$Y)) {
    stop("plsr_refit requires X and Y to have matching row counts.", call. = FALSE)
  }
  if (any(!is.finite(data$X)) || any(!is.finite(data$Y))) {
    stop("plsr_refit requires finite X and Y values.", call. = FALSE)
  }
  invisible(data)
}

.plsr_split_incremental_gain <- function(data, split, k, fit_fun) {
  train <- .predictive_subset_data(data, split$train)
  test  <- .predictive_subset_data(data, split$test)
  fit   <- fit_fun(train)

  y_mean <- matrix(fit$Ymeans, nrow = nrow(test$Y), ncol = ncol(test$Y), byrow = TRUE)
  tss <- sum((test$Y - y_mean)^2)
  if (tss <= .Machine$double.eps) {
    return(0)
  }

  pred_k <- adapter_plsr_predict(fit, test, k = k)
  pred_prev <- adapter_plsr_predict(fit, test, k = k - 1L)

  gain_k <- 1 - sum((test$Y - pred_k)^2) / tss
  gain_prev <- 1 - sum((test$Y - pred_prev)^2) / tss
  max(0, gain_k - gain_prev)
}

.plsr_checks <- function() {
  list(
    list(
      name = "predictive_x_numeric_matrix",
      detail = "`data$X` must be a numeric matrix for plsr_refit.",
      check = function(data, ...) is.list(data) &&
        is.matrix(data$X) && is.numeric(data$X)
    ),
    list(
      name = "predictive_y_numeric_matrix",
      detail = "`data$Y` must be a numeric matrix for plsr_refit.",
      check = function(data, ...) is.list(data) &&
        is.matrix(data$Y) && is.numeric(data$Y)
    ),
    list(
      name = "predictive_matching_row_count",
      detail = "`data$X` and `data$Y` must have matching row counts.",
      check = function(data, ...) is.list(data) &&
        !is.null(data$X) && !is.null(data$Y) &&
        nrow(data$X) == nrow(data$Y)
    ),
    list(
      name = "predictive_min_rows",
      detail = "plsr_refit requires at least 2 paired rows.",
      check = function(data, ...) is.list(data) &&
        is.matrix(data$X) && is.matrix(data$Y) &&
        nrow(data$X) >= 2L && nrow(data$Y) >= 2L
    ),
    list(
      name = "predictive_blocks_finite",
      detail = "`data$X` and `data$Y` must be finite matrices.",
      check = function(data, ...) is.list(data) &&
        is.matrix(data$X) && is.numeric(data$X) &&
        is.matrix(data$Y) && is.numeric(data$Y) &&
        all(is.finite(data$X)) && all(is.finite(data$Y))
    )
  )
}

#' Register the pls::plsr adapter when pls is installed
#'
#' Called from the package `.onLoad` hook. Skips registration silently if
#' `pls` is not installed.
#'
#' @keywords internal
register_plsr_refit_adapter <- function() {
  if (!requireNamespace("pls", quietly = TRUE)) {
    return(invisible(NULL))
  }
  register_infer_adapter(
    adapter_id = "plsr_refit",
    adapter    = adapter_plsr_refit(),
    overwrite  = TRUE
  )
  invisible(NULL)
}
