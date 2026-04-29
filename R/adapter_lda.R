#' multifer adapter for MASS::lda
#'
#' A Tier-1 refit adapter for the `(geneig, generalized_eigen)` family
#' using `MASS::lda()` as the fitting backend. The adapter stores the
#' fitted `lda` object together with the corresponding generalized-eigen
#' operator pair `(A, B)` where `A` is the between-class scatter and `B`
#' is the within-class scatter. Null resampling permutes class labels
#' only; rows are left in place.
#'
#' @param adapter_id Character, default `"lda_refit"`.
#' @param adapter_version Character, default `"0.0.1"`.
#'
#' @return A `multifer_adapter`.
#' @export
adapter_lda_refit <- function(adapter_id = "lda_refit",
                              adapter_version = "0.0.1") {
  require_mass <- function() {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("MASS is required by adapter_lda_refit but is not installed.",
           call. = FALSE)
    }
  }

  residualize_geneig <- function(x, k, data, ...) {
    payload <- .geneig_normalize_state(.as_geneig_payload(data))
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 0) {
      stop("`k` must be a single non-negative numeric value.", call. = FALSE)
    }
    if (k < 1L) {
      return(list(
        A = payload$A,
        B = payload$B,
        metric = payload$metric,
        state = payload$state
      ))
    }

    zero_tol <- max(1, sum(.geneig_roots(payload))) * .Machine$double.eps
    out <- payload
    for (step in seq_len(as.integer(k))) {
      out <- .geneig_deflate_state(out, zero_tol = zero_tol)
    }
    list(A = out$A, B = out$B, metric = out$metric, state = out$state)
  }

  infer_adapter(
    adapter_id      = adapter_id,
    adapter_version = adapter_version,
    shape_kinds     = "geneig",
    capabilities    = capability_matrix(
      list(
        geometry = "geneig",
        relation = "generalized_eigen",
        targets  = "component_significance"
      )
    ),

    roots = function(x, ...) x$roots,

    scores = function(x, domain = NULL, ...) x$scores,

    loadings = function(x, domain = NULL, ...) x$fit$scaling,

    truncate = function(x, k, ...) {
      k <- min(as.integer(k), length(x$roots))
      out <- x
      out$roots <- x$roots[seq_len(k)]
      out$scores <- x$scores[, seq_len(k), drop = FALSE]
      out$fit$scaling <- x$fit$scaling[, seq_len(k), drop = FALSE]
      out$fit$svd <- x$fit$svd[seq_len(k)]
      out
    },

    residualize = residualize_geneig,

    refit = function(x, new_data, ...) {
      require_mass()
      payload <- .as_lda_data(new_data)
      .fit_lda_refit(payload$X, payload$y)
    },

    null_action = function(x, data, ...) {
      payload <- .as_geneig_payload(data)
      if (is.null(payload$state$X) || is.null(payload$state$y)) {
        stop(
          "lda_refit null_action requires payload state with `X` and `y` so labels can be permuted.",
          call. = FALSE
        )
      }
      y_perm <- sample(payload$state$y)
      op <- .lda_geneig_operator(payload$state$X, y_perm)
      list(
        A = op$A,
        B = op$B,
        metric = op$metric,
        state = list(X = payload$state$X, y = y_perm)
      )
    },

    component_stat = function(x, data, k, ...) {
      payload <- .as_geneig_payload(data)
      roots <- .geneig_roots(.geneig_normalize_state(payload))
      if (k > length(roots)) {
        return(0)
      }
      .geneig_tail_ratio(roots[seq.int(k, length(roots))], zero_tol = .Machine$double.eps)
    },

    validity_level       = "conditional",
    declared_assumptions = c("labels_exchangeable"),
    checked_assumptions  = .lda_checks(),
    geneig_deflation     = "b_metric"
  )
}

.fit_lda_refit <- function(X, y) {
  fit <- MASS::lda(x = X, grouping = y)
  scores <- stats::predict(fit, X)$x
  op <- .lda_geneig_operator(X, y)
  roots <- .geneig_roots(.geneig_normalize_state(list(
    A = op$A,
    B = op$B,
    metric = op$metric
  )))

  structure(
    list(
      fit = fit,
      operator = op,
      roots = roots,
      scores = scores,
      X = X,
      y = y
    ),
    class = "multifer_lda_refit_fit"
  )
}

.lda_geneig_operator <- function(X, y, ridge = 1e-8) {
  payload <- .as_lda_data(list(X = X, y = y))
  X <- payload$X
  y <- payload$y
  grand_mean <- colMeans(X)
  classes <- levels(y)
  p <- ncol(X)
  within <- matrix(0, nrow = p, ncol = p)
  between <- matrix(0, nrow = p, ncol = p)

  for (cls in classes) {
    Xg <- X[y == cls, , drop = FALSE]
    mu_g <- colMeans(Xg)
    centered <- sweep(Xg, 2L, mu_g, "-")
    within <- within + crossprod(centered)
    mean_shift <- matrix(mu_g - grand_mean, ncol = 1L)
    between <- between + nrow(Xg) * (mean_shift %*% t(mean_shift))
  }

  within <- (within + t(within)) / 2
  between <- (between + t(between)) / 2
  within <- within + diag(ridge, nrow = p)

  geneig_operator(
    A = between,
    B = within,
    metric = "within_class"
  )
}

.as_lda_data <- function(data) {
  if (!is.list(data) || is.null(data$X) || is.null(data$y)) {
    stop("LDA data must be a list with components `X` and `y`.", call. = FALSE)
  }
  if (!is.matrix(data$X) || !is.numeric(data$X)) {
    stop("LDA data component `X` must be a numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(data$X))) {
    stop("LDA data component `X` must not contain NA, NaN, or Inf.", call. = FALSE)
  }
  y <- if (is.factor(data$y)) data$y else factor(data$y)
  if (length(y) != nrow(data$X)) {
    stop("LDA data component `y` must have one value per row of `X`.", call. = FALSE)
  }
  if (anyNA(y)) {
    stop("LDA grouping labels must not contain NA.", call. = FALSE)
  }
  list(X = data$X, y = y)
}

.as_geneig_payload <- function(data) {
  if (is_geneig_operator(data)) {
    return(list(A = data$A, B = data$B, metric = data$metric, state = NULL))
  }
  if (is.list(data) && !is.null(data$A) && !is.null(data$B)) {
    return(data)
  }
  if (is.list(data) && !is.null(data$X) && !is.null(data$y)) {
    payload <- .as_lda_data(data)
    op <- .lda_geneig_operator(payload$X, payload$y)
    return(list(
      A = op$A,
      B = op$B,
      metric = op$metric,
      state = list(X = payload$X, y = payload$y)
    ))
  }
  stop("Geneig payload must be a geneig operator, an operator list, or LDA data with `X` and `y`.",
       call. = FALSE)
}

.lda_checks <- function() {
  list(
    list(
      name = "lda_x_numeric_matrix",
      detail = "LDA requires `X` to be a numeric matrix.",
      check = function(data, ...) {
        is.list(data) && is.matrix(data$X) && is.numeric(data$X)
      }
    ),
    list(
      name = "lda_x_finite",
      detail = "LDA requires finite X values.",
      check = function(data, ...) {
        if (!is.list(data) || !is.matrix(data$X) || !is.numeric(data$X)) {
          return(TRUE)
        }
        all(is.finite(data$X))
      }
    ),
    list(
      name = "lda_paired_rows",
      detail = "LDA requires one grouping label per row of X with no missing labels.",
      check = function(data, ...) {
        if (!is.list(data) || is.null(data$X) || is.null(data$y) || !is.matrix(data$X)) {
          return(TRUE)
        }
        length(data$y) == nrow(data$X) && !anyNA(data$y)
      }
    ),
    list(
      name = "lda_grouping_factor",
      detail = "LDA requires a factor response with at least two levels.",
      check = function(data, ...) {
        if (!is.list(data) || is.null(data$y)) {
          return(TRUE)
        }
        y <- if (is.factor(data$y)) data$y else factor(data$y)
        nlevels(y) >= 2L
      }
    ),
    list(
      name = "lda_within_class_sample_size",
      detail = "Every class must have at least two observations for within-class scatter.",
      check = function(data, ...) {
        if (!is.list(data) || is.null(data$y)) {
          return(TRUE)
        }
        y <- if (is.factor(data$y)) data$y else factor(data$y)
        if (anyNA(y)) {
          return(TRUE)
        }
        tab <- table(y)
        all(tab >= 2L)
      }
    )
  )
}

#' Register the MASS::lda adapter when MASS is installed
#'
#' Called from the package `.onLoad` hook. Skips registration silently
#' if `MASS` is not installed.
#'
#' @keywords internal
register_lda_refit_adapter <- function() {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    return(invisible(NULL))
  }
  register_infer_adapter(
    adapter_id = "lda_refit",
    adapter = adapter_lda_refit(),
    overwrite = TRUE
  )
  invisible(NULL)
}
