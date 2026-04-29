# Internal executable plan compiler.
#
# The public adapter contract is declarative capabilities plus hooks. Runtime
# engines should consume only resolved callbacks. This first plan compiler is
# intentionally narrow: it covers the adapter-driven ladder path and exists so
# relation-specific hook details are absorbed before the ladder runs.

compile_infer_plan <- function(recipe,
                               adapter = NULL,
                               validate_data = NULL,
                               labels = NULL) {
  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  adapter <- adapter %||% recipe$adapter
  if (!is_infer_adapter(adapter)) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }
  if (!is.null(validate_data) && !is.function(validate_data)) {
    stop("`validate_data` must be NULL or a function.", call. = FALSE)
  }

  rel_kind <- recipe$shape$relation$kind
  component_stat <- if (identical(rel_kind, "predictive")) {
    function(fit, data, k) adapter$component_stat(fit, data, k, split = NULL)
  } else {
    function(fit, data, k) adapter$component_stat(fit, data, k)
  }

  validate_cb <- if (is.null(validate_data)) {
    function(data) invisible(data)
  } else {
    validate_data
  }

  structure(
    list(
      recipe = recipe,
      adapter = adapter,
      fit = function(previous_fit, data) adapter$refit(previous_fit, data),
      roots = function(fit) adapter$roots(fit),
      component_stat = component_stat,
      null_action = function(fit, data) adapter$null_action(fit, data),
      residualize = function(fit, k, data) adapter$residualize(fit, k, data),
      validate_data = validate_cb,
      labels = .infer_plan_labels(recipe, labels)
    ),
    class = "multifer_infer_plan"
  )
}

is_infer_plan <- function(x) inherits(x, "multifer_infer_plan")

compile_infer_execution_plan <- function(recipe,
                                         data,
                                         adapter = NULL,
                                         model = NULL) {
  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  adapter <- adapter %||% recipe$adapter
  if (!is_infer_adapter(adapter)) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }

  problem <- .recipe_problem(recipe)
  geom_kind <- problem$shape$geometry$kind
  rel_kind <- problem$shape$relation$kind
  needs_component <- "component_significance" %in% problem$targets
  engine_original_fit <- NULL

  original_fit <- function(engine_out = NULL) {
    if (!is.null(model)) {
      return(model)
    }
    if (!is.null(engine_out$original_fit)) {
      engine_original_fit <<- engine_out$original_fit
      return(engine_original_fit)
    }
    if (!is.null(engine_original_fit)) {
      return(engine_original_fit)
    }
    engine_original_fit <<- .infer_original_fit(
      adapter_obj = adapter,
      data = data,
      geom_kind = geom_kind,
      rel_kind = rel_kind
    )
    engine_original_fit
  }

  run_engine <- .compile_engine_runner(
    recipe = recipe,
    adapter = adapter,
    data = data,
    model = model,
    original_fit = original_fit,
    needs_component = needs_component,
    geom_kind = geom_kind,
    rel_kind = rel_kind
  )

  structure(
    list(
      recipe = recipe,
      adapter = adapter,
      problem = problem,
      targets = problem$targets,
      run_engine = run_engine,
      original_fit = original_fit,
      bootstrap_supported = !identical(geom_kind, "geneig"),
      bootstrap_error = if (identical(geom_kind, "geneig")) {
        paste(
          "Bootstrap stability for geneig wrappers is not wired through infer() yet.",
          "Request `targets = \"component_significance\"` for now."
        )
      } else {
        NULL
      }
    ),
    class = "multifer_execution_plan"
  )
}

is_infer_execution_plan <- function(x) inherits(x, "multifer_execution_plan")

.infer_plan_labels <- function(recipe, labels = NULL) {
  labels <- labels %||% list()
  geom_kind <- recipe$shape$geometry$kind
  rel_kind <- recipe$shape$relation$kind

  default_statistic <- if (identical(rel_kind, "predictive")) {
    sprintf("adapter predictive component_stat on %s data", geom_kind)
  } else {
    sprintf("adapter component_stat on %s data", geom_kind)
  }
  default_estimand <- if (identical(rel_kind, "predictive")) {
    "adapter predictive gain"
  } else {
    "adapter latent roots"
  }

  list(
    statistic = labels$statistic %||% default_statistic,
    null = labels$null %||% sprintf("adapter null_action on %s data", geom_kind),
    estimand = labels$estimand %||% default_estimand
  )
}

.compile_engine_runner <- function(recipe,
                                   adapter,
                                   data,
                                   model,
                                   original_fit,
                                   needs_component,
                                   geom_kind,
                                   rel_kind) {
  if (!needs_component) {
    return(function(B,
                    B_total,
                    batch_size,
                    alpha,
                    seed,
                    auto_subspace,
                    tie_threshold) {
      .zero_component_engine_out(
        adapter = adapter,
        fit = original_fit(),
        auto_subspace = auto_subspace,
        tie_threshold = tie_threshold
      )
    })
  }

  if (identical(geom_kind, "oneblock")) {
    if (identical(adapter$component_execution, "adapter")) {
      return(function(B,
                      B_total,
                      batch_size,
                      alpha,
                      seed,
                      auto_subspace,
                      tie_threshold) {
        .validate_oneblock_data(data)
        fit <- original_fit()
        run_adapter_ladder(
          recipe = recipe,
          adapter = adapter,
          data = data,
          B = B,
          B_total = B_total,
          batch_size = batch_size,
          alpha = alpha,
          seed = seed,
          auto_subspace = auto_subspace,
          tie_threshold = tie_threshold,
          original_fit = fit,
          validate_data = .validate_oneblock_data
        )
      })
    }

    return(function(B,
                    B_total,
                    batch_size,
                    alpha,
                    seed,
                    auto_subspace,
                    tie_threshold) {
      .validate_oneblock_data(data)
      run_oneblock_ladder(
        recipe = recipe,
        X = data,
        B = B,
        B_total = B_total,
        batch_size = batch_size,
        alpha = alpha,
        seed = seed,
        auto_subspace = auto_subspace,
        tie_threshold = tie_threshold
      )
    })
  }

  if (identical(geom_kind, "cross")) {
    return(function(B,
                    B_total,
                    batch_size,
                    alpha,
                    seed,
                    auto_subspace,
                    tie_threshold) {
      if (!is.list(data) || is.null(data$X) || is.null(data$Y)) {
        stop("For cross geometry, `data` must be a list with X and Y.",
             call. = FALSE)
      }
      if (identical(rel_kind, "predictive")) {
        return(run_predictive_ladder(
          recipe = recipe,
          X = data$X,
          Y = data$Y,
          B = B,
          B_total = B_total,
          batch_size = batch_size,
          alpha = alpha,
          seed = seed,
          auto_subspace = auto_subspace,
          tie_threshold = tie_threshold
        ))
      }
      run_cross_ladder(
        recipe = recipe,
        X = data$X,
        Y = data$Y,
        B = B,
        B_total = B_total,
        batch_size = batch_size,
        alpha = alpha,
        seed = seed,
        auto_subspace = auto_subspace,
        tie_threshold = tie_threshold
      )
    })
  }

  if (identical(geom_kind, "geneig")) {
    return(function(B,
                    B_total,
                    batch_size,
                    alpha,
                    seed,
                    auto_subspace,
                    tie_threshold) {
      payload <- .as_geneig_payload(data)
      operator <- geneig_operator(
        A = payload$A,
        B = payload$B,
        metric = payload$metric
      )
      run_geneig_ladder(
        recipe = recipe,
        operator = operator,
        state = payload$state,
        B = B,
        B_total = B_total,
        batch_size = batch_size,
        alpha = alpha,
        seed = seed,
        auto_subspace = auto_subspace,
        tie_threshold = tie_threshold
      )
    })
  }

  if (identical(geom_kind, "multiblock")) {
    return(function(B,
                    B_total,
                    batch_size,
                    alpha,
                    seed,
                    auto_subspace,
                    tie_threshold) {
      .validate_multiblock_data(data)
      run_multiblock_ladder(
        recipe = recipe,
        adapter = adapter,
        data = data,
        B = B,
        B_total = B_total,
        batch_size = batch_size,
        alpha = alpha,
        seed = seed,
        auto_subspace = auto_subspace,
        tie_threshold = tie_threshold,
        original_fit = model
      )
    })
  }

  if (identical(geom_kind, "adapter")) {
    return(function(B,
                    B_total,
                    batch_size,
                    alpha,
                    seed,
                    auto_subspace,
                    tie_threshold) {
      fit <- original_fit()
      run_adapter_ladder(
        recipe = recipe,
        adapter = adapter,
        data = data,
        B = B,
        B_total = B_total,
        batch_size = batch_size,
        alpha = alpha,
        seed = seed,
        auto_subspace = auto_subspace,
        tie_threshold = tie_threshold,
        original_fit = fit
      )
    })
  }

  stop(sprintf(
    "infer() supports only 'oneblock', 'cross', 'multiblock', 'geneig', and 'adapter'; got '%s'.",
    geom_kind
  ), call. = FALSE)
}

.zero_component_engine_out <- function(adapter,
                                       fit,
                                       auto_subspace = TRUE,
                                       tie_threshold = 0.01) {
  roots_observed <- adapter$roots(fit)
  if (!is.numeric(roots_observed) || length(roots_observed) == 0L ||
      any(!is.finite(roots_observed))) {
    stop("`adapter$roots()` must return a non-empty finite numeric vector.",
         call. = FALSE)
  }

  units <- form_units(
    roots_observed,
    selected = rep(FALSE, length(roots_observed)),
    group_near_ties = isTRUE(auto_subspace),
    tie_threshold = tie_threshold
  )

  list(
    units = units,
    component_tests = data.frame(),
    roots_observed = roots_observed,
    ladder_result = list(
      step_results = list(),
      rejected_through = 0L,
      total_draws_used = 0L,
      batch_schedule = integer(0L),
      stopping_boundary = "not_run"
    ),
    labels = list(
      statistic = NA_character_,
      null = NA_character_,
      estimand = "adapter latent roots"
    ),
    original_fit = fit
  )
}

compile_oneblock_ladder_plan <- function(recipe, X, max_steps = NULL) {
  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  if (recipe$shape$geometry$kind != "oneblock") {
    stop(
      paste0("`recipe` geometry must be \"oneblock\"; got \"",
             recipe$shape$geometry$kind, "\"."),
      call. = FALSE
    )
  }
  if (recipe$shape$relation$kind != "variance") {
    stop(
      paste0("`recipe` relation must be \"variance\"; got \"",
             recipe$shape$relation$kind, "\"."),
      call. = FALSE
    )
  }

  .validate_oneblock_data(X)
  Xc <- sweep(X, 2L, colMeans(X), "-")
  roots_observed <- cached_svd(Xc)$d^2
  zero_tol <- max(1, sum(Xc^2)) * .Machine$double.eps

  if (is.null(max_steps)) {
    max_steps <- min(min(dim(Xc)) - 1L, 50L)
  }
  max_steps <- as.integer(max(1L, max_steps))

  observed_stat_fn <- function(step, data) {
    s_top <- top_singular_values(data, 1L)[1L]
    total <- sum(data * data)
    if (total <= zero_tol) return(0)
    (s_top * s_top) / total
  }

  null_stat_fn <- function(step, data) {
    perm <- base::apply(data, 2L, base::sample)
    s <- top_singular_values(perm, step)
    if (length(s) < step) return(0)
    s2 <- s * s
    total_F2 <- sum(data * data)
    tail_total <- total_F2 - sum(s2[seq_len(step - 1L)])
    if (tail_total <= zero_tol) return(0)
    s2[step] / tail_total
  }

  deflate_fn <- function(step, data) {
    sv <- top_svd(data, 1L)
    resid <- data - sv$d[1L] *
      (sv$u[, 1L, drop = FALSE] %*% t(sv$v[, 1L, drop = FALSE]))
    if (sum(resid^2) <= zero_tol) {
      return(matrix(0, nrow = nrow(data), ncol = ncol(data)))
    }
    resid
  }

  structure(
    list(
      initial_data = Xc,
      max_steps = max_steps,
      roots_observed = roots_observed,
      observed_stat_fn = observed_stat_fn,
      null_stat_fn = null_stat_fn,
      deflate_fn = deflate_fn,
      labels = list(
        statistic = "Vitale P3 tail-ratio on latent variance roots",
        null = "column permutation of residual matrix",
        estimand = "latent variance roots"
      )
    ),
    class = "multifer_ladder_plan"
  )
}

compile_cross_ladder_plan <- function(recipe,
                                      X,
                                      Y,
                                      max_steps = NULL) {
  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  if (recipe$shape$geometry$kind != "cross") {
    stop(paste0("`recipe` geometry must be \"cross\"; got \"",
                recipe$shape$geometry$kind, "\"."), call. = FALSE)
  }
  rel_kind <- recipe$shape$relation$kind
  if (!(rel_kind %in% c("covariance", "correlation"))) {
    stop(paste0("`recipe` relation must be \"covariance\" or \"correlation\"; got \"",
                rel_kind, "\"."), call. = FALSE)
  }

  .validate_cross_blocks(X, Y)

  Xc <- sweep(X, 2L, colMeans(X), "-")
  Yc <- sweep(Y, 2L, colMeans(Y), "-")
  design_kind <- recipe$shape$design$kind
  allow_multiroot_correlation <- rel_kind == "correlation" &&
    design_kind %in% c("paired_rows", "nuisance_adjusted", "blocked_rows")

  cross_stat <- if (identical(rel_kind, "covariance")) {
    function(Xr, Yr) cached_svd(crossprod(Xr, Yr))$d
  } else {
    function(Qx, Qy) {
      if (ncol(Qx) == 0L || ncol(Qy) == 0L) {
        return(0)
      }
      cached_svd(crossprod(Qx, Qy))$d
    }
  }

  qx_full <- NULL
  qy_full <- NULL
  corr_groups <- NULL
  if (identical(rel_kind, "covariance")) {
    s_full <- cross_stat(Xc, Yc)
  } else {
    if (identical(design_kind, "nuisance_adjusted")) {
      resid_basis <- .nuisance_residual_basis(
        recipe$shape$design$Z,
        n = nrow(Xc),
        groups = recipe$shape$design$groups
      )
      X_corr <- crossprod(resid_basis$Q, Xc)
      Y_corr <- crossprod(resid_basis$Q, Yc)
      corr_groups <- resid_basis$groups
    } else {
      X_corr <- Xc
      Y_corr <- Yc
      if (identical(design_kind, "blocked_rows")) {
        corr_groups <- recipe$shape$design$groups
      }
    }
    qx_full <- .orth_basis(X_corr)
    qy_full <- .orth_basis(Y_corr)
    s_full <- cross_stat(qx_full, qy_full)
  }

  roots_observed <- s_full^2
  zero_tol <- max(1, sum(roots_observed)) * .Machine$double.eps

  if (is.null(max_steps)) {
    if (identical(rel_kind, "correlation") && allow_multiroot_correlation) {
      max_steps <- min(length(roots_observed), 50L)
    } else {
      max_steps <- min(min(nrow(Xc), ncol(Xc), ncol(Yc)) - 1L, 50L)
    }
  }
  max_steps <- as.integer(max(1L, max_steps))
  if (identical(rel_kind, "correlation") && !allow_multiroot_correlation) {
    max_steps <- min(max_steps, 1L)
  }

  cross_matrix <- if (identical(rel_kind, "covariance")) {
    function(Xr, Yr) base::crossprod(Xr, Yr)
  } else {
    function(Qx, Qy) {
      if (ncol(Qx) == 0L || ncol(Qy) == 0L) {
        return(matrix(0, nrow = 1L, ncol = 1L))
      }
      base::crossprod(Qx, Qy)
    }
  }

  cov_step_cache <- new.env(parent = emptyenv())
  get_cov_step_core <- function(step, data) {
    key <- as.character(step)
    hit <- cov_step_cache[[key]]
    if (!is.null(hit)) return(hit)
    core <- .cross_covariance_core(data$X, data$Y)
    cov_step_cache[[key]] <- core
    core
  }

  observed_stat_fn <- function(step, data) {
    if (identical(rel_kind, "covariance")) {
      core <- get_cov_step_core(step, data)
      if (isTRUE(core$use_core)) {
        return(.cross_covariance_core_observed_sv2(core))
      }
      M <- cross_matrix(data$X, data$Y)
    } else {
      M <- cross_matrix(data$Qx, data$Qy)
    }
    s1 <- top_singular_values(M, 1L)[1L]
    s1 * s1
  }

  null_stat_fn <- function(step, data) {
    if (identical(rel_kind, "covariance")) {
      core <- get_cov_step_core(step, data)
      if (isTRUE(core$use_core)) {
        perm <- base::sample.int(nrow(core$Uy))
        return(.cross_covariance_core_null_sv2(core, perm))
      }
      perm <- base::sample.int(nrow(data$Y))
      Yp <- data$Y[perm, , drop = FALSE]
      M <- cross_matrix(data$X, Yp)
    } else {
      perm <- if (is.null(data$groups)) {
        base::sample.int(nrow(data$Qy))
      } else {
        .restricted_row_permutation(data$groups)
      }
      Qyp <- data$Qy[perm, , drop = FALSE]
      M <- cross_matrix(data$Qx, Qyp)
    }
    s1 <- top_singular_values(M, 1L)[1L]
    s1 * s1
  }

  deflate_fn <- function(step, data) {
    if (identical(rel_kind, "covariance")) {
      core <- get_cov_step_core(step, data)
      if (isTRUE(core$use_core)) {
        next_xy <- .cross_covariance_core_deflate(core, data$X, data$Y)
        Xn <- next_xy$X
        Yn <- next_xy$Y
      } else {
        sv <- top_svd(crossprod(data$X, data$Y), 1L)
        u1 <- sv$u[, 1L, drop = FALSE]
        v1 <- sv$v[, 1L, drop = FALSE]
        Xn <- data$X - data$X %*% u1 %*% t(u1)
        Yn <- data$Y - data$Y %*% v1 %*% t(v1)
      }
      if (sum(Xn^2) + sum(Yn^2) <= zero_tol) {
        return(list(
          X = matrix(0, nrow = nrow(data$X), ncol = ncol(data$X)),
          Y = matrix(0, nrow = nrow(data$Y), ncol = ncol(data$Y))
        ))
      }
      return(list(X = Xn, Y = Yn))
    }

    if (ncol(data$Qx) == 0L || ncol(data$Qy) == 0L) {
      return(list(Qx = data$Qx, Qy = data$Qy, groups = data$groups))
    }
    sv <- top_svd(crossprod(data$Qx, data$Qy), 1L)
    tx <- data$Qx %*% sv$u[, 1L, drop = FALSE]
    ty <- data$Qy %*% sv$v[, 1L, drop = FALSE]
    Qx_next <- .orth_basis(data$Qx - tx %*% crossprod(tx, data$Qx))
    Qy_next <- .orth_basis(data$Qy - ty %*% crossprod(ty, data$Qy))
    list(Qx = Qx_next, Qy = Qy_next, groups = data$groups)
  }

  labels <- .cross_ladder_labels(rel_kind, design_kind, recipe$shape$design)
  structure(
    list(
      initial_data = if (identical(rel_kind, "covariance")) {
        list(X = Xc, Y = Yc)
      } else {
        list(Qx = qx_full, Qy = qy_full, groups = corr_groups)
      },
      max_steps = max_steps,
      roots_observed = roots_observed,
      observed_stat_fn = observed_stat_fn,
      null_stat_fn = null_stat_fn,
      deflate_fn = deflate_fn,
      labels = labels
    ),
    class = "multifer_ladder_plan"
  )
}

is_ladder_plan <- function(x) inherits(x, "multifer_ladder_plan")

.ladder_plan_result <- function(plan,
                                ladder_result,
                                auto_subspace = TRUE,
                                tie_threshold = 0.01) {
  if (!is_ladder_plan(plan)) {
    stop("`plan` must be a multifer_ladder_plan.", call. = FALSE)
  }
  rejected_through <- ladder_result$rejected_through
  roots_observed <- plan$roots_observed
  selected <- logical(length(roots_observed))
  if (rejected_through >= 1L) {
    selected[seq_len(min(rejected_through, length(selected)))] <- TRUE
  }

  units <- form_units(
    roots_observed,
    selected = selected,
    group_near_ties = isTRUE(auto_subspace),
    tie_threshold = tie_threshold
  )

  sr <- ladder_result$step_results
  component_tests <- data.frame(
    step = vapply(sr, function(x) x$step, integer(1L)),
    observed_stat = vapply(sr, function(x) x$observed_stat, double(1L)),
    p_value = vapply(sr, function(x) x$p_value, double(1L)),
    mc_se = vapply(sr, function(x) x$mc_se, double(1L)),
    r = vapply(sr, function(x) x$r, integer(1L)),
    B = vapply(sr, function(x) x$B, integer(1L)),
    selected = vapply(sr, function(x) x$selected, logical(1L)),
    stringsAsFactors = FALSE
  )

  list(
    units = units,
    component_tests = component_tests,
    roots_observed = roots_observed,
    ladder_result = ladder_result,
    labels = plan$labels
  )
}

.validate_cross_blocks <- function(X, Y) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("`Y` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(X) != nrow(Y)) {
    stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  }
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("`X` and `Y` must not contain NA, NaN, or Inf.", call. = FALSE)
  }
  if (nrow(X) < 2L || ncol(X) < 1L || ncol(Y) < 1L) {
    stop("`X` must have at least 2 rows; `X` and `Y` need at least 1 column each.",
         call. = FALSE)
  }
  invisible(NULL)
}

.orth_basis <- function(M, tol = 1e-10) {
  q <- qr(M, tol = tol)
  r <- q$rank
  if (r <= 0L) {
    return(matrix(0, nrow = nrow(M), ncol = 0L))
  }
  qr.Q(q, complete = FALSE)[, seq_len(r), drop = FALSE]
}

.cross_ladder_labels <- function(rel_kind, design_kind, design) {
  estimand_txt <- if (identical(rel_kind, "covariance")) {
    "cross-covariance singular roots"
  } else {
    "canonical correlations"
  }
  statistic_txt <- if (identical(rel_kind, "covariance")) {
    "Vitale P3 tail-ratio on cross-covariance roots"
  } else {
    "Vitale P3 tail-ratio on canonical correlations"
  }
  null_txt <- if (identical(rel_kind, "covariance")) {
    "row permutation of Y"
  } else if (identical(design_kind, "nuisance_adjusted") &&
             is.null(design$groups)) {
    "row permutation of Y in nuisance residual basis"
  } else if (identical(design_kind, "nuisance_adjusted")) {
    "row permutation of Y in nuisance residual basis, within block"
  } else if (identical(design_kind, "blocked_rows")) {
    "row permutation of Y within block"
  } else if (identical(design_kind, "paired_rows")) {
    "row permutation of Y"
  } else {
    "row permutation of Y (first-root only, conservative fallback)"
  }

  list(
    statistic = statistic_txt,
    null = null_txt,
    estimand = estimand_txt
  )
}
