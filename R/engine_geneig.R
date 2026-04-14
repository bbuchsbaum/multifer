#' Run the generalized-eigen sequential deflation engine
#'
#' Implements the test ladder for `(geneig, generalized_eigen)` recipes
#' via the shared [ladder_driver()]. At each rung the engine whitens the
#' current operator pair to
#'
#' `C = B^{-1/2} A B^{-1/2}`
#'
#' computes the ordered generalized roots of `C`, evaluates the
#' tail-ratio statistic on the leading remaining root, and deflates in
#' the `B`-metric by removing the top eigendirection in the whitened
#' space. Null resampling is delegated to the adapter's `null_action()`
#' hook because valid nulls are sub-family-specific.
#'
#' @param recipe A compiled `multifer_infer_recipe`. Must have geometry
#'   `"geneig"` and relation `"generalized_eigen"`.
#' @param operator A [geneig_operator()] object containing the operator
#'   pair `(A, B)`.
#' @param state Optional list carried through the ladder and handed to
#'   `recipe$adapter$null_action()`. This is primarily for adapter-side
#'   bookkeeping; when `state$factor_matrix` is supplied and `B = I`,
#'   the engine preserves that factor matrix through deflation so the
#'   oneblock regression can be checked exactly.
#' @param B Per-rung cap on Monte Carlo draws. Default `1000L`.
#' @param B_total Optional integer global Monte Carlo budget shared
#'   across ladder rungs. Defaults to `B * max_steps`.
#' @param batch_size Positive integer, Besag-Clifford batch size within
#'   each rung. Default `32L`.
#' @param alpha Significance threshold. Default `0.05`.
#' @param max_steps Maximum number of ladder rungs. Defaults to
#'   `min(length(roots_observed), 50L)`.
#' @param seed Integer or `NULL`.
#' @param auto_subspace Logical. When `TRUE` (default), near-tied roots
#'   are bundled via [form_units()] `group_near_ties = TRUE`.
#' @param tie_threshold Positive numeric. Relative-gap threshold used
#'   when `auto_subspace = TRUE`. Default `0.01`.
#'
#' @return A list with `units`, `component_tests`, `roots_observed`, and
#'   `ladder_result`, mirroring [run_oneblock_ladder()].
#'
#' @export
run_geneig_ladder <- function(recipe,
                              operator,
                              state         = NULL,
                              B             = 1000L,
                              B_total       = NULL,
                              batch_size    = 32L,
                              alpha         = 0.05,
                              max_steps     = NULL,
                              seed          = NULL,
                              auto_subspace = TRUE,
                              tie_threshold = 0.01) {

  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  if (recipe$shape$geometry$kind != "geneig") {
    stop(
      paste0("`recipe` geometry must be \"geneig\"; got \"",
             recipe$shape$geometry$kind, "\"."),
      call. = FALSE
    )
  }
  if (recipe$shape$relation$kind != "generalized_eigen") {
    stop(
      paste0("`recipe` relation must be \"generalized_eigen\"; got \"",
             recipe$shape$relation$kind, "\"."),
      call. = FALSE
    )
  }
  if (!is_geneig_operator(operator)) {
    stop("`operator` must be a multifer_geneig_operator (see geneig_operator()).",
         call. = FALSE)
  }
  if (!is.null(state) && !is.list(state)) {
    stop("`state` must be NULL or a list.", call. = FALSE)
  }

  initial_data <- .geneig_normalize_state(list(
    A      = operator$A,
    B      = operator$B,
    metric = operator$metric,
    state  = state
  ))

  roots_observed <- .geneig_roots(initial_data)
  zero_tol <- max(1, sum(roots_observed)) * .Machine$double.eps

  if (is.null(max_steps)) {
    max_steps <- min(length(roots_observed), 50L)
  }
  max_steps <- as.integer(max(1L, max_steps))

  observed_stat_fn <- function(step, data) {
    .geneig_observed_stat(data, zero_tol = zero_tol)
  }

  null_stat_fn <- function(step, data) {
    null_payload <- recipe$adapter$null_action(recipe$shape, data)
    null_data <- .geneig_normalize_state(null_payload, template = data)
    .geneig_observed_stat(null_data, zero_tol = zero_tol)
  }

  deflate_fn <- function(step, data) {
    .geneig_deflate_state(data, zero_tol = zero_tol)
  }

  ladder_result <- ladder_driver(
    observed_stat_fn = observed_stat_fn,
    null_stat_fn     = null_stat_fn,
    deflate_fn       = deflate_fn,
    initial_data     = initial_data,
    max_steps        = max_steps,
    B                = B,
    B_total          = B_total,
    batch_size       = batch_size,
    alpha            = alpha,
    seed             = seed
  )

  rejected_through <- ladder_result$rejected_through
  selected <- logical(length(roots_observed))
  if (rejected_through >= 1L) {
    selected[seq_len(rejected_through)] <- TRUE
  }

  units <- form_units(
    roots_observed,
    selected        = selected,
    group_near_ties = isTRUE(auto_subspace),
    tie_threshold   = tie_threshold
  )

  sr <- ladder_result$step_results
  component_tests <- data.frame(
    step          = vapply(sr, function(x) x$step,          integer(1L)),
    observed_stat = vapply(sr, function(x) x$observed_stat, double(1L)),
    p_value       = vapply(sr, function(x) x$p_value,       double(1L)),
    mc_se         = vapply(sr, function(x) x$mc_se,         double(1L)),
    r             = vapply(sr, function(x) x$r,             integer(1L)),
    B             = vapply(sr, function(x) x$B,             integer(1L)),
    selected      = vapply(sr, function(x) x$selected,      logical(1L)),
    stringsAsFactors = FALSE
  )

  list(
    units           = units,
    component_tests = component_tests,
    roots_observed  = roots_observed,
    ladder_result   = ladder_result
  )
}

.geneig_normalize_state <- function(payload, template = NULL) {
  if (is_geneig_operator(payload)) {
    payload <- list(A = payload$A, B = payload$B, metric = payload$metric)
  }
  if (!is.list(payload)) {
    stop("Geneig engine payloads must be a geneig_operator or a list.", call. = FALSE)
  }

  A <- payload$A
  B <- payload$B
  metric <- payload$metric
  state <- payload$state

  if (is.null(B) && !is.null(template)) {
    B <- template$B
  }
  if (is.null(metric) && !is.null(template)) {
    metric <- template$metric
  }
  if (is.null(state) && !is.null(template)) {
    state <- template$state
  }

  op <- geneig_operator(A = A, B = B, metric = metric)
  b_dec <- .geneig_metric_decomposition(op$B)

  list(
    A             = op$A,
    B             = op$B,
    metric        = op$metric,
    state         = state,
    B_half        = b_dec$B_half,
    B_inv_half    = b_dec$B_inv_half,
    B_is_identity = .geneig_is_identity(op$B)
  )
}

.geneig_metric_decomposition <- function(B) {
  eig <- eigen((B + t(B)) / 2, symmetric = TRUE)
  vals <- eig$values
  min_eig <- min(vals)
  tol <- max(1, max(abs(B))) * sqrt(.Machine$double.eps)
  if (!is.finite(min_eig) || min_eig <= tol) {
    stop(sprintf(
      "B must be symmetric positive definite; smallest eigenvalue = %.6g.",
      min_eig
    ), call. = FALSE)
  }

  vecs <- eig$vectors
  B_half <- vecs %*% diag(sqrt(vals), nrow = length(vals)) %*% t(vecs)
  B_inv_half <- vecs %*% diag(1 / sqrt(vals), nrow = length(vals)) %*% t(vecs)
  list(B_half = B_half, B_inv_half = B_inv_half)
}

.geneig_whiten <- function(data) {
  C <- data$B_inv_half %*% data$A %*% data$B_inv_half
  (C + t(C)) / 2
}

.geneig_roots <- function(data) {
  if (isTRUE(data$B_is_identity) &&
      is.list(data$state) &&
      !is.null(data$state$factor_matrix)) {
    return(cached_svd(data$state$factor_matrix)$d^2)
  }

  vals <- eigen(.geneig_whiten(data), symmetric = TRUE, only.values = TRUE)$values
  pmax(vals, 0)
}

.geneig_tail_ratio <- function(roots, zero_tol) {
  if (length(roots) == 0L) {
    return(0)
  }
  tail_total <- sum(roots)
  if (!is.finite(tail_total) || tail_total <= zero_tol) {
    return(0)
  }
  roots[1L] / tail_total
}

.geneig_observed_stat <- function(data, zero_tol) {
  if (isTRUE(data$B_is_identity) &&
      is.list(data$state) &&
      !is.null(data$state$factor_matrix)) {
    X <- data$state$factor_matrix
    s1 <- top_singular_values(X, 1L)[1L]
    total <- sum(X * X)
    if (!is.finite(total) || total <= zero_tol) {
      return(0)
    }
    return((s1 * s1) / total)
  }

  .geneig_tail_ratio(.geneig_roots(data), zero_tol = zero_tol)
}

.geneig_deflate_state <- function(data, zero_tol) {
  if (isTRUE(data$B_is_identity) &&
      is.list(data$state) &&
      is.matrix(data$state$factor_matrix) &&
      ncol(data$state$factor_matrix) == nrow(data$A)) {
    sv <- top_svd(data$state$factor_matrix, 1L)
    X_next <- data$state$factor_matrix -
      sv$d[1L] * (sv$u[, 1L, drop = FALSE] %*% t(sv$v[, 1L, drop = FALSE]))
    A_next <- crossprod(X_next)
    state_next <- data$state
    state_next$factor_matrix <- X_next
    return(.geneig_normalize_state(list(
      A = (A_next + t(A_next)) / 2,
      B = data$B,
      metric = data$metric,
      state = state_next
    )))
  }

  C <- .geneig_whiten(data)
  eig <- eigen(C, symmetric = TRUE)
  if (length(eig$values) == 0L || eig$values[1L] <= zero_tol) {
    A_zero <- matrix(0, nrow = nrow(data$A), ncol = ncol(data$A))
    return(.geneig_normalize_state(list(
      A = A_zero,
      B = data$B,
      metric = data$metric,
      state = data$state
    )))
  }

  u1 <- eig$vectors[, 1L, drop = FALSE]
  proj <- diag(nrow(C)) - u1 %*% t(u1)
  C_next <- proj %*% C %*% proj
  A_next <- data$B_half %*% C_next %*% data$B_half
  A_next <- (A_next + t(A_next)) / 2
  if (sum(abs(A_next)) <= zero_tol) {
    A_next[,] <- 0
  }

  .geneig_normalize_state(list(
    A = A_next,
    B = data$B,
    metric = data$metric,
    state = data$state
  ))
}

.geneig_is_identity <- function(B) {
  tol <- max(1, max(abs(B))) * sqrt(.Machine$double.eps)
  isTRUE(all.equal(B, diag(nrow(B)), tolerance = tol))
}
