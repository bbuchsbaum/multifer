make_fi_artifact <- function(loadings_by_rep, domains = names(loadings_by_rep[[1L]])) {
  structure(
    list(
      reps = lapply(loadings_by_rep, function(L) list(aligned_loadings = L)),
      R = length(loadings_by_rep),
      method_align = "sign",
      domains = domains,
      seed = NULL
    ),
    class = "multifer_bootstrap_artifact"
  )
}

slow_importance_reference <- function(L, units, unit_idx, alpha) {
  members <- attr(units, "members")
  q <- matrix(NA_real_, nrow = nrow(L), ncol = length(unit_idx))
  for (u in seq_along(unit_idx)) {
    q[, u] <- rowSums(L[, members[[unit_idx[u]]], drop = FALSE]^2)
  }
  as.numeric(q %*% alpha / sum(alpha))
}

test_that("feature_importance_from_bootstrap computes aggregate squared-loading importance", {
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "component"),
    members = list(1L, 2L),
    identifiable = c(TRUE, TRUE),
    selected = c(TRUE, FALSE)
  )
  reps <- list(
    list(aligned_loadings = list(X = matrix(c(1, 0, 0, 1), nrow = 2))),
    list(aligned_loadings = list(X = matrix(c(1, 0, 0, 1), nrow = 2)))
  )
  art <- make_fi_artifact(lapply(reps, `[[`, "aligned_loadings"))

  out <- feature_importance_from_bootstrap(
    art, units, k = "selected", normalize = "none"
  )

  expect_s3_class(out, "multifer_feature_importance")
  expect_equal(out$scope, c("aggregate", "aggregate"))
  expect_equal(out$estimate, c(1, 0))
  expect_false("p_value" %in% names(out))
})

test_that("feature_importance_from_bootstrap handles subspace units and selected semantics", {
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("subspace", "component"),
    members = list(c(1L, 2L), 3L),
    identifiable = c(FALSE, TRUE),
    selected = c(TRUE, FALSE)
  )
  L <- matrix(
    c(1, 2, 0,
      0, 0, 3),
    nrow = 2,
    byrow = TRUE
  )
  reps <- list(
    list(aligned_loadings = list(X = L)),
    list(aligned_loadings = list(X = L))
  )
  art <- make_fi_artifact(lapply(reps, `[[`, "aligned_loadings"))

  out <- feature_importance_from_bootstrap(
    art, units, k = "selected", scope = "both", normalize = "none"
  )

  agg <- out[out$scope == "aggregate", ]
  unit <- out[out$scope == "unit", ]
  expect_equal(agg$estimate, c(5, 0))
  expect_equal(unit$unit_id, c("u1", "u1"))
  expect_equal(unit$estimate, c(5, 0))
})

test_that("feature_importance_from_bootstrap supports root weights, ranks, and top-m frequency", {
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "component"),
    members = list(1L, 2L),
    identifiable = c(TRUE, TRUE),
    selected = c(TRUE, TRUE)
  )
  L1 <- matrix(c(1, 0, 0, 1), nrow = 2)
  L2 <- matrix(c(0, 1, 1, 0), nrow = 2)
  reps <- list(
    list(aligned_loadings = list(X = L1)),
    list(aligned_loadings = list(X = L2))
  )
  art <- make_fi_artifact(lapply(reps, `[[`, "aligned_loadings"))

  out <- feature_importance_from_bootstrap(
    art, units, k = "all", weights = "root", roots = c(3, 1),
    normalize = "block", top_m = 1L
  )

  expect_equal(sum(out$estimate), 1)
  expect_true(all(out$rank_mean >= 1))
  expect_true(all(out$top_m_frequency >= 0 & out$top_m_frequency <= 1))
})

test_that("feature_importance_from_bootstrap matches an independent reference calculation", {
  set.seed(81)
  units <- infer_units(
    unit_id = c("u1", "u2", "u3"),
    unit_type = c("component", "subspace", "component"),
    members = list(1L, c(2L, 3L), 4L),
    identifiable = c(TRUE, FALSE, TRUE),
    selected = c(TRUE, TRUE, FALSE)
  )
  L1 <- matrix(stats::rnorm(24), nrow = 6)
  L2 <- matrix(stats::rnorm(24), nrow = 6)
  rownames(L1) <- rownames(L2) <- paste0("v", seq_len(6))
  roots <- c(4, 3, 2, 1)
  unit_idx <- c(1L, 2L)
  alpha <- c(roots[1], sum(roots[2:3]))
  ref <- rowMeans(cbind(
    slow_importance_reference(L1, units, unit_idx, alpha),
    slow_importance_reference(L2, units, unit_idx, alpha)
  ))

  art <- make_fi_artifact(list(list(X = L1), list(X = L2)))
  out <- feature_importance_from_bootstrap(
    art, units, k = "selected", weights = "root", roots = roots,
    normalize = "none"
  )

  expect_equal(out$variable, paste0("v", seq_len(6)))
  expect_equal(out$estimate, ref, tolerance = 1e-14)
  expect_true(all(is.finite(out$estimate)))
  expect_true(all(out$estimate >= 0))
})

test_that("feature importance is invariant to component sign flips", {
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "component"),
    members = list(1L, 2L),
    identifiable = c(TRUE, TRUE),
    selected = c(TRUE, TRUE)
  )
  L <- matrix(c(1, -2, 3, 4, -5, 6), nrow = 3, ncol = 2)
  art <- make_fi_artifact(list(list(X = L)))
  art_flipped <- make_fi_artifact(list(list(X = sweep(L, 2L, c(-1, 1), `*`))))

  out <- feature_importance_from_bootstrap(art, units, k = "all", normalize = "none")
  out_flipped <- feature_importance_from_bootstrap(art_flipped, units, k = "all", normalize = "none")

  expect_equal(out$estimate, out_flipped$estimate, tolerance = 1e-14)
  expect_equal(out$rank_mean, out_flipped$rank_mean, tolerance = 1e-14)
})

test_that("feature importance obeys raw scale and block-normalization metamorphisms", {
  units <- infer_units(
    unit_id = "u1",
    unit_type = "component",
    members = list(1L),
    identifiable = TRUE,
    selected = TRUE
  )
  L <- matrix(c(1, 2, 3), ncol = 1)
  art <- make_fi_artifact(list(list(X = L)))
  art_scaled <- make_fi_artifact(list(list(X = 10 * L)))

  raw <- feature_importance_from_bootstrap(art, units, normalize = "none")
  raw_scaled <- feature_importance_from_bootstrap(art_scaled, units, normalize = "none")
  norm <- feature_importance_from_bootstrap(art, units, normalize = "block")
  norm_scaled <- feature_importance_from_bootstrap(art_scaled, units, normalize = "block")

  expect_equal(raw_scaled$estimate, 100 * raw$estimate, tolerance = 1e-14)
  expect_equal(norm_scaled$estimate, norm$estimate, tolerance = 1e-14)
  expect_equal(sum(norm$estimate), 1, tolerance = 1e-14)
})

test_that("feature importance is equivariant to variable row permutations", {
  units <- infer_units(
    unit_id = "u1",
    unit_type = "component",
    members = list(1L),
    identifiable = TRUE,
    selected = TRUE
  )
  L <- matrix(c(1, 2, 3, 4), ncol = 1)
  rownames(L) <- c("a", "b", "c", "d")
  perm <- c(3, 1, 4, 2)
  art <- make_fi_artifact(list(list(X = L)))
  art_perm <- make_fi_artifact(list(list(X = L[perm, , drop = FALSE])))

  out <- feature_importance_from_bootstrap(art, units, normalize = "none")
  out_perm <- feature_importance_from_bootstrap(art_perm, units, normalize = "none")

  expect_equal(out_perm$estimate, out$estimate[match(out_perm$variable, out$variable)])
})

test_that("feature importance rejects unavailable component members", {
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "component"),
    members = list(1L, 3L),
    identifiable = c(TRUE, TRUE),
    selected = c(TRUE, TRUE)
  )
  art <- make_fi_artifact(list(list(X = matrix(1, nrow = 2, ncol = 2))))

  expect_error(
    feature_importance_from_bootstrap(art, units, k = "selected"),
    "outside the available loading columns"
  )
})

test_that("feature importance handles empty selected set and validates weights", {
  units <- infer_units(
    unit_id = "u1",
    unit_type = "component",
    members = list(1L),
    identifiable = TRUE,
    selected = FALSE
  )
  art <- make_fi_artifact(list(list(X = matrix(1, nrow = 2, ncol = 1))))

  empty <- feature_importance_from_bootstrap(art, units, k = "selected")
  expect_s3_class(empty, "multifer_feature_importance")
  expect_equal(nrow(empty), 0L)

  expect_error(
    feature_importance_from_bootstrap(art, units, k = "all", weights = "root"),
    "`roots` must be supplied"
  )
  expect_error(
    feature_importance_from_bootstrap(art, units, k = "all", weights = 0),
    "positive sum"
  )
})

test_that("feature_importance_pvalues uses null_action without a bootstrap artifact", {
  clear_adapter_registry()

  set.seed(71)
  X <- matrix(stats::rnorm(80), 20, 4)
  adapter <- adapter_svd()
  register_infer_adapter("svd_oneblock", adapter, overwrite = TRUE)
  recipe <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = adapter
  )
  fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(fit), selected = c(TRUE, FALSE, FALSE, FALSE))

  out <- feature_importance_pvalues(
    recipe = recipe,
    adapter = adapter,
    data = X,
    units = units,
    original_fit = fit,
    B = 5L,
    k = "selected",
    adjust = "maxT",
    seed = 72L
  )

  expect_s3_class(out, "multifer_feature_importance_pvalues")
  expect_equal(nrow(out), ncol(X))
  expect_true(all(out$p_value >= 1 / 6 & out$p_value <= 1))
  expect_true(all(out$p_adjusted >= out$p_value))
})

test_that("feature_importance_pvalues calls null_action B times and uses maxT", {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  L_obs <- matrix(c(2, 1), ncol = 1)
  L_null <- matrix(c(1, 3), ncol = 1)
  adapter <- structure(
    list(
      refit = function(x, new_data, ...) {
        if (isTRUE(attr(new_data, "null"))) list(L = L_null) else list(L = L_obs)
      },
      loadings = function(x, domain = NULL, ...) x$L,
      roots = function(x, ...) 1,
      null_action = function(x, data, ...) {
        calls$n <- calls$n + 1L
        structure(data, null = TRUE)
      }
    ),
    class = "multifer_adapter"
  )
  recipe <- structure(
    list(
      shape = typed_shape(geometry("oneblock"), relation("variance"), exchangeable_rows()),
      validity_level = "conditional"
    ),
    class = "multifer_infer_recipe"
  )
  units <- infer_units(
    unit_id = "u1",
    unit_type = "component",
    members = list(1L),
    identifiable = TRUE,
    selected = TRUE
  )

  out <- feature_importance_pvalues(
    recipe = recipe,
    adapter = adapter,
    data = matrix(0, 2, 1),
    original_fit = list(L = L_obs),
    units = units,
    B = 4L,
    adjust = "maxT",
    normalize = "none"
  )

  expect_equal(calls$n, 4L)
  expect_equal(out$p_value, c(1 / 5, 1))
  expect_equal(out$p_adjusted, c(1, 1))
})

test_that("feature_importance_pvalues drives a geneig (LDA) adapter end-to-end", {
  skip_if_not_installed("MASS")

  adapter <- adapter_lda_refit()
  data <- list(X = as.matrix(iris[, 1:4]), y = iris$Species)
  recipe <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    adapter = adapter,
    targets = "component_significance"
  )
  fit <- adapter$refit(NULL, data)
  rank <- ncol(adapter$loadings(fit))
  units <- form_units(adapter$roots(fit)[seq_len(rank)])

  out <- feature_importance_pvalues(
    recipe = recipe,
    adapter = adapter,
    data = data,
    units = units,
    original_fit = fit,
    B = 25L,
    adjust = "BH",
    seed = 901L
  )

  expect_s3_class(out, "multifer_feature_importance_pvalues")
  expect_equal(nrow(out), ncol(data$X))
  expect_setequal(out$variable, colnames(data$X))
  expect_equal(unique(out$domain), "X")
  expect_true(all(out$p_value >= 1 / 26 & out$p_value <= 1))
  expect_true(all(out$mc_uncertainty >= 0))
})

test_that("feature_importance_pvalues drives a multiblock adapter end-to-end", {
  adapter <- make_stub_multiblock_adapter()
  data <- make_stub_multiblock_data()
  recipe <- infer_recipe(
    geometry = "multiblock",
    relation = "variance",
    adapter = adapter,
    targets = "component_significance"
  )
  fit <- adapter$refit(NULL, data)
  units <- form_units(adapter$roots(fit))

  out <- feature_importance_pvalues(
    recipe = recipe,
    adapter = adapter,
    data = data,
    units = units,
    original_fit = fit,
    B = 20L,
    adjust = "none",
    seed = 902L
  )

  expect_s3_class(out, "multifer_feature_importance_pvalues")
  expect_setequal(unique(out$domain), names(data))
  expected_rows <- sum(vapply(data, ncol, integer(1L)))
  expect_equal(nrow(out), expected_rows)
  expect_true(all(out$p_value >= 1 / 21 & out$p_value <= 1))
})
