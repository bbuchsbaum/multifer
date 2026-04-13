test_that("subspace_stability_from_bootstrap returns one row per unit", {
  ensure_default_adapters()
  set.seed(601)
  signal <- matrix(rnorm(40 * 3), 40, 3) %*% diag(c(8, 6, 4)) %*%
    matrix(rnorm(3 * 6), 3, 6)
  X <- signal + 0.1 * matrix(rnorm(240), 40, 6)
  rec <- infer_recipe(geometry = "oneblock", relation = "variance",
                      adapter = "prcomp_oneblock")
  adapter <- get_infer_adapter("prcomp_oneblock")
  fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(fit))
  art <- bootstrap_fits(rec, adapter, X, fit, units, R = 12, seed = 7)
  tbl <- subspace_stability_from_bootstrap(art, fit, adapter, units)
  expect_s3_class(tbl, "multifer_subspace_stability")
  expect_equal(nrow(tbl), nrow(units))
  expect_true(all(tbl$principal_angle_max >= tbl$principal_angle_mean))
})

test_that("strong-signal units are labeled stable", {
  ensure_default_adapters()
  set.seed(602)
  signal <- matrix(rnorm(40 * 3), 40, 3) %*% diag(c(20, 15, 10)) %*%
    matrix(rnorm(3 * 5), 3, 5)
  X <- signal + 0.1 * matrix(rnorm(200), 40, 5)
  rec <- infer_recipe(geometry = "oneblock", relation = "variance",
                      adapter = "prcomp_oneblock")
  adapter <- get_infer_adapter("prcomp_oneblock")
  fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(fit))
  art <- bootstrap_fits(rec, adapter, X, fit, units, R = 15, seed = 1)
  tbl <- subspace_stability_from_bootstrap(art, fit, adapter, units)
  # The leading unit on a strong fit should be stable.
  expect_true(tbl$stability_label[tbl$unit_id == "u1"] %in%
              c("stable", "tied"))
})

test_that("subspace_stability handles cross geometry across two domains", {
  ensure_default_adapters()
  set.seed(603)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  rec <- infer_recipe(geometry = "cross", relation = "covariance",
                      adapter = "cross_svd")
  adapter <- get_infer_adapter("cross_svd")
  fit <- adapter$refit(NULL,
                       list(X = X, Y = Y, relation = "covariance"))
  units <- form_units(adapter$roots(fit))
  art <- bootstrap_fits(rec, adapter, list(X = X, Y = Y), fit, units,
                        R = 8, seed = 9)
  tbl <- subspace_stability_from_bootstrap(art, fit, adapter, units)
  expect_equal(nrow(tbl), nrow(units))
  expect_true(all(is.finite(tbl$principal_angle_mean) |
                  is.na(tbl$principal_angle_mean)))
})

test_that("subspace_stability marks units beyond a truncated fit rank unknown", {
  skip_if_not_installed("multivarious")
  ensure_default_adapters()

  set.seed(604)
  X <- matrix(rnorm(40 * 6), 40, 6)
  rec <- infer_recipe(geometry = "oneblock", relation = "variance",
                      adapter = "multivarious_pca")
  adapter <- get_infer_adapter("multivarious_pca")
  fit <- multivarious::pca(X, ncomp = 5)
  units <- infer_units(
    unit_id = sprintf("u%d", 1:6),
    unit_type = rep("component", 6),
    members = lapply(1:6, function(i) as.integer(i)),
    identifiable = rep(TRUE, 6),
    selected = rep(FALSE, 6)
  )
  art <- bootstrap_fits(rec, adapter, X, fit, units, R = 4, seed = 13)
  tbl <- subspace_stability_from_bootstrap(art, fit, adapter, units)

  expect_true(is.na(tbl$principal_angle_mean[tbl$unit_id == "u6"]))
  expect_true(is.na(tbl$principal_angle_max[tbl$unit_id == "u6"]))
  expect_equal(tbl$stability_label[tbl$unit_id == "u6"], "unknown")
})

test_that("subspace_stability rejects malformed inputs", {
  expect_error(subspace_stability_from_bootstrap("nope", NULL, NULL,
                                                  form_units(c(1))),
               "multifer_bootstrap_artifact")
})
