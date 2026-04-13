test_that("variable_stability_from_bootstrap returns a valid table on oneblock", {
  ensure_default_adapters()
  set.seed(401)
  X <- matrix(rnorm(200), 40, 5)
  rec <- infer_recipe(geometry = "oneblock", relation = "variance",
                      adapter = "prcomp_oneblock")
  adapter <- get_infer_adapter("prcomp_oneblock")
  fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(fit))
  art <- bootstrap_fits(rec, adapter, X, fit, units, R = 10, seed = 7)
  tbl <- variable_stability_from_bootstrap(art, units)
  expect_s3_class(tbl, "multifer_variable_stability")
  expect_true(nrow(tbl) == nrow(units) * ncol(X))
  expect_true(all(tbl$stability >= 0 & tbl$stability <= 1))
  expect_true(all(tbl$domain == "X"))
})

test_that("variable_stability_from_bootstrap covers both domains for cross", {
  ensure_default_adapters()
  set.seed(402)
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
  tbl <- variable_stability_from_bootstrap(art, units)
  expect_setequal(unique(tbl$domain), c("X", "Y"))
  expect_true(nrow(tbl) > 0L)
})

test_that("stability is closer to 1 on a strong-signal fit than on noise", {
  ensure_default_adapters()
  set.seed(403)
  signal <- matrix(rnorm(40 * 3), 40, 3) %*% diag(c(8, 6, 4)) %*%
    matrix(rnorm(3 * 5), 3, 5)
  X_signal <- signal + 0.1 * matrix(rnorm(200), 40, 5)
  X_noise  <- matrix(rnorm(200), 40, 5)
  rec <- infer_recipe(geometry = "oneblock", relation = "variance",
                      adapter = "prcomp_oneblock")
  adapter <- get_infer_adapter("prcomp_oneblock")
  fit_signal <- adapter$refit(NULL, X_signal)
  fit_noise  <- adapter$refit(NULL, X_noise)
  u_signal <- form_units(adapter$roots(fit_signal))
  u_noise  <- form_units(adapter$roots(fit_noise))
  art_signal <- bootstrap_fits(rec, adapter, X_signal, fit_signal,
                               u_signal, R = 15, seed = 1)
  art_noise  <- bootstrap_fits(rec, adapter, X_noise, fit_noise,
                               u_noise, R = 15, seed = 1)
  tbl_signal <- variable_stability_from_bootstrap(art_signal, u_signal)
  tbl_noise  <- variable_stability_from_bootstrap(art_noise, u_noise)
  # Mean stability of the leading unit should be higher under signal.
  s1_signal <- mean(tbl_signal$stability[tbl_signal$unit_id == "u1"])
  s1_noise  <- mean(tbl_noise$stability[tbl_noise$unit_id == "u1"])
  expect_true(s1_signal > s1_noise)
})

test_that("variable_stability rejects malformed inputs", {
  expect_error(variable_stability_from_bootstrap("nope", form_units(c(1))),
               "multifer_bootstrap_artifact")
  art <- structure(list(reps = list(), R = 0L, domains = "X",
                        method_align = "sign", seed = NULL),
                   class = "multifer_bootstrap_artifact")
  expect_error(variable_stability_from_bootstrap(art, "nope"),
               "multifer_units")
})
