test_that("score_stability_from_bootstrap returns a valid table on oneblock", {
  ensure_default_adapters()
  set.seed(501)
  X <- matrix(rnorm(200), 40, 5)
  rec <- infer_recipe(geometry = "oneblock", relation = "variance",
                      adapter = "prcomp_oneblock")
  adapter <- get_infer_adapter("prcomp_oneblock")
  fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(fit))
  art <- bootstrap_fits(rec, adapter, X, fit, units, R = 10, seed = 7)
  tbl <- score_stability_from_bootstrap(art, X, units)
  expect_s3_class(tbl, "multifer_score_stability")
  # n observations per (unit, domain) row count.
  expect_equal(nrow(tbl), nrow(units) * nrow(X))
  expect_true(all(tbl$lower <= tbl$estimate))
  expect_true(all(tbl$upper >= tbl$estimate))
})

test_that("score_stability_from_bootstrap projects ORIGINAL observations", {
  # The defining design lock: even if a particular observation is
  # OMITTED from a given bootstrap rep's resample, it still appears in
  # the score_stability output, because we project the ORIGINAL data
  # through each rep's fit.
  ensure_default_adapters()
  set.seed(502)
  X <- matrix(rnorm(200), 40, 5)
  rec <- infer_recipe(geometry = "oneblock", relation = "variance",
                      adapter = "prcomp_oneblock")
  adapter <- get_infer_adapter("prcomp_oneblock")
  fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(fit))
  art <- bootstrap_fits(rec, adapter, X, fit, units, R = 5, seed = 11)

  # Under R = 5 it is essentially guaranteed at least one observation is
  # missing from at least one rep. But every observation still appears
  # in the output.
  tbl <- score_stability_from_bootstrap(art, X, units)
  expect_equal(sort(unique(tbl$observation[tbl$unit_id == "u1"])),
               seq_len(nrow(X)))
})

test_that("score_stability_from_bootstrap handles cross with both domains", {
  ensure_default_adapters()
  set.seed(503)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  rec <- infer_recipe(geometry = "cross", relation = "covariance",
                      adapter = "cross_svd")
  adapter <- get_infer_adapter("cross_svd")
  fit <- adapter$refit(NULL,
                       list(X = X, Y = Y, relation = "covariance"))
  units <- form_units(adapter$roots(fit))
  art <- bootstrap_fits(rec, adapter, list(X = X, Y = Y), fit, units,
                        R = 8, seed = 13)
  tbl <- score_stability_from_bootstrap(art, list(X = X, Y = Y), units)
  expect_setequal(unique(tbl$domain), c("X", "Y"))
})

test_that("score_stability_from_bootstrap does not require stored aligned_scores", {
  ensure_default_adapters()
  set.seed(504)
  X <- matrix(rnorm(200), 40, 5)
  rec <- infer_recipe(geometry = "oneblock", relation = "variance",
                      adapter = "prcomp_oneblock")
  adapter <- get_infer_adapter("prcomp_oneblock")
  fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(fit))

  art_full <- bootstrap_fits(rec, adapter, X, fit, units,
                             R = 8, seed = 17, store_aligned_scores = TRUE)
  art_light <- bootstrap_fits(rec, adapter, X, fit, units,
                              R = 8, seed = 17, store_aligned_scores = FALSE)

  tbl_full <- score_stability_from_bootstrap(art_full, X, units)
  tbl_light <- score_stability_from_bootstrap(art_light, X, units)

  expect_equal(tbl_light, tbl_full, tolerance = 1e-12)
})

test_that("score_stability rejects malformed inputs", {
  expect_error(score_stability_from_bootstrap("nope", matrix(0, 2, 2),
                                              form_units(c(1))),
               "multifer_bootstrap_artifact")
  art <- structure(list(reps = list(), R = 0L, domains = "X",
                        method_align = "sign", seed = NULL),
                   class = "multifer_bootstrap_artifact")
  expect_error(score_stability_from_bootstrap(art, matrix(0, 2, 2), "nope"),
               "multifer_units")
})
