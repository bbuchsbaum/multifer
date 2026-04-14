default_adapter_negative_fixture <- function(adapter_id) {
  switch(
    adapter_id,
    "svd_oneblock" = list(
      data = cbind(matrix(stats::rnorm(24), nrow = 8L), 1),
      geometry = "oneblock",
      relation = "variance",
      pattern = "oneblock_columns_have_variance"
    ),
    "prcomp_oneblock" = list(
      data = cbind(matrix(stats::rnorm(24), nrow = 8L), 1),
      geometry = "oneblock",
      relation = "variance",
      pattern = "oneblock_columns_have_variance"
    ),
    "multivarious_pca" = list(
      data = cbind(matrix(stats::rnorm(24), nrow = 8L), 1),
      geometry = "oneblock",
      relation = "variance",
      pattern = "oneblock_columns_have_variance"
    ),
    "cross_svd" = list(
      data = list(
        X = cbind(matrix(stats::rnorm(40), nrow = 10L), 1),
        Y = matrix(stats::rnorm(30), nrow = 10L)
      ),
      geometry = "cross",
      relation = "covariance",
      pattern = "cross_columns_have_variance"
    ),
    "cancor_cross" = list(
      data = list(
        X = cbind(matrix(stats::rnorm(40), nrow = 10L), 1),
        Y = matrix(stats::rnorm(30), nrow = 10L)
      ),
      geometry = "cross",
      relation = "correlation",
      pattern = "cross_columns_have_variance"
    ),
    "multivarious_plsc" = list(
      data = list(
        X = cbind(matrix(stats::rnorm(40), nrow = 10L), 1),
        Y = matrix(stats::rnorm(30), nrow = 10L)
      ),
      geometry = "cross",
      relation = "covariance",
      pattern = "cross_columns_have_variance"
    ),
    "lda_refit" = list(
      data = list(
        X = matrix(stats::rnorm(24), nrow = 8L),
        y = factor(rep("a", 8L))
      ),
      geometry = "geneig",
      relation = "generalized_eigen",
      pattern = "lda_grouping_factor"
    ),
    "plsr_refit" = list(
      data = list(
        X = matrix(stats::rnorm(24), nrow = 8L),
        Y = matrix("bad", nrow = 8L, ncol = 1L)
      ),
      geometry = "cross",
      relation = "predictive",
      pattern = "predictive_y_numeric_matrix"
    ),
    stop(sprintf("No negative fixture registered for adapter '%s'.", adapter_id),
         call. = FALSE)
  )
}

test_that("every default adapter ships concrete checked assumptions", {
  ensure_default_adapters()

  ids <- list_infer_adapters()
  expect_true(length(ids) >= 4L)

  for (adapter_id in ids) {
    adapter <- get_infer_adapter(adapter_id)
    expect_true(length(adapter$checked_assumptions) > 0L, info = adapter_id)

    for (check in adapter$checked_assumptions) {
      expect_true(is.list(check), info = adapter_id)
      expect_true(is.character(check$name) && length(check$name) == 1L && nzchar(check$name),
                  info = adapter_id)
      expect_true(is.character(check$detail) && length(check$detail) == 1L && nzchar(check$detail),
                  info = adapter_id)
      expect_true(is.function(check$check), info = paste(adapter_id, check$name))
    }
  }
})

test_that("every default adapter fails strict infer() on a named negative fixture", {
  ensure_default_adapters()

  for (adapter_id in list_infer_adapters()) {
    fixture <- default_adapter_negative_fixture(adapter_id)
    expect_error(
      infer(
        adapter = adapter_id,
        data = fixture$data,
        geometry = fixture$geometry,
        relation = fixture$relation,
        strict = TRUE,
        B = 9L,
        R = 0L,
        seed = 1L
      ),
      fixture$pattern,
      info = adapter_id
    )
  }
})

test_that("every default adapter records a failed named check in permissive mode", {
  ensure_default_adapters()

  for (adapter_id in list_infer_adapters()) {
    fixture <- default_adapter_negative_fixture(adapter_id)
    adapter <- get_infer_adapter(adapter_id)
    recipe <- infer_recipe(
      geometry = fixture$geometry,
      relation = fixture$relation,
      adapter = adapter
    )
    results <- run_adapter_checks(adapter, fixture$data, recipe = recipe, strict = FALSE)

    expect_true(fixture$pattern %in% names(results), info = adapter_id)
    expect_false(isTRUE(results[[fixture$pattern]]$passed), info = adapter_id)
  }
})
