test_that("geometry accepts the scalar valid kinds", {
  for (kind in c("oneblock", "cross", "multiblock", "adapter")) {
    g <- geometry(kind)
    expect_true(is_geometry(g))
    expect_equal(g$kind, kind)
    expect_s3_class(g, paste0("multifer_geometry_", kind))
  }
})

test_that("geometry constructs a validated geneig geometry", {
  A <- matrix(c(2, 0.5, 0.5, 1), nrow = 2)
  B <- matrix(c(1.5, 0.2, 0.2, 1), nrow = 2)

  g <- geometry("geneig", A = A, B = B, metric = "within_covariance")

  expect_true(is_geometry(g))
  expect_equal(g$kind, "geneig")
  expect_s3_class(g, "multifer_geometry_geneig")
  expect_equal(g$A, A)
  expect_equal(g$B, B)
  expect_equal(g$metric, "within_covariance")
})

test_that("geometry rejects unknown kinds", {
  expect_error(geometry("pca"), "Unknown geometry")
  expect_error(geometry(""), "Unknown geometry")
})

test_that("geometry rejects non-scalar / non-character input", {
  expect_error(geometry(c("oneblock", "cross")), "single non-NA")
  expect_error(geometry(NA_character_), "single non-NA")
  expect_error(geometry(1), "single non-NA")
  expect_error(geometry(NULL), "single non-NA")
})

test_that("geneig geometry rejects invalid operators with specific errors", {
  A_bad <- matrix(c(1, 2, 0, 1), nrow = 2)
  A <- diag(c(2, 1))
  B_nonsym <- matrix(c(1, 1, 0, 1), nrow = 2)
  B_not_pd <- matrix(c(1, 0, 0, 0), nrow = 2)

  expect_error(
    geometry("geneig", A = A_bad, B = diag(2), metric = "B"),
    "`A` must be symmetric"
  )
  expect_error(
    geometry("geneig", A = A, B = B_nonsym, metric = "B"),
    "`B` must be symmetric"
  )
  expect_error(
    geometry("geneig", A = A, B = B_not_pd, metric = "B"),
    "smallest eigenvalue = 0"
  )
  expect_error(
    geometry("geneig", A = A, B = diag(3), metric = "B"),
    "identical dimensions"
  )
  expect_error(
    geometry("geneig", A = A, B = diag(2), metric = ""),
    "single non-empty character string"
  )
})

test_that("geometry prints", {
  g <- geometry("cross")
  expect_output(print(g), "multifer_geometry: cross")
})

# geneig_operator validation (multifer-093.2) ----------------------------------

test_that("geneig_operator accepts valid symmetric A + SPD B", {
  set.seed(1)
  A <- crossprod(matrix(rnorm(30), 10, 3))
  B <- crossprod(matrix(rnorm(30), 10, 3)) + diag(3)
  op <- geneig_operator(A, B)
  expect_true(is_geneig_operator(op))
  expect_equal(op$dim, 3L)
  expect_equal(op$metric, "b_metric")
  expect_s3_class(op, "multifer_geneig_operator")
})

test_that("geneig_operator rejects non-square A", {
  expect_error(
    geneig_operator(matrix(1, 2, 3), diag(3)),
    "must be square"
  )
})

test_that("geneig_operator rejects mismatched dimensions", {
  expect_error(
    geneig_operator(diag(3), diag(4)),
    "matching dimensions"
  )
})

test_that("geneig_operator rejects asymmetric A", {
  A <- matrix(c(1, 2, 0, 1), 2, 2)
  B <- diag(2)
  expect_error(geneig_operator(A, B), "A.*must be symmetric")
})

test_that("geneig_operator rejects asymmetric B", {
  A <- diag(2)
  B <- matrix(c(1, 2, 0, 1), 2, 2)
  expect_error(geneig_operator(A, B), "B.*must be symmetric")
})

test_that("geneig_operator rejects non-positive-definite B via Cholesky", {
  A <- diag(2)
  B <- matrix(c(1, 2, 2, 1), 2, 2)  # symmetric but eigenvalues = {3, -1}
  expect_error(
    geneig_operator(A, B),
    "symmetric positive definite"
  )
})

test_that("geneig_operator rejects NA / Inf", {
  A <- diag(2); A[1, 1] <- NA
  B <- diag(2)
  expect_error(geneig_operator(A, B), "NA, NaN, or Inf")
})

test_that("geneig_operator rejects non-matrix inputs", {
  expect_error(geneig_operator("not a matrix", diag(2)),
               "`A` must be a numeric matrix")
  expect_error(geneig_operator(diag(2), "not a matrix"),
               "`B` must be a numeric matrix")
})

test_that("geneig_operator accepts custom metric label", {
  op <- geneig_operator(diag(2), diag(2), metric = "within_class")
  expect_equal(op$metric, "within_class")
})

test_that("geneig_operator rejects invalid metric arg", {
  expect_error(geneig_operator(diag(2), diag(2), metric = c("a", "b")),
               "single non-NA")
})

test_that("geneig_operator prints with dim and metric", {
  op <- geneig_operator(diag(3), diag(3), metric = "within_class")
  expect_output(print(op),
                "multifer_geneig_operator: dim = 3, metric = within_class")
})

test_that("is_geneig_operator only TRUE for operator objects", {
  op <- geneig_operator(diag(2), diag(2))
  expect_true(is_geneig_operator(op))
  expect_false(is_geneig_operator(diag(2)))
  expect_false(is_geneig_operator(list(A = diag(2), B = diag(2))))
})

# geneig type-system plumbing (multifer-093.2) ---------------------------------

test_that("typed_shape(geneig, generalized_eigen, exchangeable_rows) composes", {
  ts <- typed_shape(
    geometry("geneig"),
    relation("generalized_eigen"),
    exchangeable_rows()
  )
  expect_true(is_typed_shape(ts))
  expect_equal(ts$geometry$kind, "geneig")
  expect_equal(ts$relation$kind, "generalized_eigen")
})

test_that("infer_recipe accepts a mock geneig adapter and refuses other relations", {
  residualize_geneig <- function(x, k, data, ...) data

  mock <- infer_adapter(
    adapter_id = "mock_geneig_093_2",
    adapter_version = "0.0.1",
    shape_kinds = "geneig",
    capabilities = capability_matrix(
      list(geometry = "geneig", relation = "generalized_eigen",
           targets = "component_significance")
    ),
    roots = function(x, ...) x$lambda,
    null_action = function(x, data, ...) data,
    component_stat = function(x, data, k, ...) 1,
    residualize = residualize_geneig,
    validity_level = "conditional",
    declared_assumptions = character(0),
    checked_assumptions = list(),
    geneig_deflation = "b_metric"
  )

  rec <- infer_recipe(
    geometry = "geneig", relation = "generalized_eigen",
    adapter = mock, targets = "component_significance"
  )
  expect_true(is_infer_recipe(rec))
  expect_equal(rec$shape$design$kind, "exchangeable_rows")

  expect_error(
    infer_recipe(
      geometry = "geneig", relation = "variance",
      adapter = mock, targets = "component_significance"
    ),
    "requires relation 'generalized_eigen'"
  )
})

test_that("infer_recipe refuses geneig geometry on non-geneig adapters", {
  ensure_default_adapters()
  pca_adapter <- get_infer_adapter("svd_oneblock")
  expect_error(
    infer_recipe(
      geometry = "geneig", relation = "generalized_eigen",
      adapter = pca_adapter
    ),
    "shape_kinds"
  )
})
