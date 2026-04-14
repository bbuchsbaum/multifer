#' Geometry of a latent operator
#'
#' The first component of a typed shape. Declares the geometric family
#' of the latent operator being inferred on, independent of the relation
#' (variance / covariance / correlation / generalized eigenvalue) and
#' independent of the experimental design.
#'
#' @param kind One of `"oneblock"`, `"cross"`, `"multiblock"`, `"geneig"`.
#' @param A Optional for `kind = "geneig"`: the left operator in the
#'   generalized eigenvalue problem `A v = lambda B v`. When supplied,
#'   must be a symmetric square matrix.
#' @param B Optional for `kind = "geneig"`: the metric operator in the
#'   generalized eigenvalue problem. When supplied, must be a symmetric
#'   positive definite square matrix with the same dimensions as `A`.
#' @param metric Optional for `kind = "geneig"`: a scalar label
#'   describing the metric or whitening convention carried by `B`.
#'
#' @return An object of class `multifer_geometry` (and a `kind`-specific
#'   subclass).
#'
#' @details
#' - `oneblock` — a single data block, e.g. PCA, kernel PCA.
#' - `cross` — two paired blocks, e.g. PLS, CCA, reduced-rank regression.
#' - `multiblock` — three or more blocks with shared structure,
#'   e.g. JIVE, AJIVE, multiblock projectors.
#' - `geneig` — generalized-eigen / contrastive / metric-weighted fits,
#'   e.g. cPCA++, discriminant analysis.
#'
#' @export
geometry <- function(kind, A = NULL, B = NULL, metric = NULL) {
  valid <- c("oneblock", "cross", "multiblock", "geneig")
  if (!is.character(kind) || length(kind) != 1L || is.na(kind)) {
    stop("`kind` must be a single non-NA character string.", call. = FALSE)
  }
  if (!(kind %in% valid)) {
    stop(sprintf(
      "Unknown geometry `%s`. Must be one of: %s.",
      kind, paste(valid, collapse = ", ")
    ), call. = FALSE)
  }

  payload <- list(kind = kind)
  if (identical(kind, "geneig")) {
    if (!all(vapply(list(A, B, metric), is.null, logical(1)))) {
      payload <- c(payload, .validate_geneig_geometry(A = A, B = B, metric = metric))
    }
  }

  structure(
    payload,
    class = c(paste0("multifer_geometry_", kind), "multifer_geometry")
  )
}

#' @export
print.multifer_geometry <- function(x, ...) {
  cat("<multifer_geometry: ", x$kind, ">\n", sep = "")
  invisible(x)
}

#' Test whether an object is a multifer geometry
#' @param x Object to test.
#' @export
is_geometry <- function(x) inherits(x, "multifer_geometry")

.validate_geneig_geometry <- function(A, B, metric) {
  if (any(vapply(list(A, B, metric), is.null, logical(1)))) {
    stop(
      "For `kind = \"geneig\"`, `A`, `B`, and `metric` must be supplied together.",
      call. = FALSE
    )
  }
  .validate_square_symmetric_matrix(A, name = "A")
  .validate_square_symmetric_matrix(B, name = "B")

  if (!identical(dim(A), dim(B))) {
    stop("`A` and `B` must have identical dimensions for `kind = \"geneig\"`.",
         call. = FALSE)
  }
  if (!is.character(metric) || length(metric) != 1L || is.na(metric) ||
      nchar(metric) == 0L) {
    stop("`metric` must be a single non-empty character string for `kind = \"geneig\"`.",
         call. = FALSE)
  }

  eigvals <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eigvals)
  tol <- .matrix_symmetry_tolerance(B)
  if (!is.finite(min_eig) || min_eig <= tol) {
    stop(sprintf(
      "B must be symmetric positive definite; smallest eigenvalue = %.6g.",
      min_eig
    ), call. = FALSE)
  }

  list(A = A, B = B, metric = metric)
}

.validate_square_symmetric_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || length(dim(x)) != 2L ||
      nrow(x) != ncol(x)) {
    stop(sprintf(
      "`%s` must be a numeric square matrix for `kind = \"geneig\"`.",
      name
    ), call. = FALSE)
  }

  tol <- .matrix_symmetry_tolerance(x)
  if (max(abs(x - t(x))) > tol) {
    stop(sprintf(
      "`%s` must be symmetric for `kind = \"geneig\"`.",
      name
    ), call. = FALSE)
  }

  invisible(NULL)
}

.matrix_symmetry_tolerance <- function(x) {
  max(1, max(abs(x))) * sqrt(.Machine$double.eps)
}

#' Validated operator pair for a `geneig` problem
#'
#' Constructs a validated `(A, B)` operator pair for a generalized
#' eigenvalue problem `A v = lambda B v`. This is a data-side helper,
#' not a shape: it carries the symmetric operator `A` and the
#' symmetric positive-definite metric `B` plus optional metadata, and
#' validates them at construction time so the engine does not have to
#' re-check on every rung.
#'
#' The matching typed shape is
#' `typed_shape(geometry("geneig"), relation("generalized_eigen"), design)`.
#'
#' @param A A square, symmetric numeric matrix. Symmetry is checked to
#'   a tolerance proportional to `max(abs(A)) * sqrt(.Machine$double.eps)`.
#' @param B A square, symmetric positive-definite numeric matrix of the
#'   same dimension as `A`. Positive-definiteness is verified via a
#'   Cholesky attempt; failure produces a specific error message
#'   naming the failing pivot.
#' @param metric Character scalar labelling the metric that `B`
#'   encodes (e.g. `"within_class"` for LDA, `"background_cov"` for
#'   contrastive PCA). Purely descriptive — not consumed by the
#'   engine. Default `"b_metric"`.
#'
#' @details
#' The constructor exists so that geneig adapters can accept a single
#' data argument and know statically that `A` and `B` have been
#' validated. See `notes/engine_geneig_spec.md` for the full engine
#' spec, including why the deflation is B-orthogonal rather than
#' Euclidean.
#'
#' @return An object of class `multifer_geneig_operator`, structurally
#'   a list with fields `A`, `B`, `metric`, and `dim`.
#'
#' @examples
#' A <- crossprod(matrix(stats::rnorm(30), 10, 3))
#' B <- crossprod(matrix(stats::rnorm(30), 10, 3)) + diag(3)
#' op <- geneig_operator(A, B)
#'
#' @export
geneig_operator <- function(A, B, metric = "b_metric") {
  if (!is.matrix(A) || !is.numeric(A)) {
    stop("`A` must be a numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(B) || !is.numeric(B)) {
    stop("`B` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(A) != ncol(A)) {
    stop(sprintf("`A` must be square; got %d x %d.", nrow(A), ncol(A)),
         call. = FALSE)
  }
  if (nrow(B) != ncol(B)) {
    stop(sprintf("`B` must be square; got %d x %d.", nrow(B), ncol(B)),
         call. = FALSE)
  }
  if (nrow(A) != nrow(B)) {
    stop(sprintf("`A` and `B` must have matching dimensions; got %d and %d.",
                 nrow(A), nrow(B)), call. = FALSE)
  }
  if (!is.character(metric) || length(metric) != 1L || is.na(metric)) {
    stop("`metric` must be a single non-NA character string.", call. = FALSE)
  }
  if (any(!is.finite(A)) || any(!is.finite(B))) {
    stop("`A` and `B` must not contain NA, NaN, or Inf.", call. = FALSE)
  }

  tol_A <- max(1, max(abs(A))) * sqrt(.Machine$double.eps)
  if (max(abs(A - t(A))) > tol_A) {
    stop(sprintf(
      "`A` must be symmetric; max |A - t(A)| = %.3g exceeds tolerance %.3g.",
      max(abs(A - t(A))), tol_A
    ), call. = FALSE)
  }

  tol_B <- max(1, max(abs(B))) * sqrt(.Machine$double.eps)
  if (max(abs(B - t(B))) > tol_B) {
    stop(sprintf(
      "`B` must be symmetric; max |B - t(B)| = %.3g exceeds tolerance %.3g.",
      max(abs(B - t(B))), tol_B
    ), call. = FALSE)
  }

  chol_try <- tryCatch(chol((B + t(B)) / 2),
                       error = function(e) e)
  if (inherits(chol_try, "error")) {
    stop(sprintf(
      "`B` must be symmetric positive definite; Cholesky failed: %s",
      conditionMessage(chol_try)
    ), call. = FALSE)
  }

  structure(
    list(
      A      = A,
      B      = B,
      metric = metric,
      dim    = nrow(A)
    ),
    class = "multifer_geneig_operator"
  )
}

#' @export
print.multifer_geneig_operator <- function(x, ...) {
  cat("<multifer_geneig_operator: dim = ", x$dim,
      ", metric = ", x$metric, ">\n", sep = "")
  invisible(x)
}

#' Test whether an object is a validated geneig operator pair
#' @param x Object to test.
#' @export
is_geneig_operator <- function(x) inherits(x, "multifer_geneig_operator")
