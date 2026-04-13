#' Linear operator interface
#'
#' Minimum surface (from Part 4 section 33) that the Phase 1.5 fast path
#' and core-space updates consume. A `multifer_linear_operator` exposes
#' `apply`, `apply_t`, `dim`, and optionally `gram` and `core_update`.
#' The interface is deliberately tiny so that dense matrices, sparse
#' matrices, implicit residualizers, kernel operators, and (later)
#' randomized sketches all satisfy it.
#'
#' This is an interface, not an optimization in itself. It exists so
#' that `core()` / `update_core()` hooks in Phase 1.5 adapters can
#' accept "operator-shaped" originals without forcing a dense materialization.
#'
#' Contract:
#' \describe{
#'   \item{`apply(v)`}{`n x m` operator times `p x k` matrix. Returns `n x k`.}
#'   \item{`apply_t(v)`}{Transposed apply: `p x n` times `n x k`. Returns `p x k`.}
#'   \item{`dim()`}{Integer vector `c(n, p)`.}
#'   \item{`gram()`}{Optional. Returns `p x p` Gram matrix `t(op) %*% op`.
#'     If not supplied, callers compute it via `apply_t(apply(I))`.}
#'   \item{`core_update(indices, ...)`}{Optional. Advisory hook used by
#'     Phase 1.5 adapters to accelerate row-indexed perturbations.}
#' }
#'
#' @param apply Function. Applies the operator to a matrix.
#' @param apply_t Function. Applies the transposed operator.
#' @param dim Integer vector of length 2: `c(nrow, ncol)`.
#' @param gram Optional function returning `t(op) %*% op`.
#' @param core_update Optional function for Phase 1.5 core-space updates.
#' @param repr Optional list of internal data the operator wraps (useful
#'   for debugging and for `as.matrix()` materialization).
#'
#' @return An object of class `multifer_linear_operator`.
#' @export
linear_operator <- function(apply,
                            apply_t,
                            dim,
                            gram = NULL,
                            core_update = NULL,
                            repr = NULL) {

  if (!is.function(apply)) {
    stop("`apply` must be a function.", call. = FALSE)
  }
  if (!is.function(apply_t)) {
    stop("`apply_t` must be a function.", call. = FALSE)
  }
  if (!is.numeric(dim) || length(dim) != 2L || any(!is.finite(dim)) ||
      any(dim != as.integer(dim)) || any(dim < 0L)) {
    stop("`dim` must be an integer vector of length 2 with non-negative values.",
         call. = FALSE)
  }
  dim <- as.integer(dim)
  if (!is.null(gram) && !is.function(gram)) {
    stop("`gram` must be NULL or a function.", call. = FALSE)
  }
  if (!is.null(core_update) && !is.function(core_update)) {
    stop("`core_update` must be NULL or a function.", call. = FALSE)
  }

  structure(
    list(
      apply       = apply,
      apply_t     = apply_t,
      dim         = dim,
      gram        = gram,
      core_update = core_update,
      repr        = repr
    ),
    class = "multifer_linear_operator"
  )
}

#' Test whether `x` is a `multifer_linear_operator`.
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_linear_operator <- function(x) {
  inherits(x, "multifer_linear_operator")
}

#' Coerce an object to a `multifer_linear_operator`.
#'
#' @param x Object to coerce.
#' @param ... Method-specific arguments.
#' @return A `multifer_linear_operator`.
#' @export
as_linear_operator <- function(x, ...) UseMethod("as_linear_operator")

#' @method as_linear_operator multifer_linear_operator
#' @rdname as_linear_operator
#' @export
as_linear_operator.multifer_linear_operator <- function(x, ...) x

#' @method as_linear_operator matrix
#' @rdname as_linear_operator
#' @export
as_linear_operator.matrix <- function(x, ...) {
  if (!is.numeric(x)) {
    stop("`as_linear_operator.matrix` requires a numeric matrix.", call. = FALSE)
  }
  n <- nrow(x); p <- ncol(x)

  apply_fn <- function(v) {
    if (is.vector(v)) v <- matrix(v, ncol = 1L)
    if (nrow(v) != p) {
      stop(sprintf("apply: input has %d rows; operator expects %d.", nrow(v), p),
           call. = FALSE)
    }
    x %*% v
  }
  apply_t_fn <- function(v) {
    if (is.vector(v)) v <- matrix(v, ncol = 1L)
    if (nrow(v) != n) {
      stop(sprintf("apply_t: input has %d rows; operator expects %d.", nrow(v), n),
           call. = FALSE)
    }
    base::crossprod(x, v)
  }
  gram_fn <- function() base::crossprod(x)

  linear_operator(
    apply      = apply_fn,
    apply_t    = apply_t_fn,
    dim        = c(n, p),
    gram       = gram_fn,
    repr       = list(matrix = x)
  )
}

#' @method dim multifer_linear_operator
#' @export
dim.multifer_linear_operator <- function(x) x$dim

#' Apply a linear operator to a matrix.
#'
#' @param op A `multifer_linear_operator`.
#' @param v A numeric matrix or vector compatible with `op`.
#' @return A numeric matrix.
#' @export
lo_apply <- function(op, v) {
  if (!is_linear_operator(op)) {
    stop("`op` must be a multifer_linear_operator.", call. = FALSE)
  }
  op$apply(v)
}

#' Apply the transposed linear operator.
#'
#' @param op A `multifer_linear_operator`.
#' @param v A numeric matrix or vector.
#' @return A numeric matrix.
#' @export
lo_apply_t <- function(op, v) {
  if (!is_linear_operator(op)) {
    stop("`op` must be a multifer_linear_operator.", call. = FALSE)
  }
  op$apply_t(v)
}

#' Gram matrix of a linear operator.
#'
#' Returns `t(op) %*% op` as a `p x p` numeric matrix. Uses the operator's
#' `gram` hook if present; otherwise falls back to
#' `apply_t(apply(I_p))`.
#'
#' @param op A `multifer_linear_operator`.
#' @return A `p x p` numeric matrix.
#' @export
lo_gram <- function(op) {
  if (!is_linear_operator(op)) {
    stop("`op` must be a multifer_linear_operator.", call. = FALSE)
  }
  if (!is.null(op$gram)) {
    return(op$gram())
  }
  p <- op$dim[2L]
  op$apply_t(op$apply(diag(1, p)))
}

#' Materialize a linear operator as a dense matrix.
#'
#' Applies the operator to the `p x p` identity. Intended for tests and
#' small operators only.
#'
#' @param x A `multifer_linear_operator`.
#' @param ... Ignored.
#' @method as.matrix multifer_linear_operator
#' @export
as.matrix.multifer_linear_operator <- function(x, ...) {
  if (!is.null(x$repr) && !is.null(x$repr$matrix)) {
    return(x$repr$matrix)
  }
  p <- x$dim[2L]
  x$apply(diag(1, p))
}

#' @method print multifer_linear_operator
#' @export
print.multifer_linear_operator <- function(x, ...) {
  cat("<multifer_linear_operator>\n")
  cat(sprintf("  dim:         %d x %d\n", x$dim[1L], x$dim[2L]))
  cat(sprintf("  gram hook:   %s\n", if (is.null(x$gram)) "no" else "yes"))
  cat(sprintf("  core_update: %s\n", if (is.null(x$core_update)) "no" else "yes"))
  invisible(x)
}
