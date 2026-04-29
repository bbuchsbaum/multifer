make_stub_multiblock_data <- function(n = 28L) {
  z <- rnorm(n)
  list(
    alpha = cbind(z + rnorm(n, sd = 0.2), rnorm(n)),
    beta  = cbind(0.8 * z + rnorm(n, sd = 0.2), rnorm(n), rnorm(n)),
    gamma = cbind(-0.6 * z + rnorm(n, sd = 0.2), rnorm(n))
  )
}

make_stub_multiblock_adapter <- function(ncomp = 4L) {
  fit_blocks <- function(blocks) {
    nm <- names(blocks)
    if (is.null(nm) || any(!nzchar(nm))) {
      nm <- paste0("block", seq_along(blocks))
      names(blocks) <- nm
    }
    X <- do.call(cbind, blocks)
    Xc <- sweep(X, 2L, colMeans(X), "-")
    sv <- svd(Xc)
    k <- min(ncomp, length(sv$d), ncol(sv$v))
    lens <- vapply(blocks, ncol, integer(1L))
    ends <- cumsum(lens)
    starts <- ends - lens + 1L
    block_indices <- Map(seq, starts, ends)
    names(block_indices) <- nm
    structure(
      list(
        d = sv$d[seq_len(k)],
        v = sv$v[, seq_len(k), drop = FALSE],
        s = sv$u[, seq_len(k), drop = FALSE] %*% diag(sv$d[seq_len(k)], nrow = k),
        domains = nm,
        block_indices = block_indices
      ),
      class = "stub_multiblock_fit"
    )
  }

  split_blocks <- function(X, template) {
    lens <- vapply(template, ncol, integer(1L))
    ends <- cumsum(lens)
    starts <- ends - lens + 1L
    out <- Map(function(i, j) X[, seq.int(i, j), drop = FALSE], starts, ends)
    names(out) <- names(template)
    out
  }

  infer_adapter(
    adapter_id = "stub_multiblock",
    adapter_version = "0.0.1",
    shape_kinds = "multiblock",
    capabilities = capability_matrix(
      list(
        geometry = "multiblock",
        relation = "variance",
        targets = c(
          "component_significance",
          "variable_stability",
          "score_stability",
          "subspace_stability"
        )
      )
    ),
    domains = function(x, data = NULL, ...) {
      if (!is.null(x) && !is.null(x$domains)) return(x$domains)
      if (!is.null(names(data)) && all(nzchar(names(data)))) return(names(data))
      paste0("block", seq_along(data))
    },
    roots = function(x, ...) x$d^2,
    scores = function(x, domain = NULL, ...) x$s,
    loadings = function(x, domain = NULL, ...) {
      if (is.null(domain) || identical(domain, "global")) {
        return(x$v)
      }
      x$v[x$block_indices[[domain]], , drop = FALSE]
    },
    residualize = function(x, k, data, ...) {
      if (!identical(k, 1L)) {
        stop("stub multiblock residualize expects k = 1 on current residual.",
             call. = FALSE)
      }
      X <- do.call(cbind, data)
      Xc <- sweep(X, 2L, colMeans(X), "-")
      V <- x$v[, 1L, drop = FALSE]
      split_blocks(Xc - Xc %*% V %*% t(V), data)
    },
    refit = function(x, new_data, ...) fit_blocks(new_data),
    null_action = function(x, data, ...) {
      out <- lapply(data, function(block) {
        block[sample.int(nrow(block)), , drop = FALSE]
      })
      names(out) <- names(data)
      out
    },
    component_stat = function(x, data, k, ...) {
      if (!identical(k, 1L)) {
        stop("stub multiblock component_stat expects k = 1 on current residual.",
             call. = FALSE)
      }
      X <- do.call(cbind, data)
      Xc <- sweep(X, 2L, colMeans(X), "-")
      s2 <- svd(Xc)$d^2
      denom <- sum(s2)
      if (denom <= .Machine$double.eps) return(0)
      s2[1L] / denom
    },
    validity_level = "conditional",
    declared_assumptions = c("aligned_rows", "rows_exchangeable"),
    checked_assumptions = list(
      list(
        name = "aligned_multiblock_rows",
        check = function(data, ...) {
          is.list(data) &&
            all(vapply(data, is.matrix, logical(1L))) &&
            length(unique(vapply(data, nrow, integer(1L)))) == 1L
        },
        detail = "all blocks must be matrices with the same row count"
      )
    )
  )
}
