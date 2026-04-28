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
      X <- do.call(cbind, data)
      Xc <- sweep(X, 2L, colMeans(X), "-")
      kk <- min(k, ncol(x$v))
      V <- x$v[, seq_len(kk), drop = FALSE]
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
      X <- do.call(cbind, data)
      Xc <- sweep(X, 2L, colMeans(X), "-")
      s2 <- svd(Xc)$d^2
      if (k > length(s2)) return(0)
      denom <- sum(s2[k:length(s2)])
      if (denom <= .Machine$double.eps) return(0)
      s2[k] / denom
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

test_that("infer dispatches multiblock data through adapter hooks", {
  set.seed(2001)
  adapter <- make_stub_multiblock_adapter()
  dat <- make_stub_multiblock_data()

  res <- infer(
    adapter = adapter,
    data = dat,
    geometry = "multiblock",
    relation = "variance",
    B = 9L,
    R = 3L,
    seed = 19L
  )

  expect_s3_class(res, "infer_result")
  expect_equal(res$provenance$adapter_id, "stub_multiblock")
  expect_true(all(names(dat) %in% unique(res$variable_stability$domain)))
  expect_true(all(names(dat) %in% unique(res$score_stability$domain)))
  expect_gt(nrow(res$component_tests), 0L)
})

test_that("bootstrap_fits uses arbitrary multiblock domains", {
  set.seed(2002)
  adapter <- make_stub_multiblock_adapter()
  dat <- make_stub_multiblock_data()
  rec <- infer_recipe(
    geometry = "multiblock",
    relation = "variance",
    adapter = adapter,
    targets = c("variable_stability", "score_stability", "subspace_stability")
  )
  fit <- adapter$refit(NULL, dat)
  units <- form_units(adapter$roots(fit))

  art <- bootstrap_fits(
    recipe = rec,
    adapter = adapter,
    data = dat,
    original_fit = fit,
    units = units,
    R = 2L,
    seed = 33L
  )

  expect_s3_class(art, "multifer_bootstrap_artifact")
  expect_equal(art$domains, names(dat))
  expect_equal(names(art$reps[[1L]]$aligned_loadings), names(dat))
})

test_that("multiblock geometry requires aligned matrix blocks", {
  adapter <- make_stub_multiblock_adapter()
  bad <- list(alpha = matrix(1, 3, 2), beta = matrix(1, 4, 2))

  expect_error(
    infer(
      adapter = adapter,
      data = bad,
      geometry = "multiblock",
      relation = "variance",
      targets = "component_significance",
      B = 3L,
      R = 0L
    ),
    "same row count"
  )
})
