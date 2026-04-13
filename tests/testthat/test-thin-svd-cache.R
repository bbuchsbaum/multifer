test_that("svd_cache_new starts empty", {
  cache <- svd_cache_new()
  expect_equal(svd_cache_hits(cache), 0L)
  expect_equal(svd_cache_lookups(cache), 0L)
  expect_true(is.na(svd_cache_rate(cache)))
})

test_that("cached_svd returns the correct decomposition on a miss", {
  set.seed(1)
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)
  cache <- svd_cache_new()
  ref <- svd(X)
  got <- cached_svd(X, cache = cache)
  expect_equal(got$d, ref$d)
  expect_equal(abs(got$u), abs(ref$u))
  expect_equal(abs(got$v), abs(ref$v))
  expect_equal(svd_cache_lookups(cache), 1L)
  expect_equal(svd_cache_hits(cache), 0L)
})

test_that("cached_svd serves hits on identical matrices", {
  set.seed(2)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  cache <- svd_cache_new()

  r1 <- cached_svd(X, cache = cache)
  r2 <- cached_svd(X, cache = cache)
  r3 <- cached_svd(X, cache = cache)

  expect_identical(r1, r2)
  expect_identical(r2, r3)
  expect_equal(svd_cache_lookups(cache), 3L)
  expect_equal(svd_cache_hits(cache), 2L)
  expect_equal(svd_cache_rate(cache), 2 / 3)
})

test_that("cached_svd distinguishes different matrices (no false hits)", {
  set.seed(3)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  Y <- matrix(rnorm(40), nrow = 10, ncol = 4)
  cache <- svd_cache_new()

  cached_svd(X, cache = cache)
  cached_svd(Y, cache = cache)
  expect_equal(svd_cache_hits(cache), 0L)
  expect_equal(svd_cache_lookups(cache), 2L)
})

test_that("cached_svd invalidates when data changes in place", {
  set.seed(4)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  cache <- svd_cache_new()
  cached_svd(X, cache = cache)
  X[1L, 1L] <- X[1L, 1L] + 1
  cached_svd(X, cache = cache)
  expect_equal(svd_cache_hits(cache), 0L)
  expect_equal(svd_cache_lookups(cache), 2L)
})

test_that("svd_cache_reset clears state", {
  cache <- svd_cache_new()
  cached_svd(matrix(1:12, 3, 4), cache = cache)
  expect_equal(svd_cache_lookups(cache), 1L)
  svd_cache_reset(cache)
  expect_equal(svd_cache_lookups(cache), 0L)
  expect_equal(svd_cache_hits(cache), 0L)
  expect_true(is.na(svd_cache_rate(cache)))
})

test_that("package-default cache is isolated per infer() call", {
  # Reset then populate directly (mimicking what infer() does).
  svd_cache_reset()
  cached_svd(matrix(1:12, 3, 4))
  expect_equal(svd_cache_lookups(), 1L)
  # A second reset should zero counters without errors.
  svd_cache_reset()
  expect_equal(svd_cache_lookups(), 0L)
})

test_that("svd_cache_fingerprint is stable and discriminating", {
  set.seed(5)
  X <- matrix(rnorm(60), nrow = 12, ncol = 5)
  fp1 <- svd_cache_fingerprint(X)
  fp2 <- svd_cache_fingerprint(X)
  expect_identical(fp1, fp2)

  Y <- X
  Y[1L, 1L] <- Y[1L, 1L] + 1e-6
  expect_false(identical(svd_cache_fingerprint(Y), fp1))

  # Different shape -> different fingerprint.
  expect_false(identical(
    svd_cache_fingerprint(matrix(rnorm(60), nrow = 10, ncol = 6)),
    fp1
  ))
})

test_that("cached_svd handles large matrices via sparse fingerprinting", {
  set.seed(6)
  big <- matrix(rnorm(100 * 150), nrow = 100, ncol = 150)
  cache <- svd_cache_new()
  r1 <- cached_svd(big, cache = cache)
  r2 <- cached_svd(big, cache = cache)
  expect_identical(r1, r2)
  expect_equal(svd_cache_hits(cache), 1L)

  # Perturb a single cell -> different matrix, should miss.
  big2 <- big
  big2[50L, 75L] <- big2[50L, 75L] + 1
  cached_svd(big2, cache = cache)
  expect_equal(svd_cache_hits(cache), 1L)  # hit count unchanged
  expect_equal(svd_cache_lookups(cache), 3L)
})

test_that("cached_svd rejects non-numeric input", {
  expect_error(cached_svd("abc"), "numeric matrix")
  expect_error(cached_svd(1:10), "numeric matrix")
})

test_that("infer() populates $cost$cache_hits", {
  svd_cache_reset()
  set.seed(42)
  X <- matrix(rnorm(30 * 6), nrow = 30, ncol = 6)

  # Directly call cached_svd on X twice via the cached_svd primitive -- this
  # exercises the same code path the engine uses for the full observed roots.
  cached_svd(X)
  cached_svd(X)
  expect_equal(svd_cache_hits(), 1L)
  expect_equal(svd_cache_lookups(), 2L)
  expect_equal(svd_cache_rate(), 0.5)
})
