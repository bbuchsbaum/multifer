#' Match subspace components by maximum absolute inner product
#'
#' Finds an integer permutation \code{p} such that \code{Vb[, p]} is the
#' best match to \code{Vref} in terms of summed absolute inner products.
#'
#' @details
#' `method = "auto"` uses \pkg{clue}'s Hungarian implementation when
#' installed. Without \pkg{clue}, it uses an exact exhaustive assignment for
#' small component counts and falls back to the historical greedy matcher for
#' larger problems. The return value carries diagnostics as attributes:
#' `match_score`, `match_margin`, `ambiguous_match`, and `match_method`.
#'
#' @param Vref Numeric matrix. Reference loading matrix, n x k.
#' @param Vb Numeric matrix. Bootstrap/replicate loading matrix, n x k.
#'   Must have the same dimensions as \code{Vref}.
#' @param method Character scalar. One of `"auto"`, `"greedy"`, or
#'   `"hungarian"`.
#' @param metric Character scalar. Currently `"abs_inner_product"`.
#' @param ambiguity_tol Numeric tolerance for declaring an ambiguous match.
#' @param diagnostics Logical; when `TRUE`, attach matching diagnostics as
#'   attributes to the returned permutation.
#'
#' @return An integer vector of length \code{ncol(Vref)} giving the
#'   permutation \code{p} such that \code{Vb[, p]} best matches \code{Vref}.
#'
#' @export
match_components <- function(Vref,
                             Vb,
                             method = c("auto", "greedy", "hungarian"),
                             metric = c("abs_inner_product", "cosine"),
                             ambiguity_tol = 1e-8,
                             diagnostics = FALSE) {
  method <- match.arg(method)
  metric <- match.arg(metric)
  if (!is.matrix(Vref) || !is.matrix(Vb)) {
    stop("`Vref` and `Vb` must both be matrices.", call. = FALSE)
  }
  if (!identical(dim(Vref), dim(Vb))) {
    stop("`Vref` and `Vb` must have identical dimensions.", call. = FALSE)
  }
  k <- ncol(Vref)
  if (k == 0L) {
    return(integer(0L))
  }

  # absolute inner product matrix: entry [i, j] = |Vref[,i] . Vb[,j]|
  ip <- abs(t(Vref) %*% Vb)   # k x k
  if (metric == "cosine") {
    ref_norm <- sqrt(colSums(Vref * Vref))
    rep_norm <- sqrt(colSums(Vb * Vb))
    denom <- outer(ref_norm, rep_norm, `*`)
    denom[denom == 0] <- NA_real_
    ip <- ip / denom
    ip[!is.finite(ip)] <- 0
  }

  if (method == "greedy") {
    perm <- .match_components_greedy(ip)
    actual_method <- "greedy"
  } else {
    optimal <- .match_components_optimal(ip)
    if (is.null(optimal)) {
      if (method == "hungarian") {
        stop(
          "`method = \"hungarian\"` requires package `clue` or k <= 8 for the base-R exact fallback.",
          call. = FALSE
        )
      }
      perm <- .match_components_greedy(ip)
      actual_method <- "greedy"
    } else {
      perm <- optimal$perm
      actual_method <- optimal$method
    }
  }

  if (isTRUE(diagnostics)) {
    diag <- .match_components_diagnostics(ip, perm, actual_method,
                                          ambiguity_tol)
    attr(perm, "match_score") <- diag$match_score
    attr(perm, "match_margin") <- diag$match_margin
    attr(perm, "ambiguous_match") <- diag$ambiguous_match
    attr(perm, "match_method") <- diag$match_method
  }
  perm
}

.match_components_greedy <- function(score) {
  k <- nrow(score)
  perm <- integer(k)
  used <- logical(k)
  for (i in seq_len(k)) {
    row_scores <- score[i, ]
    row_scores[used] <- -Inf
    j <- which.max(row_scores)
    perm[i] <- j
    used[j] <- TRUE
  }
  perm
}

.match_components_optimal <- function(score) {
  k <- nrow(score)
  if (requireNamespace("clue", quietly = TRUE)) {
    return(list(
      perm = as.integer(clue::solve_LSAP(score, maximum = TRUE)),
      method = "hungarian"
    ))
  }
  if (k > 8L) {
    return(NULL)
  }
  best_score <- -Inf
  best_perm <- seq_len(k)
  used <- logical(k)
  current <- integer(k)
  search <- function(i, total) {
    if (i > k) {
      if (total > best_score) {
        best_score <<- total
        best_perm <<- current
      }
      return(invisible(NULL))
    }
    for (j in seq_len(k)) {
      if (!used[j]) {
        used[j] <<- TRUE
        current[i] <<- j
        search(i + 1L, total + score[i, j])
        used[j] <<- FALSE
      }
    }
    invisible(NULL)
  }
  search(1L, 0)
  list(perm = best_perm, method = "exhaustive")
}

.match_components_diagnostics <- function(score, perm, method, ambiguity_tol) {
  k <- nrow(score)
  assigned <- score[cbind(seq_len(k), perm)]
  margin <- rep(Inf, k)
  for (i in seq_len(k)) {
    alternatives <- score[i, -perm[i], drop = TRUE]
    if (length(alternatives) > 0L) {
      margin[i] <- assigned[i] - max(alternatives)
    }
  }
  list(
    match_score = sum(assigned),
    match_margin = min(margin),
    ambiguous_match = is.finite(min(margin)) && min(margin) <= ambiguity_tol,
    match_method = method
  )
}


#' Sign-align columns of a matched loading matrix
#'
#' For each column \code{j}, flips the sign of \code{Vb_matched[, j]} if
#' the inner product with \code{Vref[, j]} is negative, so that
#' \code{diag(t(Vref) \%*\% aligned)} is everywhere non-negative.
#'
#' @param Vref Numeric matrix. Reference loading matrix, n x k.
#' @param Vb_matched Numeric matrix. Loading matrix already permuted to match
#'   \code{Vref} (e.g. output of \code{Vb[, match_components(Vref, Vb)]}).
#'   Must have the same dimensions as \code{Vref}.
#'
#' @return A numeric matrix of the same dimensions as \code{Vb_matched} with
#'   sign-corrected columns.
#'
#' @export
align_sign <- function(Vref, Vb_matched) {
  if (!is.matrix(Vref) || !is.matrix(Vb_matched)) {
    stop("`Vref` and `Vb_matched` must both be matrices.", call. = FALSE)
  }
  if (!identical(dim(Vref), dim(Vb_matched))) {
    stop("`Vref` and `Vb_matched` must have identical dimensions.", call. = FALSE)
  }

  signs <- sign(colSums(Vref * Vb_matched))
  # colSums(Vref * Vb_matched) == diag(t(Vref) %*% Vb_matched), cheaper
  # replace any zero inner product with +1 (arbitrary, no flip needed)
  signs[signs == 0] <- 1L

  Vb_matched %*% diag(signs, nrow = length(signs))
}


#' Procrustes-align a loading matrix to a reference
#'
#' Rotates \code{Vb_matched} by the closest orthogonal matrix Q to
#' \code{Vref}. Given the cross-product \code{M = t(Vref) \%*\% Vb_matched}
#' with SVD \code{M = U D V^t}, the optimal Q is \code{V \%*\% t(U)} and
#' the aligned matrix is \code{Vb_matched \%*\% V \%*\% t(U)}.
#'
#' This helper is retained for backward compatibility with legacy
#' bootstrap workflows. It is not the recommended default for inferential
#' use in \pkg{multifer}; prefer permutation matching plus sign correction
#' for separated components, and subspace summaries when roots are tied or
#' nearly tied.
#'
#' @param Vref Numeric matrix. Reference loading matrix, n x k.
#' @param Vb_matched Numeric matrix. Loading matrix already permuted to match
#'   \code{Vref}. Must have the same dimensions as \code{Vref}.
#'
#' @return A numeric matrix of the same dimensions as \code{Vb_matched}
#'   rotated to be as close as possible to \code{Vref} under orthogonal
#'   Procrustes.
#'
#' @export
align_procrustes <- function(Vref, Vb_matched) {
  if (!is.matrix(Vref) || !is.matrix(Vb_matched)) {
    stop("`Vref` and `Vb_matched` must both be matrices.", call. = FALSE)
  }
  if (!identical(dim(Vref), dim(Vb_matched))) {
    stop("`Vref` and `Vb_matched` must have identical dimensions.", call. = FALSE)
  }

  M   <- t(Vref) %*% Vb_matched
  sv  <- svd(M)
  # optimal orthogonal Q = V %*% t(U)
  Q   <- sv$v %*% t(sv$u)
  Vb_matched %*% Q
}


#' Principal angles between two subspaces
#'
#' Computes the principal (canonical) angles between the subspaces spanned by
#' the columns of \code{A} and \code{B}.
#'
#' @details
#' When \pkg{multivarious} is installed, this function delegates to
#' \code{multivarious::prinang} (which returns the same quantity). Otherwise
#' the angles are computed from first principles: orthonormalise both bases
#' via \code{qr.Q}, form the cross-Gram matrix, and take \code{acos} of the
#' clamped singular values.
#'
#' Both \code{A} and \code{B} should have full column rank. If not, the result
#' is numerically unreliable (a warning is not issued here; callers are
#' responsible for checking rank).
#'
#' @param A Numeric matrix whose columns span the first subspace.
#' @param B Numeric matrix whose columns span the second subspace.
#'
#' @return A numeric vector of principal angles in radians, sorted ascending,
#'   of length \code{min(ncol(A), ncol(B))}.
#'
#' @export
principal_angles <- function(A, B) {
  if (!is.matrix(A) || !is.matrix(B)) {
    stop("`A` and `B` must both be matrices.", call. = FALSE)
  }
  if (nrow(A) != nrow(B)) {
    stop("`A` and `B` must have the same number of rows.", call. = FALSE)
  }
  if (ncol(A) == 0L || ncol(B) == 0L) {
    return(numeric(0L))
  }

  if (requireNamespace("multivarious", quietly = TRUE) &&
      !is.null(getExportedValue <- tryCatch(
        { getExportedValue("multivarious", "prinang"); TRUE },
        error = function(e) FALSE
      )) &&
      isTRUE(getExportedValue)) {
    return(sort(multivarious::prinang(A, B)))
  }

  # fallback: compute from first principles
  Qa <- qr.Q(qr(A))
  Qb <- qr.Q(qr(B))

  # singular values of the cross-Gram matrix are the cosines of the angles
  cosines <- svd(t(Qa) %*% Qb)$d

  # clamp to [-1, 1] to guard against numerical noise before acos
  cosines <- pmin(pmax(cosines, -1), 1)

  sort(acos(cosines))
}


#' Align loadings: match components then apply sign or legacy Procrustes alignment
#'
#' Convenience pipeline. First calls \code{\link{match_components}} to find
#' the best column permutation, then applies either sign alignment
#' (\code{\link{align_sign}}) or orthogonal Procrustes alignment
#' (\code{\link{align_procrustes}}). The recommended default is
#' \code{method = "sign"}; \code{"procrustes"} is retained only for
#' backward compatibility with legacy workflows.
#'
#' CRITICAL ordering: \code{match_components} is always called before any
#' sign or rotation step. Bootstrap replicates that swap neighboring
#' components would otherwise poison downstream stability summaries if sign
#' alignment were applied without first resolving the permutation.
#'
#' @param Vref Numeric matrix. Reference loading matrix, n x k.
#' @param Vb Numeric matrix. Bootstrap/replicate loading matrix, n x k.
#' @param method Character scalar. One of \code{"sign"} (default) or
#'   \code{"procrustes"} (legacy-only).
#'
#' @return A numeric matrix of the same dimensions as \code{Vb} with columns
#'   permuted and then sign- or legacy-Procrustes-aligned to \code{Vref}.
#'
#' @export
align_loadings <- function(Vref, Vb, method = c("sign", "procrustes")) {
  method <- match.arg(method)

  if (!is.matrix(Vref) || !is.matrix(Vb)) {
    stop("`Vref` and `Vb` must both be matrices.", call. = FALSE)
  }
  if (!identical(dim(Vref), dim(Vb))) {
    stop("`Vref` and `Vb` must have identical dimensions.", call. = FALSE)
  }

  # Step 1: resolve column permutation
  p         <- match_components(Vref, Vb)
  Vb_perm   <- Vb[, p, drop = FALSE]

  # Step 2: sign or Procrustes alignment
  if (method == "sign") {
    align_sign(Vref, Vb_perm)
  } else {
    align_procrustes(Vref, Vb_perm)
  }
}
