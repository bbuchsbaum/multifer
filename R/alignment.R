#' Match subspace components by maximum absolute inner product
#'
#' Finds an integer permutation \code{p} such that \code{Vb[, p]} is the
#' best match to \code{Vref} in terms of summed absolute inner products.
#'
#' @details
#' Implementation note: the Hungarian algorithm (optimal assignment) is not
#' available in base R. This function uses a greedy fallback: for each column
#' of \code{Vref} in order, the unused column of \code{Vb} with the largest
#' absolute inner product is assigned. The greedy solution is suboptimal in
#' pathological cases (e.g. two nearly equal scores) but is acceptable for
#' v1 where components are assumed to be reasonably well separated. A future
#' version may add optional Hungarian via \pkg{clue} or \pkg{RcppHungarian}.
#'
#' @param Vref Numeric matrix. Reference loading matrix, n x k.
#' @param Vb Numeric matrix. Bootstrap/replicate loading matrix, n x k.
#'   Must have the same dimensions as \code{Vref}.
#'
#' @return An integer vector of length \code{ncol(Vref)} giving the
#'   permutation \code{p} such that \code{Vb[, p]} best matches \code{Vref}.
#'
#' @export
match_components <- function(Vref, Vb) {
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

  perm    <- integer(k)
  used    <- logical(k)

  for (i in seq_len(k)) {
    row_scores <- ip[i, ]
    row_scores[used] <- -Inf
    j        <- which.max(row_scores)
    perm[i]  <- j
    used[j]  <- TRUE
  }

  perm
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
