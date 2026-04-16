.fold_restore_seed <- function(old_seed) {
  if (!is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}

.fold_resolve_ids <- function(n, folds = NULL, n_folds = 5L, seed = NULL) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      n != as.integer(n) || n < 2L) {
    stop("`n` must be an integer >= 2.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.null(folds)) {
    if (length(folds) != n) {
      stop("`folds` must have length `n`.", call. = FALSE)
    }
    if (anyNA(folds)) {
      stop("`folds` must not contain missing values.", call. = FALSE)
    }
    fold_ids <- as.integer(as.factor(folds))
  } else {
    if (!is.numeric(n_folds) || length(n_folds) != 1L || is.na(n_folds) ||
        n_folds != as.integer(n_folds) || n_folds < 2L) {
      stop("`n_folds` must be an integer >= 2 when `folds` is NULL.",
           call. = FALSE)
    }
    n_folds <- as.integer(min(n_folds, n))

    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                             inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv)
      } else {
        NULL
      }
      on.exit(.fold_restore_seed(old_seed), add = TRUE)
      set.seed(as.integer(seed))
    }

    order <- sample.int(n)
    fold_ids <- integer(n)
    fold_ids[order] <- rep(seq_len(n_folds), length.out = n)
  }

  counts <- tabulate(fold_ids)
  if (length(counts) < 2L || any(counts == 0L)) {
    stop("`folds` must define at least 2 non-empty folds.", call. = FALSE)
  }
  if (any(counts >= n)) {
    stop("Each fold must leave at least one observation for training.",
         call. = FALSE)
  }

  fold_ids
}

.fold_payload_n <- function(data, what = "data") {
  if (is.matrix(data) || is.data.frame(data)) {
    return(nrow(data))
  }

  if (is.atomic(data) && length(dim(data)) <= 1L) {
    return(length(data))
  }

  if (is.list(data)) {
    counts <- vapply(
      seq_along(data),
      function(i) .fold_payload_n(data[[i]], what = sprintf("%s$%s",
                                                            what,
                                                            names(data)[i] %||% i)),
      integer(1L)
    )

    if (length(counts) == 0L) {
      stop(sprintf("`%s` must not be an empty list.", what), call. = FALSE)
    }

    if (length(unique(counts)) != 1L) {
      stop(sprintf("`%s` must be row-aligned across all elements.", what),
           call. = FALSE)
    }

    return(counts[[1L]])
  }

  stop(
    sprintf(
      "`%s` must be a row-aligned matrix, data frame, vector, or nested list.",
      what
    ),
    call. = FALSE
  )
}

.fold_subset_row_aligned <- function(data, idx) {
  payload_n <- function(x, what = "data") {
    if (is.matrix(x) || is.data.frame(x)) {
      return(nrow(x))
    }

    if (is.atomic(x) && length(dim(x)) <= 1L) {
      return(length(x))
    }

    if (is.list(x)) {
      counts <- vapply(
        seq_along(x),
        function(i) payload_n(x[[i]], what = sprintf("%s$%s",
                                                     what,
                                                     names(x)[i] %||% i)),
        integer(1L)
      )

      if (length(counts) == 0L) {
        stop(sprintf("`%s` must not be an empty list.", what), call. = FALSE)
      }

      if (length(unique(counts)) != 1L) {
        stop(sprintf("`%s` must be row-aligned across all elements.", what),
             call. = FALSE)
      }

      return(counts[[1L]])
    }

    stop(
      sprintf(
        "`%s` must be a row-aligned matrix, data frame, vector, or nested list.",
        what
      ),
      call. = FALSE
    )
  }

  subset_impl <- function(x) {
    payload_n(x)

    if (is.matrix(x) || is.data.frame(x)) {
      return(x[idx, , drop = FALSE])
    }

    if (is.atomic(x) && length(dim(x)) <= 1L) {
      return(x[idx])
    }

    if (is.list(x)) {
      out <- lapply(x, subset_impl)
      names(out) <- names(x)
      return(out)
    }

    stop("Unsupported row-aligned payload.", call. = FALSE)
  }

  subset_impl(data)
}

.fold_make_splits <- function(fold_ids) {
  if (length(fold_ids) < 2L || anyNA(fold_ids)) {
    stop("`fold_ids` must not contain missing values and must have length >= 2.",
         call. = FALSE)
  }
  counts <- tabulate(fold_ids)
  n <- length(fold_ids)

  if (length(counts) < 2L || any(counts == 0L)) {
    stop("`fold_ids` must define at least 2 non-empty folds.", call. = FALSE)
  }
  if (any(counts >= n)) {
    stop("Each fold must leave at least one observation for training.",
         call. = FALSE)
  }

  lapply(unique(fold_ids), function(fold) {
    list(
      fold = fold,
      train = which(fold_ids != fold),
      test = which(fold_ids == fold)
    )
  })
}

.fold_map <- function(data,
                      fold_ids,
                      fold_fun,
                      ...,
                      subset_data = .fold_subset_row_aligned,
                      backend = c("sequential", "mirai", "auto"),
                      seeds = NULL) {
  backend <- match.arg(backend)
  splits <- .fold_make_splits(fold_ids)
  extra_args <- list(...)

  worker <- function(split, data, subset_data, fold_fun, extra_args) {
    train_data <- subset_data(data, split$train)
    test_data <- subset_data(data, split$test)
    do.call(
      fold_fun,
      c(
        list(train_data = train_data, test_data = test_data, split = split),
        extra_args
      )
    )
  }

  multifer_parallel_lapply(
    X = splits,
    FUN = worker,
    data = data,
    subset_data = subset_data,
    fold_fun = fold_fun,
    extra_args = extra_args,
    seeds = seeds,
    backend = backend
  )
}
