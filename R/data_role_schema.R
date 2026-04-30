#' Declare a role in an adapter data schema
#'
#' A role describes how one named element of an adapter-owned data payload
#' aligns to rows, columns, both axes, or neither. Role declarations are
#' compile-time metadata; they do not resample data by themselves.
#'
#' @param kind Character scalar. One of `"primary"`, `"constraint"`,
#'   `"metric"`, `"design_index"`, or `"static"`.
#' @param axes Character vector. Each element must be `"row"` or `"col"`.
#'   Use `character(0)` for static roles.
#' @param axis Optional singular alias for `axes`.
#' @param policy Optional character scalar for metric roles. One of
#'   `"identity"`, `"diagonal_observation_weights"`,
#'   `"diagonal_position_weights"`, `"full_spd_no_replacement_only"`, or
#'   `"adapter_owned"`.
#'
#' @return A `multifer_data_role` object.
#' @export
role <- function(kind,
                 axes = character(0),
                 axis = NULL,
                 policy = NULL) {
  valid_kinds <- c("primary", "constraint", "metric", "design_index", "static")
  if (!is.character(kind) || length(kind) != 1L || is.na(kind) ||
      !(kind %in% valid_kinds)) {
    stop(sprintf("`kind` must be one of: %s.",
                 paste(valid_kinds, collapse = ", ")), call. = FALSE)
  }
  if (!is.null(axis)) {
    if (!identical(axes, character(0))) {
      stop("Use either `axes` or `axis`, not both.", call. = FALSE)
    }
    axes <- axis
  }
  if (!is.character(axes) || anyNA(axes)) {
    stop("`axes` must be a character vector.", call. = FALSE)
  }
  axes <- unique(axes)
  bad_axes <- setdiff(axes, c("row", "col"))
  if (length(bad_axes) > 0L) {
    stop(sprintf("unknown role axis/axes: %s.",
                 paste(bad_axes, collapse = ", ")), call. = FALSE)
  }
  if (identical(kind, "static") && length(axes) > 0L) {
    stop("Static roles must use `axes = character(0)`.", call. = FALSE)
  }
  valid_policies <- c("identity", "diagonal_observation_weights",
                      "diagonal_position_weights",
                      "full_spd_no_replacement_only", "adapter_owned")
  if (identical(kind, "metric")) {
    if (is.null(policy)) {
      policy <- "identity"
    }
    if (!is.character(policy) || length(policy) != 1L || is.na(policy) ||
        !(policy %in% valid_policies)) {
      stop(sprintf("metric `policy` must be one of: %s.",
                   paste(valid_policies, collapse = ", ")), call. = FALSE)
    }
  } else if (!is.null(policy)) {
    stop("`policy` is only valid for metric roles.", call. = FALSE)
  }

  structure(
    list(kind = kind, axes = axes, policy = policy),
    class = "multifer_data_role"
  )
}

#' Test whether an object is a data role
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_data_role <- function(x) inherits(x, "multifer_data_role")

#' Declare an adapter data-role schema
#'
#' `data_role_schema()` names the roles present in an adapter's data payload.
#' It is the declarative layer that later bootstrap/null planners can resolve
#' against actual runtime data.
#'
#' @param ... Named `role()` objects.
#'
#' @return A named list with class `multifer_data_role_schema`.
#' @export
data_role_schema <- function(...) {
  roles <- list(...)
  if (length(roles) == 0L) {
    stop("`data_role_schema()` requires at least one named role.",
         call. = FALSE)
  }
  nms <- names(roles)
  if (is.null(nms) || any(!nzchar(nms)) || anyNA(nms)) {
    stop("all roles in `data_role_schema()` must be named.", call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("role names in `data_role_schema()` must be unique.", call. = FALSE)
  }
  bad <- nms[!vapply(roles, is_data_role, logical(1L))]
  if (length(bad) > 0L) {
    stop(sprintf("schema entries must be built with `role()`: %s.",
                 paste(bad, collapse = ", ")), call. = FALSE)
  }

  structure(roles, class = "multifer_data_role_schema")
}

#' Test whether an object is a data-role schema
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_data_role_schema <- function(x) inherits(x, "multifer_data_role_schema")

#' Declare null-action semantics
#'
#' `null_spec()` is descriptive provenance metadata for an adapter's
#' operational `null_action()`. It does not replace the null callback.
#'
#' @param kind Character scalar. One of `"within_role_randomization"`,
#'   `"between_role_decoupling"`, or `"parametric_draw"`.
#' @param randomized_role Character scalar naming the role randomized by the
#'   null, or `NA` when not role-specific.
#' @param reference_role Optional character scalar naming the reference role.
#' @param axis Optional character scalar, `"row"` or `"col"`.
#' @param preserves Character vector describing invariants preserved by the
#'   null action.
#'
#' @return A `multifer_null_spec` object.
#' @export
null_spec <- function(kind,
                      randomized_role = NA_character_,
                      reference_role = NA_character_,
                      axis = NA_character_,
                      preserves = character(0)) {
  valid_kinds <- c("within_role_randomization", "between_role_decoupling",
                   "parametric_draw")
  if (!is.character(kind) || length(kind) != 1L || is.na(kind) ||
      !(kind %in% valid_kinds)) {
    stop(sprintf("`kind` must be one of: %s.",
                 paste(valid_kinds, collapse = ", ")), call. = FALSE)
  }
  scalar_chr <- function(x, name, allow_na = TRUE) {
    if (!is.character(x) || length(x) != 1L ||
        (!allow_na && is.na(x))) {
      stop(sprintf("`%s` must be a character scalar.", name), call. = FALSE)
    }
    x
  }
  randomized_role <- scalar_chr(randomized_role, "randomized_role")
  reference_role <- scalar_chr(reference_role, "reference_role")
  axis <- scalar_chr(axis, "axis")
  if (!is.na(axis) && !(axis %in% c("row", "col"))) {
    stop("`axis` must be \"row\", \"col\", or NA.", call. = FALSE)
  }
  if (!is.character(preserves) || anyNA(preserves)) {
    stop("`preserves` must be a character vector.", call. = FALSE)
  }

  structure(
    list(
      kind = kind,
      randomized_role = randomized_role,
      reference_role = reference_role,
      axis = axis,
      preserves = preserves
    ),
    class = "multifer_null_spec"
  )
}

#' Test whether an object is a null specification
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_null_spec <- function(x) inherits(x, "multifer_null_spec")
