# v1 boundary drift-guard.
#
# This test fails CI if the package's v1 core boundary -- what the
# notes, wrappers, and evidence all agree on -- starts drifting.
# Concretely it asserts that DESCRIPTION dependencies, the README
# adapter table, the live adapter registry, and the infer wrapper
# defaults all describe the same package. Adding a new mature adapter,
# demoting an adapter to narrow, or moving multivarious out of Imports
# should all fail here until every downstream artifact is updated in
# the same commit.

pkg_root <- function() {
  # testthat::test_path() is relative to tests/testthat/, so going
  # up two levels lands at the package root.
  testthat::test_path("..", "..")
}

parse_readme_adapter_table <- function(readme_path) {
  lines <- readLines(readme_path, warn = FALSE)
  # Find the start of the adapter table: the header row contains
  # "Adapter" and "Maturity" columns.
  header_idx <- which(
    grepl("\\| Adapter", lines) & grepl("Maturity", lines)
  )
  if (length(header_idx) == 0L) {
    stop("README adapter table header not found.")
  }
  header <- header_idx[[1L]]
  body_start <- header + 2L  # skip separator row
  body <- lines[seq(body_start, length(lines))]
  body <- body[grepl("^\\|", body)]
  body <- body[!grepl("^\\| ?---", body)]
  # Cut off where the table ends (first line not starting with |)
  body <- body[seq_len(min(which(!nzchar(body)),
                           length(body) + 1L) - 1L)]

  rows <- lapply(body, function(ln) {
    cells <- strsplit(ln, "\\|", fixed = FALSE)[[1L]]
    cells <- trimws(cells)
    cells <- cells[nzchar(cells)]
    if (length(cells) < 4L) return(NULL)
    list(
      adapter_id = gsub("`", "", gsub("\\s*\\(.*\\)", "", cells[[1L]])),
      geometry   = cells[[2L]],
      relation   = cells[[3L]],
      maturity   = gsub("\\*", "", cells[[4L]])
    )
  })
  rows <- Filter(Negate(is.null), rows)
  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}

test_that("DESCRIPTION declares multivarious as a hard dependency", {
  desc_path <- file.path(pkg_root(), "DESCRIPTION")
  skip_if_not(file.exists(desc_path), "DESCRIPTION not found from test cwd")
  desc <- read.dcf(desc_path)
  imports <- if ("Imports" %in% colnames(desc)) desc[1L, "Imports"] else ""
  suggests <- if ("Suggests" %in% colnames(desc)) desc[1L, "Suggests"] else ""
  expect_true(grepl("multivarious", imports, fixed = TRUE))
  expect_false(grepl("multivarious", suggests, fixed = TRUE))
})

test_that("README adapter table is consistent with the live registry", {
  readme_path <- file.path(pkg_root(), "README.md")
  skip_if_not(file.exists(readme_path), "README.md not found from test cwd")

  ensure_default_adapters()
  readme <- parse_readme_adapter_table(readme_path)
  registry <- list_infer_adapters(details = TRUE)

  # Every registered adapter appears in README with a matching
  # maturity label (adapter id match modulo relation suffix like
  # "cross_svd (cov)").
  expect_true(nrow(readme) > 0L)
  expect_true(nrow(registry) > 0L)

  readme_adapters <- unique(readme$adapter_id)
  registry_adapters <- unique(registry$adapter_id)

  missing_from_readme <- setdiff(registry_adapters, readme_adapters)
  expect_length(missing_from_readme, 0L)

  missing_from_registry <- setdiff(readme_adapters, registry_adapters)
  expect_length(missing_from_registry, 0L)

  for (id in registry_adapters) {
    reg_maturity <- unique(registry$maturity[registry$adapter_id == id])
    readme_maturity <- unique(readme$maturity[readme$adapter_id == id])
    expect_length(reg_maturity, 1L)
    expect_true(all(readme_maturity == reg_maturity),
                info = sprintf("README maturity for %s: %s, registry: %s",
                               id,
                               paste(readme_maturity, collapse = ","),
                               reg_maturity))
  }
})

test_that("README opener is byte-equal to the package_vision positioning paragraph", {
  readme_path <- file.path(pkg_root(), "README.md")
  vision_path <- file.path(pkg_root(), "notes", "package_vision.md")
  skip_if_not(file.exists(readme_path))
  skip_if_not(file.exists(vision_path))

  readme_lines <- readLines(readme_path, warn = FALSE)
  vision_lines <- readLines(vision_path, warn = FALSE)

  # Find the vision paragraph under '## Recommended positioning'.
  vision_header <- which(vision_lines == "## Recommended positioning")
  skip_if(length(vision_header) == 0L,
          "'Recommended positioning' header not found in package_vision.md")
  vision_tail <- vision_lines[seq(vision_header[[1L]] + 1L,
                                  length(vision_lines))]
  vision_para_lines <- vision_tail[grepl("^\\*\\*.+\\*\\*$", vision_tail)]
  skip_if(length(vision_para_lines) == 0L,
          "positioning paragraph not found in package_vision.md")
  vision_para <- vision_para_lines[[1L]]

  # The README opener must contain that exact paragraph verbatim.
  readme_has_it <- any(readme_lines == vision_para)
  expect_true(
    readme_has_it,
    info = paste0(
      "README opener does not contain the verbatim positioning ",
      "paragraph from notes/package_vision.md."
    )
  )
})

test_that("README adapter table uses only canonical maturity labels", {
  readme_path <- file.path(pkg_root(), "README.md")
  skip_if_not(file.exists(readme_path))
  readme <- parse_readme_adapter_table(readme_path)
  expect_true(all(readme$maturity %in% multifer_maturity_levels()))
})

test_that("infer_pca default adapter is explicit and multivarious-backed", {
  ensure_default_adapters()

  set.seed(9001)
  X <- matrix(rnorm(120), 30, 4)
  res <- infer_pca(X, B = 29L, R = 4L, seed = 5L)
  expect_equal(res$provenance$adapter_id, "multivarious_pca")

  explicit <- infer_pca(
    X,
    adapter = "prcomp_oneblock",
    B = 29L,
    R = 4L,
    seed = 6L
  )
  expect_equal(explicit$provenance$adapter_id, "prcomp_oneblock")
})

test_that("infer_plsc default adapter is explicit and multivarious-backed", {
  ensure_default_adapters()

  set.seed(9003)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  res <- infer_plsc(X, Y, B = 29L, R = 4L, seed = 7L)
  expect_equal(res$provenance$adapter_id, "multivarious_plsc")

  explicit <- infer_plsc(
    X,
    Y,
    adapter = "cross_svd",
    B = 29L,
    R = 4L,
    seed = 9L
  )
  expect_equal(explicit$provenance$adapter_id, "cross_svd")
})

test_that("infer_cca still defaults to cancor_cross", {
  ensure_default_adapters()
  set.seed(9004)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  res <- infer_cca(X, Y, B = 29L, R = 4L, seed = 8L)
  expect_equal(res$provenance$adapter_id, "cancor_cross")
})
