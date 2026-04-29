# Test helper: re-register the adapters that .onLoad normally provides.
# Use this at the top of any test that needs to look up an adapter by id
# AFTER another test or test file has called clear_adapter_registry().
ensure_default_adapters <- function() {
  clear_adapter_registry()
  register_infer_adapter("svd_oneblock",      adapter_svd(),         overwrite = TRUE)
  register_infer_adapter("prcomp_oneblock",   adapter_prcomp(),      overwrite = TRUE)
  register_infer_adapter("cross_svd",         adapter_cross_svd(),   overwrite = TRUE)
  register_infer_adapter("cancor_cross",      adapter_cancor(),      overwrite = TRUE)
  register_infer_adapter("multivarious_pca",
                         adapter_multivarious_pca(),  overwrite = TRUE)
  register_infer_adapter("multivarious_plsc",
                         adapter_multivarious_plsc(), overwrite = TRUE)
  register_infer_adapter("multivarious_cca",
                         adapter_multivarious_cca(),  overwrite = TRUE)
  if (requireNamespace("MASS", quietly = TRUE)) {
    register_infer_adapter("lda_refit",
                           adapter_lda_refit(), overwrite = TRUE)
  }
  if (requireNamespace("pls", quietly = TRUE)) {
    register_infer_adapter("plsr_refit",
                           adapter_plsr_refit(), overwrite = TRUE)
  }
  invisible(NULL)
}
