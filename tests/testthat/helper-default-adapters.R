# Test helper: re-register the adapters that .onLoad normally provides.
# Use this at the top of any test that needs to look up an adapter by id
# AFTER another test or test file has called clear_adapter_registry().
ensure_default_adapters <- function() {
  clear_adapter_registry()
  register_infer_adapter("svd_oneblock",      adapter_svd(),         overwrite = TRUE)
  register_infer_adapter("prcomp_oneblock",   adapter_prcomp(),      overwrite = TRUE)
  register_infer_adapter("cross_svd",         adapter_cross_svd(),   overwrite = TRUE)
  register_infer_adapter("cancor_cross",      adapter_cancor(),      overwrite = TRUE)
  if (requireNamespace("multivarious", quietly = TRUE)) {
    register_infer_adapter("multivarious_pca",
                           adapter_multivarious_pca(),  overwrite = TRUE)
    register_infer_adapter("multivarious_plsc",
                           adapter_multivarious_plsc(), overwrite = TRUE)
  }
  invisible(NULL)
}
