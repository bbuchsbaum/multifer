.multifer_bench_registry <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  assign("oneblock_null", list(
    name      = "oneblock_null",
    geometry  = "oneblock",
    purpose   = "Null calibration for one-block shapes: over-selection rate and p-value calibration",
    locked    = TRUE,
    generator = bench_oneblock_null,
    targets   = bench_oneblock_null_targets
  ), envir = .multifer_bench_registry)

  assign("oneblock_shadowing", list(
    name      = "oneblock_shadowing",
    geometry  = "oneblock",
    purpose   = "Power under the Vitale shadowing case: strong early roots hide weaker late roots",
    locked    = TRUE,
    generator = bench_oneblock_shadowing,
    targets   = bench_oneblock_shadowing_targets
  ), envir = .multifer_bench_registry)

  assign("cross_null", list(
    name      = "cross_null",
    geometry  = "cross",
    purpose   = "Null calibration for cross-block shapes: zero cross-association with within-block structure",
    locked    = TRUE,
    generator = bench_cross_null,
    targets   = bench_cross_null_targets
  ), envir = .multifer_bench_registry)

  assign("speed_agreement", list(
    name      = "speed_agreement",
    geometry  = "cross",
    purpose   = "Wall-clock speedup and decision agreement: core-update vs. refit at medium and large scale",
    locked    = TRUE,
    generator = bench_speed_agreement,
    targets   = bench_speed_agreement_targets
  ), envir = .multifer_bench_registry)

  register_oneblock_baser_adapters()
  register_cross_baser_adapters()
  register_multivarious_pca_adapter()
  register_multivarious_plsc_adapter()
  register_lda_refit_adapter()
  register_plsr_refit_adapter()
}
