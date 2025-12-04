# zzz_param_helpers.R -- internal stubs for legacy helpers

#' @keywords internal
get_param_df <- function(...) {
  stop("get_param_df() is not implemented in this release of samOptiPro. ",
       "This helper was only used by internal plotting code. ",
       "Please refactor your workflow to use diagnostics_by_target().")
}

#' @keywords internal
derive_sampler_params <- function(...) {
  stop("derive_sampler_params() is not implemented in this release of samOptiPro. ",
       "This helper was only used by internal plotting code and has been ",
       "removed from the public API.")
}
