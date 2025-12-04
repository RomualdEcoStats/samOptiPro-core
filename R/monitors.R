# monitors.R -- automatic monitor selection

#' @importFrom utils capture.output sessionInfo
NULL
#' @importFrom data.table :=
NULL

#' Internal helper to sanitize monitor root names
#'
#' @param x Character vector of candidate monitor roots.
#'
#' @return A character vector of cleaned, unique roots.
#' @keywords internal
.default_sanitize_roots <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[nzchar(x)]
  # vire tokens parasites courants
  bad <- grepl("^\\d+$", x) | x %in% c("thin", "=", ":", ",")
  unique(x[!bad])
}


#' Ensure monitors exist on the model (variable-aware)
#'
#' @param model nimbleModel
#' @param monitors character() de noms de variables (ex: "z","p","sd_proc","logLik")
#' @return sous-ensemble valide de `monitors`, avec warning soft si certains manquent
#' @export
#' @keywords internal
ensure_monitors_exist <- function(model, monitors) {
  monitors <- .default_sanitize_roots(monitors)
  if (!length(monitors)) return(character(0))

  all_nodes <- model$getNodeNames(stochOnly = FALSE, includeData = FALSE)
  all_vars  <- unique(sub("\\[.*\\]$", "", all_nodes))

  keep <- monitors %in% all_vars
  if (!all(keep)) {
    bad <- monitors[!keep]
    warning("Dropping non-existing monitor variables: ", paste(bad, collapse = ","))
  }
  unique(monitors[keep])
}


default_monitors <- function(model, opts = samOptiPro_options()) {
  include_data   <- isTRUE(opts$include_data)
  include_logLik <- isTRUE(opts$include_logLik)

  stoch_nodes <- model$getNodeNames(stochOnly = TRUE, includeData = include_data)
  stoch_vars  <- unique(sub("\\[.*\\]$", "", stoch_nodes))
  mons <- .default_sanitize_roots(stoch_vars)

  if (include_logLik) {
    all_nodes <- model$getNodeNames(stochOnly = FALSE, includeData = TRUE)
    all_vars  <- unique(sub("\\[.*\\]$", "", all_nodes))
    if ("logLik" %in% all_vars) mons <- unique(c(mons, "logLik"))
  }

  if (!is.null(opts$extra_monitors)) {
    mons <- unique(c(mons, .default_sanitize_roots(opts$extra_monitors)))
  }
  mons
}


