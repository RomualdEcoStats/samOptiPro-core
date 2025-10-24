# monitors.R -- automatic monitor selection

#' Discover default monitors from a nimble model (variable-level)
#'
#' Retourne des NOMS DE VARIABLES (sans indices) pour laisser le runner
#' poser ensuite les monitors correctement via conf$setMonitors / setMonitors2.
#'
#' @param model nimbleModel
#' @param opts list from samOptiPro_options(); champs supportes :
#'   - include_data      : logical, inclure les data nodes dans l'exploration (rarement utile)
#'   - include_logLik    : logical, inclure "logLik" dans la selection par defaut
#'   - extra_monitors    : character(), variables supplementaires a monitorer (samples)
#'   - extra_monitors2   : character(), variables supplementaires a monitorer en monitors2 (samples2)
#' @return character() de noms de variables sans indices (ex: "sd_proc","z","p","logLik")
#' @export
.default_sanitize_roots <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[nzchar(x)]
  # vire tokens parasites courants
  bad <- grepl("^\\d+$", x) | x %in% c("thin", "=", ":", ",")
  unique(x[!bad])
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


#' Ensure monitors exist on the model (variable-aware)
#'
#' Valide une liste de monitors donnee au NIVEAU VARIABLE (sans indices).
#' Ne "coupe" plus selon les noeuds indexes ; on verifie juste que la VARIABLE existe.
#'
#' @param model nimbleModel
#' @param monitors character() de noms de variables (ex: "z","p","sd_proc","logLik")
#' @return sous-ensemble valide de `monitors`, avec warning soft si certains manquent
#' @export
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


