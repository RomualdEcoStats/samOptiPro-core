
#' Robustly check initial values against a compiled NIMBLE model
#'
#' @param model a NIMBLE model object (or function that builds one)
#' @param inits list of initial values
#' @param silent logical
#' @return logical TRUE if valid; otherwise throws informative error
#' @export
#' @keywords internal
checkInits <- function(model, inits, silent = FALSE) {
  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("Package 'nimble' is required for checkInits().")
  }
  mod <- model
  if (is.function(model)) mod <- model()
  ok <- TRUE
  err <- NULL
  tryCatch({
# TODO(samOptiPro): remplacer nimble::setInits(...) par <modele>$setInits(...) (appel méthode d'objet)
    nimble::setInits(mod, inits)
  }, error = function(e) {
    ok <<- FALSE; err <<- e
  })
  if (!ok) {
    if (!silent) {
      msg("Invalid inits: ", conditionMessage(err))
    }
    stop(err)
  }
  invisible(TRUE)
}

#' Check inits then run MCMC
#'
#' @param run_fn function with signature function(inits) -> list(samples=samples, runtime_s=seconds)
#' @param inits list (or list of lists for chains)
#' @param ... forwarded to run_fn
#' @return list(samples, runtime_s)
#' @export
#' @keywords internal
checkInitsAndRun <- function(run_fn, inits, ...) {
  if (is.list(inits) && length(inits) > 0 && is.list(inits[[1]])) {
    # multiple chains inits
    for (i in seq_along(inits)) {
      run_fn(inits[[i]], ...)
    }
  }
  run_fn(inits, ...)
}
