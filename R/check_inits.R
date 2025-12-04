
#' Robustly check initial values against a compiled NIMBLE model
#'
#' This helper tries to apply a list of initial values to a NIMBLE model and
#' fails early with an informative error if something is incompatible.
#'
#' @param model A NIMBLE model object, or a function that returns either
#'   a NIMBLE model or a list containing a `$model` component.
#' @param inits List of initial values, typically matching node/variable names
#'   used in the model.
#' @param silent Logical; if FALSE (default), a short diagnostic message is
#'   emitted before the error is thrown.
#'
#' @return Invisibly returns `TRUE` if the initial values are valid; otherwise
#'   throws an error describing the offending field(s).
#' @export
#' @keywords internal
checkInits <- function(model, inits, silent = FALSE) {
  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("Package 'nimble' is required for checkInits().")
  }

  ## Normalize 'model' input:
  ##  - if it's a builder, call it;
  ##  - if it returns a list with $model, use that;
  ##  - otherwise assume it is already a NIMBLE model.
  mod <- model
  if (is.function(model)) {
    mod <- model()
  }
  if (is.list(mod) && !is.null(mod$model)) {
    mod <- mod$model
  }

  ## Safety check: we expect an object with a setInits() method
  if (!("setInits" %in% names(mod))) {
    stop("checkInits(): provided object does not expose a setInits() method; ",
         "did you pass a valid nimbleModel (or its builder)?")
  }

  ok  <- TRUE
  err <- NULL

  tryCatch({
    ## IMPORTANT: use the model method, not nimble::setInits()
    mod$setInits(inits)
  }, error = function(e) {
    ok  <<- FALSE
    err <<- e
  })

  if (!ok) {
    if (!silent && exists("msg", mode = "function")) {
      ## msg() est ton helper interne; on ne l'appelle que s'il existe
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


