# runners.R -- robust generic runners
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

## ---- Exported functions (4) ----

#' Build a fresh MCMC configuration with automatic monitors if missing
#' @param model nimbleModel (ou compile; on accepte les deux)
#' @param monitors optional character vector; if NULL, auto-discover
#' @param opts samOptiPro_options()
#' @export
#' @keywords internal
build_conf_with_monitors <- function(model, monitors = NULL, opts = samOptiPro_options()) {
  .configure_with_monitors(model, monitors = monitors, opts = opts)
}


#' Run baseline et retourne un mcmc.list fusionne (samples + samples2)
#' @export
#' @keywords internal
run_baseline_coda <- function(build_fn, niter,
                              nburnin = floor(0.25 * niter),
                              thin = 1L,
                              monitors = NULL,
                              nchains = 1L,
                              drop_loglik = FALSE,
                              opts = samOptiPro_options()) {
  rb <- run_baseline_config(build_fn, niter = niter, nburnin = nburnin, thin = thin,
                            monitors = monitors, nchains = nchains, opts = opts)
  ml <- as_mcmc_list_sop(rb$samples, rb$samples2, drop_loglik = drop_loglik, thin = thin)
  list(samples = ml, runtime_s = rb$runtime_s, conf = rb$conf)
}


#' Run a baseline RW/slice configuration (sequential)
#'
#' @param build_fn Function with no arguments that returns a list with at
#'   least \code{model}, \code{cmodel} and \code{conf}, or a builder that
#'   is passed to an internal fresh-build wrapper.
#' @param niter Total number of MCMC iterations.
#' @param nburnin Number of initial iterations discarded as burn-in.
#' @param thin Thinning interval (keep 1 draw every \code{thin} iterations).
#' @param monitors Optional character vector of monitor roots; if \code{NULL},
#'   default monitors are inferred.
#' @param nchains Integer; number of chains to run in \code{runMCMC}.
#' @param opts List of options as returned by \code{samOptiPro_options()}.
#'
#' @return A list with components \code{samples}, \code{samples2},
#'   \code{runtime_s} and \code{conf}.
#' @export
run_baseline_config <- function(build_fn, niter,
                                nburnin = floor(0.25 * niter),
                                thin = 1L,
                                monitors = NULL,
                                nchains = 1L,
                                opts = samOptiPro_options()) {
  built <- .fresh_build(build_fn, monitors = monitors, thin = thin, opts = opts)


  # Build MCMC en R puis compile avec 'project = built$cmodel'
  mcmc  <- nimble::buildMCMC(built$conf)
  cmcmc <- nimble::compileNimble(mcmc, project = built$cmodel, resetFunctions = TRUE)

  out   <- .run_and_collect(cmcmc, niter = niter, nburnin = nburnin, thin = thin, nchains = nchains)

  list(
    samples   = out$samples,
    samples2  = out$samples2,
    runtime_s = out$runtime_s,
    conf      = built$conf
  )
}


#' Run HMC/NUTS on all eligible nodes with slice fallback
#'
#' @param build_fn Function with no arguments returning a list with
#'   \code{model}, optional \code{cmodel}, and optional \code{monitors}.
#' @param niter Total number of MCMC iterations.
#' @param nburnin Number of burn-in iterations discarded.
#' @param thin Thinning interval (keep 1 draw every \code{thin} iterations).
#' @param monitors Optional character vector of monitor roots; if \code{NULL},
#'   defaults are inferred.
#' @param nchains Integer; number of chains for \code{runMCMC}.
#' @param opts List of options as returned by \code{samOptiPro_options()}.
#'
#' @return A list with \code{samples}, \code{samples2}, \code{runtime_s},
#'   \code{conf}, and a logical \code{hmc_applied} flag.
#' @export
run_hmc_all_nodes <- function(build_fn, niter,
                              nburnin = floor(0.25 * niter),
                              thin = 1L,
                              monitors = NULL,
                              nchains = 1L,
                              opts = samOptiPro_options()) {
  built <- .fresh_build(build_fn, monitors = monitors, thin = thin, opts = opts)


  hmc_ok <- TRUE
  tryCatch({
    if (!requireNamespace("nimbleHMC", quietly = TRUE)) stop("nimbleHMC not installed")
    nimbleHMC::configureHMC(built$conf, model = built$model)
  }, error = function(e) {
    message("configureHMC failed: ", e$message)
    hmc_ok <<- FALSE
  })

  # Rattrapage des noeuds non couverts
  uns <- try(built$conf$getUnsampledNodes(), silent = TRUE)
  if (!inherits(uns, "try-error") && length(uns)) {
    for (u in uns) built$conf$addSampler(u, type = "slice")
  }

  mcmc  <- nimble::buildMCMC(built$conf)
  cmcmc <- nimble::compileNimble(mcmc, project = built$cmodel, resetFunctions = TRUE)

  out <- .run_and_collect(cmcmc, niter = niter, nburnin = nburnin, thin = thin, nchains = nchains)
  list(
    samples     = out$samples,
    samples2    = out$samples2,
    runtime_s   = out$runtime_s,
    conf        = built$conf,
    hmc_applied = hmc_ok
  )
}


## ---- Internal functions (5) ----

#' @keywords internal
.fresh_build <- function(build_fn, monitors = NULL, thin = 1L, opts = samOptiPro_options()) {
  stopifnot(is.function(build_fn))

  # (1) Construire la "build" (supporte le pattern fabrique)
  built <- build_fn()
  if (is.function(built)) built <- built()

  if (!is.list(built) || is.null(built$model)) {
    stop("samOptiPro: build_fn() must return a list with at least $model.")
  }

  # (2) Verifier la signature NIMBLE et recuperer le modele R-level
  mdl <- built$model
  if (!.sop_is_model(mdl)) {
    stop("samOptiPro: build_fn()$model does not resemble a NIMBLE model (classes: ",
         paste(class(mdl), collapse = "/"), ")")
  }
  mdl <- .sop_get_uncompiled_model(mdl)

  # (3) Compiler le modele si necessaire
  cmod <- built$cmodel %||% try(nimble::compileNimble(mdl), silent = TRUE)
  if (inherits(cmod, "try-error") || is.null(cmod)) {
    stop("samOptiPro: model compilation failure.")
  }

  # (4) Monitors & configuration MCMC
  mons_in <- monitors %||% built$monitors
  # Optionnel : check des monitors
  # if (!is.null(mons_in)) {
  #   missing <- setdiff(mons_in, nimble::getNodeNames(mdl, stochOnly = FALSE, includeData = TRUE))
  #   if (length(missing)) warning("Monitors not found in model and will be ignored: ", paste(missing, collapse = ", "))
  # }
  conf <- .configure_with_monitors(mdl, monitors = mons_in, thin = thin, thin2 = thin, opts = opts)

  # (5) Sortie standardisee
  list(
    model           = mdl,
    cmodel          = cmod,
    conf            = conf,
    monitors        = attr(conf, "._sop_monitors_roots")  %||% character(0),
    monitors2       = attr(conf, "._sop_monitors2_roots") %||% character(0),
    monitors_nodes  = attr(conf, "._sop_monitors_nodes")  %||% character(0),
    monitors2_nodes = attr(conf, "._sop_monitors2_nodes") %||% character(0)
    # + eventuellement d'autres champs utiles de `built` si tu veux les propager
  )
}


#' @keywords internal
.merge_mcmc <- function(samples, samples2) {
  grab_mcpar <- function(obj) {
    if (inherits(obj, "mcmc.list") && length(obj)) return(attr(obj[[1]], "mcpar"))
    if (inherits(obj, "mcmc"))                     return(attr(obj, "mcpar"))
    NULL
  }
  mp <- grab_mcpar(samples); if (is.null(mp)) mp <- grab_mcpar(samples2)
  start <- if (!is.null(mp)) mp[1] else 1L
  thin  <- if (!is.null(mp)) mp[3] else 1L

  s1 <- .safe_as_mcmc_list(samples,  start = start, thin = thin)
  s2 <- .safe_as_mcmc_list(samples2, start = start, thin = thin)

  if (is.null(s1) && is.null(s2))
    stop("No samples (samples/samples2) were produced.")

  if (is.null(s2)) return(s1)
  if (is.null(s1)) return(s2)

  stopifnot(length(s1) == length(s2))
  coda::mcmc.list(lapply(seq_along(s1), function(i) {
    m1 <- as.matrix(s1[[i]]); m2 <- as.matrix(s2[[i]])
    if (NCOL(m1) == 0L) return(s2[[i]])
    if (NCOL(m2) == 0L) return(s1[[i]])
    mp <- attr(s1[[i]], "mcpar"); if (is.null(mp)) mp <- attr(s2[[i]], "mcpar")
    start <- if (!is.null(mp)) mp[1] else 1L
    thin  <- if (!is.null(mp)) mp[3] else 1L
    end   <- if (!is.null(mp)) mp[2] else (start + (NROW(m1) - 1L) * thin)
    coda::mcmc(cbind(m1, m2), start = start, end = end, thin = thin)
  }))
}


#' @keywords internal
.run_and_collect <- function(cmcmc, niter, nburnin, thin, nchains = 1L) {
  t0 <- proc.time()
  res <- nimble::runMCMC(
    cmcmc,
    niter    = as.integer(niter),
    nburnin  = as.integer(nburnin),
    thin     = as.integer(thin),
    nchains  = as.integer(nchains),
    samplesAsCodaMCMC = FALSE,  # <- pour garantir un retour {samples, samples2}
    summary  = FALSE,
    WAIC     = FALSE
  )
  t1 <- proc.time()
  runtime_s <- unname((t1 - t0)[["elapsed"]])
  if (!is.finite(runtime_s)) runtime_s <- NA_real_

  # 'res' est une list avec $samples et eventuellement $samples2
  list(
    samples   = res$samples,
    samples2  = res$samples2 %||% NULL,
    runtime_s = as.numeric(runtime_s)
  )
}


#' @keywords internal
.safe_as_mcmc_list <- function(x, start = 1L, thin = 1L) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "mcmc.list")) return(x)
  if (inherits(x, "mcmc"))      return(coda::mcmc.list(x))
  if (is.matrix(x) || is.data.frame(x)) {
    if (NCOL(x) == 0L) return(NULL)
    return(coda::mcmc.list(coda::mcmc(as.matrix(x), start = start, thin = thin)))
  }
  if (is.list(x)) {
    lst <- lapply(x, function(mat) {
      if (is.null(mat)) return(NULL)
      if (is.matrix(mat) || is.data.frame(mat)) {
        if (NCOL(mat) == 0L) return(NULL)
        return(coda::mcmc(as.matrix(mat), start = start, thin = thin))
      }
      if (inherits(mat, "mcmc")) return(mat)
      NULL
    })
    lst <- Filter(Negate(is.null), lst)
    if (!length(lst)) return(NULL)
    return(coda::mcmc.list(lst))
  }
  NULL
}


  grab_mcpar <- function(obj) {
    if (inherits(obj, "mcmc.list") && length(obj)) return(attr(obj[[1]], "mcpar"))
    if (inherits(obj, "mcmc"))                     return(attr(obj, "mcpar"))
    NULL
  }


