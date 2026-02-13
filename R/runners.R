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
#' Run a slice-only baseline MCMC configuration in parallel
#'
#' Runs multiple independent MCMC chains using a **strict slice-only sampling
#' strategy** (`onlySlice` semantics), in parallel, on either a PSOCK or FORK
#' backend. All stochastic nodes are updated exclusively with slice samplers,
#' with no use of random-walk, block, conjugate, or gradient-based samplers.
#'
#' This function is primarily intended as a **robust convergence-oriented
#' baseline** for complex or non-differentiable models (e.g. models involving
#' truncated distributions, constrained parameters, or piecewise likelihoods),
#' for which default sampler configurations or HMC/NUTS are inappropriate.
#'
#' Implementation notes:
#' - Each worker builds and configures its own model and MCMC independently.
#' - All existing samplers are removed and replaced by slice samplers on all
#'   stochastic nodes.
#' - Parallel execution can rely on a PSOCK cluster (portable) or a FORK backend
#'   (Unix-alike systems only).
#' - The reported runtime measures the wall-clock time spent in `runMCMC()`,
#'   excluding model construction and compilation.
#'
#' Builder contract:
#' - `build_fn()` must work either as `build_fn()` or
#'   `build_fn(chain_id = <int>, export_global = <logical>)`.
#' - It must return a list containing at least:
#'   - `model`  : an uncompiled `nimbleModel`
#'   - `cmodel` : the compiled model (`compileNimble(model)`)
#'   - `conf`   : an MCMC configuration object
#'
#' @param build_fn Function building per-chain model objects (see Builder contract).
#' @param niter Integer; total number of MCMC iterations.
#' @param nburnin Integer; number of burn-in iterations.
#' @param thin Integer; thinning interval.
#' @param monitors Optional character vector of monitored nodes.
#' @param nchains Integer; number of chains to run.
#' @param seed Integer; base seed (chain `k` uses `seed + k - 1`).
#' @param n_cores Integer; number of parallel workers (default = `nchains`).
#' @param extra_export Optional character vector of additional objects to export
#'   to workers.
#' @param parallel_backend Character; parallel backend to use (`"PSOCK"`,
#'   `"FORK"`, or `"AUTO"`).
#' @param worker_log Logical; whether to emit minimal per-worker diagnostic logs.
#' @param useConjugacy Logical; passed for API consistency but ignored (slice-only).
#'
#' @return A list with components:
#' - `samples`    : posterior samples (combined across chains when possible)
#' - `runtime_s`  : wall-clock runtime (seconds) for the parallel `runMCMC()` phase
#' - `conf`       : `NULL` in parallel mode; per-chain configurations are attached
#'                  as attribute `conf_by_chain`
#'
#' @export
run_only_slice_parallel <- function(build_fn,
                                    niter,
                                    nburnin  = floor(0.25 * niter),
                                    thin     = 1L,
                                    monitors = NULL,
                                    nchains  = 1L,
                                    seed     = 1L,
                                    n_cores  = nchains,
                                    extra_export = character(0),
                                    parallel_backend = c("PSOCK", "FORK", "AUTO"),
                                    worker_log = TRUE,
                                    useConjugacy = FALSE) {
  stopifnot(is.function(build_fn))
  niter    <- as.integer(niter)
  nburnin  <- as.integer(nburnin)
  thin     <- as.integer(thin)
  nchains  <- as.integer(nchains)
  n_cores  <- as.integer(n_cores)
  seed     <- as.integer(seed)

  stopifnot(niter >= 1L, nburnin >= 0L, thin >= 1L, nchains >= 1L, n_cores >= 1L)

  supports_chain_id <- ("chain_id" %in% names(formals(build_fn)))
  backend <- match.arg(parallel_backend)
  if (backend == "AUTO") backend <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"

  ## --- builder
  .build_one <- function(chain_id) {
    if (supports_chain_id) {
      built <- build_fn(chain_id = chain_id, export_global = FALSE)
    } else {
      built <- build_fn()
    }
    if (is.null(built$model)) stop("build_fn must return at least $model.")
    built
  }

  ## --- run one chain (historical-compatible)
  .run_one <- function(chain_id) {
    chain_id <- as.integer(chain_id)
    set.seed(seed + chain_id)

    built <- .build_one(chain_id)
    m <- built$model

    ## Reconfigure baseline MCMC exactly like your script
    conf <- nimble::configureMCMC(
      m,
      monitors     = monitors,
      useConjugacy = isTRUE(useConjugacy),
      onlySlice    = TRUE
    )

    mcmc  <- nimble::buildMCMC(conf)

    ## Compile model first
    Cm <- nimble::compileNimble(m, showCompilerOutput = FALSE)

    ## Compile MCMC using project = *R model* (as in your script)
    Cmc <- nimble::compileNimble(mcmc, project = m)

    t0 <- proc.time()
    samples <- nimble::runMCMC(
      Cmc,
      niter             = niter,
      nburnin           = nburnin,
      thin              = thin,
      setSeed           = seed + chain_id,
      samplesAsCodaMCMC = TRUE,
      summary           = FALSE,
      WAIC              = FALSE
    )
    rt <- as.numeric((proc.time() - t0)[["elapsed"]])

    if (!inherits(samples, "mcmc")) {
      tmp <- try(coda::as.mcmc(samples), silent = TRUE)
      if (!inherits(tmp, "try-error")) samples <- tmp
    }

    list(ok = TRUE, chain_id = chain_id, samples = samples, runtime_s = rt, conf = conf, model = m)
  }

  ## --- sequential
  if (nchains == 1L || !supports_chain_id) {
    out <- .run_one(1L)
    return(list(samples = out$samples, runtime_s = out$runtime_s, conf = out$conf, slice_only = TRUE))
  }

  ## --- FORK
  if (backend == "FORK") {
    if (.Platform$OS.type != "unix") stop("FORK backend is only available on Unix-alikes.")
    if (!requireNamespace("parallel", quietly = TRUE)) stop("parallel required")

    ans <- parallel::mclapply(
      X = seq_len(nchains),
      FUN = function(ch) {
        tryCatch(.run_one(ch),
                 error = function(e) list(ok = FALSE, chain_id = as.integer(ch), message = conditionMessage(e)))
      },
      mc.cores       = n_cores,
      mc.preschedule = FALSE
    )

    if (any(vapply(ans, is.null, logical(1)))) stop("FORK hard-crash detected (NULL worker result).")

    errs <- Filter(function(x) is.list(x) && isFALSE(x$ok), ans)
    if (length(errs)) {
      msg <- paste0(vapply(errs, function(x) sprintf("Chain %d: %s", x$chain_id, x$message), character(1)),
                    collapse = "\n")
      stop("FORK worker error(s):\n", msg)
    }

    samples_by_chain <- lapply(ans, `[[`, "samples")
    runtimes         <- vapply(ans, `[[`, numeric(1), "runtime_s")
    conf_by_chain    <- lapply(ans, `[[`, "conf")

    samples_out <- try(coda::mcmc.list(samples_by_chain), silent = TRUE)
    if (inherits(samples_out, "try-error")) samples_out <- samples_by_chain

    res <- list(samples = samples_out, runtime_s = max(runtimes), conf = NULL, slice_only = TRUE)
    attr(res, "conf_by_chain") <- conf_by_chain
    res$conf_by_chain <- conf_by_chain
    return(res)
  }

  ## --- PSOCK
  if (!requireNamespace("parallel", quietly = TRUE)) stop("parallel required")

  log_file <- NULL
  if (isTRUE(worker_log)) log_file <- tempfile(pattern = "psock_workers_", fileext = ".log")

  cl <- parallel::makeCluster(
    n_cores,
    type    = "PSOCK",
    outfile = if (isTRUE(worker_log)) log_file else NULL
  )
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(nimble))
    suppressPackageStartupMessages(library(coda))
    options(warn = 1)
    NULL
  })

  export_env <- environment()

  ## optional extra exports (avoid big objects)
  if (length(extra_export)) {
    banned <- c("mydata", "myconstants", "Data_nimble", "Const_nimble", "data_", "const_")
    if (any(extra_export %in% banned)) {
      stop("extra_export contains large objects that should not be PSOCK-exported: ",
           paste(intersect(extra_export, banned), collapse = ", "),
           ". Capture them inside build_fn (closure) instead.")
    }
    for (nm in extra_export) {
      if (!exists(nm, envir = export_env, inherits = FALSE) &&
          exists(nm, envir = .GlobalEnv, inherits = TRUE)) {
        assign(nm, get(nm, envir = .GlobalEnv, inherits = TRUE), envir = export_env)
      }
    }
  }

  parallel::clusterExport(
    cl,
    varlist = unique(c(
      "build_fn", "niter", "nburnin", "thin", "monitors", "seed",
      "supports_chain_id", "useConjugacy",
      ".build_one", ".run_one",
      extra_export
    )),
    envir = export_env
  )

  parallel::clusterSetRNGStream(cl, iseed = seed)

  ok <- TRUE
  err_msg <- NULL
  ans <- tryCatch(
    parallel::parLapply(cl, X = seq_len(nchains), fun = function(ch) {
      tryCatch(.run_one(ch),
               error = function(e) list(ok = FALSE, chain_id = as.integer(ch), message = conditionMessage(e)))
    }),
    error = function(e) {
      ok <<- FALSE
      err_msg <<- conditionMessage(e)
      NULL
    }
  )

  if (!ok || is.null(ans)) {
    if (isTRUE(worker_log) && !is.null(log_file) && file.exists(log_file)) {
      ll <- readLines(log_file, warn = FALSE)
      tail_txt <- paste(tail(ll, 160), collapse = "\n")
      stop("PSOCK crashed (", err_msg, "). Worker log tail:\n\n", tail_txt)
    }
    stop("PSOCK crashed (", err_msg, ").")
  }

  errs <- Filter(function(x) is.list(x) && isFALSE(x$ok), ans)
  if (length(errs)) {
    msg <- paste0(vapply(errs, function(x) sprintf("Chain %d: %s", x$chain_id, x$message), character(1)),
                  collapse = "\n")
    if (isTRUE(worker_log) && !is.null(log_file) && file.exists(log_file)) {
      ll <- readLines(log_file, warn = FALSE)
      tail_txt <- paste(tail(ll, 160), collapse = "\n")
      stop("PSOCK worker error(s):\n", msg, "\n\nWorker log tail:\n\n", tail_txt)
    }
    stop("PSOCK worker error(s):\n", msg)
  }

  samples_by_chain <- lapply(ans, `[[`, "samples")
  runtimes         <- vapply(ans, `[[`, numeric(1), "runtime_s")
  conf_by_chain    <- lapply(ans, `[[`, "conf")

  samples_out <- try(coda::mcmc.list(samples_by_chain), silent = TRUE)
  if (inherits(samples_out, "try-error")) samples_out <- samples_by_chain

  res <- list(samples = samples_out, runtime_s = max(runtimes), conf = NULL, slice_only = TRUE)
  attr(res, "conf_by_chain") <- conf_by_chain
  res$conf_by_chain <- conf_by_chain

  res
}

#' Run a baseline NIMBLE MCMC configuration in parallel (timing = runMCMC only)
#'
#' Runs multiple MCMC chains on a PSOCK cluster using a *baseline* NIMBLE
#' configuration. The reported runtime follows the historic semantics used in
#' samOptiPro: **`runtime_s` measures only the wall-clock time spent inside
#' `nimble::runMCMC()` across all chains**, excluding model building, MCMC setup,
#' and compilation.
#'
#' Implementation notes:
#' - Each worker calls `build_fn()` to build (and possibly compile) the model and
#'   MCMC configuration **outside** the timed section.
#' - Each worker then builds + compiles the MCMC and runs `nimble::runMCMC()`; only
#'   the `runMCMC()` call is timed per chain.
#' - The master records a global wall-clock time around the parallel `runMCMC()`
#'   phase (start just before `parLapply()`, stop after all chains return).
#'
#' Builder contract:
#' - `build_fn()` must work as either `build_fn()` or
#'   `build_fn(chain_id = <int>, export_global = <logical>)`.
#' - It must return a list containing at least:
#'   - `model`  : an uncompiled `nimbleModel` (recommended; used for diagnostics)
#'   - `cmodel` : the compiled model (`compileNimble(model)`)
#'   - `conf`   : an MCMC configuration (e.g. from `nimble::configureMCMC()`)
#'   Optionally it may include `monitors`.
#'
#' @param build_fn Function that builds per-chain model objects (see Builder contract).
#' @param niter Integer, total MCMC iterations.
#' @param nburnin Integer, burn-in iterations.
#' @param thin Integer, thinning interval.
#' @param monitors Optional character vector of additional monitors to add.
#' @param nchains Integer, number of chains to run.
#' @param seed Integer, base seed; chain `k` uses `seed + k - 1`.
#' @param n_cores Integer, number of PSOCK workers (defaults to `nchains`).
#' @param extra_export Character vector of additional objects to export to workers.
#' @param samplesAsCodaMCMC Logical; if TRUE, `runMCMC()` returns coda objects.
#' @param export_global Logical; if TRUE export full `.GlobalEnv` (use with care).
#' @param opts Optional samOptiPro options object (reserved for future use).
#'
#' @return A list with:
#' - `samples`          : list of samples (one element per chain; type depends on `samplesAsCodaMCMC`)
#' - `runtime_s`        : global wall time (seconds) for the parallel `runMCMC()` phase only
#' - `runtime_by_chain` : per-chain wall time (seconds) for `runMCMC()` only
#'
#' @export
run_baseline_config_parallel <- function(build_fn,
                                         niter,
                                         nburnin,
                                         thin,
                                         monitors = NULL,
                                         nchains  = 3L,
                                         seed     = 1L,
                                         n_cores  = nchains,
                                         extra_export = character(0),
                                         samplesAsCodaMCMC = FALSE,
                                         export_global = FALSE,
                                         opts = NULL) {

  stopifnot(is.function(build_fn))
  nchains <- as.integer(nchains)
  n_cores <- as.integer(n_cores)
  stopifnot(nchains >= 1L, n_cores >= 1L)

  ## ---- Probe the builder (generic contract) ----
  parts0 <- try(build_fn(), silent = TRUE)
  if (inherits(parts0, "try-error") || !is.list(parts0)) {
    parts0 <- try(build_fn(chain_id = 1L, export_global = isTRUE(export_global)), silent = TRUE)
  }
  if (inherits(parts0, "try-error") || !is.list(parts0)) {
    stop("build_fn must return a list. It should work as build_fn() or build_fn(chain_id=, export_global=).")
  }

  need <- c("model", "cmodel", "conf")
  miss <- setdiff(need, names(parts0))
  if (length(miss)) stop("build_fn() is missing required fields: ", paste(miss, collapse = ", "))

  ## If monitors not provided, try builder-provided monitors
  if (is.null(monitors) || !length(monitors)) {
    if (!is.null(parts0$monitors) && length(parts0$monitors)) {
      monitors <- parts0$monitors
    }
  }

  ## ---- Cluster ----
  if (!requireNamespace("parallel", quietly = TRUE)) stop("parallel required")
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  ## ---- Load packages + define custom nimbleFunctions on each worker ----
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(nimble))
    suppressPackageStartupMessages(library(coda))

    ## nlchoose (log binomial coefficient)
    nlchoose <- nimble::nimbleFunction(
      run = function(n = nimble::double(0), k = nimble::double(0)) {
        nimble::returnType(nimble::double(0))
        k <- round(k); n <- round(n)
        if (n == k | k == 0) return(0)
        if (n < k) return(-Inf)
        return(sum(log((n - k + 1):n)) - sum(log(1:k)))
      }
    )

    ## hypergeometric log-density
    dhyperg <- nimble::nimbleFunction(
      run = function(x  = nimble::double(0),
                     nx = nimble::double(0),
                     N1 = nimble::double(0),
                     N0 = nimble::double(0),
                     log = nimble::integer(0, default = 1)) {
        nimble::returnType(nimble::double(0))
        logProb <- nlchoose(N1, x) + nlchoose(N0, nx - x) - nlchoose(N0 + N1, nx)
        if (log) return(logProb)
        return(exp(logProb))
      }
    )

    ## hypergeometric RNG (n=1)
    rhyperg <- nimble::nimbleFunction(
      run = function(n  = nimble::integer(0, default = 1),
                     nx = nimble::double(0),
                     N1 = nimble::double(0),
                     N0 = nimble::double(0)) {
        nimble::returnType(nimble::double(0))
        if (n != 1) print("rhyperg only allows n=1; using n = 1")

        ## Use NIMBLE RNG (supported in nimbleFunctions)
        dev <- nimble::runif(1, 0, 1)

        minx <- max(nx - N0, 0)
        maxx <- min(nx, N1)

        xs <- numeric(length = maxx - minx + 1)
        ps <- xs

        xs[1] <- minx
        ps[1] <- dhyperg(minx, nx, N1, N0, log = 0)

        if (maxx >= (minx + 1)) {
          for (i in (minx + 1):maxx) {
            xs[i - minx + 1] <- i
            ps[i - minx + 1] <- dhyperg(i, nx, N1, N0, log = 0) + ps[i - minx]
          }
        }

        return(xs[which(!dev > ps)[1]])
      }
    )

    assign("nlchoose", nlchoose, .GlobalEnv)
    assign("dhyperg",  dhyperg,  .GlobalEnv)
    assign("rhyperg",  rhyperg,  .GlobalEnv)

    NULL
  })

  ## ---- Export strategy ----
  parallel::clusterExport(cl, varlist = c("build_fn"), envir = environment())

  if (isTRUE(export_global)) {
    all_names <- ls(envir = .GlobalEnv, all.names = TRUE)
    to_export <- unique(c(all_names, extra_export))
    parallel::clusterExport(cl, varlist = to_export, envir = .GlobalEnv)
  } else {
    to_export <- unique(extra_export)
    to_export <- to_export[to_export %in% ls(envir = .GlobalEnv, all.names = TRUE)]
    if (length(to_export)) parallel::clusterExport(cl, varlist = to_export, envir = .GlobalEnv)
  }

  ## ---- Worker (compile OUTSIDE timing; time only runMCMC) ----
  worker_fun <- function(chain_id,
                         build_fn,
                         niter, nburnin, thin,
                         monitors,
                         seed,
                         samplesAsCodaMCMC,
                         export_global) {

    chain_id <- as.integer(chain_id)
    set.seed(as.integer(seed) + chain_id - 1L)

    parts <- try(build_fn(chain_id = chain_id, export_global = isTRUE(export_global)), silent = TRUE)
    if (inherits(parts, "try-error") || !is.list(parts)) {
      parts <- build_fn()
      parts$chain_id <- chain_id
    }

    if (!is.list(parts) || is.null(parts$conf) || is.null(parts$cmodel)) {
      stop("build_fn must return a list with at least $conf and $cmodel (and ideally $model).")
    }

    conf <- parts$conf

    if (!is.null(monitors) && length(monitors)) {
      try(conf$addMonitors(monitors), silent = TRUE)
    } else if (!is.null(parts$monitors) && length(parts$monitors)) {
      try(conf$addMonitors(parts$monitors), silent = TRUE)
    }

    mcmc  <- nimble::buildMCMC(conf)
    cmcmc <- nimble::compileNimble(mcmc, project = parts$cmodel, resetFunctions = TRUE)

    t0  <- Sys.time()
    smp <- nimble::runMCMC(
      cmcmc,
      niter   = as.integer(niter),
      nburnin = as.integer(nburnin),
      thin    = as.integer(thin),
      nchains = 1L,
      samplesAsCodaMCMC = isTRUE(samplesAsCodaMCMC),
      summary = FALSE,
      WAIC    = FALSE
    )
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    list(samples = smp, runtime_s = dt)
  }

  ## ---- Run (historic timing: wall time of the parallel runMCMC phase) ----
  time_start <- Sys.time()

  ans <- parallel::parLapply(
    cl,
    X = seq_len(nchains),
    fun = function(ch, build_fn, niter, nburnin, thin, monitors, seed, samplesAsCodaMCMC, export_global) {
      worker_fun(ch, build_fn, niter, nburnin, thin, monitors, seed, samplesAsCodaMCMC, export_global)
    },
    build_fn = build_fn,
    niter    = niter,
    nburnin  = nburnin,
    thin     = thin,
    monitors = monitors,
    seed     = seed,
    samplesAsCodaMCMC = samplesAsCodaMCMC,
    export_global = export_global
  )

  time_end  <- Sys.time()
  runtime_s <- as.numeric(difftime(time_end, time_start, units = "secs"))

  samples  <- lapply(ans, `[[`, "samples")
  runtimes <- vapply(ans, `[[`, numeric(1), "runtime_s")

  list(
    samples          = samples,
    runtime_s        = runtime_s,
    runtime_by_chain = runtimes
  )
}
#' Run a baseline NIMBLE MCMC configuration in parallel (robust “large-model” variant)
#'
#' @description
#' `run_baseline_config_parallel_bis()` is the *large-model* counterpart of
#' `run_baseline_config_parallel()`. It is designed to remain stable on very large
#' NIMBLE models (e.g., >80,000 nodes; WGNAS-scale) and on HPC/cluster environments,
#' where naïve PSOCK parallelization can trigger memory blow-ups or fail due to
#' missing closure bindings on workers.
#'
#' Compared to the non-`_bis` function, this version:
#' \itemize{
#'   \item supports resolving `build_fn` by *name* (character) on workers, including optional namespace lookup;
#'   \item supports injecting `build_fn` as text (`opts$build_fn_text`) and exporting the *captured environment*
#'         of the builder closure (`opts$build_fn_object`, `opts$builder_env_names`) to prevent
#'         errors like “objet 'hindcast_' introuvable”;
#'   \item avoids mandatory disk I/O (no log files required), suitable for clusters with restricted write access;
#'   \item times only the MCMC execution phase (optionally using `$run(time=TRUE)`), keeping compilation outside
#'         the per-chain runtime by construction;
#'   \item says memory: optional worker-side `gc()` clean-up after each chain to reduce peak usage on large graphs;
#'   \item can optionally return the *actual* `nimbleMCMCconf` used on each worker (`opts$return_conf = TRUE`)
#'         so downstream diagnostics can reuse the exact sampler/monitor configuration.
#' }
#'
#' The function assumes that `build_fn` returns at least a list with `model` and `conf`
#' (plus optional `monitors`). The returned object includes `samples`, `runtime_s`,
#' `runtime_by_chain`, and `sampler_times`, and optionally `conf_by_chain`.
#'
#' @param build_fn A builder function, or a single character string giving the builder name.
#'        The builder must return a list with at least `model` and `conf`.
#' @param niter,nburnin,thin MCMC settings passed to NIMBLE.
#' @param monitors Optional monitor vector (if `NULL`, may fall back to `parts$monitors`).
#' @param nchains Number of parallel chains to run.
#' @param seed Base RNG seed (each chain uses `seed + chain_id - 1`).
#' @param n_cores Number of PSOCK workers.
#' @param extra_export Optional additional symbols to export to workers.
#' @param samplesAsCodaMCMC Passed to `nimble::runMCMC()` fallback (when not using `$run()`).
#' @param export_global If `TRUE`, exports the whole master `.GlobalEnv` (not recommended unless needed).
#' @param opts List of advanced options (cluster type, packages, builder injection, timings).
#'
#' @return A list with:
#' \itemize{
#'   \item `samples`: list of per-chain samples (matrix or coda, depending on options),
#'   \item `runtime_s`: wall-clock time of the full parallel phase,
#'   \item `runtime_by_chain`: numeric vector of per-chain MCMC runtimes,
#'   \item `sampler_times`: list of `Cmcmc$getTimes()` outputs (if enabled),
#'   \item `conf_by_chain`: optional list of `nimbleMCMCconf` (if `opts$return_conf = TRUE`),
#'   \item `opts_effective`: the effective options used.
#' }
#'
#' @export
run_baseline_config_parallel_bis <- function(build_fn,
                                             niter,
                                             nburnin,
                                             thin,
                                             monitors = NULL,
                                             nchains  = 3L,
                                             seed     = 1L,
                                             n_cores  = nchains,
                                             extra_export = character(0),
                                             samplesAsCodaMCMC = FALSE,
                                             export_global = FALSE,
                                             opts = NULL) {

  stopifnot(is.function(build_fn) || (is.character(build_fn) && length(build_fn) == 1L))
  stopifnot(is.numeric(niter), length(niter) == 1L, niter >= 1)
  stopifnot(is.numeric(nburnin), length(nburnin) == 1L, nburnin >= 0)
  stopifnot(is.numeric(thin), length(thin) == 1L, thin >= 1)

  nchains <- as.integer(nchains)
  n_cores <- as.integer(n_cores)
  stopifnot(nchains >= 1L, n_cores >= 1L)

  if (!requireNamespace("parallel", quietly = TRUE)) stop("parallel required")

  `%||%` <- function(x, y) if (is.null(x)) y else x
  opts <- if (is.null(opts)) list() else opts

  ## ---- Builder resolution strategy ----
  build_fn_name <- NULL
  build_fn_ns   <- NULL
  build_fn_text <- opts$build_fn_text %||% NULL

  ## ---- Solution B (added, non-invasive) ----
  ## If TRUE: return the *actual* conf used on each worker (post-monitors)
  return_conf <- isTRUE(opts$return_conf %||% FALSE)

  # NEW: optionally export the *environment bindings* of the builder closure
  # (needed for closures like build_M <- local({hindcast_ <- hindcast; ...; function(...) ...})
  # when you inject code via build_fn_text).
  #
  # Provide opts$build_fn_object = build_M (the function object from master).
  # Optionally restrict what to export with opts$builder_env_names (character vector).
  build_fn_object <- opts$build_fn_object %||% NULL
  builder_env_names <- opts$builder_env_names %||% NULL
  builder_env_list <- NULL

  if (is.character(build_fn)) {
    build_fn_name <- build_fn
    build_fn <- NULL
  } else if (is.function(build_fn)) {
    if (!is.null(opts$build_fn_name)) build_fn_name <- as.character(opts$build_fn_name)
  } else {
    stop("build_fn must be a function or a single character string (function name).")
  }

  if (!is.null(opts$build_fn_ns))   build_fn_ns   <- as.character(opts$build_fn_ns)
  if (!is.null(opts$build_fn_name)) build_fn_name <- as.character(opts$build_fn_name)

  if (!is.null(build_fn_text)) {
    stopifnot(is.character(build_fn_text), length(build_fn_text) >= 1L)
    if (is.null(build_fn_name) || !nzchar(build_fn_name)) {
      stop("opts$build_fn_text provided but build_fn_name is NULL. Provide build_fn as a name or set opts$build_fn_name.")
    }
    if (!is.null(build_fn_object) && is.function(build_fn_object)) {
      e <- environment(build_fn_object)
      nm <- ls(envir = e, all.names = TRUE)
      if (!is.null(builder_env_names)) {
        stopifnot(is.character(builder_env_names))
        nm <- intersect(nm, builder_env_names)
      }
      builder_env_list <- if (length(nm)) mget(nm, envir = e, inherits = FALSE) else list()
    }
  }

  ## ---- Cluster options ----
  cluster_type <- opts$cluster_type %||% "PSOCK"
  outfile <- opts$outfile %||% ""

  prefer_run_object <- isTRUE(opts$prefer_run_object %||% TRUE)
  reset_on_run      <- isTRUE(opts$reset_on_run %||% TRUE)
  mcmc_time         <- isTRUE(opts$mcmc_time %||% TRUE)
  progressBar       <- isTRUE(opts$progressBar %||% FALSE)

  resetFunctions     <- isTRUE(opts$resetFunctions %||% TRUE)
  showCompilerOutput <- isTRUE(opts$showCompilerOutput %||% FALSE)
  gc_worker          <- isTRUE(opts$gc_worker %||% TRUE)

  worker_libpaths <- opts$worker_libpaths %||% NULL
  worker_packages <- opts$worker_packages %||% character(0)

  cl <- parallel::makeCluster(n_cores, type = cluster_type, outfile = outfile)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  ## ---- Load required packages on workers ----
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(nimble))
    if (requireNamespace("coda", quietly = TRUE)) NULL
    NULL
  })

  ## ---- Set .libPaths + load user packages on workers ----
  parallel::clusterExport(
    cl,
    varlist = c("worker_libpaths", "worker_packages"),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    if (!is.null(worker_libpaths) && length(worker_libpaths)) {
      .libPaths(unique(c(worker_libpaths, .libPaths())))
    }
    if (length(worker_packages)) {
      for (p in worker_packages) {
        suppressPackageStartupMessages(library(p, character.only = TRUE))
      }
    }
    NULL
  })

  ## ---- Inject builder closure environment bindings (if provided) ----
  parallel::clusterExport(
    cl,
    varlist = c("builder_env_list"),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    if (!is.null(builder_env_list) && length(builder_env_list)) {
      for (nm in names(builder_env_list)) {
        assign(nm, builder_env_list[[nm]], envir = .GlobalEnv)
      }
    }
    NULL
  })

  ## ---- Inject builder code on workers if requested ----
  parallel::clusterExport(
    cl,
    varlist = c("build_fn_text", "build_fn_name"),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    if (!is.null(build_fn_text)) {
      val <- try(eval(parse(text = build_fn_text), envir = .GlobalEnv), silent = TRUE)
      if (!inherits(val, "try-error")) {
        if (!exists(build_fn_name, envir = .GlobalEnv, inherits = FALSE) && is.function(val)) {
          assign(build_fn_name, val, envir = .GlobalEnv)
        }
      }
      if (!exists(build_fn_name, envir = .GlobalEnv, inherits = FALSE)) {
        stop("build_fn_text evaluated but did not create a function named '", build_fn_name, "' on worker.")
      }
    }
    NULL
  })

  ## ---- Export light scalars/options + (optional) the function itself ----
  parallel::clusterExport(
    cl,
    varlist = c(
      "niter", "nburnin", "thin", "monitors", "seed",
      "samplesAsCodaMCMC",
      "build_fn_name", "build_fn_ns",
      "prefer_run_object", "reset_on_run", "mcmc_time", "progressBar",
      "resetFunctions", "showCompilerOutput", "gc_worker",
      "return_conf"   # <- added (Solution B)
    ),
    envir = environment()
  )

  if (is.null(build_fn_name)) {
    parallel::clusterExport(cl, varlist = "build_fn", envir = environment())
  }

  ## ---- Optional exports ----
  if (isTRUE(export_global)) {
    all_names <- ls(envir = .GlobalEnv, all.names = TRUE)
    to_export <- unique(c(all_names, extra_export))
    parallel::clusterExport(cl, varlist = to_export, envir = .GlobalEnv)
  } else {
    to_export <- unique(extra_export)
    to_export <- to_export[to_export %in% ls(envir = .GlobalEnv, all.names = TRUE)]
    if (length(to_export)) parallel::clusterExport(cl, varlist = to_export, envir = .GlobalEnv)
  }

  ## ---- Worker function ----
  worker_fun <- function(chain_id) {
    chain_id <- as.integer(chain_id)
    set.seed(as.integer(seed) + chain_id - 1L)

    bf <- NULL

    # 1) Namespace resolution (exported or internal)
    if (!is.null(build_fn_name) && !is.null(build_fn_ns)) {
      ns <- try(asNamespace(build_fn_ns), silent = TRUE)
      if (!inherits(ns, "try-error") && exists(build_fn_name, envir = ns, inherits = TRUE)) {
        bf <- get(build_fn_name, envir = ns, inherits = TRUE)
      }
    }

    # 2) GlobalEnv resolution (injected)
    if (is.null(bf) && !is.null(build_fn_name) &&
        exists(build_fn_name, envir = .GlobalEnv, inherits = FALSE)) {
      bf <- get(build_fn_name, envir = .GlobalEnv, inherits = FALSE)
    }

    # 3) Exported function object fallback
    if (is.null(bf) && is.null(build_fn_name)) bf <- build_fn

    if (is.null(bf) || !is.function(bf)) {
      stop("Cannot resolve build_fn on worker. Tried: ",
           if (!is.null(build_fn_ns)) paste0(build_fn_ns, " namespace, then ") else "",
           ".GlobalEnv, then exported function object.")
    }

    parts <- bf(chain_id = chain_id, export_global = FALSE)
    if (!is.list(parts)) stop("build_fn must return a list")

    need <- c("model", "conf")
    miss <- setdiff(need, names(parts))
    if (length(miss)) stop("build_fn() is missing required fields: ", paste(miss, collapse = ", "))

    mon_use <- monitors
    if (is.null(mon_use) || !length(mon_use)) mon_use <- parts$monitors %||% NULL
    if (!is.null(mon_use) && length(mon_use)) {
      try(parts$conf$addMonitors(mon_use), silent = TRUE)
    }

    ## ---- Solution B (added): capture the actual conf used (post-monitors)
    conf_used <- NULL
    if (isTRUE(return_conf)) conf_used <- parts$conf

    Rmcmc <- nimble::buildMCMC(parts$conf)

    Cmcmc <- nimble::compileNimble(
      Rmcmc,
      project = parts$model,
      resetFunctions = isTRUE(resetFunctions),
      showCompilerOutput = isTRUE(showCompilerOutput)
    )

    t0 <- Sys.time()

    if (prefer_run_object && is.function(Cmcmc$run) && !isTRUE(samplesAsCodaMCMC)) {
      Cmcmc$run(
        niter       = as.integer(niter),
        nburnin     = as.integer(nburnin),
        thin        = as.integer(thin),
        time        = isTRUE(mcmc_time),
        progressBar = isTRUE(progressBar),
        reset       = isTRUE(reset_on_run)
      )
      smp <- as.matrix(Cmcmc$mvSamples)
    } else {
      smp <- nimble::runMCMC(
        Cmcmc,
        niter   = as.integer(niter),
        nburnin = as.integer(nburnin),
        thin    = as.integer(thin),
        nchains = 1L,
        samplesAsCodaMCMC = isTRUE(samplesAsCodaMCMC),
        summary = FALSE,
        WAIC    = FALSE
      )
    }

    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    st <- NULL
    if (isTRUE(mcmc_time) && is.function(Cmcmc$getTimes)) {
      st <- try(Cmcmc$getTimes(), silent = TRUE)
      if (inherits(st, "try-error")) st <- NULL
    }

    if (gc_worker) {
      ## do not rm conf_used if requested (it must survive in return object)
      rm(parts, Rmcmc, Cmcmc)
      gc(); gc()
    }

    out <- list(samples = smp, runtime_s = dt, sampler_times = st)
    if (isTRUE(return_conf)) out$conf <- conf_used
    out
  }

  ## ---- Run ----
  time_start <- Sys.time()
  ans <- parallel::parLapply(cl, X = seq_len(nchains), fun = worker_fun)
  time_end  <- Sys.time()

  #runtime_s <- as.numeric(difftime(time_end, time_start, units = "secs")) #Overall runtime_s
  runtime_s <- max(vapply(ans, `[[`, numeric(1), "runtime_s"), na.rm = TRUE)
  out <- list(
    samples          = lapply(ans, `[[`, "samples"),
    runtime_s        = runtime_s,
    runtime_by_chain = vapply(ans, `[[`, numeric(1), "runtime_s"),
    sampler_times    = lapply(ans, `[[`, "sampler_times"),
    opts_effective   = list(
      cluster_type = cluster_type,
      build_fn_name = build_fn_name,
      build_fn_ns   = build_fn_ns,
      worker_libpaths = worker_libpaths,
      worker_packages = worker_packages,
      build_fn_text = !is.null(build_fn_text),
      builder_env_exported = !is.null(builder_env_list) && length(builder_env_list) > 0,
      prefer_run_object = prefer_run_object,
      reset_on_run = reset_on_run,
      mcmc_time = mcmc_time,
      progressBar = progressBar,
      resetFunctions = resetFunctions,
      showCompilerOutput = showCompilerOutput,
      gc_worker = gc_worker,
      return_conf = return_conf  # <- added (Solution B)
    )
  )

  ## ---- Solution B (added): expose conf_by_chain when requested
  if (isTRUE(return_conf)) {
    out$conf_by_chain <- lapply(ans, `[[`, "conf")
  }

  out
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

  # Catching up with uncovered nodes
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
#' Run HMC/NUTS on all eligible nodes (with slice fallback), optionally in parallel
#'
#' @description
#' Configures gradient-based MCMC using \pkg{nimbleHMC} (HMC/NUTS) for all nodes
#' that can be handled by HMC. Any remaining unsampled nodes are automatically
#' assigned a slice sampler as a robust fallback. If \code{nchains > 1} and
#' \code{build_fn} supports a \code{chain_id} argument, chains are run in true
#' parallel (PSOCK). Otherwise, execution falls back to a safe sequential mode.
#'
#' @param build_fn Builder function that creates a fresh NIMBLE model and returns a list
#'   with at least \code{model}, \code{cmodel} and \code{conf}. For parallel execution,
#'   the builder should accept \code{chain_id} and build/compile chain-specific objects.
#' @param niter Total number of MCMC iterations.
#' @param nburnin Number of burn-in iterations discarded.
#' @param thin Thinning interval (keep 1 draw every \code{thin} iterations).
#' @param monitors Optional character vector of monitor roots. If \code{NULL}, monitors
#'   may be inferred by \code{.fresh_build()} or taken from the builder output.
#' @param nchains Integer; number of chains. If \code{nchains > 1} and the builder
#'   supports \code{chain_id}, chains are run in parallel (PSOCK).
#' @param opts List of options as returned by \code{samOptiPro_options()} (e.g. seed,
#'   export strategy handled by \code{.fresh_build()}).
#'
#' @details
#' Runtime is measured for the sampling phase only (around \code{runMCMC()}), not for
#' model/MCMC compilation. In parallel mode, \code{runtime_s} is the maximum runtime
#' across chains (slowest chain). Because each chain has its own MCMC configuration
#' (after HMC + slice fallback), \code{conf} is returned in sequential mode only.
#'
#' @return A list with components:
#' \describe{
#'   \item{samples}{Posterior draws (typically a \code{coda::mcmc.list} in parallel mode).}
#'   \item{samples2}{Optional secondary draws (as produced by \code{.run_and_collect()}).}
#'   \item{runtime_s}{Runtime (seconds) of the sampling phase; slowest chain in parallel mode.}
#'   \item{conf}{MCMC configuration (\code{nimbleMCMCconf}) in sequential mode; \code{NULL} in parallel.}
#'   \item{hmc_applied}{Logical; \code{TRUE} if HMC/NUTS configuration succeeded (all chains in parallel).}
#'   \item{conf_by_chain}{List of per-chain MCMC configurations (parallel mode); also attached as an attribute.}
#' }
#'
#' @importFrom stats setNames
#' @importFrom utils capture.output
#'
#' @examples
#' \dontrun{
#' res <- run_hmc_all_nodes(
#'   build_fn  = build_M,
#'   niter     = 2000,
#'   nburnin   = 500,
#'   thin      = 2,
#'   monitors  = c("N", "logit_theta"),
#'   nchains   = 4
#' )
#' }
#'
#' @export
run_hmc_all_nodes_parallel <- function(build_fn, niter,
                              nburnin = floor(0.25 * niter),
                              thin = 1L,
                              monitors = NULL,
                              nchains = 1L,
                              opts = samOptiPro_options()) {

  stopifnot(is.function(build_fn))
  nchains <- as.integer(nchains)
  stopifnot(nchains >= 1L)

  supports_chain_id <- ("chain_id" %in% names(formals(build_fn)))

  ## ---- Sequential fallback (safe) ----
  if (!supports_chain_id || nchains == 1L) {

    built <- .fresh_build(build_fn, monitors = monitors, thin = thin, opts = opts)

    hmc_ok <- TRUE
    tryCatch({
      if (!requireNamespace("nimbleHMC", quietly = TRUE)) stop("nimbleHMC not installed")
      nimbleHMC::configureHMC(built$conf, model = built$model)
    }, error = function(e) {
      message("configureHMC failed: ", e$message)
      hmc_ok <<- FALSE
    })

    uns <- try(built$conf$getUnsampledNodes(), silent = TRUE)
    if (!inherits(uns, "try-error") && length(uns)) {
      for (u in uns) built$conf$addSampler(u, type = "slice")
    }

    mcmc  <- nimble::buildMCMC(built$conf)
    cmcmc <- nimble::compileNimble(mcmc, project = built$cmodel, resetFunctions = TRUE)

    out <- .run_and_collect(cmcmc, niter = niter, nburnin = nburnin, thin = thin, nchains = nchains)

    return(list(
      samples     = out$samples,
      samples2    = out$samples2,
      runtime_s   = out$runtime_s,
      conf        = built$conf,     # unambiguous in sequential mode
      hmc_applied = hmc_ok
    ))
  }

  ## ---- True parallel over chains (PSOCK) ----
  if (!requireNamespace("parallel", quietly = TRUE)) stop("parallel required")
  if (!requireNamespace("nimbleHMC", quietly = TRUE)) stop("nimbleHMC not installed")

  n_cores <- nchains
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(nimble))
    suppressPackageStartupMessages(library(coda))
    suppressPackageStartupMessages(library(nimbleHMC))
    NULL
  })

  parallel::clusterExport(
    cl,
    varlist = c("build_fn", "monitors", "thin", "opts", ".fresh_build", ".run_and_collect"),
    envir = environment()
  )

  ## Seed base
  seed0 <- if (!is.null(opts$seed)) as.integer(opts$seed) else 1L
  parallel::clusterSetRNGStream(cl, iseed = seed0)

  worker <- function(chain_id, build_fn, niter, nburnin, thin, monitors, opts, seed0) {

    chain_id <- as.integer(chain_id)
    set.seed(seed0 + chain_id)

    built <- .fresh_build(build_fn, chain_id = chain_id, monitors = monitors, thin = thin, opts = opts)

    hmc_ok <- TRUE
    tryCatch({
      nimbleHMC::configureHMC(built$conf, model = built$model)
    }, error = function(e) {
      message("configureHMC failed (chain ", chain_id, "): ", e$message)
      hmc_ok <<- FALSE
    })

    uns <- try(built$conf$getUnsampledNodes(), silent = TRUE)
    if (!inherits(uns, "try-error") && length(uns)) {
      for (u in uns) built$conf$addSampler(u, type = "slice")
    }

    mcmc  <- nimble::buildMCMC(built$conf)
    cmcmc <- nimble::compileNimble(mcmc, project = built$cmodel, resetFunctions = TRUE)

    out <- .run_and_collect(cmcmc, niter = niter, nburnin = nburnin, thin = thin, nchains = 1L)

    list(
      samples     = out$samples,
      samples2    = out$samples2,
      runtime_s   = out$runtime_s,
      conf        = built$conf,  # per-chain conf
      hmc_applied = hmc_ok
    )
  }

  ans <- parallel::parLapply(
    cl,
    X = seq_len(nchains),
    fun = function(ch, build_fn, niter, nburnin, thin, monitors, opts, seed0) {
      worker(ch, build_fn, niter, nburnin, thin, monitors, opts, seed0)
    },
    build_fn  = build_fn,
    niter     = niter,
    nburnin   = nburnin,
    thin      = thin,
    monitors  = monitors,
    opts      = opts,
    seed0     = seed0
  )

  ## Aggregate samples
  samples_by_chain  <- lapply(ans, `[[`, "samples")
  samples2_by_chain <- lapply(ans, `[[`, "samples2")
  runtimes          <- vapply(ans, `[[`, numeric(1), "runtime_s")
  hmc_flags         <- vapply(ans, `[[`, logical(1), "hmc_applied")
  conf_by_chain     <- lapply(ans, `[[`, "conf")

  samples_out <- try(coda::mcmc.list(samples_by_chain), silent = TRUE)
  if (inherits(samples_out, "try-error")) samples_out <- samples_by_chain

  samples2_out <- try(coda::mcmc.list(samples2_by_chain), silent = TRUE)
  if (inherits(samples2_out, "try-error")) samples2_out <- samples2_by_chain

  ## Return: keep API (`conf`) but make it unambiguous:
  ## - conf = NULL in parallel
  ## - attach conf_by_chain as an attribute (non-breaking) + explicit element
  res <- list(
    samples     = samples_out,
    samples2    = samples2_out,
    runtime_s   = max(runtimes),
    conf        = NULL,                 # unambiguous in parallel mode
    hmc_applied = all(hmc_flags)
  )

  attr(res, "conf_by_chain") <- conf_by_chain
  res$conf_by_chain <- conf_by_chain   # explicit (addition is backward-compatible)

  res
}

