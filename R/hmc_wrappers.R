#' Configure and run HMC/NUTS safely
#'
#' Sets up an HMC/NUTS configuration via \pkg{nimbleHMC}, compiles the MCMC,
#' runs it, and returns basic diagnostics.
#'
#' @param build_fn A function that builds the model and returns at least
#'   a list with components \code{model} and \code{conf}.
#' @param niter,nburnin,thin Integers; MCMC iterations, burn-in, and thinning.
#' @param monitors Character vector of node names to monitor.
#' @param nchains Integer; number of chains to run.
#' @param out_dir Character or \code{NULL}; if not \code{NULL}, figures are saved there.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{conf}{The configured MCMC configuration.}
#'     \item{res}{Raw run output (samples, runtime, …).}
#'     \item{diag_tbl}{A data frame of diagnostics (e.g., R-hat, ESS/s).}
#'   }
#'
#' @examples
#' \dontrun{
#'   out <- configure_hmc_safely(
#'     build_fn = my_build_fn, niter = 2000, nburnin = 1000, thin = 1,
#'     monitors = c("beta[1]", "sigma"), nchains = 4, out_dir = "outputs"
#'   )
#'   out$diag_tbl
#' }
#' @export
configure_hmc_safely <- function(build_fn,
                                 niter, nburnin, thin,
                                 monitors, nchains,
                                 out_dir = NULL) {
  b <- .fresh_build(build_fn, monitors = monitors, thin = thin)
  conf <- b$conf

  # Configure NUTS/HMC
  nimbleHMC::configureHMC(conf, model = b$model)

  # Catch unsampled non-likelihood nodes (add slice)
  uns <- try(conf$getUnsampledNodes(), silent = TRUE)
  if (!inherits(uns, "try-error") && length(uns)) {
    uns <- uns[!grepl("^logLik(\\[.*\\])?$|log_?lik|logdens|lpdf", uns, perl = TRUE, ignore.case = TRUE)]
    for (u in uns) conf$addSampler(u, type = "slice")
  }

  # Compile & run
  cmcmc <- .compile_mcmc_with_build(conf, b, reset = TRUE, show = FALSE)
  out   <- .run_and_collect(cmcmc, niter = niter, nburnin = nburnin, thin = thin, nchains = nchains)

  ml    <- as_mcmc_list_sop(out$samples, out$samples2, drop_loglik = FALSE, thin = thin)
  dg    <- compute_diag_from_mcmc(ml, runtime_s = out$runtime_s)

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    .plot_rhat_bar(dg, nodes = colnames(as.matrix(ml[[1L]])),
                   out_file = file.path(out_dir, "rhat_bar.png"))
  }

  list(conf = conf, res = out, diag_tbl = dg)
}
