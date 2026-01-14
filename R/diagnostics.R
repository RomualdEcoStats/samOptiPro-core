# diagnostics.R -- target-level diagnostics (time + step proxy)

#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x


## ---- Exported functions (11) ----

#' Build a standard diagnostics table from MCMC samples
#'
#' Constructs a basic diagnostics table from an \code{mcmc.list}, including:
#' \itemize{
#'   \item \code{target}: parameter name,
#'   \item \code{ESS}: minimum effective sample size across chains,
#'   \item \code{AE_ESS_per_it}: algorithmic efficiency (ESS / iterations),
#'   \item \code{ESS_per_sec}: computational efficiency (ESS / runtime),
#'   \item \code{time_s_per_ESS}: seconds per effective sample,
#'   \item \code{Rhat}: Gelman-Rubin convergence diagnostic,
#'   \item \code{Family}: top-level family (prefix before the first \code{[}).
#' }
#'
#' @param samples An \code{mcmc.list} (or a single \code{mcmc} / \code{matrix}).
#' @param runtime_s Numeric scalar; wall-clock runtime in seconds.
#'
#' @return A \code{data.frame} named \code{diag_tbl}.
#' @export
compute_diag_from_mcmc <- function(samples, runtime_s) {
  stopifnot(is.numeric(runtime_s), length(runtime_s) == 1L, is.finite(runtime_s))

  ## Input normalization
  if (inherits(samples, "mcmc"))  samples <- coda::mcmc.list(samples)
  if (is.matrix(samples))         samples <- coda::mcmc.list(coda::mcmc(samples))
  if (!inherits(samples, "mcmc.list")) {
    stop("compute_diag_from_mcmc: 'samples' must be mcmc.list, mcmc or matrix.")
  }

  ## mcpar / number of kept iterations (first chain as reference)
  mcpar <- attr(samples[[1]], "mcpar")
  n_keep <- if (!is.null(mcpar) && length(mcpar) >= 3L) {
    as.integer((mcpar[2] - mcpar[1]) / mcpar[3] + 1L)
  } else {
    nrow(samples[[1]])
  }

  mats <- lapply(samples, function(m) as.matrix(m))
  params <- Reduce(intersect, lapply(mats, colnames))
  if (is.null(params) || !length(params)) {
    stop("compute_diag_from_mcmc: no common parameters between chains.")
  }
  params <- as.character(params)

  ess_min <- numeric(length(params)); names(ess_min) <- params
  rhat    <- rep(NA_real_, length(params)); names(rhat) <- params

  ## ESS (min over chains) and Rhat per parameter
  for (p in params) {
    ess_chain <- vapply(
      mats,
      function(M) {
        v <- M[, p]
        as.numeric(tryCatch(coda::effectiveSize(v), error = function(e) NA_real_))
      },
      numeric(1)
    )
    ess_min[p] <- suppressWarnings(min(ess_chain, na.rm = TRUE))

    ml_param <- coda::mcmc.list(lapply(seq_along(samples), function(k) {
      m <- as.matrix(samples[[k]])
      coda::mcmc(
        m[, p, drop = FALSE],
        start = if (!is.null(mcpar)) mcpar[1] else 1,
        end   = if (!is.null(mcpar)) mcpar[2] else (nrow(m) + (mcpar[1] %||% 1) - 1),
        thin  = if (!is.null(mcpar)) mcpar[3] else 1
      )
    }))
    rhat[p] <- tryCatch({
      gd <- coda::gelman.diag(ml_param, autoburnin = FALSE)$psrf
      as.numeric(gd[1])
    }, error = function(e) NA_real_)
  }

  ess_min[!is.finite(ess_min)] <- NA_real_

  ## AE = ESS / iterations; ESS/s; time per ESS
  AE              <- ess_min / max(n_keep, 1L)
  ESS_per_sec     <- if (is.finite(runtime_s) && runtime_s > 0) ess_min / runtime_s else NA_real_
  time_s_per_ESS  <- runtime_s / pmax(ess_min, .Machine$double.eps)

  out <- data.frame(
    target         = params,
    ESS            = as.numeric(ess_min),
    AE_ESS_per_it  = as.numeric(AE),
    ESS_per_sec    = as.numeric(ESS_per_sec),
    time_s_per_ESS = as.numeric(time_s_per_ESS),
    Rhat           = as.numeric(rhat),
    stringsAsFactors = FALSE
  )
  out$Family <- sub("\\[.*", "", out$target)

  ## Filtre vectoriel: ESS finie strictement positive
  keep <- is.finite(out$ESS) & (out$ESS > 0)
  out  <- out[keep, , drop = FALSE]

  rownames(out) <- NULL
  out
}


#' Scalable diagnostics from MCMC samples (block-wise)
#'
#' Efficiently computes convergence and efficiency diagnostics (ESS, R-hat,
#' algorithmic efficiency, and computational efficiency) from large
#' \code{mcmc.list} objects, designed for very high-dimensional hierarchical models
#' (for example, more than 90,000 parameters). The function operates in column
#' blocks to control memory usage and can optionally leverage the \pkg{posterior}
#' package for faster and rank-normalized diagnostics.
#'
#' @details
#' This implementation is optimized for large Bayesian stock assessment or
#' life-cycle models (for example, WGNAS, GEREM, Scorff LCM) where standard
#' \pkg{coda} routines become memory bound. It avoids unnecessary matrix copies,
#' truncates all chains to the shortest length, and computes diagnostics in
#' column blocks to maintain stability under limited RAM.
#'
#' \strong{Formulas:}
#' \deqn{ESS_{worst} = \min_c ESS_c(\theta_i)}{}
#' \deqn{ESS_{total} = ESS(\text{pooled chains})}{}
#' Algorithmic efficiency: \eqn{AE = ESS / n_{iter}}.\\
#' Computational efficiency: \eqn{CE = ESS / t_{run}}.
#'
#' @param samples An \code{mcmc.list}, \code{mcmc}, or numeric matrix of MCMC samples.
#' @param runtime_s Numeric scalar; runtime in seconds (total or per chain).
#' @param compute_rhat Character, one of \code{"none"}, \code{"classic"},
#'   \code{"split"}, or \code{"both"} (default \code{"none"}).
#' @param ess_for Character, one of \code{"both"}, \code{"worst"}, or \code{"total"}.
#' @param ignore_patterns Character vector of regular expressions used to drop
#'   parameters by name (for example \code{"^lifted_"}, \code{"^logProb_"}).
#' @param cols_by Integer; number of columns per processing block
#'   (typically \eqn{\ge 1000}). Can be tuned automatically via
#'   \code{target_block_ram_gb}.
#' @param warn_mem_gb Numeric; memory threshold (GiB) above which a warning
#'   is issued.
#' @param step_timeout_s Numeric; per-block time limit in seconds.
#' @param runtime_is_total Logical; if \code{FALSE}, \code{runtime_s} is assumed
#'   to be per chain and is multiplied by the number of chains.
#' @param use_posterior Character; \code{"never"} or \code{"if_available"} to use
#'   the \pkg{posterior} package for ESS and R-hat.
#' @param target_block_ram_gb Numeric; if non-NA, automatically computes
#'   \code{cols_by} to target this RAM usage (GiB) per block.
#' @param verbose Logical; if \code{TRUE}, display progress messages.
#'
#' @return
#' A \code{data.frame} with one row per parameter and the following columns:
#' \describe{
#'   \item{target}{Parameter name.}
#'   \item{ESS_worst}{Minimum ESS across chains.}
#'   \item{ESS_total}{Combined ESS across all chains.}
#'   \item{AE_worst, AE_total}{Algorithmic efficiencies.}
#'   \item{ESS_per_sec_worst, ESS_per_sec_total}{Computational efficiencies.}
#'   \item{time_s_per_ESS_worst, time_s_per_ESS_total}{Seconds per effective sample.}
#'   \item{Rhat_classic, Rhat_split}{Gelman-Rubin diagnostics (classic / split).}
#'   \item{Family}{Top-level node family (extracted from \code{target}).}
#' }
#'
#' @note
#' Stability empirically validated for models with at least
#' \eqn{90,000} parameters and at most \eqn{10} chains.
#'
#' @seealso
#' \code{\link[coda]{effectiveSize}}, \code{\link[posterior]{ess_bulk}},
#' \code{\link[posterior]{rhat}}, Vehtari et al. (2021),
#' *Bayesian Analysis* 16(2):667--718.
#'
#' @examples
#' \dontrun{
#' res_diag <- compute_diag_from_mcmc_vect(
#'   samples   = my_mcmc,
#'   runtime_s = 5400,
#'   compute_rhat         = "both",
#'   ess_for              = "both",
#'   target_block_ram_gb  = 2
#' )
#' head(res_diag)
#' }
#'
#' @keywords internal
#' @noRd
compute_diag_from_mcmc_vect <- function(
    samples, runtime_s,
    compute_rhat        = c("none", "classic", "split", "both"),
    ess_for             = c("both", "worst", "total"),
    ignore_patterns     = c("^lifted_", "^logProb_"),
    cols_by             = 5000L,
    warn_mem_gb         = 10,
    step_timeout_s      = Inf,
    runtime_is_total    = TRUE,
    use_posterior       = c("never", "if_available"),
    target_block_ram_gb = NA_real_,
    verbose             = TRUE
) {
  t0 <- proc.time()[3]
  stopifnot(is.numeric(runtime_s), length(runtime_s) == 1L, is.finite(runtime_s))
  compute_rhat   <- match.arg(compute_rhat)
  ess_for        <- match.arg(ess_for)
  use_posterior  <- match.arg(use_posterior)

  vmsg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  ## Input normalization
  if (inherits(samples, "mcmc"))  samples <- coda::mcmc.list(samples)
  if (is.matrix(samples))         samples <- coda::mcmc.list(coda::mcmc(samples))
  if (!inherits(samples, "mcmc.list")) {
    stop("compute_diag_from_mcmc_vect: 'samples' must be mcmc.list, mcmc or matrix.")
  }

  m <- length(samples)
  if (m < 1L) stop("compute_diag_from_mcmc_vect: empty mcmc.list.")

  mats0   <- lapply(samples, function(x) as.matrix(x))
  params  <- Reduce(intersect, lapply(mats0, colnames))
  if (!length(params)) stop("compute_diag_from_mcmc_vect: no common parameters between chains.")
  mats0   <- lapply(mats0, function(M) M[, params, drop = FALSE])

  ## Filter by ignore_patterns
  if (!is.null(ignore_patterns) && length(ignore_patterns)) {
    bad  <- Reduce(`|`, lapply(ignore_patterns, function(p) grepl(p, params, perl = TRUE)))
    keep <- !bad
    if (!any(keep)) stop("compute_diag_from_mcmc_vect: all parameters removed by 'ignore_patterns'.")
    params <- params[keep]
    mats0  <- lapply(mats0, function(M) M[, keep, drop = FALSE])
  }

  P     <- length(params)
  n_vec <- vapply(mats0, nrow, integer(1))
  n_min <- min(n_vec)

  ## Single truncation to the shortest chain length
  mats <- lapply(mats0, function(M) {
    if (nrow(M) > n_min) M[seq_len(n_min), , drop = FALSE] else M
  })
  rm(mats0); gc(FALSE)

  ## Memory estimation (very rough upper bound)
  approx_gb <- (sum(n_vec) * P * 8) / (1024^3)
  if (is.finite(warn_mem_gb) && approx_gb > warn_mem_gb) {
    vmsg("[compute_diag_from_mcmc_vect] Warning: approximate size ~%.1f GiB; consider reducing monitors/thin.", approx_gb)
  }

  ## Block size choice
  if (is.finite(target_block_ram_gb) &&
      !is.na(target_block_ram_gb) &&
      target_block_ram_gb > 0) {
    safety <- 1.3
    bytes_target <- target_block_ram_gb * (1024^3)
    ## matrices used in ESS_total: m chains * n_min rows * 8 bytes per double
    ## cols_by ~= bytes_target / (8 * n_min * m * safety)
    cols_by_auto <- floor(bytes_target / (8 * n_min * max(1, m) * safety))
    cols_by_auto <- max(256L, cols_by_auto)
    if (isTRUE(verbose)) {
      vmsg("Automatic blocking: cols_by ~ %d (target %.1f GiB).", cols_by_auto, target_block_ram_gb)
    }
    cols_by <- cols_by_auto
  } else {
    cols_by <- as.integer(max(3000L, cols_by))
  }

  nb <- as.integer(ceiling(P / cols_by))
  block_idx <- split(seq_len(P), ceiling(seq_len(P) / cols_by))

  ## Helpers
  safe_time_check <- function(t_start, limit_s, label) {
    if (is.finite(limit_s) && (proc.time()[3] - t_start) > limit_s) {
      stop(sprintf("%s: timeout after %.1fs.", label, limit_s))
    }
  }

  ## Optional robust diagnostics via posterior
  have_posterior <- (use_posterior == "if_available" &&
                       requireNamespace("posterior", quietly = TRUE))

  ## ---------- ESS_worst ----------
  do_worst  <- ess_for %in% c("both", "worst")
  ESS_worst <- rep(NA_real_, P)
  if (do_worst) {
    vmsg("ESS_worst: %d parameters in %d blocks.", P, nb)
    t_w <- proc.time()[3]
    for (b in seq_len(nb)) {
      ix <- block_idx[[b]]
      if (have_posterior && m >= 1L) {
        ## compute ESS chain-by-chain using posterior::ess_bulk
        ess_block <- matrix(NA_real_, nrow = length(ix), ncol = m)
        for (k in seq_len(m)) {
          Mk <- mats[[k]][, ix, drop = FALSE]
          dm <- posterior::as_draws_matrix(Mk)
          ess_block[, k] <- posterior::ess_bulk(dm)
        }
      } else {
        ess_block <- matrix(NA_real_, nrow = length(ix), ncol = m)
        for (k in seq_len(m)) {
          Mk <- mats[[k]][, ix, drop = FALSE]
          ess_block[, k] <- coda::effectiveSize(coda::mcmc(Mk))
        }
      }
      ESS_worst[ix] <- apply(
        ess_block, 1L,
        function(v) suppressWarnings(min(v, na.rm = TRUE))
      )
      if (b %% 20 == 0) gc(FALSE)
      if (b %% 10 == 0) vmsg("  Block %d/%d done.", b, nb)
      safe_time_check(t_w, step_timeout_s, "ESS_worst")
    }
  }

  ## ---------- ESS_total ----------
  do_total  <- (m >= 2L) && ess_for %in% c("both", "total")
  ESS_total <- rep(NA_real_, P)
  if (do_total) {
    vmsg("ESS_total: combining %d chains in %d blocks.", m, nb)
    t_t <- proc.time()[3]
    for (b in seq_len(nb)) {
      ix <- block_idx[[b]]
      if (have_posterior) {
        ## bind chains column-wise, then convert to draws_array
        Mbind <- do.call(cbind, lapply(mats, function(Mk) Mk[, ix, drop = FALSE]))
        it <- n_min; ch <- m; vv <- length(ix)
        arr <- array(NA_real_, dim = c(it, ch, vv))
        for (j in seq_len(vv)) {
          colidx <- j + (0:(ch - 1)) * vv
          arr[, , j] <- as.matrix(Mbind[, colidx, drop = FALSE])
        }
        da <- posterior::as_draws_array(arr)
        ESS_total[ix] <- posterior::ess_bulk(da)
      } else {
        trunc_list <- coda::mcmc.list(
          lapply(seq_len(m), function(k) coda::mcmc(mats[[k]][, ix, drop = FALSE]))
        )
        ESS_total[ix] <- coda::effectiveSize(trunc_list)
      }
      if (b %% 20 == 0) gc(FALSE)
      if (b %% 10 == 0) vmsg("  Block %d/%d done.", b, nb)
      safe_time_check(t_t, step_timeout_s, "ESS_total")
    }
  }

  ## ---------- R-hat ----------
  Rhat_classic <- Rhat_split <- rep(NA_real_, P)
  if (m >= 2L && compute_rhat %in% c("classic", "both", "split")) {
    vmsg("R-hat (%s): %d blocks.", compute_rhat, nb)
    t_r <- proc.time()[3]
    use_rank <- have_posterior
    for (b in seq_len(nb)) {
      ix <- block_idx[[b]]

      ## classic R-hat
      if (compute_rhat %in% c("classic", "both")) {
        if (use_rank) {
          Mbind <- do.call(cbind, lapply(mats, function(Mk) Mk[, ix, drop = FALSE]))
          it <- n_min; ch <- m; vv <- length(ix)
          arr <- array(NA_real_, dim = c(it, ch, vv))
          for (j in seq_len(vv)) {
            colidx <- j + (0:(ch - 1)) * vv
            arr[, , j] <- as.matrix(Mbind[, colidx, drop = FALSE])
          }
          da <- posterior::as_draws_array(arr)
          Rhat_classic[ix] <- posterior::rhat(da)
        } else {
          means_mat <- matrix(NA_real_, nrow = m, ncol = length(ix))
          vars_mat  <- matrix(NA_real_, nrow = m, ncol = length(ix))
          for (k in seq_len(m)) {
            Mk <- mats[[k]][, ix, drop = FALSE]
            muk  <- colMeans(Mk)
            ssq  <- colSums(Mk * Mk)
            vark <- (ssq - n_min * muk * muk) / max(1, (n_min - 1))
            means_mat[k, ] <- muk
            vars_mat[k, ]  <- pmax(vark, 0)
          }
          W      <- colMeans(vars_mat)
          mu_bar <- colMeans(means_mat)
          sumsq  <- colSums(means_mat * means_mat)
          var_mu <- (sumsq - m * mu_bar * mu_bar) / max(1, (m - 1))
          B      <- n_min * pmax(var_mu, 0)
          var_hat <- ((n_min - 1) / max(1, n_min)) * W + (B / max(1, n_min))
          ok <- (W > 0) & is.finite(W) & is.finite(var_hat)
          tmp <- rep(NA_real_, length(ix))
          tmp[ok] <- sqrt(pmax(var_hat[ok] / W[ok], 0))
          Rhat_classic[ix] <- tmp
        }
      }

      ## split R-hat
      if (compute_rhat %in% c("split", "both")) {
        h <- n_min %/% 2L
        if (h >= 2L) {
          if (use_rank) {
            M1 <- lapply(mats, function(Mk) Mk[1:h, ix, drop = FALSE])
            M2 <- lapply(mats, function(Mk) Mk[(n_min - h + 1):n_min, ix, drop = FALSE])
            mats_split <- c(M1, M2)
            Mbind <- do.call(cbind, lapply(mats_split, function(Mk) Mk))
            it <- h; ch <- 2 * m; vv <- length(ix)
            arr <- array(NA_real_, dim = c(it, ch, vv))
            for (j in seq_len(vv)) {
              colidx <- j + (0:(ch - 1)) * vv
              arr[, , j] <- as.matrix(Mbind[, colidx, drop = FALSE])
            }
            da <- posterior::as_draws_array(arr)
            Rhat_split[ix] <- posterior::rhat(da)
          } else {
            halves_means <- list()
            halves_vars  <- list()
            for (k in seq_len(m)) {
              M1 <- mats[[k]][1:h, ix, drop = FALSE]
              M2 <- mats[[k]][(n_min - h + 1):n_min, ix, drop = FALSE]
              for (Mm in list(M1, M2)) {
                mu  <- colMeans(Mm)
                ssq <- colSums(Mm * Mm)
                v   <- (ssq - nrow(Mm) * mu * mu) / max(1, (nrow(Mm) - 1))
                halves_means[[length(halves_means) + 1L]] <- mu
                halves_vars[[length(halves_vars) + 1L]]   <- pmax(v, 0)
              }
            }
            means_mat <- do.call(rbind, halves_means)
            vars_mat  <- do.call(rbind, halves_vars)
            W      <- colMeans(vars_mat)
            mu_bar <- colMeans(means_mat)
            sumsq  <- colSums(means_mat * means_mat)
            var_mu <- (sumsq - nrow(means_mat) * mu_bar * mu_bar) /
              max(1, (nrow(means_mat) - 1))
            n_eff_block <- h
            B <- n_eff_block * pmax(var_mu, 0)
            var_hat <- ((n_eff_block - 1) / max(1, n_eff_block)) * W +
              (B / max(1, n_eff_block))
            ok <- (W > 0) & is.finite(W) & is.finite(var_hat)
            tmp <- rep(NA_real_, length(ix))
            tmp[ok] <- sqrt(pmax(var_hat[ok] / W[ok], 0))
            Rhat_split[ix] <- tmp
          }
        } else {
          Rhat_split[ix] <- NA_real_
        }
      }

      if (b %% 20 == 0) gc(FALSE)
      if (b %% 10 == 0) vmsg("  Block %d/%d done.", b, nb)
      safe_time_check(t_r, step_timeout_s, "R-hat")
    }
  }

  ## ---------- AE / CE ----------
  total_runtime <- if (runtime_is_total) runtime_s else runtime_s * max(1L, m)

  AE_worst <- ESS_worst / max(1L, n_min)
  AE_total <- if (m >= 2L && do_total) {
    ESS_total / max(1L, n_min * m)
  } else {
    rep(NA_real_, P)
  }

  ESSps_worst <- if (total_runtime > 0) ESS_worst / total_runtime else rep(NA_real_, P)
  ESSps_total <- if (total_runtime > 0 && m >= 2L && do_total) {
    ESS_total / total_runtime
  } else {
    rep(NA_real_, P)
  }

  t_per_ESS_worst <- total_runtime / pmax(ESS_worst, .Machine$double.eps)
  t_per_ESS_total <- if (m >= 2L && do_total) {
    total_runtime / pmax(ESS_total, .Machine$double.eps)
  } else {
    rep(NA_real_, P)
  }

  out <- data.frame(
    target               = params,
    ESS_worst            = as.numeric(ESS_worst),
    ESS_total            = as.numeric(ESS_total),
    AE_worst             = as.numeric(AE_worst),
    AE_total             = as.numeric(AE_total),
    ESS_per_sec_worst    = as.numeric(ESSps_worst),
    ESS_per_sec_total    = as.numeric(ESSps_total),
    time_s_per_ESS_worst = as.numeric(t_per_ESS_worst),
    time_s_per_ESS_total = as.numeric(t_per_ESS_total),
    Rhat_classic         = as.numeric(
      if (compute_rhat %in% c("classic", "both")) Rhat_classic else rep(NA_real_, P)
    ),
    Rhat_split           = as.numeric(
      if (compute_rhat %in% c("split", "both")) Rhat_split else rep(NA_real_, P)
    ),
    stringsAsFactors     = FALSE
  )
  out$Family <- sub("\\[.*", "", out$target)

  rownames(out) <- NULL
  vmsg("Finished in %.1fs.", proc.time()[3] - t0)
  out
}


#' Configure and run HMC/NUTS safely (global / subset / auto)
#'
#' Minimal, strict wrapper around NIMBLE + nimbleHMC that:
#' \enumerate{
#'   \item builds the model,
#'   \item configures NUTS either globally (model signature),
#'         on a subset of nodes (conf + model + nodes signature),
#'         or in auto mode,
#'   \item compiles the MCMC,
#'   \item runs MCMC.
#' }
#'
#' No per-node fallback such as \code{addSampler("NUTS")} is attempted:
#' on error, the function stops.
#'
#' @param build_fn Function; model builder (e.g. \code{build_M}) used to
#'   create the NIMBLE model before configuring HMC/NUTS.
#'   If \code{.fresh_build()} exists, it will be called as
#'   \code{.fresh_build(build_fn, monitors, thin)}.
#'   Otherwise, \code{build_fn()} must return a list with at least
#'   \code{model} and \code{conf} (NIMBLE objects).
#' @param niter Integer. Total number of iterations.
#' @param nburnin Integer (>= 0). Number of burn-in iterations (0 allowed).
#' @param thin Integer (>= 1). Thinning interval.
#' @param nchains Integer (>= 1). Number of chains.
#' @param monitors Character vector or \code{NULL}. Passed through to
#'   building and/or global HMC configuration.
#' @param nuts_mode Character. One of \code{c("all","subset","auto","none")}:
#'   \describe{
#'     \item{"all"}{Global NUTS via \code{nimbleHMC::configureHMC(model, monitors = ..., enableWAIC = ..., buildDerivs = ...)}.}
#'     \item{"subset"}{NUTS only on \code{nuts_nodes} via
#'       \code{nimbleHMC::configureHMC(conf, model, nodes = ..., ...)}.}
#'     \item{"auto"}{Let \code{nimbleHMC::configureHMC(conf, model, ...)} choose targets automatically.}
#'     \item{"none"}{Do not modify the existing configuration.}
#'   }
#' @param nuts_nodes Character vector of node names (can include indexed forms
#'   such as \code{"beta[]"}) when \code{nuts_mode = "subset"}.
#' @param enable_WAIC Logical. Only relevant for \code{"all"} (model-signature path).
#' @param buildDerivs Logical. Passed to \code{nimbleHMC::configureHMC()}.
#' @param compiler_cores Integer. Passed to
#'   \code{nimbleOptions(numCompilerCores = ...)}.
#' @param show_compiler_output Logical. Passed to
#'   \code{nimbleOptions(showCompilerOutput = ...)}.
#' @param project_name Character or \code{NULL}. If non-\code{NULL}, sets
#'   \code{options(nimbleProjectName = ...)} to sanitize the C++ cache.
#' @param out_dir Character or \code{NULL}. If non-\code{NULL}, saves a
#'   \code{sessionInfo()} log and, optionally, samples as \code{.rds}.
#' @param inits \code{NULL} or a list of per-chain initial values to pass to
#'   \code{nimble::runMCMC()} (length must match \code{nchains} if provided).
#' @param save_samples Logical. If \code{TRUE} (default) and \code{out_dir} is
#'   non-\code{NULL}, save \code{samples} as \code{.rds}. If \code{FALSE},
#'   skip saving samples.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{conf}: final MCMC configuration (with HMC/NUTS if configured),
#'   \item \code{cmcmc}: compiled MCMC object,
#'   \item \code{samples}: \code{coda::mcmc.list} of posterior samples,
#'   \item \code{runtime_s}: wall time in seconds,
#'   \item \code{build}: the build object returned by \code{build_fn} or
#'         \code{.fresh_build()}.
#' }
#'
#' @export
configure_hmc_safely_bis <- function(build_fn,
                                     niter, nburnin, thin, nchains,
                                     monitors = NULL,
                                     nuts_mode = c("all", "subset", "auto", "none"),
                                     nuts_nodes = NULL,
                                     enable_WAIC = FALSE,
                                     buildDerivs = TRUE,
                                     compiler_cores = max(1L, parallel::detectCores(TRUE) %/% 2L),
                                     show_compiler_output = TRUE,
                                     project_name = NULL,
                                     out_dir = NULL,
                                     inits = NULL,
                                     save_samples = TRUE) {

  ## ---------- Local helpers ----------
  .as_char_vec <- function(x) {
    if (is.null(x)) return(character(0))
    if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
    x <- suppressWarnings(try(as.character(x), silent = TRUE))
    if (inherits(x, "try-error")) return(character(0))
    attributes(x) <- NULL
    unname(x[!is.na(x) & nzchar(x)])
  }

  .is_ll <- function(x) {
    grepl("^logLik(\\[.*\\])?$|log_?lik|logdens|lpdf", x,
          perl = TRUE, ignore.case = TRUE)
  }

  .safe_scalar_targets <- function(model, nodes) {
    nodes <- .as_char_vec(nodes)
    if (!length(nodes)) return(character(0))
    nodes <- nodes[!grepl("^lifted_|^logProb_", nodes)]
    nodes <- nodes[!.is_ll(nodes)]
    expd  <- suppressWarnings(
      try(model$expandNodeNames(nodes, returnScalarComponents = TRUE),
          silent = TRUE)
    )
    .as_char_vec(expd)
  }

  .mkdir <- function(p) {
    if (!is.null(p) && nzchar(p) && !dir.exists(p)) {
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
    }
  }

  .log <- function(...) {
    cat(sprintf("%s %s\n",
                format(Sys.time(), "%H:%M:%S"),
                paste0(..., collapse = "")))
  }

  .head6 <- function(x) paste(utils::head(x, 6), collapse = ", ")

  ## ---------- Checks / options ----------
  nuts_mode <- match.arg(nuts_mode)
  stopifnot(is.function(build_fn))

  stopifnot(all(c(niter, thin, nchains) >= 1L))
  stopifnot(nburnin >= 0L)
  stopifnot(niter > nburnin)

  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("Package 'nimble' is required.")
  }
  if (!requireNamespace("nimbleHMC", quietly = TRUE) && nuts_mode != "none") {
    stop("Package 'nimbleHMC' is required to configure NUTS/HMC.")
  }

  nimble::nimbleOptions(
    showCompilerOutput = isTRUE(show_compiler_output),
    numCompilerCores   = as.integer(max(1L, compiler_cores))
  )
  if (!is.null(project_name) && nzchar(project_name)) {
    options(nimbleProjectName = project_name)
  }

  if (!is.null(out_dir)) {
    .mkdir(out_dir)
    .log("[IO] out_dir = ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
  }

  ## ---------- Build ----------
  if (exists(".fresh_build", mode = "function")) {
    b <- .fresh_build(build_fn, monitors = monitors, thin = thin)
    if (!all(c("model", "conf") %in% names(b))) {
      stop("`.fresh_build()` must return at least `model` and `conf`.")
    }
  } else {
    base <- build_fn()
    if (!(is.list(base) && all(c("model", "conf") %in% names(base)))) {
      stop("Without `.fresh_build`, `build_fn()` must return a list with `model` and `conf`.")
    }
    b <- base
  }
  model <- b$model
  conf  <- b$conf

  ## ---------- HMC/NUTS configuration ----------
  if (identical(nuts_mode, "all")) {
    conf <- tryCatch({
      nimbleHMC::configureHMC(
        model,
        monitors    = monitors,
        enableWAIC  = isTRUE(enable_WAIC),
        buildDerivs = isTRUE(buildDerivs)
      )
    }, error = function(e) {
      stop(sprintf("[HMC/all] configureHMC(model, ...) failed: %s",
                   conditionMessage(e)))
    })
  } else if (identical(nuts_mode, "subset")) {
    nodes <- .safe_scalar_targets(model, nuts_nodes)
    if (!length(nodes)) {
      stop("[HMC/subset] No valid `nuts_nodes` after sanitization/expansion (empty or non-scalar).")
    }
    .log(sprintf("[HMC] subset: %d targets (head: %s)",
                 length(nodes), .head6(nodes)))
    conf <- tryCatch({
      nimbleHMC::configureHMC(
        conf, model,
        nodes       = nodes,
        buildDerivs = isTRUE(buildDerivs)
      )
    }, error = function(e) {
      stop(sprintf("[HMC/subset] configureHMC(conf, model, nodes = ...) failed: %s",
                   conditionMessage(e)))
    })
  } else if (identical(nuts_mode, "auto")) {
    conf <- tryCatch({
      nimbleHMC::configureHMC(
        conf, model,
        buildDerivs = isTRUE(buildDerivs)
      )
    }, error = function(e) {
      stop(sprintf("[HMC/auto] configureHMC(conf, model) failed: %s",
                   conditionMessage(e)))
    })
  } else {
    .log("[HMC] mode = 'none': leaving MCMC configuration unchanged.")
  }

  ## ---------- Safety net: unsampled nodes ----------
  uns <- try(conf$getUnsampledNodes(), silent = TRUE)
  if (!inherits(uns, "try-error") && length(uns)) {
    uns <- uns[!grepl("^logLik(\\[.*\\])?$|log_?lik|logdens|lpdf",
                      uns, perl = TRUE, ignore.case = TRUE)]
    if (length(uns)) {
      .log(sprintf("[safety] add slice on %d unsampled nodes (head: %s)",
                   length(uns), .head6(uns)))
      for (u in uns) {
        conf$addSampler(u, type = "slice")
      }
    }
  }

  ## ---------- Log HMC/NUTS targets ----------
  invisible(try({
    sams   <- conf$getSamplers()
    is_hmc <- vapply(sams, function(s) any(grepl("HMC|NUTS", class(s))), logical(1))
    if (any(is_hmc)) {
      targs <- unique(unlist(lapply(
        sams[is_hmc],
        function(s) model$expandNodeNames(s$target, returnScalarComponents = TRUE)
      )))
      targs <- .as_char_vec(targs)
      .log(sprintf("[HMC] %d HMC/NUTS samplers; ~%d scalar targets (head: %s)",
                   sum(is_hmc), length(targs), .head6(targs)))
    } else {
      .log("[HMC] no HMC/NUTS sampler detected in `conf`.")
    }
  }, silent = TRUE))

  ## ---------- Compile ----------
  if (exists(".compile_mcmc_with_build", mode = "function")) {
    cmcmc <- .compile_mcmc_with_build(conf, b, reset = TRUE, show = FALSE)
  } else {
    mcmc   <- nimble::buildMCMC(conf)
    cmodel <- nimble::compileNimble(model, resetFunctions = TRUE)
    cmcmc  <- nimble::compileNimble(mcmc, project = cmodel)
  }

  ## ---------- Run MCMC ----------
  t0 <- proc.time()
  samples <- nimble::runMCMC(
    cmcmc,
    niter            = as.integer(niter),
    nburnin          = as.integer(nburnin),
    thin             = as.integer(thin),
    nchains          = as.integer(nchains),
    inits            = inits,
    samplesAsCodaMCMC = TRUE
  )
  runtime_s <- as.numeric((proc.time() - t0)[["elapsed"]])
  .log(sprintf("[run] done in %.1f s  (niter = %d, burnin = %d, thin = %d, chains = %d)",
               runtime_s, niter, nburnin, thin, nchains))

  ## ---------- IO ----------
  if (!is.null(out_dir)) {
    ts    <- format(Sys.time(), "%Y%m%d_%H%M%S")
    f_rds <- file.path(out_dir, sprintf("samples_nuts_%s.rds", ts))
    f_log <- file.path(out_dir, sprintf("runlog_%s.txt",      ts))

    if (isTRUE(save_samples)) {
      suppressWarnings(try(saveRDS(samples, f_rds), silent = TRUE))
      .log("[IO] saved RDS: ", f_rds)
    } else {
      .log("[IO] save_samples = FALSE: skipping RDS samples.")
    }

    suppressWarnings(try(
      cat(capture.output(sessionInfo()), sep = "\n", file = f_log),
      silent = TRUE
    ))
    .log("[IO] saved log: ", f_log)
  }

  invisible(list(
    conf       = conf,
    cmcmc      = cmcmc,
    samples    = samples,
    runtime_s  = runtime_s,
    build      = b
  ))
}
#' Diagnose model structure, dependencies, and sampler time (parameter and family levels)
#'
#' Inspect a NIMBLE model to extract the universe of nodes, classify stochastic
#' vs deterministic nodes, compute downstream dependencies per node, map
#' configured samplers to their target nodes, and (optionally) profile
#' per-sampler run time. Produces tidy tables and publication-quality figures
#' at both the \emph{parameter} level and the \emph{family} level (where
#' "family" = base variable name before any indices, replacing indices by
#' square brackets, for example \code{"beta[]"}).
#'
#' Key features:
#' \itemize{
#'   \item Robust filtering of nodes (ignored patterns, removal list).
#'   \item Downstream dependency counts per node and per family.
#'   \item Per-sampler times aggregated to parameters and families.
#'   \item Optional auto-profiling of samplers via a short MCMC run.
#'   \item Non-regressive: original per-parameter plots are preserved by default.
#' }
#'
#' @param model A compiled or uncompiled \code{nimbleModel} (required).
#' @param include_data Logical; include data nodes when enumerating nodes
#'   (default \code{FALSE}).
#' @param removed_nodes Character vector of nodes to exclude explicitly
#'   (default \code{NULL}).
#' @param ignore_patterns Character vector of regular expressions used to
#'   exclude nodes (for example \code{"^lifted_"}, \code{"^logProb_"}).
#' @param make_plots Logical; if \code{TRUE}, generate ggplot objects and
#'   optionally save them (default \code{TRUE}).
#' @param output_dir Directory where figures/CSVs will be saved; if \code{NULL},
#'   nothing is written to disk (default \code{NULL}).
#' @param save_csv Logical; if \code{TRUE}, write CSV exports (dependencies per
#'   node, family tables) (default \code{FALSE}).
#' @param node_of_interest Optional character scalar naming a node to highlight
#'   or subset (reserved for user logic) (default \code{NULL}).
#' @param sampler_times Optional numeric vector of per-sampler times aligned
#'   with \code{nimble::configureMCMC(model)$getSamplers()}.
#' @param sampler_times_unit Character label for time axis
#'   (for example \code{"seconds"}, \code{"ms"}) (default \code{"seconds"}).
#' @param auto_profile Logical; if \code{TRUE} and \code{sampler_times} is
#'   \code{NULL}, profile sampler times automatically via a short MCMC run
#'   (default \code{TRUE}).
#' @param profile_niter Integer; iterations used by the auto-profiler
#'   (default \code{3000L}).
#' @param profile_burnin Integer; burn-in iterations for auto-profiler
#'   (kept for API compatibility; not used internally) (default \code{100L}).
#' @param profile_thin Integer; thinning for auto-profiler (kept for API
#'   compatibility; not used internally) (default \code{1L}).
#' @param profile_seed Optional integer seed for reproducibility
#'   (default \code{NULL}).
#' @param np Proportion in \code{(0, 1]} used elsewhere for "worst sampler"
#'   selection (kept for API compatibility) (default \code{0.10}).
#' @param by_family Logical; if \code{TRUE}, compute and plot family-level
#'   summaries in addition to parameter-level summaries (default \code{TRUE}).
#' @param family_stat One of \code{c("median", "mean", "sum")}; summary
#'   statistic for family-level aggregation (default \code{"median"}).
#' @param time_normalize One of \code{c("none", "per_node")}; if
#'   \code{"per_node"}, divide family time by the number of distinct nodes
#'   in the family (default \code{"none"}).
#' @param only_family_plots Logical; if \code{TRUE}, only family-level figures
#'   are exported (parameter-level plots are not written) (default \code{FALSE}).
#' @param ... Additional arguments forwarded to \code{nimble::configureMCMC()}.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{dependencies_df}: data.frame of (node, dependency) pairs.
#'   \item \code{dep_counts}: data.frame with per-parameter downstream
#'         dependency counts.
#'   \item \code{samplers_df}: data.frame listing samplers and their target
#'         nodes (list-column).
#'   \item \code{per_param_times}: data.frame with per-parameter aggregated
#'         sampler time.
#'   \item \code{deps_df}: tidy data.frame used for parameter-level dependency
#'         plotting.
#'   \item \code{sampler_df}: tidy data.frame used for parameter-level time
#'         plotting.
#'   \item \code{fam_deps_df}: family-level dependency summary
#'         (one statistic per family).
#'   \item \code{fam_time_df}: family-level time summary (optionally normalized).
#'   \item \code{plots}: list of ggplot objects (some may be \code{NULL}):
#'         \code{plot_dependencies}, \code{plot_sampler_time},
#'         \code{plot_combined}, \code{plot_dependencies_family},
#'         \code{plot_sampler_time_family}, \code{plot_combined_family}.
#' }
#'
#' @examples
#' \dontrun{
#' res <- diagnose_model_structure(
#'   model            = my_nimble_model,
#'   make_plots       = TRUE,
#'   output_dir       = "outputs/diagnostics",
#'   save_csv         = TRUE,
#'   by_family        = TRUE,
#'   family_stat      = "median",
#'   time_normalize   = "per_node",
#'   only_family_plots = FALSE
#' )
#' }
#'
#' @export
diagnose_model_structure <- function(model,
                                     include_data        = FALSE,
                                     removed_nodes       = NULL,
                                     ignore_patterns     = c("^lifted_", "^logProb_"),
                                     make_plots          = TRUE,
                                     output_dir          = NULL,
                                     save_csv            = FALSE,
                                     node_of_interest    = NULL,
                                     sampler_times       = NULL,
                                     sampler_times_unit  = "seconds",
                                     auto_profile        = TRUE,
                                     profile_niter       = 3000L,
                                     profile_burnin      = 100L,  # kept for compat (unused)
                                     profile_thin        = 1L,    # kept for compat (unused)
                                     profile_seed        = NULL,
                                     np                  = 0.10,
                                     by_family           = TRUE,
                                     family_stat         = c("median", "mean", "sum"),
                                     time_normalize      = c("none", "per_node"),
                                     only_family_plots   = FALSE,
                                     ...) {
  stopifnot(!missing(model))

  `%||%` <- function(a, b) if (is.null(a)) b else a
  .opt_logical <- function(opt, default) isTRUE(getOption(opt, default = default))
  .opt_integer <- function(opt, default) {
    val <- suppressWarnings(as.integer(getOption(opt, default = default)))
    if (is.na(val)) default else val
  }
  .to_family <- function(x) sub("\\[.*\\]", "", x)

  quiet <- FALSE
  .timed <- function(label, expr) {
    if (!quiet) cat(sprintf("\n[TIME] %s ...\n", label))
    tm <- system.time(res <- eval.parent(substitute(expr)))
    if (!quiet) {
      cat(sprintf(
        "[TIME] %s: user=%.3fs, sys=%.3fs, elapsed=%.3fs\n",
        label, tm["user.self"], tm["sys.self"], tm["elapsed"]
      ))
    }
    res
  }

  .diag_fast            <- .opt_logical("sop.diag.fast", FALSE)
  .profile_max_samplers <- .opt_integer("sop.diag.profile.max_samplers", Inf)
  .max_examples         <- .opt_integer("sop.diag.max_examples", 5L)

  family_stat    <- match.arg(family_stat)
  time_normalize <- match.arg(time_normalize)
  if (!isTRUE(by_family)) only_family_plots <- FALSE
  if (length(np) != 1 || !is.finite(np) || np <= 0 || np > 1) {
    warning("Invalid 'np'; using 0.10")
    np <- 0.10
  }

  .ignore_re <- if (length(ignore_patterns)) paste(ignore_patterns, collapse = "|") else NULL
  .is_ignored <- function(x) {
    if (is.null(.ignore_re)) rep(FALSE, length(x)) else grepl(.ignore_re, x, perl = TRUE)
  }

  ## 1) Nodes
  all_nodes_raw <- .timed("Get all node names", model$getNodeNames(includeData = include_data))
  removed_nodes <- unique(c(removed_nodes %||% character(0)))
  base_nodes    <- .timed("Filter ignore patterns and removed_nodes", {
    z <- all_nodes_raw[!.is_ignored(all_nodes_raw)]
    setdiff(z, removed_nodes)
  })

  ## 2) Stochastic / deterministic
  stoch_all <- .timed("Get stochastic node names",
                      model$getNodeNames(stochOnly = TRUE, includeData = include_data))
  stochastic_nodes    <- intersect(stoch_all, base_nodes)
  deterministic_nodes <- setdiff(base_nodes, stochastic_nodes)
  if (!quiet) {
    cat(sprintf(
      "\n[STRUCTURE]\n- # stochastic nodes   : %d\n- # deterministic nodes: %d\n",
      length(stochastic_nodes), length(deterministic_nodes)
    ))
  }

  ## 3) Dimensions
  dims <- .timed("getVarInfo() over base variables", {
    base_vars <- unique(gsub("\\[.*\\]", "", base_nodes))
    dv <- vector("list", length(base_vars))
    names(dv) <- base_vars
    for (i in seq_along(base_vars)) {
      info <- try(model$getVarInfo(base_vars[i]), silent = TRUE)
      dv[[i]] <- if (inherits(info, "try-error") || is.null(info)) {
        NA_integer_
      } else {
        info$nDim
      }
    }
    dv
  })

  ## 4) Dependencies
  deps_env <- new.env(parent = emptyenv())
  get_deps <- function(nd) {
    if (exists(nd, envir = deps_env, inherits = FALSE)) {
      return(get(nd, envir = deps_env))
    }
    d <- try(
      model$getDependencies(nodes = nd, self = FALSE, downstream = TRUE),
      silent = TRUE
    )
    if (inherits(d, "try-error") || is.null(d) || !length(d)) {
      res <- character(0)
    } else {
      d   <- d[!.is_ignored(d)]
      res <- intersect(d, base_nodes)
    }
    assign(nd, res, envir = deps_env)
    res
  }

  dep_counts <- dependencies_df <- NULL
  deps_pack <- .timed("Compute dependencies for all base nodes (memoized)", {
    node_dependencies  <- lapply(base_nodes, get_deps)
    names(node_dependencies) <- base_nodes
    total_dependencies <- vapply(node_dependencies, length, integer(1))

    dep_counts <- data.frame(
      parameter          = base_nodes,
      total_dependencies = as.numeric(total_dependencies),
      stringsAsFactors   = FALSE
    )

    any_non_empty <- total_dependencies > 0L
    dependencies_df <- if (any(any_non_empty)) {
      nodes_rep <- rep(base_nodes[any_non_empty], total_dependencies[any_non_empty])
      deps_cat  <- unlist(node_dependencies[any_non_empty], use.names = FALSE)
      data.frame(node = nodes_rep, dependency = deps_cat, stringsAsFactors = FALSE)
    } else {
      data.frame(node = character(0), dependency = character(0), stringsAsFactors = FALSE)
    }

    list(dep_counts = dep_counts, dependencies_df = dependencies_df)
  })
  dep_counts      <- deps_pack$dep_counts
  dependencies_df <- deps_pack$dependencies_df

  ## 5) Raw CSV
  if (isTRUE(save_csv)) {
    .timed("Write CSV dependencies_per_node", {
      if (is.null(output_dir)) {
        warning("save_csv=TRUE but output_dir is NULL: CSV not written.")
      } else {
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }
        utils::write.csv(
          dependencies_df,
          file      = file.path(output_dir, "dependencies_per_node.csv"),
          row.names = FALSE
        )
      }
    })
  }

  ## 6) configureMCMC + samplers
  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("Package 'nimble' is required for nimble::configureMCMC().")
  }

  mcmc_conf <- .timed("nimble::configureMCMC(model, ...)", {
    conf <- try(nimble::configureMCMC(model, ...), silent = TRUE)
    if (inherits(conf, "try-error")) {
      stop("Failed to run nimble::configureMCMC(model, ...).")
    }
    conf
  })

  samplers_df <- .timed("Extract samplers and targets", {
    smp_list <- mcmc_conf$getSamplers()
    if (length(smp_list) == 0L) {
      return(data.frame(
        Type        = character(0),
        TargetNodes = I(list()),
        stringsAsFactors = FALSE
      ))
    }
    Type        <- character(length(smp_list))
    TargetNodes <- vector("list", length(smp_list))
    for (i in seq_along(smp_list)) {
      s   <- smp_list[[i]]
      tgt <- s$target %||% character(0)
      if (length(tgt)) {
        tgt <- tgt[!.is_ignored(tgt)]
        tgt <- intersect(tgt, base_nodes)
        tgt <- intersect(tgt, stochastic_nodes)
        if (length(tgt) == 0L) tgt <- NA_character_
      } else {
        tgt <- NA_character_
      }
      Type[i]          <- paste(class(s), collapse = "/")
      TargetNodes[[i]] <- tgt
    }
    data.frame(Type = Type, TargetNodes = I(TargetNodes), stringsAsFactors = FALSE)
  })

  n_samplers <- nrow(samplers_df)
  if (length(np) != 1 || !is.finite(np) || np <= 0 || np > 1) {
    warning("Invalid 'np'; using 0.10")
    np <- 0.10
  }
  prop_worst <- if (n_samplers > 0L) {
    max(1L, min(n_samplers, ceiling(n_samplers * np)))
  } else {
    0L
  }

  has_sampler_nodes <- unique(unlist(samplers_df$TargetNodes))
  has_sampler_nodes <- has_sampler_nodes[!is.na(has_sampler_nodes)]

  ## 7) Sampler profiling
  per_param_times <- data.frame(parameter = character(0), sampler_time = numeric(0))

  if (isTRUE(auto_profile) && is.null(sampler_times)) {
    .timed("Auto-profile samplers (internal MCMC run)", {
      if (.diag_fast && is.finite(.profile_max_samplers) &&
          n_samplers > .profile_max_samplers) {
        warning("Too many samplers for fast diag; skipping profiling.")
      } else if (n_samplers == 0L) {
        warning("No samplers to profile; sampler_times remains NULL.")
      } else {
        st <- try({
          mcmc_obj <- nimble::buildMCMC(mcmc_conf)

          cmodel_local <- if ("RmodelBaseClass" %in% class(model)) {
            nimble::compileNimble(model)
          } else {
            model
          }

          cm <- nimble::compileNimble(
            mcmc_obj,
            project        = cmodel_local,
            resetFunctions = TRUE
          )

          old_opt <- nimble::nimbleOptions(MCMCprogressBar = FALSE)
          on.exit(nimble::nimbleOptions(MCMCprogressBar = old_opt), add = TRUE)

          if ("resetTimes" %in% ls(cm)) {
            try(cm$resetTimes(), silent = TRUE)
          }
          if (!is.null(profile_seed)) set.seed(profile_seed)

          cm$run(as.integer(profile_niter), time = TRUE)
          tm <- cm$getTimes()
          tm
        }, silent = TRUE)

        if (inherits(st, "try-error") || is.null(st)) {
          warning("Internal profiling failed; sampler_times remains NULL.")
        } else {
          st <- as.numeric(st)
          if (length(st) == 0L) {
            warning("Internal profiling returned length 0; sampler_times remains NULL.")
          } else {
            sampler_times <- st
            if (!quiet) {
              cat(sprintf("[PROFILE] Raw sampler_times length (internal) = %d.\n",
                          length(sampler_times)))
            }
          }
        }
      }
    })
  } else if (!is.null(sampler_times) && !quiet) {
    cat(sprintf("[PROFILE] Using user-supplied sampler_times (length = %d).\n",
                length(sampler_times)))
  }

  ## Harmonize sampler_times length
  if (!is.null(sampler_times)) {
    if (!is.numeric(sampler_times)) {
      warning("sampler_times is not numeric; ignoring sampler times.")
      sampler_times <- NULL
    } else if (length(sampler_times) != n_samplers) {
      if (length(sampler_times) == 0L) {
        warning(sprintf(
          "sampler_times has length 0 but there are %d samplers; using NA for all samplers.",
          n_samplers
        ))
        sampler_times <- rep(NA_real_, n_samplers)
      } else if (length(sampler_times) > n_samplers) {
        warning(sprintf(
          "sampler_times has length %d > n_samplers = %d; truncating.",
          length(sampler_times), n_samplers
        ))
        sampler_times <- sampler_times[seq_len(n_samplers)]
      } else {
        warning(sprintf(
          "sampler_times has length %d < n_samplers = %d; padding with NA.",
          length(sampler_times), n_samplers
        ))
        sampler_times <- c(
          sampler_times,
          rep(NA_real_, n_samplers - length(sampler_times))
        )
      }
    }
  }

  if (!is.null(sampler_times) && n_samplers > 0L && length(has_sampler_nodes) > 0L) {
    .timed("Aggregate sampler times per parameter", {
      idx_ok <- which(!vapply(
        samplers_df$TargetNodes,
        function(x) length(x) == 1 && is.na(x),
        logical(1)
      ))
      if (length(idx_ok)) {
        param_vec <- unlist(samplers_df$TargetNodes[idx_ok], use.names = FALSE)
        time_vec  <- rep(
          sampler_times[idx_ok],
          times = vapply(samplers_df$TargetNodes[idx_ok], length, integer(1))
        )
        keep <- param_vec %in% has_sampler_nodes
        if (any(keep)) {
          long_tbl <- data.frame(
            parameter    = param_vec[keep],
            sampler_time = time_vec[keep],
            stringsAsFactors = FALSE
          )
          per_param_times <- stats::aggregate(
            sampler_time ~ parameter,
            data = long_tbl,
            FUN  = sum
          )
        }
      }
    })
  }

  if (n_samplers > 0L && isTRUE(make_plots) && is.null(sampler_times)) {
    message("[diagnose_model_structure] 'sampler_times' is NULL -> sampler-time plots will be empty.")
  }

  ## 8) Tidy tables
  deps_df <- .timed("Assemble deps_df (node order by total_dependencies)", {
    if (is.null(dep_counts) || !nrow(dep_counts)) {
      data.frame(parameter = factor(character(0)), total_dependencies = numeric(0))
    } else {
      ord_nodes <- dep_counts$parameter[order(-dep_counts$total_dependencies)]
      dep_map   <- stats::setNames(dep_counts$total_dependencies, dep_counts$parameter)
      data.frame(
        parameter          = factor(ord_nodes, levels = ord_nodes),
        total_dependencies = as.numeric(dep_map[ord_nodes]),
        stringsAsFactors   = FALSE
      )
    }
  })

  sampler_df <- .timed("Assemble sampler_df (restricted to nodes with sampler)", {
    st_named    <- stats::setNames(per_param_times$sampler_time, per_param_times$parameter)
    time_nodes  <- intersect(as.character(deps_df$parameter), has_sampler_nodes)
    sampler_vec <- unname(st_named[time_nodes])
    ok          <- !is.na(sampler_vec)
    time_nodes  <- time_nodes[ok]
    sampler_vec <- sampler_vec[ok]
    if (length(time_nodes)) {
      data.frame(
        parameter    = factor(time_nodes, levels = time_nodes),
        sampler_time = as.numeric(sampler_vec),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(parameter = factor(character(0)), sampler_time = numeric(0))
    }
  })

  ## 9) node_of_interest
  node_of_interest_deps <- .timed("Compute node_of_interest downstream", {
    if (is.null(node_of_interest)) return(NULL)
    if (!(node_of_interest %in% base_nodes)) {
      warning(sprintf(
        "node_of_interest = '%s' not present after filtering; skipping.",
        node_of_interest
      ))
      return(NULL)
    }
    deps <- try(
      model$getDependencies(
        nodes = node_of_interest,
        self  = FALSE,
        downstream = TRUE
      ),
      silent = TRUE
    )
    if (inherits(deps, "try-error") || is.null(deps) || !length(deps)) return(character(0))
    deps <- deps[!.is_ignored(deps)]
    intersect(deps, base_nodes)
  })

  ## 10) Aggregation by family
  fam_deps_df <- .timed("Aggregate dependencies by family", {
    if (!nrow(deps_df)) {
      data.frame(family = character(0), deps_stat = numeric(0))
    } else {
      tmp <- deps_df
      tmp$fam <- .to_family(as.character(tmp$parameter))
      out <- stats::aggregate(
        total_dependencies ~ fam,
        data = tmp,
        FUN  = switch(
          family_stat,
          median = stats::median,
          mean   = base::mean,
          sum    = base::sum
        )
      )
      names(out) <- c("family", "deps_stat")
      out
    }
  })

  fam_time_df <- .timed("Aggregate sampler time by family", {
    if (!nrow(sampler_df)) {
      data.frame(family = character(0), time_stat = numeric(0))
    } else {
      tmp <- sampler_df
      tmp$fam <- .to_family(as.character(tmp$parameter))
      out <- stats::aggregate(
        sampler_time ~ fam,
        data = tmp,
        FUN  = switch(
          family_stat,
          median = stats::median,
          mean   = base::mean,
          sum    = base::sum
        )
      )
      names(out) <- c("family", "time_stat")
      if (time_normalize == "per_node") {
        n_by_fam <- stats::aggregate(
          parameter ~ fam,
          data = tmp,
          FUN  = function(z) length(unique(z))
        )
        names(n_by_fam) <- c("family", "n_nodes")
        out <- merge(out, n_by_fam, by = "family", all.x = TRUE)
        out$time_stat <- out$time_stat / pmax(out$n_nodes, 1L)
        out$n_nodes   <- NULL
      }
      out
    }
  })

  if (nrow(fam_deps_df) > 0L) {
    ord_fam_all <- fam_deps_df$family[order(-fam_deps_df$deps_stat)]
  } else {
    ord_fam_all <- character(0)
  }
  if (nrow(fam_deps_df) > 0L) {
    fam_deps_df$family <- factor(fam_deps_df$family, levels = ord_fam_all)
  }

  if (nrow(fam_time_df) > 0L) {
    ord_fam_time <- if (length(ord_fam_all)) {
      intersect(ord_fam_all, fam_time_df$family)
    } else {
      fam_time_df$family
    }
    fam_time_df$family <- factor(fam_time_df$family, levels = ord_fam_time)
  }

  ## 11) Plots
  plot_dependencies <- plot_sampler_time <- plot_combined <- NULL
  plot_dependencies_family <- plot_sampler_time_family <- plot_combined_family <- NULL

  if (isTRUE(make_plots)) {
    .timed("Make plots", {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("ggplot2 not installed: plots will not be produced.")
      } else {
        base_theme <- ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(
            plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
            axis.title.x     = ggplot2::element_text(face = "bold"),
            axis.title.y     = ggplot2::element_text(face = "bold"),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position  = "top",
            legend.title     = ggplot2::element_text(face = "bold"),
            axis.text.x      = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
          )

        if (!isTRUE(only_family_plots)) {
          if (nrow(deps_df) > 0L) {
            plot_dependencies <- ggplot2::ggplot(
              deps_df,
              ggplot2::aes(x = parameter, y = total_dependencies, fill = "Total dependencies")
            ) +
              ggplot2::geom_col(width = 1) +
              ggplot2::scale_x_discrete(expand = c(0, 0)) +
              ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
              ggplot2::scale_fill_manual(
                name   = "",
                values = c("Total dependencies" = "orange")
              ) +
              ggplot2::labs(
                title = "Number of dependencies per node",
                x     = "Node",
                y     = "Number of dependencies"
              ) +
              base_theme
          }
          if (nrow(sampler_df) > 0L) {
            time_label <- sprintf("Time in samplers (%s)", sampler_times_unit)
            plot_sampler_time <- ggplot2::ggplot(
              sampler_df,
              ggplot2::aes(x = parameter, y = sampler_time, fill = time_label)
            ) +
              ggplot2::geom_col(width = 1) +
              ggplot2::scale_x_discrete(expand = c(0, 0)) +
              ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
              ggplot2::scale_fill_manual(
                name   = "",
                values = stats::setNames("blue", time_label)
              ) +
              ggplot2::labs(
                title = "Time spent in samplers per node",
                x     = "Node",
                y     = time_label
              ) +
              base_theme
          }
          if (!is.null(plot_dependencies) && !is.null(plot_sampler_time)) {
            if (requireNamespace("patchwork", quietly = TRUE)) {
              plot_combined <- plot_dependencies +
                plot_sampler_time +
                patchwork::plot_layout(ncol = 2)
            } else if (requireNamespace("gridExtra", quietly = TRUE)) {
              plot_combined <- gridExtra::arrangeGrob(
                plot_dependencies, plot_sampler_time, ncol = 2
              )
            }
          }
        }

        if (isTRUE(by_family)) {
          if (nrow(fam_deps_df) > 0L) {
            lab_dep <- sprintf("Number of dependencies (%s by family)", family_stat)
            plot_dependencies_family <- ggplot2::ggplot(
              fam_deps_df,
              ggplot2::aes(x = family, y = deps_stat, fill = "Deps by family")
            ) +
              ggplot2::geom_col(width = 1) +
              ggplot2::scale_x_discrete(expand = c(0, 0)) +
              ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
              ggplot2::scale_fill_manual(
                name   = "",
                values = c("Deps by family" = "orange")
              ) +
              ggplot2::labs(
                title = lab_dep,
                x     = "Family",
                y     = family_stat
              ) +
              base_theme +
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
              )
          }
          if (nrow(fam_time_df) > 0L) {
            lab_tim <- sprintf(
              "Time in samplers (%s by family, %s)",
              family_stat,
              if (time_normalize == "per_node") "per node" else "raw"
            )
            plot_sampler_time_family <- ggplot2::ggplot(
              fam_time_df,
              ggplot2::aes(x = family, y = time_stat, fill = "Time by family")
            ) +
              ggplot2::geom_col(width = 1) +
              ggplot2::scale_x_discrete(expand = c(0, 0)) +
              ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
              ggplot2::scale_fill_manual(
                name   = "",
                values = c("Time by family" = "blue")
              ) +
              ggplot2::labs(
                title = lab_tim,
                x     = "Family",
                y     = sprintf("Time (%s)", sampler_times_unit)
              ) +
              base_theme +
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
              )
          }
          if (!is.null(plot_dependencies_family) && !is.null(plot_sampler_time_family)) {
            if (requireNamespace("patchwork", quietly = TRUE)) {
              plot_combined_family <- plot_dependencies_family +
                plot_sampler_time_family +
                patchwork::plot_layout(ncol = 2)
            } else if (requireNamespace("gridExtra", quietly = TRUE)) {
              plot_combined_family <- gridExtra::arrangeGrob(
                plot_dependencies_family,
                plot_sampler_time_family,
                ncol = 2
              )
            }
          }
        }

        ## Save to disk
        if (!is.null(output_dir)) {
          if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
          }

          if (isTRUE(save_csv)) {
            if (nrow(fam_deps_df) > 0L) {
              utils::write.csv(
                fam_deps_df,
                file      = file.path(output_dir,
                                      sprintf("deps_by_family_%s.csv", family_stat)),
                row.names = FALSE
              )
            }
            if (nrow(fam_time_df) > 0L) {
              utils::write.csv(
                fam_time_df,
                file = file.path(
                  output_dir,
                  sprintf("sampler_time_by_family_%s_%s.csv",
                          family_stat, time_normalize)
                ),
                row.names = FALSE
              )
            }
          }

          if (!isTRUE(only_family_plots)) {
            if (!is.null(plot_dependencies)) {
              ggplot2::ggsave(
                file.path(output_dir, "dependencies_per_parameter.png"),
                plot   = plot_dependencies,
                width  = 12, height = 6, dpi = 300
              )
            }
            if (!is.null(plot_sampler_time)) {
              ggplot2::ggsave(
                file.path(output_dir, "sampler_time_per_parameter.png"),
                plot   = plot_sampler_time,
                width  = 12, height = 6, dpi = 300
              )
            }
            if (!is.null(plot_combined)) {
              if (inherits(plot_combined, "ggplot")) {
                ggplot2::ggsave(
                  file.path(output_dir, "deps_and_sampler_time_side_by_side.png"),
                  plot   = plot_combined,
                  width  = 14, height = 6, dpi = 300
                )
              } else {
                grDevices::png(
                  file.path(output_dir, "deps_and_sampler_time_side_by_side.png"),
                  width = 14, height = 6, units = "in", res = 300
                )
                grid::grid.newpage()
                grid::grid.draw(plot_combined)
                grDevices::dev.off()
              }
            }
          }

          if (nrow(fam_deps_df) > 0L && !is.null(plot_dependencies_family)) {
            ggplot2::ggsave(
              file.path(output_dir, sprintf("deps_by_family_%s.png", family_stat)),
              plot   = plot_dependencies_family,
              width  = 12, height = 6, dpi = 300
            )
          }
          if (nrow(fam_time_df) > 0L && !is.null(plot_sampler_time_family)) {
            ggplot2::ggsave(
              file.path(
                output_dir,
                sprintf("sampler_time_by_family_%s_%s.png",
                        family_stat, time_normalize)
              ),
              plot   = plot_sampler_time_family,
              width  = 12, height = 6, dpi = 300
            )
          }
          if (!is.null(plot_combined_family)) {
            if (inherits(plot_combined_family, "ggplot")) {
              ggplot2::ggsave(
                file.path(
                  output_dir,
                  sprintf("family_deps_and_time_side_by_side_%s_%s.png",
                          family_stat, time_normalize)
                ),
                plot   = plot_combined_family,
                width  = 14, height = 6, dpi = 300
              )
            } else {
              grDevices::png(
                file.path(
                  output_dir,
                  sprintf("family_deps_and_time_side_by_side_%s_%s.png",
                          family_stat, time_normalize)
                ),
                width = 14, height = 6, units = "in", res = 300
              )
              grid::grid.newpage()
              grid::grid.draw(plot_combined_family)
              grDevices::dev.off()
            }
          }
        } else {
          message("[diagnose_model_structure] 'output_dir' is NULL -> plots will not be saved to disk.")
        }

        ## Print (without duplication)
        if (isTRUE(make_plots)) {
          if (!isTRUE(only_family_plots)) {
            if (!is.null(plot_combined)) {
              if (inherits(plot_combined, "ggplot")) {
                print(plot_combined)
              } else {
                grid::grid.newpage()
                grid::grid.draw(plot_combined)
              }
            } else {
              if (!is.null(plot_dependencies)) print(plot_dependencies)
              if (!is.null(plot_sampler_time)) print(plot_sampler_time)
            }
          }

          if (isTRUE(by_family)) {
            if (!is.null(plot_combined_family)) {
              if (inherits(plot_combined_family, "ggplot")) {
                print(plot_combined_family)
              } else {
                grid::grid.newpage()
                grid::grid.draw(plot_combined_family)
              }
            } else if (!is.null(plot_dependencies_family)) {
              print(plot_dependencies_family)
            } else if (!is.null(plot_sampler_time_family)) {
              print(plot_sampler_time_family)
            }
          }
        }
      }
    })
  }

  ## 12) Return
  .timed("Assemble return object", {
    invisible(list(
      stochastic_nodes      = stochastic_nodes,
      deterministic_nodes   = deterministic_nodes,
      dims                  = dims,
      nodes                 = NULL,
      dependencies_df       = dependencies_df,
      dep_counts            = dep_counts,
      node_of_interest_deps = node_of_interest_deps,
      mcmc_conf             = mcmc_conf,
      samplers_df           = samplers_df,
      n_samplers            = n_samplers,
      prop_worst            = prop_worst,
      per_param_times       = per_param_times,
      deps_df               = deps_df,
      sampler_df            = sampler_df,
      fam_deps_df           = fam_deps_df,
      fam_time_df           = fam_time_df,
      plots = list(
        plot_dependencies        = plot_dependencies,
        plot_sampler_time        = plot_sampler_time,
        plot_combined            = plot_combined,
        plot_dependencies_family = plot_dependencies_family,
        plot_sampler_time_family = plot_sampler_time_family,
        plot_combined_family     = plot_combined_family
      )
    ))
  })
}

#' Target-level diagnostics (time + optional step-size proxy)
#'
#' Builds a robust MCMC configuration, profiles per-sampler time with
#' \code{\link{profile_sampler_times}}, maps samplers to their targets, and
#' optionally computes a step-size proxy per target via \code{\link{proxy_step_sd}}
#' if \code{samples} are provided.
#'
#' Note: uses a tolerant default monitor set followed by a robust configuration
#' (no \code{monitors = } passed directly to \code{configureMCMC}).
#'
#' @param build_fn A model builder function returning at least a list with
#'   components \code{model} (nimbleModel) and \code{cmodel} (compiled model).
#'   It may also return additional elements (e.g. \code{monitors}) which are
#'   ignored here.
#' @param opts     A list of options (for example from \code{samOptiPro_options()}).
#' @param niter    Integer; iterations for time profiling (defaults to
#'   \code{opts$time_profile_n}).
#' @param samples  Optional \code{coda::mcmc.list} used to compute the step-size
#'   proxy via \code{proxy_step_sd()}.
#'
#' @return A \code{data.frame} with columns \code{target}, \code{type},
#'   \code{time_s}, \code{step_sd}.
#' @export
#' @keywords internal
diagnostics_by_target <- function(build_fn,
                                  opts    = samOptiPro_options(),
                                  niter   = opts$time_profile_n,
                                  samples = NULL) {

  stopifnot(is.function(build_fn))

  # Build model and compiled model (builder is responsible for both)
  built <- build_fn()

  # Tolerant monitors -> then robust configuration
  mons <- default_monitors(built$model, opts)
  mons <- ensure_monitors_exist(built$model, mons)
  conf <- .configure_with_monitors(built$model, monitors = mons)

  # Profile per-sampler time
  time_s <- profile_sampler_times(conf, built$cmodel, niter = niter)

  # Map samplers to targets
  df_map <- sampler_df_from_conf(conf)  # columns: name, type, target, scale (fallback)
  df_map$time_s <- time_s

  # Optional step proxy from MCMC samples
  if (!is.null(samples)) {
    pats <- unique(gsub("\\[.*$", "[", df_map$target))
    step_df <- do.call(rbind, lapply(pats, function(p) {
      vals <- proxy_step_sd(samples, p)
      if (!length(vals)) return(NULL)
      data.frame(
        target  = names(vals),
        step_sd = as.numeric(vals),
        stringsAsFactors = FALSE
      )
    }))
    if (!is.null(step_df) && nrow(step_df) > 0L) {
      df_map <- merge(df_map, step_df, by = "target", all.x = TRUE, sort = FALSE)
    } else {
      df_map$step_sd <- NA_real_
    }
  } else {
    df_map$step_sd <- NA_real_
  }

  df_map[, c("target", "type", "time_s", "step_sd")]
}
#' Profile sampler times per MCMC sampler
#'
#' Compile (once) and run a NIMBLE MCMC configured in \code{conf}, with
#' per-sampler timing enabled, then return the vector of times (seconds) for
#' each sampler defined in \code{conf}. A small in-session cache can be added
#' by the caller if desired, but this function itself performs no caching.
#'
#' @param conf A NIMBLE MCMC configuration object (for example returned by
#'   \code{nimble::configureMCMC(model)}).
#' @param cmodel A compiled \code{nimbleModel} used as \code{project = } for
#'   \code{nimble::compileNimble()}.
#' @param niter Integer scalar; number of MCMC iterations to run for profiling.
#'   Default \code{5e4}.
#'
#' @return A numeric vector of length equal to the number of samplers in
#'   \code{conf}, giving the per-sampler elapsed time (in seconds). If NIMBLE
#'   reports fewer or more entries than samplers, the vector is padded or
#'   truncated accordingly with \code{NA}s.
#'
#' @examples
#' \dontrun{
#'   model  <- nimbleModel(code, constants = Const, data = Data, inits = Inits)
#'   cmodel <- compileNimble(model)
#'   conf   <- configureMCMC(model)
#'   ts <- profile_sampler_times(conf, cmodel, niter = 1e4)
#'   print(ts)
#' }
#'
#' @importFrom nimble buildMCMC compileNimble nimbleOptions
#' @export
profile_sampler_times <- function(conf, cmodel, niter = 5e4) {
  # Samplers
  samplers <- conf$getSamplers()
  ns <- length(samplers)

  if (ns == 0L) {
    warning("profile_sampler_times(): no samplers; returning numeric(0).")
    return(numeric(0))
  }

  # Build + compile MCMC
  m  <- nimble::buildMCMC(conf)
  cm <- nimble::compileNimble(m, project = cmodel, resetFunctions = TRUE)

  old_opt <- nimble::nimbleOptions(MCMCprogressBar = FALSE)
  on.exit(nimble::nimbleOptions(MCMCprogressBar = old_opt), add = TRUE)

  if ("resetTimes" %in% ls(cm)) {
    try(cm$resetTimes(), silent = TRUE)
  }

  cm$run(as.integer(niter), time = TRUE)

  tm <- cm$getTimes()

  if (is.null(tm) || length(tm) == 0L) {
    warning("profile_sampler_times(): cm$getTimes() returned NULL or length 0; using NA.")
    return(rep(NA_real_, ns))
  }

  tm <- as.numeric(tm)
  if (length(tm) > ns) tm <- tm[seq_len(ns)]
  if (length(tm) < ns) tm <- c(tm, rep(NA_real_, ns - length(tm)))

  tm
}
#' Step-size proxy: sd(diff(chain)) on columns matching a pattern
#'
#' Computes, for each parameter whose name matches \code{pattern}, the standard
#' deviation of successive differences \code{sd(diff(.))}, as a simple
#' step-size proxy.
#'
#' @param samples A \code{coda::mcmc.list}.
#' @param pattern A string; pattern to match against column names (fixed match).
#'
#' @return A named numeric vector giving \code{sd(diff(.))} per parameter.
#' @export
#' @keywords internal
proxy_step_sd <- function(samples, pattern) {
  stopifnot(inherits(samples, "mcmc.list"))
  mx <- as.matrix(samples[[1]])
  cols <- grep(pattern, colnames(mx), fixed = TRUE)
  if (!length(cols)) return(numeric(0))
  apply(mx[, cols, drop = FALSE], 2, function(v) stats::sd(diff(v), na.rm = TRUE))
}
#' Run structural diagnostics and (optional) HMC/NUTS smoke test
#'
#' Inspects a NIMBLE model produced by \code{build_fn()} to:
#' \itemize{
#'   \item count stochastic vs. deterministic nodes (optionally including data),
#'   \item parse the model code to detect non-differentiable functions and BUGS-style truncation,
#'   \item rebuild a per-node table with distribution, support, and bounds,
#'   \item tag \strong{HMC showstoppers} (for example discrete latents, simplex constraints,
#'         truncation, non-differentiable deterministic operations feeding latents),
#'   \item optionally attempt a short HMC/NUTS run (via \pkg{nimbleHMC}) to check practical feasibility.
#' }
#'
#' Compared to \code{diagnose_model_structure()}, this function does not filter nodes by
#' ignore patterns or explicit removals and does not compute dependency fan-out.
#' Its emphasis is \emph{HMC suitability} rather than structural load.
#'
#' @param build_fn A zero-arg function returning a list with at least
#'   \code{model} (\code{nimbleModel}); optional elements: \code{cmodel},
#'   \code{monitors}, \code{code_text} (or \code{code}) used to strengthen
#'   non-differentiability detection.
#' @param include_data Logical; include data nodes in counts and scans.
#'   Default \code{FALSE}.
#' @param try_hmc Logical; if \code{TRUE}, attempt a brief HMC/NUTS run.
#'   Default \code{TRUE}.
#' @param niter Integer; total iterations for the HMC test. Default \code{50L}.
#' @param nburnin Integer; burn-in for the HMC test. Default \code{10L}.
#' @param seed Integer; RNG seed for the HMC run. Default \code{1L}.
#'
#' @return An invisible list with components:
#' \itemize{
#'   \item \code{diag}: list with \code{nodes} (data.frame of node metadata and HMC showstopper tags),
#'         \code{nondiff_signals}, \code{code_scan}, \code{hmc_globally_ok}.
#'   \item \code{hmc}: list describing the HMC trial (\code{ok}, \code{error}, \code{details}).
#' }
#'
#' @examples
#' \dontrun{
#' out <- run_structure_and_hmc_test(my_builder, include_data = FALSE, try_hmc = TRUE)
#' if (!out$hmc$ok) message("HMC not feasible: ", out$hmc$error)
#' }
#'
#' @export
run_structure_and_hmc_test <- function(build_fn,
                                       include_data = FALSE,
                                       try_hmc     = TRUE,
                                       niter       = 50L,
                                       nburnin     = 10L,
                                       seed        = 1L) {
  # ---- Compat: internal util
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # ---- Internal options (no API impact) ----
  fast_scan         <- isTRUE(getOption("sop.fast_scan", FALSE))
  max_print_nodes   <- as.integer(getOption("sop.max_print_nodes", 5000L))
  max_show_examples <- as.integer(getOption("sop.max_show_examples", 5L))

  # ---- Small utilities ----
  to_logical <- function(x) {
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    if (is.factor(x))  x <- as.character(x)
    if (is.character(x)) {
      y <- tolower(trimws(x))
      return(y %in% c("true","t","1","yes","y"))
    }
    rep(FALSE, length(x))
  }
  strip_index <- function(x) sub("\\[.*\\]$", "", x)

  normalize_dist <- function(d) {
    if (is.null(d) || length(d) == 0 || all(is.na(d))) return(NA_character_)
    d <- tolower(trimws(d))
    if (!nzchar(d)) return(NA_character_)
    d <- sub("^dbernoulli$", "dbern", d)
    d <- sub("^dbin$",       "dbinom", d)
    d <- sub("^dmulti$",     "dmultinom", d)
    d <- sub("^dnbinom$",    "dnegbin", d)
    d <- sub("^ddirich(let)?$", "ddirichlet", d)
    d <- sub("^ddirch(let)?$",  "ddirichlet", d)
    d <- sub("^dhyperg(eometric)?$", "dhypergeom", d)
    d <- sub("^dhyper$",             "dhypergeom", d)
    d <- sub("^dpois_rate$", "dpois", d)
    d <- sub("^dtri(angle)?$",   "dtriangle",  d)
    d <- sub("^dhalfnorm(al)?$", "dhalfnorm",  d)
    d <- sub("^dhalfcauchy$",    "dhalfcauchy", d)
    d
  }

  discrete_dists <- c(
    "dbern","dbinom","dcat","dmultinom",
    "dgeom","dnegbin","dpois","dhypergeom","dinterval"
  )
  simplex_dists  <- c("ddirichlet")

  # log-normal has natural support (0, +Inf) but is not a user truncation
  bounded_support_dists <- c(
    "dbeta","dunif","dtriangle","dhalfnorm","dhalfcauchy","dlnorm"
  )

  support_of <- function(dist_name, lb, ub) {
    d <- normalize_dist(dist_name)
    if (is.na(d)) return("unknown")
    if (d %in% discrete_dists)        return("discrete")
    if (d %in% simplex_dists)         return("simplex")
    if (d %in% bounded_support_dists) return("bounded-continuous")
    "continuous"
  }

  .sop_supports_derivs <- function(model) {
    bd1 <- try(model$modelDef$buildDerivs, silent = TRUE)
    if (!inherits(bd1, "try-error") && !is.null(bd1)) return(isTRUE(bd1))
    bd2 <- try(model$buildDerivs, silent = TRUE)
    if (!inherits(bd2, "try-error") && !is.null(bd2)) return(isTRUE(bd2))
    if (requireNamespace("nimbleHMC", quietly = TRUE)) {
      conf <- try(nimble::configureMCMC(model), silent = TRUE)
      if (!inherits(conf, "try-error")) {
        ok <- TRUE
        tryCatch({
          nimbleHMC::configureHMC(conf, model = model)
        }, error = function(e) {
          msg <- tolower(conditionMessage(e))
          if (grepl("deriv", msg) || grepl("buildderiv", msg)) ok <<- FALSE
        })
        return(ok)
      }
    }
    FALSE
  }

  # ==== Build ====
  stopifnot(is.function(build_fn))
  parts <- build_fn()
  m     <- parts$model
  cm    <- parts$cmodel %||% NULL
  mons  <- parts$monitors %||% character(0)

  # Prefer explicit code text for scanning (if provided by the builder)
  code_override <- parts$code_text %||% parts$code %||% NULL
  if (!is.null(code_override) && !is.character(code_override)) {
    code_override <- try(paste(deparse(code_override), collapse = "\n"), silent = TRUE)
    if (inherits(code_override, "try-error")) code_override <- NULL
  }

  # ---- Structure counts (once) ----
  n_stoch <- length(m$getNodeNames(stochOnly = TRUE, includeData = include_data))
  n_all   <- length(m$getNodeNames(includeData = include_data))
  n_det   <- n_all - n_stoch

  cat("\n[STRUCTURE]\n")
  cat(sprintf("- # stochastic nodes   : %d\n", n_stoch))
  cat(sprintf("- # deterministic nodes: %d\n", n_det))

  # ---- Non-diff function detection ----
  .sop_detect_nondiff_functions <- function(model,
                                            nondiff_candidates = c(
                                              "round","floor","ceiling","trunc",
                                              "abs","max","min","step","ifElse","ifelse","equals"
                                            ),
                                            code_text_override = NULL,
                                            include_data = FALSE,
                                            fast_scan = FALSE,
                                            max_print_nodes = 5000L) {
    out <- character(0)
    code_txt <- ""

    if (is.character(code_text_override) && length(code_text_override)) {
      code_txt <- paste(code_text_override, collapse = "\n")
    } else {
      code_try <- try({
        code_txt <- paste(deparse(model$modelDef$code, width.cutoff = 500), collapse = "\n")
      }, silent = TRUE)

      if (inherits(code_try, "try-error") || is.null(code_txt) || !nzchar(code_txt)) {
        n_all <- length(try(model$getNodeNames(includeData = include_data), silent = TRUE))
        if (!fast_scan && is.finite(n_all) && n_all <= max_print_nodes) {
          code_txt <- tryCatch(
            paste(utils::capture.output(print(model)), collapse = "\n"),
            error = function(e) ""
          )
        } else {
          code_txt <- ""
        }
      }
    }

    if (nzchar(code_txt)) {
      pat <- paste0("\\b(", paste(nondiff_candidates, collapse="|"), ")\\s*\\(")
      hit <- gregexpr(pat, code_txt, perl = TRUE, ignore.case = TRUE)
      if (length(hit) && hit[[1]][1] != -1) {
        out <- nondiff_candidates
      }
    }

    list(
      hits     = intersect(tolower(out), tolower(nondiff_candidates)),
      code_txt = code_txt
    )
  }

  nd <- .sop_detect_nondiff_functions(
    m,
    code_text_override = code_override,
    include_data = include_data,
    fast_scan = fast_scan,
    max_print_nodes = max_print_nodes
  )
  nondiff_hits     <- nd$hits
  code_txt_scanned <- nd$code_txt

  # ---- BUGS-style truncation "T(a,b)" via text scan ----
  truncation_detected_bugst <- nzchar(code_txt_scanned) &&
    grepl("\\bT\\s*\\(", code_txt_scanned, perl = TRUE)

  # ---- Memoization for dependencies (for non-diff propagation) ----
  dep_cache <- new.env(parent = emptyenv())
  get_deps_memo <- function(target, upstream, downstream, includeData, stochOnly) {
    key <- paste(target, upstream, downstream, includeData, stochOnly, sep = "|")
    if (exists(key, envir = dep_cache, inherits = FALSE)) return(get(key, envir = dep_cache))
    val <- tryCatch(
      m$getDependencies(
        target = target,
        upstream = upstream,
        downstream = downstream,
        includeData = includeData,
        stochOnly = stochOnly
      ),
      error = function(e) character(0)
    )
    assign(key, val, envir = dep_cache)
    val
  }

  varinfo_cache <- new.env(parent = emptyenv())
  get_varinfo_dim <- function(vname) {
    if (exists(vname, envir = varinfo_cache, inherits = FALSE)) return(get(vname, envir = varinfo_cache))
    vi <- tryCatch(m$getVarInfo(vname), error = function(e) NULL)
    res <- if (is.null(vi) || is.null(vi$nDim)) NA_integer_ else as.integer(vi$nDim)
    assign(vname, res, envir = varinfo_cache)
    res
  }

  # ---- Rebuild node table ----
  rebuild_nodes_df <- function(model, include_data) {
    all_nodes <- model$getNodeNames(includeData = include_data)
    stoch     <- model$getNodeNames(stochOnly = TRUE, includeData = include_data)

    is_stoch <- all_nodes %in% stoch
    is_data  <- if (isTRUE(include_data)) {
      out <- logical(length(all_nodes))
      for (i in seq_along(all_nodes)) {
        out[i] <- tryCatch(model$isData(all_nodes[i]), error = function(e) FALSE)
      }
      out
    } else {
      rep(FALSE, length(all_nodes))
    }

    get_dist  <- function(n) tryCatch(tolower(model$getDistribution(n)), error = function(e) NA_character_)
    get_bound <- function(n, side) suppressWarnings(
      tryCatch(model$getBound(n, side), error = function(e) NA_real_)
    )

    dist_raw <- rep(NA_character_, length(all_nodes))
    lower    <- rep(NA_real_,      length(all_nodes))
    upper    <- rep(NA_real_,      length(all_nodes))

    idx_st <- which(is_stoch)
    if (length(idx_st)) {
      nn_st <- all_nodes[idx_st]
      dist_raw[idx_st] <- vapply(nn_st, get_dist,  character(1))
      lower[idx_st]    <- vapply(nn_st, get_bound, numeric(1), side="lower")
      upper[idx_st]    <- vapply(nn_st, get_bound, numeric(1), side="upper")
    }

    dist <- ifelse(is_stoch, vapply(dist_raw, normalize_dist, character(1)), NA_character_)
    var  <- sub("\\[.*\\]$", "", all_nodes)

    support <- ifelse(
      is_stoch,
      mapply(support_of, dist, lower, upper, SIMPLIFY = TRUE),
      NA_character_
    )

    is_trunc <- is_stoch & support == "continuous" &
      (is.finite(lower) | is.finite(upper))

    out <- data.frame(
      node         = all_nodes,
      var          = var,
      is_stoch     = is_stoch,
      is_data      = is_data,
      dist_raw     = dist_raw,
      dist         = dist,
      support      = support,
      lower        = lower,
      upper        = upper,
      is_truncated = is_trunc,
      stringsAsFactors = FALSE
    )

    uvars <- unique(out$var)
    dims_map <- vapply(uvars, get_varinfo_dim, integer(1))
    names(dims_map) <- uvars
    out$dims <- unname(dims_map[out$var])

    out
  }

  nodes <- rebuild_nodes_df(m, include_data)
  if ("node" %in% names(nodes)) {
    nodes$node <- as.character(nodes$node)
  }

  # Half-Cauchy detection via T(dt(0,1,1), 0, Inf)
  nodes$is_halfcauchy_like <- FALSE
  idx_dt <- which(nodes$is_stoch & !nodes$is_data & nodes$dist == "dt")
  if (length(idx_dt)) {
    for (k in idx_dt) {
      lb <- suppressWarnings(as.numeric(nodes$lower[k]))
      ub <- suppressWarnings(as.numeric(nodes$upper[k]))

      df_k <- suppressWarnings(tryCatch(
        m$getParam(nodes$node[k], "df"),
        error = function(e) NA_real_
      ))
      mu_k <- suppressWarnings(tryCatch(
        m$getParam(nodes$node[k], "mu"),
        error = function(e) NA_real_
      ))

      is_hc <- !is.na(df_k) && abs(df_k - 1) < 1e-8 &&
        !is.na(mu_k) && abs(mu_k) < 1e-8 &&
        !is.na(lb) && lb == 0 &&
        !is.na(ub) && !is.finite(ub)

      if (isTRUE(is_hc)) {
        nodes$is_halfcauchy_like[k] <- TRUE
      }
    }
  }

  cat("\n[NON-DIFF INDICATORS]\n")
  cat(sprintf(
    "- Non-diff functions detected: %s\n",
    if (length(nondiff_hits)) paste(sort(unique(nondiff_hits)), collapse = ",") else "None"
  ))

  if (nzchar(code_txt_scanned)) {
    dists_in_code_raw <- unique(tolower(unlist(regmatches(
      code_txt_scanned,
      gregexpr("\\bd[a-zA-Z_]+\\b", code_txt_scanned, perl = TRUE)
    ))))
    dists_in_code <- sort(unique(
      vapply(dists_in_code_raw, normalize_dist, character(1))
    ))
  } else {
    dists_in_code <- sort(unique(stats::na.omit(nodes$dist[nodes$is_stoch])))
  }
  cat(sprintf(
    "- Distributions found in code : %s\n",
    if (length(dists_in_code)) paste(dists_in_code, collapse = ", ") else "None"
  ))

  bounded_latent_trunc <- any(
    nodes$is_stoch & !nodes$is_data &
      (nodes$support == "continuous") &
      (is.finite(nodes$lower) | is.finite(nodes$upper)),
    na.rm = TRUE
  )

  cat(sprintf(
    "- BUGS-style truncation 'T(a,b)' spotted: %s\n",
    if (isTRUE(truncation_detected_bugst)) "Yes" else "No"
  ))
  cat(sprintf(
    "- Bounded latent nodes (implicit truncation via finite bounds on continuous support): %s\n",
    if (bounded_latent_trunc) "Yes" else "No"
  ))

  # ---- Tagging non-diff-deterministic-op ----
  if (length(nondiff_hits)) {
    nodes$hmc_showstopper_reason <- nodes$hmc_showstopper_reason %||% NA_character_

    if (nzchar(code_txt_scanned) && any(nondiff_hits == "round")) {
      code_lines  <- strsplit(code_txt_scanned, "\n", fixed = TRUE)[[1]]
      round_lines <- grep("\\bround\\s*\\(", code_lines, value = TRUE, perl = TRUE)

      vars_in_lines <- unique(gsub("\\[.*?\\]$", "", unlist(regmatches(
        round_lines,
        gregexpr("\\b[A-Za-z_][A-Za-z0-9_]*(\\[[^\\]]+\\])?", round_lines, perl = TRUE)
      ))))

      vars_in_model <- unique(nodes$var)
      seed_vars     <- intersect(vars_in_lines, vars_in_model)

      if (length(seed_vars)) {
        affected_latents <- unique(unlist(lapply(
          seed_vars,
          function(v) get_deps_memo(
            v,
            upstream    = FALSE,
            downstream  = TRUE,
            includeData = FALSE,
            stochOnly   = TRUE
          )
        )))

        idx_aff <- match(affected_latents, nodes$node, nomatch = 0L)
        idx_aff <- idx_aff[idx_aff > 0L]
        if (length(idx_aff)) {
          nodes$hmc_showstopper_reason[idx_aff] <- ifelse(
            is.na(nodes$hmc_showstopper_reason[idx_aff]),
            "non-diff-deterministic-op",
            nodes$hmc_showstopper_reason[idx_aff]
          )
        }
      }
    }
  }

  # ---- Generic showstoppers ----
  if (!"hmc_showstopper_reason" %in% names(nodes)) nodes$hmc_showstopper_reason <- NA_character_
  latent_mask <- nodes$is_stoch & !nodes$is_data

  nodes$hmc_showstopper_reason[
    latent_mask & nodes$support == "discrete" &
      is.na(nodes$hmc_showstopper_reason)
  ] <- "discrete-latent"

  nodes$hmc_showstopper_reason[
    latent_mask & nodes$support == "simplex" &
      is.na(nodes$hmc_showstopper_reason)
  ] <- "simplex-constraint"

  nodes$hmc_showstopper_reason[
    latent_mask & nodes$is_truncated &
      is.na(nodes$hmc_showstopper_reason)
  ] <- "truncation"

  nodes$hmc_showstopper_reason[
    latent_mask &
      (nodes$dist == "dhalfcauchy" | nodes$is_halfcauchy_like) &
      is.na(nodes$hmc_showstopper_reason)
  ] <- "half-cauchy-latent"

  summarise_showstoppers <- function(nodes_df) {
    ss <- nodes_df
    if (!"hmc_showstopper_reason" %in% names(ss)) ss$hmc_showstopper_reason <- NA_character_
    if (!"is_data" %in% names(ss))               ss$is_data               <- FALSE
    ss$is_data <- to_logical(ss$is_data)
    ss <- ss[!is.na(ss$hmc_showstopper_reason) & !ss$is_data, , drop = FALSE]

    if (!nrow(ss)) {
      return(list(
        summary  = data.frame(reason = character(0), n_nodes = integer(0)),
        examples = data.frame(reason = character(0), examples = character(0))
      ))
    }

    tt <- as.data.frame(table(ss$hmc_showstopper_reason), stringsAsFactors = FALSE)
    names(tt) <- c("reason","n_nodes")
    tt$n_nodes <- as.integer(tt$n_nodes)
    tt <- tt[order(-tt$n_nodes, tt$reason), , drop = FALSE]

    ex <- do.call(rbind, lapply(split(ss, ss$hmc_showstopper_reason), function(d) {
      data.frame(
        reason   = as.character(d$hmc_showstopper_reason[1]),
        examples = paste(utils::head(unique(as.character(d$node)), max_show_examples), collapse = ","),
        stringsAsFactors = FALSE
      )
    }))
    rownames(ex) <- NULL
    list(summary = tt, examples = ex)
  }

  ss <- summarise_showstoppers(nodes)
  cat("\n[HMC/NUTS SHOWSTOPPERS]\n")
  if (nrow(ss$summary)) {
    print(ss$summary, row.names = FALSE)
    cat("\nExamples by reason:\n")
    print(ss$examples, row.names = FALSE)
  } else {
    reasons <- unique(stats::na.omit(nodes$hmc_showstopper_reason[!nodes$is_data]))
    if (length(reasons)) {
      cat("(no table -- fallback summary)\n")
      cat("Reasons: ", paste(sort(reasons), collapse = ", "), "\n", sep = "")
    } else {
      cat("- No explicit showstoppers detected at latent nodes.\n")
    }
  }

  hsk <- nodes$hmc_showstopper_reason %||% rep(NA_character_, nrow(nodes))
  hmc_globally_ok <- all(is.na(hsk[!nodes$is_data]))

  ## --- STRICT NON-DIFF RULE (global veto for HMC/NUTS) ---
  ## Any occurrence of these functions anywhere in the model code
  ## will force HMC/NUTS to be considered structurally unsuitable.
  strict_nondiff <- c("round", "floor", "ceiling", "trunc", "step", "equals")

  if (length(nondiff_hits)) {
    bad_funs <- intersect(tolower(nondiff_hits), strict_nondiff)
    if (length(bad_funs) > 0L) {
      hmc_globally_ok <- FALSE

      ## Optionally tag all latent nodes with a generic reason if none exists,
      ## so the message in the smoke test has something explicit to report.
      if (all(is.na(nodes$hmc_showstopper_reason[latent_mask]))) {
        nodes$hmc_showstopper_reason[latent_mask] <- "global-nondiff-function"
        ss <- summarise_showstoppers(nodes)  # refresh summary if needed later
      }

      cat("\n[STRICT NON-DIFF RULE]\n")
      cat(
        "- Found non-differentiable function(s) considered fatal for HMC/NUTS: ",
        paste(sort(unique(bad_funs)), collapse = ", "),
        "\n- Overriding structural decision: HMC/NUTS globally disabled.\n",
        sep = ""
      )
    }
  }

  cat(sprintf(
    "\n- Structurally suitable for HMC/NUTS? %s\n",
    if (isTRUE(hmc_globally_ok)) "Yes" else "No"
  ))

  # ------ Optional HMC/NUTS smoke test ------
  hmc_res <- list(ok = FALSE, error = "HMC test not requested.", details = NULL)
  if (isTRUE(try_hmc)) {
    cat("\n[HMC/NUTS SMOKE TEST]\n")
    has_hmc <- suppressWarnings(requireNamespace("nimbleHMC", quietly = TRUE))
    if (!has_hmc) {
      msg <- "Package 'nimbleHMC' is not available (not installed or not loaded)."
      cat("- HMC/NUTS empirical test: NOT RUN\n  reason: ", msg, "\n", sep = "")
      hmc_res <- list(ok = FALSE, error = msg, details = NULL)
    } else if (!isTRUE(hmc_globally_ok)) {
      reasons <- unique(stats::na.omit(nodes$hmc_showstopper_reason[!nodes$is_data]))
      if (!length(reasons)) {
        reasons <- "global-nondiff-function"
      }
      msg <- sprintf(
        "Model is not structurally suitable for HMC/NUTS: %s",
        paste(reasons, collapse = ", ")
      )
      cat("- HMC/NUTS empirical test: NOT RUN\n  reason: ", msg, "\n", sep = "")
      hmc_res <- list(ok = FALSE, error = msg, details = NULL)
    } else {
      model_has_derivs <- .sop_supports_derivs(m)
      if (!isTRUE(model_has_derivs)) {
        msg <- "NUTS not attempted: the model was not built with buildDerivs = TRUE."
        cat("- HMC/NUTS empirical test: NOT RUN\n  reason: ", msg, "\n", sep = "")
        hmc_res <- list(ok = FALSE, error = msg, details = NULL)
      } else {
        base_conf <- try(nimble::configureMCMC(m), silent = TRUE)
        if (inherits(base_conf, "try-error")) {
          msg <- paste("configureMCMC failed:", as.character(base_conf))
          cat("- HMC/NUTS empirical test: NOT RUN\n  reason: ", msg, "\n", sep = "")
          hmc_res <- list(ok = FALSE, error = msg, details = NULL)
        } else {
          if (length(mons)) try(base_conf$addMonitors(mons), silent = TRUE)

          candidate_nodes <- nodes$node[
            latent_mask & is.na(nodes$hmc_showstopper_reason)
          ]
          candidate_nodes <- unique(as.character(candidate_nodes[!is.na(candidate_nodes)]))

          if (length(candidate_nodes)) {
            valid_nodes <- m$getNodeNames(includeData = include_data)
            candidate_nodes <- intersect(candidate_nodes, valid_nodes)
          }

          okH   <- TRUE
          err_m <- NULL
          tryCatch({
            if (length(candidate_nodes)) {
              nimbleHMC::configureHMC(base_conf, model = m, nodes = candidate_nodes)
            } else {
              nimbleHMC::configureHMC(base_conf, model = m)
            }
          }, error = function(e) {
            okH <- FALSE
            # intentionally generic message (no localisation propagation)
            err_m <- "configureHMC failed (internal error: model/sampler configuration not supported; see ?nimbleHMC::configureHMC)."
          })

          if (!okH) {
            msg <- err_m %||% "configureHMC failed (internal error; HMC configuration not supported for this model)."
            cat("- HMC/NUTS empirical test: FAILED\n  reason: ", msg, "\n", sep = "")
            hmc_res <- list(ok = FALSE, error = msg, details = NULL)
          } else {
            mcmc  <- try(nimble::buildMCMC(base_conf), silent = TRUE)
            if (inherits(mcmc, "try-error")) {
              msg <- paste("buildMCMC failed:", as.character(mcmc))
              cat("- HMC/NUTS empirical test: FAILED\n  reason: ", msg, "\n", sep = "")
              hmc_res <- list(ok = FALSE, error = msg, details = NULL)
            } else {
              cm_obj <- if (!is.null(cm)) cm else try(nimble::compileNimble(m), silent = TRUE)
              cmcmc  <- try(
                nimble::compileNimble(mcmc, project = cm_obj, resetFunctions = TRUE),
                silent = TRUE
              )
              if (inherits(cmcmc, "try-error")) {
                msg <- paste("compileNimble failed:", as.character(cmcmc))
                cat("- HMC/NUTS empirical test: FAILED\n  reason: ", msg, "\n", sep = "")
                hmc_res <- list(ok = FALSE, error = msg, details = NULL)
              } else {
                set.seed(seed)
                run <- try(
                  nimble::runMCMC(
                    cmcmc,
                    niter   = as.integer(niter),
                    nburnin = as.integer(nburnin),
                    thin    = 1L, nchains = 1L,
                    samplesAsCodaMCMC = FALSE,
                    summary = FALSE, WAIC = FALSE
                  ),
                  silent = TRUE
                )
                if (inherits(run, "try-error")) {
                  msg <- paste("runMCMC failed:", as.character(run))
                  cat("- HMC/NUTS empirical test: FAILED\n  reason: ", msg, "\n", sep = "")
                  hmc_res <- list(ok = FALSE, error = msg, details = NULL)
                } else {
                  cat("- HMC/NUTS empirical test: PASSED\n")
                  hmc_res <- list(ok = TRUE, error = NULL, details = run)
                }
              }
            }
          }
        }
      }
    }
  }

  invisible(list(
    diag = list(
      nodes = nodes,
      nondiff_signals = list(
        functions_found      = sort(unique(nondiff_hits)),
        distributions_found  = sort(unique(dists_in_code)),
        truncation_bugst     = truncation_detected_bugst,
        bounded_latent_trunc = bounded_latent_trunc
      ),
      code_scan       = code_txt_scanned,
      hmc_globally_ok = hmc_globally_ok
    ),
    hmc  = hmc_res
  ))
}

#' Test and compare MCMC strategies on selected bottleneck nodes
#'
#' Runs a small, reproducible workflow to:
#' \enumerate{
#'   \item build and run a baseline MCMC configuration,
#'   \item optionally try full-model HMC/NUTS (if \code{nimbleHMC} is available
#'         and the model supports derivatives),
#'   \item then apply alternative samplers in a \emph{strict}, user-defined order
#'         on one or two bottleneck targets (singleton and then optional block
#'         on the union).
#' }
#'
#' Results and diagnostic plots (R-hat bars and trace/density) are written under
#' \code{out_dir}. When \code{ask = TRUE}, interactive yes/no prompts allow you
#' to stop early at each step.
#'
#' @details
#' The procedure:
#' \enumerate{
#'   \item Build and run a baseline MCMC using \code{run_baseline_config()}.
#'   \item Optionally run full-model HMC/NUTS via \code{configure_hmc_safely()}.
#'   \item Select \code{nbot} bottleneck node(s) from diagnostics (or from
#'         \code{force_singletons}), then:
#'     \itemize{
#'       \item apply \code{strict_scalar_seq} in order on the first node;
#'       \item if \code{nbot >= 2}, build the union \code{\{node1, node2\}}
#'             (or \code{force_union_nodes}) and apply \code{strict_block_seq}.
#'     }
#'   \item For each step, compile, run, compute diagnostics, and save plots.
#' }
#'
#' @param build_fn Function (or prebuilt list with \code{$model}, \code{$conf})
#'   that returns a fresh build object used by samplers in this package.
#' @param monitors Optional character vector of monitors passed to the build.
#' @param try_hmc Logical. If \code{TRUE}, try a full-model HMC/NUTS run first
#'   (ignored by the “surgical” singleton/block steps).
#' @param nchains Integer number of MCMC chains for each run.
#' @param pilot_niter Integer total iterations used in baseline and tests.
#' @param pilot_burnin Integer burn-in iterations.
#' @param thin Integer thinning interval.
#' @param out_dir Directory where outputs (plots, etc.) will be written.
#' @param nbot Integer. Number of bottleneck targets to operate on (1 or 2
#'   typical; if \code{>= 2}, a block step on the union is attempted).
#' @param strict_scalar_seq Character vector of scalar samplers to try in order.
#'   Supported values include \code{"NUTS"}, \code{"slice"}, \code{"RW"}.
#' @param strict_block_seq Character vector of block samplers to try in order.
#'   Supported values include \code{"NUTS_block"}, \code{"AF_slice"},
#'   \code{"RW_block"}.
#' @param force_singletons Optional character vector of node names to force as
#'   singleton targets (first \code{nbot} valid nodes are used).
#' @param force_union_nodes Optional character vector of node names to define
#'   the union for the block phase (must contain \code{>= 2} valid nodes).
#' @param force_union Deprecated alias of \code{force_union_nodes}.
#' @param ask Logical. If \code{TRUE}, ask before moving to the next step.
#' @param ask_before_hmc Logical. If \code{TRUE}, ask before running full HMC.
#' @param block_max Integer cap on the size of the block union.
#' @param slice_control List of controls passed to \code{"slice"} samplers.
#' @param rw_control List of controls passed to \code{"RW"} samplers.
#' @param rwblock_control List of controls passed to \code{"RW_block"} samplers.
#' @param af_slice_control List of controls passed to \code{"AF_slice"} samplers.
#' @param slice_max_contractions Integer safety cap for slice contractions.
#'
#' @return A list with elements such as:
#' \describe{
#'   \item{status}{Character status string.}
#'   \item{mode}{Character mode (e.g. \code{"HMC_full"}, \code{"surgical_*"}).}
#'   \item{baseline}{Baseline run info (runtime, samples, diagnostics).}
#'   \item{targets}{Chosen target node names.}
#'   \item{steps}{List of steps; each contains nodes, sampler, results, and directory.}
#' }
#'
#' @section Side effects:
#' Creates subfolders and PNG files (R-hat bars, traces/densities) under
#' \code{out_dir}. May load/unload compiled DLLs while switching samplers.
#'
#' @seealso \code{\link{configure_hmc_safely}},
#'   \code{\link{plot_convergence_checks}}, \code{\link{plot_bottlenecks}}
#'
#' @examples
#' \dontrun{
#' res <- test_strategy(
#'   build_fn = my_build_fn,
#'   monitors = c("theta","beta"),
#'   try_hmc  = TRUE,
#'   nbot     = 2,
#'   out_dir  = "outputs/diagnostics"
#' )
#' }
#'
#' @export
#' @keywords internal
test_strategy <- function(build_fn,
                          monitors            = NULL,   # optional, pass-through
                          try_hmc             = TRUE,   # full-model HMC path only; surgical ignores
                          nchains             = 3L,
                          pilot_niter         = 2e4,
                          pilot_burnin        = 5e3,
                          thin                = 1L,
                          out_dir             = "outputs/diagnostics",
                          nbot                = 1L,
                          # strict sequences (order enforced)
                          strict_scalar_seq   = c("NUTS","slice","RW"),
                          strict_block_seq    = c("NUTS_block","AF_slice","RW_block"),
                          # forcing
                          force_singletons    = NULL,    # e.g. c("theta[1]","beta")
                          force_union_nodes   = NULL,    # e.g. c("theta[1]","beta")
                          # alias kept for your previous calls
                          force_union         = NULL,
                          # interaction
                          ask                 = TRUE,
                          ask_before_hmc      = TRUE,
                          # safety caps
                          block_max           = 20L,
                          # sampler controls
                          slice_control       = list(),
                          rw_control          = list(),
                          rwblock_control     = list(adaptScaleOnly = TRUE),
                          af_slice_control    = list(),
                          slice_max_contractions = 5000L) {

  `%||%` <- function(x, y) if (is.null(x)) y else x
  stopifnot(nbot >= 1L)
  if (is.list(build_fn) && !is.null(build_fn$model)) {
    obj <- build_fn
    build_fn <- function() obj
  }
  stopifnot(is.function(build_fn))
  if (!requireNamespace("nimble", quietly = TRUE))
    stop("samOptiPro: 'nimble' is required.")
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- utils
  say   <- function(...) {
    msg <- try(sprintf(...), silent = TRUE)
    if (!inherits(msg,"try-error")) cat(msg,"\n")
  }
  is_ll <- function(x) {
    grepl("^logLik(\\[.*\\])?$|log_?lik|logdens|lpdf",
          x, perl = TRUE, ignore.case = TRUE)
  }

  # (A) sanitize niter/burnin once, and reuse
  .sanitize_iters <- local({
    warned <- FALSE
    function(niter, nburnin) {
      niter   <- as.integer(niter)
      nburnin <- as.integer(nburnin)
      if (is.na(niter)   || niter   < 1L) niter   <- 1L
      if (is.na(nburnin) || nburnin < 0L) nburnin <- 0L
      if (nburnin >= niter) {
        new_nburnin <- max(0L, niter - 1L)
        if (!warned) {
          message(sprintf(
            "[info] Adjusting nburnin (%d) to %d because nburnin must be < niter (%d).",
            nburnin, new_nburnin, niter
          ))
          warned <<- TRUE
        }
        nburnin <- new_nburnin
      }
      list(niter = niter, nburnin = nburnin)
    }
  })

  # (B) strict yes/no with flush + loop
  .ask_yes_no_strict <- function(prompt) {
    if (!isTRUE(ask) || !interactive()) return(TRUE)
    while (TRUE) {
      cat(paste0(prompt, " (yes/no): "))
      utils::flush.console()
      ans <- try(readline(), silent = TRUE)
      if (inherits(ans, "try-error")) return(TRUE) # non-interactive fallback = proceed
      ans <- tolower(trimws(ans))
      if (ans %in% c("y","yes","o","oui"))         return(TRUE)
      if (ans %in% c("n","no","non","0"))         return(FALSE)
      cat("Please answer 'yes' or 'no'.\n")
    }
  }

  # barrier & restore for smooth scalar->block
  .barrier_before_block <- function(tag = "scalar_to_block") {
    try(nimble::clearCompiled(), silent = TRUE); gc()
    new_bd <- file.path(
      tempdir(),
      paste0(
        "samOptiPro_build_",
        gsub("[^A-Za-z0-9_]+","_", tag), "_",
        as.integer(stats::runif(1, 1e9, 9e9))
      )
    )
    dir.create(new_bd, recursive = TRUE, showWarnings = FALSE)
    assign(".sop_old_builddir", nimble::nimbleOptions("buildDir"), inherits = FALSE)
    nimble::nimbleOptions(buildDir = new_bd)
    if (.Platform$OS.type == "windows") Sys.sleep(0.3) else Sys.sleep(0.05)
    invisible(TRUE)
  }
  .restore_builddir <- function() {
    old <- try(get(".sop_old_builddir", inherits = FALSE), silent = TRUE)
    if (!inherits(old, "try-error") && !is.null(old))
      try(nimble::nimbleOptions(buildDir = old), silent = TRUE)
    invisible(TRUE)
  }

  .ensure_unsampled <- function(conf) {
    uns <- try(conf$getUnsampledNodes(), silent = TRUE)
    if (!inherits(uns, "try-error") && length(uns)) {
      uns <- uns[!is_ll(uns)]
      for (u in uns) conf$addSampler(u, type = "slice")
    }
  }

  # ---------- compile & run (aligned with family version: use compute_diag_from_mcmc) ----------
  .compile_and_run <- function(step_tag, build_obj, conf, niter, nburnin) {
    it <- .sanitize_iters(niter, nburnin)
    niter   <- it$niter
    nburnin <- it$nburnin
    attempts <- 0L
    last_err <- NULL
    repeat {
      attempts <- attempts + 1L
      try(nimble::clearCompiled(), silent = TRUE); gc()
      res <- try({
        .ensure_unsampled(conf)
        cmcmc <- .compile_mcmc_with_build(conf, build_obj, reset = TRUE, show = FALSE)
        out   <- .run_and_collect(
          cmcmc,
          niter   = niter,
          nburnin = nburnin,
          thin    = thin,
          nchains = nchains
        )
        ml <- as_mcmc_list_sop(
          out$samples, out$samples2,
          drop_loglik = FALSE, thin = thin
        )
        # diagnostics: same pipeline as baseline
        dg <- compute_diag_from_mcmc(ml, runtime_s = out$runtime_s)
        if (!is.null(cmcmc) && is.list(cmcmc) && "unloadDLL" %in% names(cmcmc))
          try(cmcmc$unloadDLL(), silent = TRUE)
        list(out = out, ml = ml, dg = dg)
      }, silent = TRUE)
      if (!inherits(res, "try-error")) return(res)
      last_err <- tryCatch(
        conditionMessage(attr(res, "condition")),
        error = function(e) as.character(res)
      )
      if (attempts == 1L)
        message(sprintf("[info] First attempt failed at %s; retrying cleanly...", step_tag))
      else
        message(sprintf("[retry %d @ %s] %s", attempts, step_tag, last_err))

      build_obj <- .fresh_build(build_fn, monitors = monitors, thin = thin)
      conf      <- build_obj$conf
      if (attempts >= 3L)
        stop(sprintf(
          "Compilation failed after %d attempts at step '%s': %s",
          attempts, step_tag, last_err
        ))
    }
  }

  # assigners
  .add_scalar <- function(conf, nodes, type, has_hmc) {
    if (identical(type, "NUTS")) {
      if (isTRUE(has_hmc)) {
        for (n in nodes) conf$addSampler(target = n, type = "NUTS")
      } else {
        cat("[Info] nimbleHMC not available -> falling back to slice (scalar).\n")
        for (n in nodes) conf$addSampler(n, "slice", slice_control)
      }
      return(invisible())
    }
    if (identical(type, "slice")) {
      for (n in nodes) conf$addSampler(n, "slice", slice_control)
      return(invisible())
    }
    if (identical(type, "RW")) {
      for (n in nodes) conf$addSampler(n, "RW", rw_control)
      return(invisible())
    }
    for (n in nodes) conf$addSampler(n, "slice", slice_control)
  }

  .add_block <- function(conf, nodes, type, has_hmc, Sc = NULL) {
    nodes <- unique(nodes)
    if (length(nodes) < 2L) {
      if (type == "NUTS_block") .add_scalar(conf, nodes, "NUTS", has_hmc) else
        if (type == "AF_slice")   for (n in nodes) conf$addSampler(n, "AF_slice", af_slice_control) else
          if (type == "RW_block") for (n in nodes) conf$addSampler(n, "RW", rw_control) else
            for (n in nodes) conf$addSampler(n, "slice", slice_control)
      return(invisible())
    }
    if (identical(type, "NUTS_block")) {
      if (isTRUE(has_hmc)) {
        conf$addSampler(target = nodes, type = "NUTS")
      } else {
        cat("[Info] nimbleHMC not available -> falling back to AF_slice (block).\n")
        conf$addSampler(target = nodes, type = "AF_slice", control = af_slice_control)
      }
      return(invisible())
    }
    if (identical(type, "AF_slice")) {
      ok <- TRUE
      tryCatch({
        conf$addSampler(target = nodes, type = "AF_slice", control = af_slice_control)
      }, error = function(e) ok <<- FALSE)
      if (!ok) for (n in nodes) conf$addSampler(n, "slice", slice_control)
      return(invisible())
    }
    if (identical(type, "RW_block")) {
      ctrl <- rwblock_control
      if (!is.null(Sc)) ctrl$propCov <- Sc
      conf$addSampler(nodes, "RW_block", ctrl)
      return(invisible())
    }
    for (n in nodes) conf$addSampler(n, "slice", slice_control)
  }

  .metrics_for <- function(dg, nodes) {
    nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
    if (is.null(dg) || !nrow(dg)) {
      return(list(AE = NA_real_, CE = NA_real_, Rhat = NA_real_))
    }
    if (!"target" %in% names(dg)) {
      if ("Parameter" %in% names(dg)) {
        dg$target <- as.character(dg$Parameter)
      } else {
        return(list(AE = NA_real_, CE = NA_real_, Rhat = NA_real_))
      }
    }
    keep <- dg$target %in% nodes & !is_ll(dg$target)
    if (!any(keep)) {
      return(list(AE = NA_real_, CE = NA_real_, Rhat = NA_real_))
    }
    list(
      AE   = stats::median(dg$AE_ESS_per_it[keep], na.rm = TRUE),
      CE   = stats::median(dg$ESS_per_sec[keep],   na.rm = TRUE),
      Rhat = if (any(keep)) suppressWarnings(max(dg$Rhat[keep], na.rm = TRUE)) else NA_real_
    )
  }

  .print_step <- function(title, label, nodes, sampler, runtime, met, base_dx, base_rt) {
    nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
    say("--- %s ---", title)
    if (!is.null(label)) say("Target(s): %s", label)
    say("Nodes: %s", paste(nodes, collapse = ", "))
    say("Sampler: %s", sampler)
    say("Runtime_s: %.2f (baseline: %.2f)", runtime %||% NA_real_, base_rt %||% NA_real_)
    say("AE median (ESS/iter): %.3g (baseline: %.3g)",
        met$AE %||% NA_real_,
        stats::median(base_dx$AE_ESS_per_it, na.rm = TRUE))
    say("CE median (ESS/s):    %.3g (baseline: %.3g)",
        met$CE %||% NA_real_,
        stats::median(base_dx$ESS_per_sec,   na.rm = TRUE))
    say("Rhat max:             %.3g (baseline max: %.3g)",
        met$Rhat %||% NA_real_,
        suppressWarnings(max(base_dx$Rhat, na.rm = TRUE)))
  }

  # ---------- 0) Build + baseline ----------
  b0  <- .fresh_build(build_fn, monitors = monitors, thin = thin)
  mdl <- b0$model
  it  <- .sanitize_iters(pilot_niter, pilot_burnin)

  rb  <- run_baseline_config(
    build_fn = build_fn,
    niter    = it$niter,
    nburnin  = it$nburnin,
    thin     = thin,
    monitors = monitors,
    nchains  = nchains
  )

  base_ml <- as_mcmc_list_sop(
    rb$samples, rb$samples2,
    drop_loglik = FALSE, thin = thin
  )
  base_dg <- compute_diag_from_mcmc(base_ml, runtime_s = rb$runtime_s)
  base_dx <- base_dg[!is_ll(base_dg$target), , drop = FALSE]

  say("Baseline runtime_s: %.2f s", rb$runtime_s %||% NA_real_)
  say("Baseline median AE(ESS/iter): %.3g",
      stats::median(base_dx$AE_ESS_per_it, na.rm = TRUE))
  say("Baseline median CE(ESS/s):    %.3g",
      stats::median(base_dx$ESS_per_sec,   na.rm = TRUE))

  # ---------- 1) Optional full-model HMC (using run_hmc_all_nodes) ----------
  dg_struct <- try(diagnose_model_structure(mdl), silent = TRUE)
  suppressWarnings(has_hmc <- requireNamespace("nimbleHMC", quietly = TRUE))
  deriv_ok <- .sop_supports_derivs(mdl)
  blockers <- character(0)
  if (!inherits(dg_struct,"try-error") && !is.null(dg_struct)) {
    if (isTRUE(dg_struct$has_truncation))   blockers <- c(blockers, "truncation")
    if (isTRUE(dg_struct$has_simplex))      blockers <- c(blockers, "simplex-constraint")
    if (isTRUE(dg_struct$has_non_diff_fun)) blockers <- c(blockers, "non-diff-function")
  }
  nuts_ok_global <- isTRUE(try_hmc) && has_hmc && deriv_ok && (length(blockers) == 0L)

  if (isTRUE(nuts_ok_global)) {
    proceed_hmc <- .ask_yes_no_strict(sprintf(
      "Baseline ready. Runtime=%.2fs; median AE=%.3g; median CE=%.3g\nProceed with full-model HMC/NUTS?",
      rb$runtime_s %||% NA_real_,
      median(base_dx$AE_ESS_per_it, na.rm = TRUE),
      median(base_dx$ESS_per_sec,   na.rm = TRUE)
    ))

    if (proceed_hmc) {
      hmc_try <- try(
        run_hmc_all_nodes(
          build_fn = build_fn,
          niter    = it$niter,
          nburnin  = it$nburnin,
          thin     = thin,
          monitors = monitors,
          nchains  = nchains,
          opts     = samOptiPro_options()
        ),
        silent = TRUE
      )

      if (inherits(hmc_try, "try-error")) {
        cat("[Warn] Full-model HMC (run_hmc_all_nodes) threw an error -> continuing with surgical singleton plan.\n")
      } else if (!is.list(hmc_try) ||
                 !all(c("samples","samples2","runtime_s") %in% names(hmc_try))) {
        cat("[Warn] run_hmc_all_nodes returned an unexpected structure -> continuing with surgical singleton plan.\n")
      } else {
        hmc_ml <- as_mcmc_list_sop(
          hmc_try$samples,
          hmc_try$samples2,
          drop_loglik = FALSE,
          thin = thin
        )
        hmc_dg <- compute_diag_from_mcmc(hmc_ml, runtime_s = hmc_try$runtime_s)

        cat(sprintf("HMC runtime_s: %.3f\n", hmc_try$runtime_s %||% NA_real_))
        cat(sprintf("HMC median AE: %.3g ; CE: %.3g ; max Rhat: %.3g\n",
                    median(hmc_dg$AE_ESS_per_it, na.rm = TRUE),
                    median(hmc_dg$ESS_per_sec,   na.rm = TRUE),
                    if (all(is.na(hmc_dg$Rhat))) NA_real_ else max(hmc_dg$Rhat, na.rm = TRUE)))

        return(list(
          status   = "completed",
          mode     = "HMC_full",
          baseline = list(
            runtime_s = rb$runtime_s,
            samples   = base_ml,
            diag_tbl  = base_dg
          ),
          hmc      = list(
            run      = hmc_try,
            mcmc     = hmc_ml,
            diag_tbl = hmc_dg
          ),
          messages = "Full HMC (run_hmc_all_nodes) completed."
        ))
      }
    } else {
      cat("User declined full-model HMC/NUTS. Switching to surgical singleton plan.\n")
    }
  }

  # ---------- 2) Choose singleton bottlenecks (robust + forcing) ----------
  stoch_nodes <- mdl$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  stoch_nodes <- stoch_nodes[!is_ll(stoch_nodes)]

  if (is.null(force_union_nodes) && !is.null(force_union))
    force_union_nodes <- force_union

  pick_nodes <- character(0)

  # From forcing
  if (!is.null(force_singletons) && length(force_singletons)) {
    pn      <- unique(force_singletons)
    present <- intersect(
      pn,
      intersect(stoch_nodes, unique(colnames(as.matrix(base_ml[[1]]))))
    )
    pick_nodes <- utils::head(present, nbot)
  }

  # From diagnostics at node-level
  if (!length(pick_nodes)) {
    node_like <- grepl("\\[", base_dx$target) # heuristic
    if (any(node_like)) {
      ord  <- order(base_dx$ESS_per_sec[node_like], decreasing = FALSE)
      cand <- base_dx$target[node_like][ord]
      pick_nodes <- intersect(utils::head(cand, nbot), stoch_nodes)
    }
  }

  # From samples ESS/s
  if (!length(pick_nodes)) {
    if (!requireNamespace("coda", quietly = TRUE))
      stop("samOptiPro: 'coda' package required for node-level ESS.")
    common_cols <- Reduce(intersect, lapply(base_ml, function(m) colnames(as.matrix(m))))
    node_cols   <- intersect(common_cols, stoch_nodes)
    if (length(node_cols)) {
      ess_per_node <- sapply(node_cols, function(nm) {
        ess_ch <- sapply(base_ml, function(m) {
          X <- as.matrix(m)
          if (nm %in% colnames(X)) coda::effectiveSize(X[, nm]) else NA_real_
        })
        sum(ess_ch, na.rm = TRUE) / (rb$runtime_s %||% 1)
      })
      ord <- order(ess_per_node, decreasing = FALSE)
      pick_nodes <- utils::head(names(ess_per_node)[ord], nbot)
      pick_nodes <- pick_nodes[!is.na(pick_nodes)]
    }
  }

  # Last resort
  if (!length(pick_nodes)) {
    present <- intersect(stoch_nodes, unique(colnames(as.matrix(base_ml[[1]]))))
    pick_nodes <- utils::head(present, nbot)
  }
  if (!length(pick_nodes))
    stop("No stochastic singleton targets to operate on.")

  say("Selected singleton bottlenecks (nbot=%d): %s",
      length(pick_nodes), paste(pick_nodes, collapse = ", "))

  steps <- list()

  # ======================= CASE nbot = 1 =======================
  if (nbot == 1L) {
    node1 <- pick_nodes[1]
    for (samp in unique(as.character(strict_scalar_seq))) {
      bS    <- .fresh_build(build_fn, monitors = monitors, thin = thin)
      confS <- bS$conf
      try(confS$removeSamplers(node1), silent = TRUE)
      .add_scalar(confS, node1, samp, has_hmc)

      itS  <- .sanitize_iters(pilot_niter, pilot_burnin)
      resS <- .compile_and_run(
        paste0("nbot1_", gsub("[^A-Za-z0-9_]+","_", node1), "_", samp),
        bS, confS, itS$niter, itS$nburnin
      )
      metS <- .metrics_for(resS$dg, node1)
      .print_step(
        "Scalar plan on singleton",
        node1, node1, samp,
        resS$out$runtime_s, metS, base_dx, rb$runtime_s
      )

      pdir <- file.path(
        out_dir,
        sprintf("singleton_%s_%s", gsub("[^A-Za-z0-9_]", "_", node1), samp)
      )
      if (!dir.exists(pdir))
        dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
      .plot_rhat_bar(resS$dg, nodes = node1,
                     out_file = file.path(pdir, "rhat_bar.png"))
      .plot_traces(resS$ml, nodes = node1,
                   out_file_prefix = file.path(pdir, "trace_"))

      steps <- c(steps, list(list(
        level   = "singleton-scalar",
        nodes   = node1,
        sampler = samp,
        res     = resS,
        dir     = pdir
      )))

      if (!.ask_yes_no_strict(sprintf("Proceed to the next sampler for '%s'?", node1)))
        break
    }
    return(list(
      status  = "completed",
      mode    = "surgical_nbot1_singleton",
      baseline = rb,
      targets  = pick_nodes,
      steps    = steps
    ))
  }

  # ======================= CASE nbot >= 2 =======================
  # Step 1 -- singleton 1 scalar strict
  node1 <- pick_nodes[1]
  for (samp in unique(as.character(strict_scalar_seq))) {
    bS    <- .fresh_build(build_fn, monitors = monitors, thin = thin)
    confS <- bS$conf
    try(confS$removeSamplers(node1), silent = TRUE)
    .add_scalar(confS, node1, samp, has_hmc)

    itS  <- .sanitize_iters(pilot_niter, pilot_burnin)
    resS <- .compile_and_run(
      paste0("nbot2_scalar_", gsub("[^A-Za-z0-9_]+","_", node1), "_", samp),
      bS, confS, itS$niter, itS$nburnin
    )
    metS <- .metrics_for(resS$dg, node1)
    .print_step(
      "Scalar plan on singleton",
      node1, node1, samp,
      resS$out$runtime_s, metS, base_dx, rb$runtime_s
    )

    pdir <- file.path(
      out_dir,
      sprintf("singleton_%s_%s", gsub("[^A-Za-z0-9_]", "_", node1), samp)
    )
    if (!dir.exists(pdir))
      dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    .plot_rhat_bar(resS$dg, nodes = node1,
                   out_file = file.path(pdir, "rhat_bar.png"))
    .plot_traces(resS$ml, nodes = node1,
                 out_file_prefix = file.path(pdir, "trace_"))

    steps <- c(steps, list(list(
      level   = "singleton-scalar",
      nodes   = node1,
      sampler = samp,
      res     = resS,
      dir     = pdir
    )))

    if (!.ask_yes_no_strict(sprintf("Proceed to the next sampler for '%s'?", node1)))
      break
  }

  # Step 2 -- Block strict on union {node1 ? node2} (or forced union)
  node2 <- pick_nodes[2] %||% NA_character_

  union_nodes <- if (!is.null(force_union_nodes) && length(force_union_nodes) >= 2L) {
    present <- intersect(
      unique(force_union_nodes),
      intersect(
        mdl$getNodeNames(stochOnly = TRUE, includeData = FALSE),
        unique(colnames(as.matrix(base_ml[[1]])))
      )
    )
    unique(present)
  } else {
    unique(c(node1, node2))
  }

  union_nodes <- union_nodes[!is.na(union_nodes) & nzchar(union_nodes)]
  union_nodes <- unique(union_nodes)
  if (length(union_nodes) > block_max)
    union_nodes <- utils::head(union_nodes, block_max)

  if (length(union_nodes) < 2L) {
    message(sprintf(
      "[info] Skipping block phase: need at least 2 valid nodes, got %d (%s).",
      length(union_nodes),
      paste(union_nodes, collapse = ", ")
    ))
    return(list(
      status  = "completed",
      mode    = "surgical_nbot2_singleton_no_block",
      baseline = rb,
      targets  = pick_nodes,
      steps    = steps
    ))
  }

  .barrier_before_block("scalar_to_block_singletons")
  on.exit(.restore_builddir(), add = TRUE)

  # optional PD propCov from baseline
  propCov <- NULL
  M_full <- do.call(rbind, lapply(base_ml, function(m) {
    X <- as.matrix(m)
    keep <- intersect(colnames(X), union_nodes)
    X[, keep, drop = FALSE]
  }))
  if (!is.null(M_full) && is.matrix(M_full) && ncol(M_full) >= 2L) {
    S <- try(stats::cov(M_full, use = "pairwise.complete.obs"), silent = TRUE)
    if (!inherits(S,"try-error")) {
      pc <- try(.sop_make_propCov_PD(S), silent = TRUE)
      if (!inherits(pc,"try-error")) propCov <- pc
    }
  }

  for (samp in unique(as.character(strict_block_seq))) {
    bB    <- .fresh_build(build_fn, monitors = monitors, thin = thin)
    confB <- bB$conf
    for (n in union_nodes) try(confB$removeSamplers(n), silent = TRUE)
    .add_block(confB, union_nodes, samp, has_hmc, Sc = propCov)

    itB  <- .sanitize_iters(pilot_niter, pilot_burnin)
    resB <- .compile_and_run(
      paste0(
        "nbot2_block_",
        paste(gsub("[^A-Za-z0-9_]+","_", union_nodes), collapse = "_"),
        "_", samp
      ),
      bB, confB, itB$niter, itB$nburnin
    )
    metB <- .metrics_for(resB$dg, union_nodes)
    .print_step(
      "Block plan on union (singletons)",
      paste(union_nodes, collapse = " + "),
      union_nodes, samp,
      resB$out$runtime_s, metB, base_dx, rb$runtime_s
    )

    pdir <- file.path(
      out_dir,
      sprintf(
        "block_union_singletons_%s",
        paste(gsub("[^A-Za-z0-9_]+","_", union_nodes), collapse = "_")
      )
    )
    if (!dir.exists(pdir))
      dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    .plot_rhat_bar(resB$dg, nodes = union_nodes,
                   out_file = file.path(pdir, "rhat_bar.png"))
    .plot_traces(resB$ml, nodes = union_nodes,
                 out_file_prefix = file.path(pdir, "trace_"))

    steps <- c(steps, list(list(
      level   = "singletons-block",
      nodes   = union_nodes,
      sampler = samp,
      res     = resB,
      dir     = pdir
    )))

    if (!.ask_yes_no_strict("Proceed to the next BLOCK sampler for the union?"))
      break
  }

  return(list(
    status  = "completed",
    mode    = "surgical_nbot2_singleton",
    baseline = rb,
    targets  = pick_nodes,
    steps    = steps
  ))
}

#' Family-based sampler strategy: full HMC if allowed, else surgical on bottlenecks
#'
#' If the model is differentiable and \code{try_hmc = TRUE}, this runs a full
#' HMC/NUTS configuration (via \code{configure_hmc_safely()}), executes the baseline,
#' prints/plots diagnostics, and returns. Otherwise it switches to a
#' \strong{surgical strategy}: it ranks families by median efficiency, extracts
#' stochastic bottleneck nodes, and iteratively applies samplers on 1 or more nodes
#' (\code{nbot}) in the order:
#' \itemize{
#'   \item scalar: \code{NUTS} -> \code{slice} -> \code{RW}
#'   \item block (when \code{nbot >= 2}): \code{NUTS_block} -> \code{AF_slice} -> \code{RW_block}
#' }
#' After each assignment, a short MCMC is run, a side-by-side comparison versus the
#' baseline is printed (runtime, AE = ESS/iter, CE = ESS/s, Rhat), clean plots are saved
#' (Rhat bars + trace plots for the touched nodes), and the user is prompted to continue
#' unless \code{ask = FALSE}.
#'
#' @param build_fn Function; returns \code{list(model, cmodel?, monitors?)} used by the runners.
#' @param monitors Character or \code{NULL}; root names to monitor (\code{NULL} = auto).
#' @param try_hmc Logical; if \code{TRUE} and structure allows, run full HMC path.
#' @param nchains Integer; number of chains.
#' @param pilot_niter Integer; iterations for pilot runs.
#' @param pilot_burnin Integer; burn-in for pilot runs.
#' @param thin Integer; thinning factor.
#' @param out_dir Character; output directory for diagnostics/plots.
#' @param nbot Integer (\eqn{\ge} 1); number of bottleneck nodes to act on simultaneously.
#' @param strict_scalar_seq Character; sampler order for scalar mode.
#' @param strict_block_seq Character; sampler order for block mode (\code{nbot >= 2}).
#' @param force_families Character or \code{NULL}; families to force.
#' @param force_nodes List or \code{NULL}; per-family forced node vectors.
#' @param force_union Character or \code{NULL}; families to union for block stage.
#' @param ask Logical; interactive confirmation after each step.
#' @param ask_before_hmc Logical; ask before attempting full-model HMC.
#' @param block_max Integer; maximum block size.
#' @param slice_control,rw_control,rwblock_control,af_slice_control Lists; sampler controls.
#' @param slice_max_contractions Integer; informational safety for AF_slice.
#'
#' @return A list describing baseline, steps, configurations, diagnostics, and plot paths.
#'
#' @export
test_strategy_family       <- function(build_fn,
                                        monitors            = NULL,   # optional, just passed through
                                        try_hmc             = TRUE,   # only used for full-model path; surgical ignores
                                        nchains             = 3L,
                                        pilot_niter         = 10000L,
                                        pilot_burnin        = 3000L,
                                        thin                = 2L,
                                        out_dir             = "outputs/diagnostics_family",
                                        nbot                = 1L,
                                        # strict sequences (user can override; order strictly enforced)
                                        strict_scalar_seq   = c("NUTS","slice","RW"),
                                        strict_block_seq    = c("NUTS_block","AF_slice","RW_block"),
                                        # forcing
                                        force_families      = NULL,   # e.g. c("logit_theta","N")
                                        force_nodes         = NULL,   # e.g. list(logit_theta=c("logit_theta[1]",...))
                                        force_union         = NULL,   # e.g. c("logit_theta","N")
                                        # interaction
                                        ask                 = TRUE,
                                        ask_before_hmc      = TRUE,
                                        # safety caps
                                        block_max           = 20L,
                                        # sampler controls
                                        slice_control       = list(),
                                        rw_control          = list(),
                                        rwblock_control     = list(adaptScaleOnly = TRUE),
                                        af_slice_control    = list(),
                                        slice_max_contractions = 5000L) {

  `%||%` <- function(x, y) if (is.null(x)) y else x
  stopifnot(nbot >= 1L)

  if (is.list(build_fn) && !is.null(build_fn$model)) {
    obj <- build_fn
    build_fn <- function() obj
  }
  stopifnot(is.function(build_fn))

  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("samOptiPro: 'nimble' is required.")
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  say     <- function(...) {
    msg <- try(sprintf(...), silent = TRUE)
    if (!inherits(msg,"try-error")) cat(msg,"\n")
  }
  is_ll   <- function(x) grepl("^logLik(\\[.*\\])?$|log_?lik|logdens|lpdf",
                               x, perl = TRUE, ignore.case = TRUE)
  root_of <- function(x) sub("\\[.*", "", x)

  ## ---------- (A) sanitize niter/burnin (copied from test_strategy) ----------
  .sanitize_iters <- local({
    warned <- FALSE
    function(niter, nburnin) {
      niter   <- as.integer(niter)
      nburnin <- as.integer(nburnin)
      if (is.na(niter)   || niter   < 1L) niter   <- 1L
      if (is.na(nburnin) || nburnin < 0L) nburnin <- 0L
      if (nburnin >= niter) {
        new_nburnin <- max(0L, niter - 1L)
        if (!warned) {
          message(sprintf("[info] Adjusting nburnin (%d) to %d because nburnin must be < niter (%d).",
                          nburnin, new_nburnin, niter))
          warned <<- TRUE
        }
        nburnin <- new_nburnin
      }
      list(niter = niter, nburnin = nburnin)
    }
  })

  ## ---------- (B) strict yes/no helper (copied from test_strategy) ----------
  .ask_yes_no_strict <- function(prompt) {
    if (!isTRUE(ask) || !interactive()) return(TRUE)
    while (TRUE) {
      cat(paste0(prompt, " (yes/no): "))
      utils::flush.console()
      ans <- try(readline(), silent = TRUE)
      if (inherits(ans, "try-error")) return(TRUE)  # non-interactive fallback = proceed
      ans <- tolower(trimws(ans))
      if (ans %in% c("y","yes","o","oui")) return(TRUE)
      if (ans %in% c("n","no","non","0"))  return(FALSE)
      cat("Please answer 'yes' or 'no'.\n")
    }
  }

  # --------- ensure unsampled nodes get slice ----------
  .ensure_unsampled <- function(conf) {
    uns <- try(conf$getUnsampledNodes(), silent = TRUE)
    if (!inherits(uns, "try-error") && length(uns)) {
      uns <- uns[!is_ll(uns)]
      for (u in uns) conf$addSampler(u, type = "slice")
    }
  }

  ## ---------- compile + run (aligned with test_strategy) ----------
  .compile_and_run <- function(step_tag, build_obj, conf, niter, nburnin) {
    it <- .sanitize_iters(niter, nburnin)
    niter   <- it$niter
    nburnin <- it$nburnin
    attempts <- 0L
    last_err <- NULL

    repeat {
      attempts <- attempts + 1L
      try(nimble::clearCompiled(), silent = TRUE)
      gc()

      res <- try({
        .ensure_unsampled(conf)
        cmcmc <- .compile_mcmc_with_build(conf, build_obj, reset = TRUE, show = FALSE)
        out   <- .run_and_collect(
          cmcmc,
          niter   = niter,
          nburnin = nburnin,
          thin    = thin,
          nchains = nchains
        )
        ml <- as_mcmc_list_sop(out$samples, out$samples2,
                               drop_loglik = FALSE, thin = thin)

        ## IMPORTANT: use same diagnostic function as test_strategy
        dg <- compute_diag_from_mcmc(ml, runtime_s = out$runtime_s)

        if (!is.null(cmcmc) && is.list(cmcmc) && "unloadDLL" %in% names(cmcmc)) {
          try(cmcmc$unloadDLL(), silent = TRUE)
        }
        list(out = out, ml = ml, dg = dg)
      }, silent = TRUE)

      if (!inherits(res, "try-error")) return(res)

      last_err <- tryCatch(
        conditionMessage(attr(res, "condition")),
        error = function(e) as.character(res)
      )
      if (attempts == 1L) {
        message(sprintf("[info] First attempt failed at %s; retrying cleanly...", step_tag))
      } else {
        message(sprintf("[retry %d @ %s] %s", attempts, step_tag, last_err))
      }

      # fresh rebuild for next attempt
      build_obj <- .fresh_build(build_fn, monitors = monitors, thin = thin)
      conf      <- build_obj$conf

      if (attempts >= 3L) {
        stop(sprintf("Compilation failed after %d attempts at step '%s': %s",
                     attempts, step_tag, last_err))
      }
    }
  }

  # --------- sampler assigners (family-level) ----------
  .add_scalar_family <- function(conf, nodes, type, has_hmc) {
    if (identical(type, "NUTS")) {
      if (isTRUE(has_hmc)) {
        for (n in nodes) conf$addSampler(target = n, type = "NUTS")
      } else {
        cat("[Info] nimbleHMC not available -> falling back to slice (scalar).\n")
        for (n in nodes) conf$addSampler(n, "slice", slice_control)
      }
      return(invisible())
    }
    if (identical(type, "slice")) {
      for (n in nodes) conf$addSampler(n, "slice", slice_control)
      return(invisible())
    }
    if (identical(type, "RW")) {
      for (n in nodes) conf$addSampler(n, "RW", rw_control)
      return(invisible())
    }
    # fallback
    for (n in nodes) conf$addSampler(n, "slice", slice_control)
  }

  .add_block_family <- function(conf, nodes, type, has_hmc, Sc = NULL) {
    nodes <- unique(nodes)
    if (length(nodes) < 2L) {
      # graceful fallback if block is degenerate
      if (type == "NUTS_block") {
        .add_scalar_family(conf, nodes, "NUTS", has_hmc)
      } else if (type == "AF_slice") {
        for (n in nodes) conf$addSampler(n, "AF_slice", af_slice_control)
      } else if (type == "RW_block") {
        for (n in nodes) conf$addSampler(n, "RW", rw_control)
      } else {
        for (n in nodes) conf$addSampler(n, "slice", slice_control)
      }
      return(invisible())
    }

    if (identical(type, "NUTS_block")) {
      if (isTRUE(has_hmc)) {
        conf$addSampler(target = nodes, type = "NUTS")
      } else {
        cat("[Info] nimbleHMC not available -> falling back to AF_slice (block).\n")
        conf$addSampler(target = nodes, type = "AF_slice", control = af_slice_control)
      }
      return(invisible())
    }

    if (identical(type, "AF_slice")) {
      ok <- TRUE
      tryCatch({
        conf$addSampler(target = nodes, type = "AF_slice", control = af_slice_control)
      }, error = function(e) ok <<- FALSE)
      if (!ok) {
        for (n in nodes) conf$addSampler(n, "slice", slice_control)
      }
      return(invisible())
    }

    if (identical(type, "RW_block")) {
      ctrl <- rwblock_control
      if (!is.null(Sc)) ctrl$propCov <- Sc
      conf$addSampler(nodes, "RW_block", ctrl)
      return(invisible())
    }

    # fallback
    for (n in nodes) conf$addSampler(n, "slice", slice_control)
  }

  ## ---------- metrics (very proche de .metrics_for) ----------
  .fam_metrics <- function(dg, nodes) {
    nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
    if (is.null(dg) || !nrow(dg)) {
      return(list(AE = NA_real_, CE = NA_real_, Rhat = NA_real_))
    }
    if (!"target" %in% names(dg)) {
      if ("Parameter" %in% names(dg)) {
        dg$target <- as.character(dg$Parameter)
      } else {
        return(list(AE = NA_real_, CE = NA_real_, Rhat = NA_real_))
      }
    }
    keep <- dg$target %in% nodes & !is_ll(dg$target)
    if (!any(keep)) {
      return(list(AE = NA_real_, CE = NA_real_, Rhat = NA_real_))
    }
    list(
      AE   = stats::median(dg$AE_ESS_per_it[keep], na.rm = TRUE),
      CE   = stats::median(dg$ESS_per_sec[keep],   na.rm = TRUE),
      Rhat = if (any(keep)) suppressWarnings(max(dg$Rhat[keep], na.rm = TRUE)) else NA_real_
    )
  }

  .print_step <- function(title, fam_label, nodes, sampler, runtime, met, base_dx, base_rt) {
    nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
    say("--- %s ---", title)
    if (!is.null(fam_label)) say("Family: %s", fam_label)
    say("Nodes: %s", paste(nodes, collapse = ", "))
    say("Sampler: %s", sampler)
    say("Runtime_s: %.2f (baseline: %.2f)", runtime %||% NA_real_, base_rt %||% NA_real_)
    say("AE median (ESS/iter): %.3g (baseline: %.3g)",
        met$AE %||% NA_real_,
        stats::median(base_dx$AE_ESS_per_it, na.rm = TRUE))
    say("CE median (ESS/s):    %.3g (baseline: %.3g)",
        met$CE %||% NA_real_,
        stats::median(base_dx$ESS_per_sec,   na.rm = TRUE))
    say("Rhat max:             %.3g (baseline max: %.3g)",
        met$Rhat %||% NA_real_,
        suppressWarnings(max(base_dx$Rhat, na.rm = TRUE)))
  }

  # ---------- 0) Build + baseline ----------
  bld0 <- .fresh_build(build_fn, monitors = monitors, thin = thin)
  mdl  <- bld0$model

  it0 <- .sanitize_iters(pilot_niter, pilot_burnin)
  rb  <- run_baseline_config(
    build_fn = build_fn,
    niter    = it0$niter,
    nburnin  = it0$nburnin,
    thin     = thin,
    monitors = monitors,
    nchains  = nchains
  )

  base_ml <- as_mcmc_list_sop(rb$samples, rb$samples2,
                              drop_loglik = FALSE, thin = thin)
  base_dg <- compute_diag_from_mcmc(base_ml, runtime_s = rb$runtime_s)
  base_dx <- base_dg[!is_ll(base_dg$target), , drop = FALSE]

  say("Baseline runtime_s: %.2f s", rb$runtime_s %||% NA_real_)
  say("Baseline median AE(ESS/iter): %.3g",
      stats::median(base_dx$AE_ESS_per_it, na.rm = TRUE))
  say("Baseline median CE(ESS/s):    %.3g",
      stats::median(base_dx$ESS_per_sec,   na.rm = TRUE))

  # ---------- 1) Full-model HMC (optional, same spirit as test_strategy) ----------
  dg_struct <- try(diagnose_model_structure(mdl), silent = TRUE)
  suppressWarnings(has_hmc <- requireNamespace("nimbleHMC", quietly = TRUE))
  deriv_ok <- .sop_supports_derivs(mdl)
  blockers <- character(0)
  if (!inherits(dg_struct, "try-error") && !is.null(dg_struct)) {
    if (isTRUE(dg_struct$has_truncation))   blockers <- c(blockers, "truncation")
    if (isTRUE(dg_struct$has_simplex))      blockers <- c(blockers, "simplex-constraint")
    if (isTRUE(dg_struct$has_non_diff_fun)) blockers <- c(blockers, "non-diff-function")
  }
  nuts_ok_global <- isTRUE(try_hmc) && has_hmc && deriv_ok && (length(blockers) == 0L)

  if (isTRUE(nuts_ok_global) && isTRUE(ask) && isTRUE(ask_before_hmc) && isTRUE(interactive())) {

    proceed_hmc <- .ask_yes_no_strict(
      sprintf(
        "Baseline ready. Runtime=%.2fs; median AE=%.3g; median CE=%.3g\nProceed with full-model HMC/NUTS?",
        rb$runtime_s %||% NA_real_,
        median(base_dx$AE_ESS_per_it, na.rm = TRUE),
        median(base_dx$ESS_per_sec,   na.rm = TRUE)
      )
    )

    if (proceed_hmc) {

      ## --- NEW: use run_hmc_all_nodes() instead of configure_hmc_safely() ---
      hmc_try <- try(
        run_hmc_all_nodes(
          build_fn = build_fn,
          niter    = it0$niter,
          nburnin  = it0$nburnin,
          thin     = thin,
          monitors = monitors,
          nchains  = nchains,
          opts     = samOptiPro_options()  # ou opts par défaut si tu en as un autre
        ),
        silent = TRUE
      )

      if (inherits(hmc_try, "try-error")) {
        cat("[Warn] Full-model HMC (run_hmc_all_nodes) threw an error -> continuing with surgical family plan.\n")
      } else if (!is.list(hmc_try) ||
                 !all(c("samples","samples2","runtime_s") %in% names(hmc_try))) {
        cat("[Warn] run_hmc_all_nodes returned an unexpected structure -> continuing with surgical family plan.\n")
      } else {

        ## Rebuild mcmc.list and diagnostics with the SAME pipeline as baseline
        hmc_ml <- as_mcmc_list_sop(
          hmc_try$samples,
          hmc_try$samples2,
          drop_loglik = FALSE,
          thin = thin
        )

        hmc_dg <- compute_diag_from_mcmc(
          hmc_ml,
          runtime_s = hmc_try$runtime_s
        )

        cat(sprintf("HMC runtime_s: %.3f\n", hmc_try$runtime_s %||% NA_real_))
        cat(sprintf("HMC median AE: %.3g ; CE: %.3g ; max Rhat: %.3g\n",
                    median(hmc_dg$AE_ESS_per_it, na.rm = TRUE),
                    median(hmc_dg$ESS_per_sec,   na.rm = TRUE),
                    if (all(is.na(hmc_dg$Rhat))) NA_real_ else max(hmc_dg$Rhat, na.rm = TRUE)))

        ## On retourne le HMC full ici, comme avant
        return(list(
          mode     = "HMC_full",
          baseline = list(
            runtime_s = rb$runtime_s,
            samples   = base_ml,
            diag_tbl  = base_dg
          ),
          hmc      = list(
            run      = hmc_try,
            mcmc     = hmc_ml,
            diag_tbl = hmc_dg
          ),
          messages = "Full HMC (run_hmc_all_nodes) completed."
        ))
      }

    } else {
      cat("User declined full-model HMC/NUTS. Switching to surgical family plan.\n")
    }
  }


  # ---------- 2) Families selection (respecting forcing) ----------
  stoch_nodes <- mdl$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  stoch_nodes <- stoch_nodes[!is_ll(stoch_nodes)]

  if (!"Family" %in% names(base_dx)) {
    base_dx$Family <- root_of(base_dx$target)
  }

  fam_med <- stats::aggregate(ESS_per_sec ~ Family, data = base_dx, median)
  fam_med <- fam_med[order(fam_med$ESS_per_sec, decreasing = FALSE), , drop = FALSE]

  pick_fams <- character(0)

  if (!is.null(force_families) && length(force_families)) {
    pf      <- unique(force_families)
    present <- pf[pf %in% unique(base_dx$Family)]
    pick_fams <- utils::head(present, nbot)
    if (!length(pick_fams)) {
      warning("force_families provided but none matched diagnostics; using automatic selection.")
    }
  }

  if (!length(pick_fams)) {
    for (fam in fam_med$Family) {
      fam_nodes_all <- intersect(
        stoch_nodes[root_of(stoch_nodes) == fam],
        unique(base_dx$target)
      )
      if (length(fam_nodes_all)) pick_fams <- unique(c(pick_fams, fam))
      if (length(pick_fams) >= nbot) break
    }
    pick_fams <- utils::head(pick_fams, nbot)
  }

  if (!length(pick_fams)) stop("No stochastic bottleneck families found.")

  say("Selected bottleneck families (nbot=%d): %s",
      nbot, paste(pick_fams, collapse = ", "))

  get_family_nodes <- function(fam) {
    allf <- intersect(
      stoch_nodes[root_of(stoch_nodes) == fam],
      unique(base_dx$target)
    )
    forced <- force_nodes[[fam]] %||% NULL
    if (!is.null(forced)) {
      forced <- intersect(forced, allf)
      if (!length(forced)) {
        warning(sprintf("force_nodes for '%s' not found in diagnostics; using detected nodes.", fam))
      }
      return(if (length(forced)) forced else allf)
    }
    allf
  }

  steps <- list()

  # ========================== CASE nbot = 1 ==========================
  if (nbot == 1L) {
    fam1       <- pick_fams[1]
    fam1_nodes <- get_family_nodes(fam1)

    for (samp in unique(as.character(strict_scalar_seq))) {
      bS    <- .fresh_build(build_fn, monitors = monitors, thin = thin)
      confS <- bS$conf
      for (n in fam1_nodes) {
        try(confS$removeSamplers(n), silent = TRUE)
      }
      .add_scalar_family(confS, fam1_nodes, samp, has_hmc)

      resS <- .compile_and_run(
        step_tag  = paste0("nbot1_", fam1, "_", samp),
        build_obj = bS,
        conf      = confS,
        niter     = pilot_niter,
        nburnin   = pilot_burnin
      )
      metS <- .fam_metrics(resS$dg, fam1_nodes)
      .print_step(
        "Scalar plan on family",
        fam1, fam1_nodes, samp,
        resS$out$runtime_s, metS, base_dx, rb$runtime_s
      )

      pdir <- file.path(
        out_dir,
        sprintf("scalar_family_%s_%s",
                gsub("[^A-Za-z0-9_]", "_", fam1), samp)
      )
      if (!dir.exists(pdir)) {
        dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
      }
      .plot_rhat_bar(resS$dg, nodes = fam1_nodes,
                     out_file = file.path(pdir, "rhat_bar.png"))
      .plot_traces(resS$ml, nodes = fam1_nodes,
                   out_file_prefix = file.path(pdir, "trace_"))

      steps <- c(steps, list(list(
        level   = "family-scalar",
        family  = fam1,
        nodes   = fam1_nodes,
        sampler = samp,
        res     = resS,
        dir     = pdir
      )))

      cont <- .ask_yes_no_strict(
        sprintf("Proceed to the next sampler for family '%s'?", fam1)
      )
      if (!isTRUE(cont)) break
    }

    return(list(
      status   = "completed",
      mode     = "surgical_nbot1",
      baseline = rb,
      families = pick_fams,
      steps    = steps
    ))
  }

  # ========================== CASE nbot >= 2 ==========================
  # Step 1 -- Family 1 scalar strict
  fam1       <- pick_fams[1]
  fam1_nodes <- get_family_nodes(fam1)

  for (samp in unique(as.character(strict_scalar_seq))) {
    bS    <- .fresh_build(build_fn, monitors = monitors, thin = thin)
    confS <- bS$conf
    for (n in fam1_nodes) {
      try(confS$removeSamplers(n), silent = TRUE)
    }
    .add_scalar_family(confS, fam1_nodes, samp, has_hmc)

    resS <- .compile_and_run(
      step_tag  = paste0("nbot2_scalar_", fam1, "_", samp),
      build_obj = bS,
      conf      = confS,
      niter     = pilot_niter,
      nburnin   = pilot_burnin
    )
    metS <- .fam_metrics(resS$dg, fam1_nodes)
    .print_step(
      "Scalar plan on family",
      fam1, fam1_nodes, samp,
      resS$out$runtime_s, metS, base_dx, rb$runtime_s
    )

    pdir <- file.path(
      out_dir,
      sprintf("scalar_family_%s_%s",
              gsub("[^A-Za-z0-9_]", "_", fam1), samp)
    )
    if (!dir.exists(pdir)) {
      dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    }
    .plot_rhat_bar(resS$dg, nodes = fam1_nodes,
                   out_file = file.path(pdir, "rhat_bar.png"))
    .plot_traces(resS$ml, nodes = fam1_nodes,
                 out_file_prefix = file.path(pdir, "trace_"))

    steps <- c(steps, list(list(
      level   = "family-scalar",
      family  = fam1,
      nodes   = fam1_nodes,
      sampler = samp,
      res     = resS,
      dir     = pdir
    )))

    cont <- .ask_yes_no_strict(
      sprintf("Proceed to the next sampler for family '%s'?", fam1)
    )
    if (!isTRUE(cont)) break
  }

  # Step 2 -- Block strict on union of families
  fam2       <- pick_fams[2]
  fam2_nodes <- get_family_nodes(fam2)

  if (!is.null(force_union) && length(force_union) >= 2L) {
    union_fams  <- unique(force_union)
    union_nodes <- unique(unlist(lapply(union_fams, get_family_nodes)))
    fam_label   <- paste(union_fams, collapse = " + ")
  } else {
    union_nodes <- unique(c(fam1_nodes, fam2_nodes))
    fam_label   <- paste(fam1, "+", fam2)
  }

  if (length(union_nodes) > block_max) {
    union_nodes <- utils::head(union_nodes, block_max)
  }

  # optional PD propCov from baseline
  propCov <- NULL
  if (length(union_nodes) >= 2L) {
    M_full <- do.call(rbind, lapply(base_ml, function(m) {
      X <- as.matrix(m)
      keep <- intersect(colnames(X), union_nodes)
      X[, keep, drop = FALSE]
    }))
    if (!is.null(M_full) && is.matrix(M_full) && ncol(M_full) >= 2L) {
      S <- try(stats::cov(M_full, use = "pairwise.complete.obs"), silent = TRUE)
      if (!inherits(S, "try-error")) {
        pc <- try(.sop_make_propCov_PD(S), silent = TRUE)
        if (!inherits(pc, "try-error")) propCov <- pc
      }
    }
  }

  for (samp in unique(as.character(strict_block_seq))) {
    bB    <- .fresh_build(build_fn, monitors = monitors, thin = thin)
    confB <- bB$conf
    for (n in union_nodes) {
      try(confB$removeSamplers(n), silent = TRUE)
    }
    .add_block_family(confB, union_nodes, samp, has_hmc, Sc = propCov)

    resB <- .compile_and_run(
      step_tag  = paste0("nbot2_block_",
                         gsub("[^A-Za-z0-9_]+", "_", fam_label), "_", samp),
      build_obj = bB,
      conf      = confB,
      niter     = pilot_niter,
      nburnin   = pilot_burnin
    )
    metB <- .fam_metrics(resB$dg, union_nodes)
    .print_step(
      "Block plan on families union",
      fam_label, union_nodes, samp,
      resB$out$runtime_s, metB, base_dx, rb$runtime_s
    )

    pdir <- file.path(
      out_dir,
      sprintf("block_union_%s_%s",
              gsub("[^A-Za-z0-9_]", "_", fam_label), samp)
    )
    if (!dir.exists(pdir)) {
      dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    }
    .plot_rhat_bar(resB$dg, nodes = union_nodes,
                   out_file = file.path(pdir, "rhat_bar.png"))
    .plot_traces(resB$ml, nodes = union_nodes,
                 out_file_prefix = file.path(pdir, "trace_"))

    steps <- c(steps, list(list(
      level    = "families-block",
      families = fam_label,
      nodes    = union_nodes,
      sampler  = samp,
      res      = resB,
      dir      = pdir
    )))

    cont <- .ask_yes_no_strict(
      sprintf("Proceed to the next BLOCK sampler for union '%s'?", fam_label)
    )
    if (!isTRUE(cont)) break
  }

  return(list(
    status   = "completed",
    mode     = if (nbot == 1L) "surgical_nbot1" else "surgical_nbot2",
    baseline = rb,
    families = pick_fams,
    steps    = steps
  ))
}


#' Unified MCMC diagnostics (vectorised when available)
#'
#' Internal helper to compute per-target MCMC diagnostics from a
#' \code{coda::mcmc.list} object. When the vectorised implementation
#' \code{compute_diag_from_mcmc_vect()} is available in the search path,
#' it is used by default. Otherwise, the legacy
#' \code{compute_diag_from_mcmc()} implementation is used as a
#' transparent fallback.
#'
#' This function is intended to keep the rest of the pipeline agnostic
#' to the exact diagnostic engine, while always returning a single
#' tidy data frame with one row per target parameter.
#'
#' @param ml A \code{coda::mcmc.list} object (one element per chain).
#' @param runtime_s Optional numeric scalar; total runtime in seconds
#'   for the MCMC run, used to compute computational efficiency
#'   metrics (e.g. ESS per second). Can be \code{NA_real_}.
#' @param ... Additional arguments forwarded to the underlying
#'   diagnostic function (if supported).
#'
#' @return A data frame with at least the following columns:
#'   \itemize{
#'     \item \code{target}: character; parameter/node name.
#'     \item \code{ESS}: numeric; effective sample size (per parameter).
#'     \item \code{AE_ESS_per_it}: numeric; algorithmic efficiency
#'       (ESS per iteration).
#'     \item \code{ESS_per_sec}: numeric; computational efficiency
#'       (ESS per second), when \code{runtime_s} is known.
#'     \item \code{Rhat}: numeric; convergence diagnostic (if available).
#'   }
#'
#' Exact column names may depend on the underlying diagnostic engine,
#' but all downstream tools in \pkg{samOptiPro} assume at least these
#' four fields when present.
#'
#' @keywords internal
#' @noRd
.compute_diag_auto <- function(ml, runtime_s = NA_real_, ...) {
  # Prefer the vectorised engine when available
  fn_vect <- get0("compute_diag_from_mcmc_vect", mode = "function", inherits = TRUE)
  fn_legacy <- get0("compute_diag_from_mcmc", mode = "function", inherits = TRUE)

  if (is.function(fn_vect)) {
    return(fn_vect(ml, runtime_s = runtime_s, ...))
  }

  if (is.function(fn_legacy)) {
    return(fn_legacy(ml, runtime_s = runtime_s, ...))
  }

  stop("samOptiPro: no diagnostic function found. ",
       "Define 'compute_diag_from_mcmc_vect()' or 'compute_diag_from_mcmc()'.")
}

#' Scalable diagnostics from MCMC samples (classic layout, block-wise)
#'
#' Efficiently computes convergence and efficiency diagnostics from large
#' \code{mcmc.list} objects using a block-wise, vectorized strategy, while
#' preserving the classic output layout of \code{\link{compute_diag_from_mcmc}}.
#' This is intended as a drop-in, scalable replacement for
#' \code{compute_diag_from_mcmc} in very high-dimensional hierarchical models.
#'
#' @details
#' The function is designed for models with tens of thousands of monitored
#' parameters (e.g. large Bayesian stock assessment or life-cycle models such as
#' WGNAS, GEREM, Scorff LCM), where standard \pkg{coda}-based loops become
#' prohibitively slow or memory bound.
#'
#' Internally, all chains are first truncated to the shortest common length
#' and transformed into matrices with a common set of parameters. Diagnostics
#' (ESS and \eqn{\hat{R}}) are then computed in column blocks to control memory
#' use. The effective sample size used in the output,
#' \code{ESS}, corresponds to a conservative "worst chain" ESS (minimum across
#' chains), and the classic quantities are derived as:
#'
#' \itemize{
#'   \item \code{AE_ESS_per_it = ESS / n_iter},
#'   \item \code{ESS_per_sec   = ESS / runtime_s},
#'   \item \code{time_s_per_ESS = runtime_s / ESS}.
#' }
#'
#' This ensures strict backward compatibility with consumers of
#' \code{compute_diag_from_mcmc}, while benefiting from the scalable
#' implementation of \code{compute_diag_from_mcmc_vect}.
#'
#' @param samples An \code{mcmc.list}, \code{mcmc}, or numeric matrix of MCMC
#'   samples. All chains must share (at least) one common parameter name.
#' @param runtime_s Numeric scalar; total runtime in seconds for the MCMC
#'   run(s) corresponding to \code{samples}. It is assumed to be the total wall
#'   clock time for the full set of chains.
#' @param compute_rhat Character, one of \code{"none"}, \code{"classic"},
#'   \code{"split"}, or \code{"both"} to control the computation of
#'   Gelman–Rubin diagnostics. In the classic layout, only a single
#'   \code{Rhat} column is returned; when both variants are requested, a
#'   priority rule is applied internally (e.g. split-\eqn{\hat{R}} preferred
#'   when available).
#' @param ess_for Character, one of \code{"both"}, \code{"worst"},
#'   or \code{"total"}. The classic layout uses the "worst-chain" ESS for the
#'   \code{ESS} and derived AE/CE columns; the total ESS is used internally
#'   only when available.
#' @param ignore_patterns Optional character vector of regular expressions used
#'   to remove parameters by name before computing diagnostics
#'   (e.g. \code{"^lifted_"}, \code{"^logProb_"}).
#' @param cols_by Integer; number of columns per processing block. Larger values
#'   reduce overhead but increase memory usage. Can be tuned automatically via
#'   \code{target_block_ram_gb}.
#' @param warn_mem_gb Numeric; approximate memory threshold (GiB) above which a
#'   warning is issued based on the size of the chains.
#' @param step_timeout_s Numeric; per-step time limit in seconds for each
#'   diagnostic phase (ESS, total ESS, R-hat). Use \code{Inf} to disable.
#' @param runtime_is_total Logical; if \code{FALSE}, \code{runtime_s} is assumed
#'   to be a per-chain runtime and is multiplied by the number of chains before
#'   computing \code{ESS_per_sec} and \code{time_s_per_ESS}.
#' @param use_posterior Character; \code{"never"} or \code{"if_available"} to
#'   optionally leverage the \pkg{posterior} package for faster and
#'   rank-normalized ESS/\eqn{\hat{R}}. When unavailable, the function falls
#'   back to \pkg{coda}-style computations.
#' @param target_block_ram_gb Numeric; if non-\code{NA}, automatically derives
#'   \code{cols_by} to target this approximate RAM usage (GiB) per processing
#'   block, using a conservative safety factor.
#' @param verbose Logical; if \code{TRUE}, print progress messages and basic
#'   memory/blocking information.
#'
#' @return
#' A \code{data.frame} with one row per parameter and the following columns:
#' \describe{
#'   \item{target}{Parameter name (column name in the MCMC object).}
#'   \item{ESS}{Effective sample size used for efficiency summaries
#'              (worst-chain ESS).}
#'   \item{AE_ESS_per_it}{Algorithmic efficiency, defined as
#'                        \code{ESS / n_iter}.}
#'   \item{ESS_per_sec}{Computational efficiency, defined as
#'                      \code{ESS / runtime_s}.}
#'   \item{time_s_per_ESS}{Seconds per effective sample, defined as
#'                         \code{runtime_s / ESS}.}
#'   \item{Rhat}{Gelman–Rubin convergence diagnostic, if requested via
#'               \code{compute_rhat}; otherwise \code{NA}.}
#'   \item{Family}{Top-level node family (string before the first \code{"["}
#'                 in \code{target}).}
#' }
#'
#' @note
#' This function is intended as a scalable, vectorized backend for
#' \code{compute_diag_from_mcmc}, preserving its column layout so that existing
#' downstream code (plots, summaries, bottleneck detectors) continues to work
#' without modification, even for MCMC outputs with \eqn{\ge 90{,}000}
#' parameters and up to \eqn{10} chains.
#'
#' @seealso
#' \code{\link{compute_diag_from_mcmc}},
#' \code{\link[coda]{effectiveSize}},
#' \code{\link[posterior]{ess_bulk}},
#' \code{\link[posterior]{rhat}},
#' Vehtari et al. (2021), \emph{Bayesian Analysis} 16(2):667–718.
#'
#' @examples
#' \dontrun{
#' diag_alt <- compute_diag_from_mcmc_vect_alt(
#'   samples             = my_mcmc,
#'   runtime_s           = 5400,
#'   compute_rhat        = "both",
#'   ess_for             = "both",
#'   target_block_ram_gb = 2
#' )
#' head(diag_alt)
#' }
#'
#' @export
compute_diag_from_mcmc_alt <- function(samples, runtime_s) {
  stopifnot(is.numeric(runtime_s), length(runtime_s) == 1L, is.finite(runtime_s))

  ## ------------ Normalize samples ------------
  if (inherits(samples, "mcmc"))
    samples <- coda::mcmc.list(samples)
  if (is.matrix(samples))
    samples <- coda::mcmc.list(coda::mcmc(samples))
  if (!inherits(samples, "mcmc.list"))
    stop("'samples' must be mcmc.list, mcmc, or matrix.")

  m <- length(samples)
  if (m < 1L) stop("Empty mcmc.list.")

  mats <- lapply(samples, \(s) as.matrix(s))
  params <- Reduce(intersect, lapply(mats, colnames))
  if (!length(params)) stop("No common parameters across chains.")
  params <- as.character(params)

  ## remove lifted_ / logProb_
  rm_pat <- grepl("^lifted_|^logProb_", params)
  params <- params[!rm_pat]
  mats   <- lapply(mats, \(M) M[, params, drop = FALSE])

  ## align chain lengths
  n_min <- min(vapply(mats, nrow, integer(1)))
  mats  <- lapply(mats, \(M) M[seq_len(n_min), , drop = FALSE])
  it <- n_min; P <- length(params)

  ## ------------ ESS per chain (vectorized) ------------
  ESS_mat <- matrix(NA_real_, nrow = P, ncol = m)
  for (k in seq_len(m)) {
    ESS_mat[, k] <- tryCatch(
      as.numeric(coda::effectiveSize(coda::mcmc(mats[[k]]))),
      error = \(e) rep(NA_real_, P)
    )
  }
  ESS <- suppressWarnings(apply(ESS_mat, 1L, min, na.rm = TRUE))
  ESS[!is.finite(ESS)] <- NA_real_

  ## ------------ Rhat (vectorized Gelman–Rubin) ------------
  Mbig <- do.call(cbind, mats)

  chain_means <- matrix(NA_real_, m, P)
  chain_vars  <- matrix(NA_real_, m, P)

  for (k in seq_len(m)) {
    Mk <- Mbig[, ((k - 1) * P + 1):(k * P), drop = FALSE]
    mk <- colMeans(Mk)
    vk <- (colSums(Mk * Mk) - it * mk^2) / max(1, it - 1)
    chain_means[k, ] <- mk
    chain_vars[k, ]  <- pmax(vk, 0)
  }

  W <- colMeans(chain_vars)
  mu_bar <- colMeans(chain_means)
  B  <- it * pmax((colSums(chain_means^2) - m * mu_bar^2) / max(1, m - 1), 0)
  var_hat <- ((it - 1) / it) * W + B / it

  Rhat <- rep(NA_real_, P)
  good <- W > 0 & is.finite(W) & is.finite(var_hat)
  Rhat[good] <- sqrt(pmax(var_hat[good] / W[good], 0))

  ## ------------ Final statistics------------
  AE  <- ESS / it
  ESSps <- if (runtime_s > 0) ESS / runtime_s else NA_real_
  tper  <- runtime_s / pmax(ESS, .Machine$double.eps)

  out <- data.frame(
    target         = params,
    ESS            = ESS,
    AE_ESS_per_it  = AE,
    ESS_per_sec    = ESSps,
    time_s_per_ESS = tper,
    Rhat           = Rhat,
    stringsAsFactors = FALSE
  )
  out$Family <- sub("\\[.*", "", out$target)

  out <- out[is.finite(out$ESS) & out$ESS > 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}
