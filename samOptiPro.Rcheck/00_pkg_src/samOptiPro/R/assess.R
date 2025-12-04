
#' Convert an object to mcmc.list
#'
#' Converts an input object (\code{mcmc.list}, \code{mcmc}, \code{data.frame},
#' or \code{matrix}) into a \code{coda}-compatible \code{mcmc.list}.
#'
#' @param x An object of class \code{mcmc.list}, \code{mcmc}, \code{data.frame},
#'   or \code{matrix}.
#'
#' @return An object of class \code{mcmc.list}.
#'
#' @keywords internal
#' @export
#' @importFrom coda mcmc mcmc.list

as_mcmc_list <- function(x){
  if (inherits(x, 'mcmc.list')) return(x)
  if (inherits(x, 'mcmc'))      return(coda::mcmc.list(x))
  x_df <- try(as.data.frame(x), silent = TRUE)
  if (inherits(x_df, "try-error"))
    stop("as_mcmc_list: unable to convert `samples` to data.frame.")
  coda::mcmc.list(coda::mcmc(x_df))
}


#' Assess MCMC performance metrics (ESS, R-hat, AE, CE)
#'
#' Compute per-parameter and global diagnostics from an MCMC sample set.
#' Metrics include Effective Sample Size (ESS), Gelman-Rubin R-hat,
#' Algorithmic Efficiency (AE = ESS / total draws), and Computational
#' Efficiency (CE = ESS / runtime in seconds).
#'
#' @details
#' The function first coerces the input to an \code{mcmc.list} object and then
#' computes the following:
#' \itemize{
#'   \item \strong{ESS}: Effective Sample Size for each parameter.
#'   \item \strong{Rhat}: Gelman-Rubin diagnostic (computed only when two
#'     or more chains are available).
#'   \item \strong{AE}: Algorithmic efficiency, defined as ESS divided by
#'     the total number of retained post-burnin draws.
#'   \item \strong{CE}: Computational efficiency, defined as ESS divided by
#'     the total runtime in seconds.
#' }
#'
#' When the \pkg{posterior} package is installed, its implementation of
#' \code{rhat()} is used; otherwise, diagnostics from \pkg{coda} are applied.
#'
#' @param samples A list of MCMC samples or any object that can be converted
#'   by \code{as_mcmc_list()}.
#' @param runtime_s Numeric scalar. Total runtime of the MCMC run (seconds).
#' @param rhat_thresh Numeric scalar. Threshold used to flag R-hat values
#'   indicating lack of convergence. Default is 1.01.
#'
#' @return
#' A named list containing:
#' \describe{
#'   \item{\code{summary}}{A one-row tibble summarizing global diagnostics.}
#'   \item{\code{per_param}}{A tibble containing ESS, Rhat, AE, and CE
#'     for each parameter.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- run_baseline_config(build_M, niter = 2000, nburnin = 500, thin = 2)
#' perf <- assess_performance(res$samples, runtime_s = res$runtime_s)
#' perf$summary
#' }
#'
#' @seealso coda::effectiveSize, coda::gelman.diag, posterior::rhat
#' @export

assess_performance <- function(samples, runtime_s, rhat_thresh = 1.01) {
  stopifnot(is.finite(runtime_s), runtime_s >= 0)

  # Convert to mcmc.list
  ml <- as_mcmc_list(samples)
  n_chains <- length(ml)
  n_iter_by_chain <- vapply(ml, nrow, integer(1))
  N_post <- sum(n_iter_by_chain)  # total retained draws

  # --- Effective sample size per parameter
  ess <- coda::effectiveSize(ml)
  ess <- ess[is.finite(ess)]
  if (!length(ess))
    stop("assess_performance: no finite ESS could be computed.")

  # --- R-hat per parameter (if ≥ 2 chains)
  rhat_vec <- rep(NA_real_, length(ess))
  names(rhat_vec) <- names(ess)
  if (n_chains >= 2) {
    rhat_in <- NULL
    # Prefer {posterior} if available
    if (requireNamespace("posterior", quietly = TRUE)) {
      dr <- try(posterior::as_draws(ml), silent = TRUE)
      if (!inherits(dr, "try-error")) {
        rr <- try(posterior::rhat(dr), silent = TRUE)
        if (!inherits(rr, "try-error")) {
          rhat_in <- as.numeric(rr)
          names(rhat_in) <- names(rr)
        }
      }
    }
    # Fallback to coda::gelman.diag
    if (is.null(rhat_in)) {
      gd <- try(coda::gelman.diag(ml,
                                  autoburnin = FALSE,
                                  transform = FALSE,
                                  multivariate = FALSE),
                silent = TRUE)
      if (!inherits(gd, "try-error") && !is.null(gd$psrf)) {
        rhat_in <- as.numeric(gd$psrf[, 1])
        names(rhat_in) <- rownames(gd$psrf)
      }
    }
    if (!is.null(rhat_in)) {
      nm <- intersect(names(ess), names(rhat_in))
      if (length(nm)) rhat_vec[nm] <- rhat_in[nm]
    }
  }

  # --- Derived metrics: AE and CE
  AE <- ess / max(1, N_post)
  CE <- ess / max(1e-9, runtime_s)

  per_param <- tibble::tibble(
    parameter = names(ess),
    ESS       = as.numeric(ess),
    Rhat      = as.numeric(rhat_vec),
    AE        = as.numeric(AE),
    CE        = as.numeric(CE)
  )

  # --- Global summary
  ESS_total  <- sum(per_param$ESS, na.rm = TRUE)
  ESS_per_s  <- ESS_total / max(1e-9, runtime_s)
  AE_mean    <- mean(per_param$AE, na.rm = TRUE)
  AE_median  <- stats::median(per_param$AE, na.rm = TRUE)
  CE_mean    <- mean(per_param$CE, na.rm = TRUE)
  CE_median  <- stats::median(per_param$CE, na.rm = TRUE)
  n_params   <- nrow(per_param)
  rhat_ok    <- if (n_chains >= 2) mean(per_param$Rhat < rhat_thresh, na.rm = TRUE) else NA_real_

  summary <- tibble::tibble(
    runtime_s     = runtime_s,
    n_chains      = n_chains,
    n_iter        = if (length(unique(n_iter_by_chain)) == 1L) n_iter_by_chain[1] else NA_integer_,
    n_draws_total = N_post,
    n_params      = n_params,
    ESS_total     = ESS_total,
    ESS_per_s     = ESS_per_s,
    AE_mean       = AE_mean,
    AE_median     = AE_median,
    CE_mean       = CE_mean,
    CE_median     = CE_median,
    prop_rhat_ok  = rhat_ok
  )

  # Optionally sort parameters by CE for quick inspection
  per_param <- per_param[order(per_param$CE, decreasing = TRUE), ]

  list(summary = summary, per_param = per_param)
}

#' Identify per-parameter bottlenecks (low ESS/s, low ESS per draw, long time-to-target)
#'
#' This function ranks parameters according to three metrics:
#' * Computational efficiency (CE), defined as ESS per second (low values indicate bottlenecks).
#' * Algorithmic efficiency (AE), defined as ESS per post-burnin draw (low values indicate bottlenecks).
#' * Time-to-target, defined as the time required to reach a given ESS threshold at the current CE
#'   (high values indicate bottlenecks).
#'
#' Parameters with non-finite or non-positive ESS, AE, or CE are collected in the
#' \code{degenerate} field of the output and excluded from the rankings.
#'
#' When \code{sampler_params} is provided, rankings are restricted to the specified
#' parameter names. This is typically used to exclude deterministic monitors and
#' retain only stochastic nodes that have associated samplers.
#'
#' @param samples An object of class \code{mcmc.list}, \code{mcmc}, \code{matrix},
#'   or \code{data.frame}.
#' @param runtime_s Numeric scalar. Total wall-clock time of the MCMC run (in seconds).
#' @param ess_threshold Numeric scalar. Target ESS per parameter used to compute
#'   time-to-target. Default is 1000.
#' @param sampler_params Optional character vector of parameter names to retain.
#'   If provided, other parameters are ignored.
#' @param rhat_threshold Numeric scalar. Provided for API symmetry but not used
#'   in ranking. Default is NULL.
#' @param ess_per_s_min Numeric scalar. Optional threshold for computational
#'   efficiency. Parameters with CE below this value are flagged. Set to 0
#'   (default) to disable.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{type}}{Character string. Either "ok" or "degenerate_only".}
#'   \item{\code{details}}{List containing CE, AE, time-to-target, and degenerate sets.}
#'   \item{\code{per_param}}{Data frame with metrics for each parameter.}
#'   \item{\code{summary}}{Data frame summarizing the worst parameters.}
#'   \item{\code{top3}}{Data frame containing the three worst parameters according to CE.}
#' }
#'
#' @keywords internal
#' @importFrom stats median

identify_bottlenecks <- function(samples, runtime_s,
                                 ess_threshold = 1000,
                                 sampler_params = NULL,
                                 model = NULL,
                                 mcmc_conf = NULL,
                                 ignore_patterns = c("^lifted_","^logProb_"),
                                 strict_sampler_only = TRUE,
                                 auto_configure = TRUE,
                                 rhat_threshold = 1.01,
                                 ess_per_s_min = 0) {

  stopifnot(is.numeric(runtime_s), length(runtime_s) == 1L, is.finite(runtime_s))

  # --- helper: columns present in samples ---
  .sample_cols <- function(smp) {
    if (is.null(smp)) return(character(0))
    x <- try(as.data.frame(smp[[1]]), silent = TRUE)
    if (inherits(x, "try-error")) {
      x <- try(as.data.frame(smp), silent = TRUE)
      if (inherits(x, "try-error")) return(character(0))
    }
    colnames(x)
  }

  # --- helper: derive sampler-attached targets via getSamplers() ---
  .derive_sampler_params_from_conf <- function(mcmc_conf, samples, ignore_patterns, model, auto_configure) {
    cfg <- mcmc_conf
    if (is.null(cfg) && isTRUE(auto_configure) && !is.null(model)) {
      cfg <- try(nimble::configureMCMC(model), silent = TRUE)
      if (inherits(cfg, "try-error")) cfg <- NULL
    }
    tgts <- character(0)
    if (!is.null(cfg) && is.function(cfg$getSamplers)) {
      sams <- try(cfg$getSamplers(), silent = TRUE)
      if (!inherits(sams, "try-error") && length(sams)) {
        tgts <- unique(unlist(lapply(sams, function(s) s$target), use.names = FALSE))
      }
    }
    # nettoie et restreint aux colonnes présentes
    cols <- .sample_cols(samples)
    if (length(ignore_patterns)) {
      re <- paste(ignore_patterns, collapse = "|")
      tgts <- tgts[!grepl(re, tgts, perl = TRUE)]
    }
    if (length(cols)) tgts <- intersect(tgts, cols)
    unique(tgts)
  }

  # --- performance per parameter ---
  ap <- assess_performance(samples, runtime_s, rhat_thresh = rhat_threshold)
  pp <- as.data.frame(ap$per_param)

  # --- ensure sampler_params come from getSamplers() if not provided ---
  if (is.null(sampler_params) || !length(sampler_params)) {
    sampler_params <- .derive_sampler_params_from_conf(
      mcmc_conf       = mcmc_conf,
      samples         = samples,
      ignore_patterns = ignore_patterns,
      model           = model,
      auto_configure  = auto_configure
    )
  }

  # --- enforce sampler-only filter if requested ---
  if (isTRUE(strict_sampler_only)) {
    if (!length(sampler_params)) {
      stop("identify_bottlenecks: strict_sampler_only=TRUE but no sampler-attached parameters could be derived. ",
           "Provide `mcmc_conf` (preferably) or `model` (with auto_configure=TRUE), ",
           "or set strict_sampler_only=FALSE (not recommended).")
    }
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  } else if (length(sampler_params)) {
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  }

  if (!nrow(pp)) {
    stop("identify_bottlenecks: no parameters remain after sampler-only filtering. ",
         "Check that your samples contain the monitored stochastic nodes.")
  }

  # --- targets and ranks (unchanged) ---
  pp$slow_node_time  <- ifelse(pp$CE > 0, ess_threshold / pp$CE, Inf)
  pp$meet_target     <- is.finite(pp$slow_node_time) & (pp$slow_node_time <= runtime_s)
  pp$below_ess_per_s <- if (isTRUE(ess_per_s_min > 0)) (pp$CE < ess_per_s_min) else FALSE

  is_degen   <- !is.finite(pp$ESS) | (pp$ESS <= 0) | !is.finite(pp$AE) | !is.finite(pp$CE)
  degenerate <- pp[is_degen, , drop = FALSE]
  valid      <- pp[!is_degen, , drop = FALSE]

  CE_rank   <- rank(valid$CE,               ties.method = "min", na.last = "keep")
  AE_rank   <- rank(valid$AE,               ties.method = "min", na.last = "keep")
  TIME_rank <- rank(-valid$slow_node_time,  ties.method = "min", na.last = "keep")

  o_ce   <- order(CE_rank,   valid$CE,              decreasing = FALSE, na.last = NA)
  o_ae   <- order(AE_rank,   valid$AE,              decreasing = FALSE, na.last = NA)
  o_time <- order(TIME_rank, -valid$slow_node_time, decreasing = FALSE, na.last = NA)

  take <- function(ord, k) {
    idx <- ord[seq_len(min(k, length(ord)))]
    if (length(idx) == 0L) return(valid[integer(0), ])
    out <- valid[idx, , drop = FALSE]
    out$CE_rank   <- CE_rank[idx]
    out$AE_rank   <- AE_rank[idx]
    out$TIME_rank <- TIME_rank[idx]
    out
  }

  top_k      <- min(20L, nrow(valid))
  worst_ce   <- take(o_ce,   top_k)
  worst_ae   <- take(o_ae,   top_k)
  worst_time <- take(o_time, top_k)

  summary <- data.frame(
    runtime_s             = runtime_s,
    n_params              = nrow(pp),
    n_degenerate          = nrow(degenerate),
    AE_min                = min(valid$AE, na.rm = TRUE),
    AE_median             = stats::median(valid$AE, na.rm = TRUE),
    AE_mean               = mean(valid$AE, na.rm = TRUE),
    CE_min                = min(valid$CE, na.rm = TRUE),
    CE_median             = stats::median(valid$CE, na.rm = TRUE),
    CE_mean               = mean(valid$CE, na.rm = TRUE),
    slow_node_time_median = stats::median(valid$slow_node_time, na.rm = TRUE),
    prop_meet_target      = mean(valid$meet_target, na.rm = TRUE),
    ess_threshold         = ess_threshold,
    stringsAsFactors = FALSE
  )

  .reorder_cols_param <- function(df) {
    cols <- c("parameter","ESS","CE","AE","Rhat","slow_node_time",
              "meet_target","below_ess_per_s","CE_rank","AE_rank","TIME_rank")
    df[, intersect(cols, names(df)), drop = FALSE]
  }
  worst_ce   <- .reorder_cols_param(worst_ce)
  worst_ae   <- .reorder_cols_param(worst_ae)
  worst_time <- .reorder_cols_param(worst_time)
  degenerate <- .reorder_cols_param(degenerate)
  valid      <- .reorder_cols_param(valid)

  mk_top3 <- function(df, axis_label) {
    if (is.null(df) || !nrow(df)) return(NULL)
    k <- min(3L, nrow(df))
    out <- df[seq_len(k), , drop = FALSE]
    out$axis <- axis_label
    out[, intersect(c("axis","parameter","ESS","CE","AE","slow_node_time",
                      "CE_rank","AE_rank","TIME_rank"), names(out)), drop = FALSE]
  }
  top3_list <- list(
    mk_top3(worst_ce,   "CE"),
    mk_top3(worst_ae,   "AE"),
    mk_top3(worst_time, "time")
  )
  top3_list <- Filter(Negate(is.null), top3_list)
  top3 <- if (length(top3_list)) do.call(rbind, top3_list) else NULL

  if (!is.null(top3) && nrow(top3)) {
    # Remove duplicates (same parameter-axis pairs)
    top3 <- top3[!duplicated(top3[, c("parameter", "axis")]), , drop = FALSE]
    # Ensure consistent ordering by axis, then CE rank
    top3 <- top3[order(top3$axis, top3$CE_rank, na.last = TRUE), , drop = FALSE]
    rownames(top3) <- NULL
  }

  list(
    type      = "ok",
    details   = list(
      ce         = worst_ce,
      ae         = worst_ae,
      time       = worst_time,
      degenerate = degenerate
    ),
    per_param = valid,
    summary   = summary,
    top3      = top3
  )
}

#' Identify bottlenecks at the parameter-family level
#'
#' Group parameters into families and rank them by median efficiency metrics.
#'
#' @description
#' Parameters are grouped into families defined by the prefix before the first
#' \code{"["} in their name (for example, \code{"beta[1]"} and \code{"beta[2]"}
#' belong to the family \code{"beta"}). For each family, the function computes
#' median efficiency metrics and derived quantities that help identify
#' bottlenecks.
#'
#' @details
#' For each family, the following median metrics are computed:
#' \itemize{
#'   \item \code{AE_med}  = median(AE)            (low values are worse),
#'   \item \code{CE_med}  = median(CE)            (low values are worse, CE is ESS per second),
#'   \item \code{ESS_med} = median(ESS),
#'   \item \code{Rhat_med} = median(Rhat, with \code{na.rm = TRUE}).
#' }
#'
#' From these, the following diagnostics are derived:
#' \itemize{
#'   \item \code{slow_node_time} = \code{ess_threshold / CE_med}
#'     (seconds needed to reach the target ESS; higher is worse),
#'   \item \code{meet_target} = logical flag, \code{TRUE} when
#'     \code{slow_node_time <= runtime_s}.
#' }
#'
#' Families with degenerate metrics (non-finite or non-positive ESS, AE,
#' or CE) are reported in the \code{degenerate} component and excluded from
#' the ranking.
#'
#' When \code{sampler_params} is provided, only parameters whose names are
#' included in \code{sampler_params} are used to form families (typically
#' stochastic nodes that are actually sampled).
#'
#' @param samples An object containing MCMC samples; typically an object of
#'   class \code{mcmc.list}, \code{mcmc}, \code{matrix}, or \code{data.frame}.
#' @param runtime_s Numeric scalar. Wall-clock runtime of the MCMC run in seconds.
#' @param ess_threshold Numeric scalar. Target ESS per family (default is 1000).
#' @param sampler_params Optional character vector of parameter names to keep
#'   when defining families. Parameters not in this vector are ignored.
#' @param rhat_threshold Numeric scalar kept for API symmetry (not used
#'   in the ranking).
#' @param ess_per_s_min Numeric scalar. Optional CE threshold (ESS per second)
#'   used to flag families below this value. Use 0 to deactivate.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{\code{type}}{Character string, either \code{"ok"} or
#'     \code{"degenerate_only"}.}
#'   \item{\code{details}}{List with components \code{ce}, \code{ae},
#'     \code{time}, and \code{degenerate} summarising the diagnostics.}
#'   \item{\code{per_family}}{Data frame (or tibble) of metrics by family.}
#'   \item{\code{summary}}{Single-row data frame (or tibble) with global
#'     summaries across families.}
#'   \item{\code{top3}}{Data frame containing the three worst families
#'     according to the main ranking criterion.}
#' }
#'
#' @export
#' @importFrom stats median
identify_bottlenecks_family <- function(samples, runtime_s,
                                        ess_threshold = 1000,
                                        sampler_params = NULL,
                                        model = NULL,
                                        mcmc_conf = NULL,
                                        ignore_patterns = c("^lifted_","^logProb_"),
                                        strict_sampler_only = TRUE,
                                        auto_configure = TRUE,
                                        rhat_threshold = 1.01,
                                        ess_per_s_min = 0) {
  stopifnot(is.numeric(runtime_s), length(runtime_s) == 1L, is.finite(runtime_s))

  # --- base metrics
  ap <- assess_performance(samples, runtime_s, rhat_thresh = rhat_threshold)
  pp <- as.data.frame(ap$per_param)

  # --- helper: columns present in samples ---
  .sample_cols <- function(smp) {
    if (is.null(smp)) return(character(0))
    x <- try(as.data.frame(smp[[1]]), silent = TRUE)
    if (inherits(x, "try-error")) {
      x <- try(as.data.frame(smp), silent = TRUE)
      if (inherits(x, "try-error")) return(character(0))
    }
    colnames(x)
  }

  # --- helper: derive sampler-attached targets via getSamplers() ---
  .derive_sampler_params_from_conf <- function(mcmc_conf, samples, ignore_patterns, model, auto_configure) {
    cfg <- mcmc_conf
    if (is.null(cfg) && isTRUE(auto_configure) && !is.null(model)) {
      cfg <- try(nimble::configureMCMC(model), silent = TRUE)
      if (inherits(cfg, "try-error")) cfg <- NULL
    }
    tgts <- character(0)
    if (!is.null(cfg) && is.function(cfg$getSamplers)) {
      sams <- try(cfg$getSamplers(), silent = TRUE)
      if (!inherits(sams, "try-error") && length(sams)) {
        tgts <- unique(unlist(lapply(sams, function(s) s$target), use.names = FALSE))
      }
    }
    cols <- .sample_cols(samples)
    if (length(ignore_patterns)) {
      re <- paste(ignore_patterns, collapse = "|")
      tgts <- tgts[!grepl(re, tgts, perl = TRUE)]
    }
    if (length(cols)) tgts <- intersect(tgts, cols)
    unique(tgts)
  }

  # --- derive sampler params if not provided ---
  if (is.null(sampler_params) || !length(sampler_params)) {
    sampler_params <- .derive_sampler_params_from_conf(
      mcmc_conf       = mcmc_conf,
      samples         = samples,
      ignore_patterns = ignore_patterns,
      model           = model,
      auto_configure  = auto_configure
    )
  }

  # --- enforce sampler-only filter ---
  if (isTRUE(strict_sampler_only)) {
    if (!length(sampler_params)) {
      stop("identify_bottlenecks_family: strict_sampler_only=TRUE but no sampler-attached parameters could be derived. ",
           "Provide `mcmc_conf` (recommended) or `model` (with auto_configure=TRUE), ",
           "or set strict_sampler_only=FALSE (not recommended).")
    }
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  } else if (length(sampler_params)) {
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  }

  if (!nrow(pp)) {
    stop("identify_bottlenecks_family: no parameters remain after sampler-only filtering. ",
         "Check that your samples contain the monitored stochastic nodes.")
  }

  # --- 1) add family
  pp$family <- .family_of(pp$parameter)

  # --- 2) family medians
  agg <- function(df) {
    data.frame(
      ESS_med   = stats::median(df$ESS, na.rm = TRUE),
      AE_med    = stats::median(df$AE,  na.rm = TRUE),
      CE_med    = stats::median(df$CE,  na.rm = TRUE),
      Rhat_med  = if (all(is.na(df$Rhat))) NA_real_ else stats::median(df$Rhat, na.rm = TRUE),
      n_members = nrow(df),
      stringsAsFactors = FALSE
    )
  }
  per_family <- do.call(rbind, lapply(split(pp, pp$family), agg))
  per_family$family <- rownames(per_family); rownames(per_family) <- NULL

  # --- 3) time-to-target
  per_family$slow_node_time <- ifelse(per_family$CE_med > 0, ess_threshold / per_family$CE_med, Inf)
  per_family$meet_target     <- is.finite(per_family$slow_node_time) & (per_family$slow_node_time <= runtime_s)
  per_family$below_ess_per_s <- if (isTRUE(ess_per_s_min > 0)) per_family$CE_med < ess_per_s_min else FALSE

  # --- 4) degenerates
  is_degen   <- !is.finite(per_family$ESS_med) | (per_family$ESS_med <= 0) |
    !is.finite(per_family$AE_med)  | !is.finite(per_family$CE_med)
  degenerate <- per_family[is_degen, , drop = FALSE]
  valid      <- per_family[!is_degen, , drop = FALSE]

  if (nrow(valid) == 0L) {
    summary <- data.frame(
      runtime_s             = runtime_s,
      n_families            = nrow(per_family),
      n_degenerate_families = nrow(degenerate),
      AE_med_min            = NA_real_,
      AE_med_median         = NA_real_,
      AE_med_mean           = NA_real_,
      CE_med_min            = NA_real_,
      CE_med_median         = NA_real_,
      CE_med_mean           = NA_real_,
      slow_node_time_median = NA_real_,
      prop_meet_target      = mean(per_family$meet_target, na.rm = TRUE),
      ess_threshold         = ess_threshold,
      stringsAsFactors = FALSE
    )
    return(list(
      type       = "degenerate_only",
      details    = list(ce = valid, ae = valid, time = valid, degenerate = degenerate),
      per_family = per_family,
      summary    = summary,
      top3       = NULL
    ))
  }

  # --- 5) ranks (CE low, AE low, slow_node_time high)
  CE_rank   <- rank(valid$CE_med,           ties.method = "min", na.last = "keep")
  AE_rank   <- rank(valid$AE_med,           ties.method = "min", na.last = "keep")
  TIME_rank <- rank(-valid$slow_node_time,  ties.method = "min", na.last = "keep")

  o_ce   <- order(CE_rank,   valid$CE_med,           decreasing = FALSE, na.last = NA)
  o_ae   <- order(AE_rank,   valid$AE_med,           decreasing = FALSE, na.last = NA)
  o_time <- order(TIME_rank, -valid$slow_node_time,  decreasing = FALSE, na.last = NA)

  take <- function(ord, k) {
    idx <- ord[seq_len(min(k, length(ord)))]
    if (length(idx) == 0L) return(valid[integer(0), ])
    out <- valid[idx, , drop = FALSE]
    out$CE_rank   <- CE_rank[idx]
    out$AE_rank   <- AE_rank[idx]
    out$TIME_rank <- TIME_rank[idx]
    out
  }

  top_k      <- min(20L, nrow(valid))
  worst_ce   <- take(o_ce,   top_k)
  worst_ae   <- take(o_ae,   top_k)
  worst_time <- take(o_time, top_k)

  summary <- data.frame(
    runtime_s             = runtime_s,
    n_families            = nrow(per_family),
    n_degenerate_families = nrow(degenerate),
    AE_med_min            = min(valid$AE_med, na.rm = TRUE),
    AE_med_median         = stats::median(valid$AE_med, na.rm = TRUE),
    AE_med_mean           = mean(valid$AE_med, na.rm = TRUE),
    CE_med_min            = min(valid$CE_med, na.rm = TRUE),
    CE_med_median         = stats::median(valid$CE_med, na.rm = TRUE),
    CE_med_mean           = mean(valid$CE_med, na.rm = TRUE),
    slow_node_time_median = stats::median(valid$slow_node_time, na.rm = TRUE),
    prop_meet_target      = mean(valid$meet_target, na.rm = TRUE),
    ess_threshold         = ess_threshold,
    stringsAsFactors = FALSE
  )

  # --- reorder columns (CE → AE → slow_node_time)
  .reorder_cols_family <- function(df) {
    cols <- c("family","n_members","ESS_med","CE_med","AE_med","Rhat_med",
              "slow_node_time","meet_target","below_ess_per_s",
              "CE_rank","AE_rank","TIME_rank")
    df[, intersect(cols, names(df)), drop = FALSE]
  }
  worst_ce   <- .reorder_cols_family(worst_ce)
  worst_ae   <- .reorder_cols_family(worst_ae)
  worst_time <- .reorder_cols_family(worst_time)
  degenerate <- .reorder_cols_family(degenerate)
  valid      <- .reorder_cols_family(valid)

  # --- top3 (deduplicated)
  mk_top3 <- function(df, axis_label) {
    if (is.null(df) || !nrow(df)) return(NULL)
    k <- min(3L, nrow(df))
    out <- df[seq_len(k), , drop = FALSE]
    out$axis <- axis_label
    out[, intersect(c("axis","family","n_members","CE_med","AE_med","slow_node_time",
                      "CE_rank","AE_rank","TIME_rank"), names(out)), drop = FALSE]
  }
  top3_list <- list(
    mk_top3(worst_ce,   "CE"),
    mk_top3(worst_ae,   "AE"),
    mk_top3(worst_time, "time")
  )
  top3_list <- Filter(Negate(is.null), top3_list)
  top3 <- if (length(top3_list)) do.call(rbind, top3_list) else NULL

  if (!is.null(top3) && nrow(top3)) {
    top3 <- top3[!duplicated(top3[, c("family", "axis")]), , drop = FALSE]
    top3 <- top3[order(top3$axis, top3$CE_rank, na.last = TRUE), , drop = FALSE]
    rownames(top3) <- NULL
  }

  list(
    type       = "ok",
    details    = list(ce = worst_ce, ae = worst_ae, time = worst_time, degenerate = degenerate),
    per_family = valid,
    summary    = summary,
    top3       = top3
  )
}


#' Merge two MCMC sample objects into a single mcmc.list
#'
#' Combine \code{res$samples} and \code{res$samples2} (or any two objects
#' coercible to \code{mcmc} objects) into a single \code{mcmc.list}. Each
#' input is converted using \code{coda::as.mcmc} before being combined.
#'
#' @details
#' This helper is intended for internal use in cases where a model produces
#' two separate sample objects (for example, when running a baseline sampler
#' and an additional HMC sampler). Both inputs are converted to MCMC objects
#' and returned as a single list of chains.
#'
#' @return
#' An object of class \code{mcmc.list} containing all chains from both
#' inputs.
#'
#' @keywords internal
#' @importFrom coda as.mcmc mcmc mcmc.list
#' @export
merge_mcmc_samples <- function(res) {
  to_mlist <- function(x) {
    if (inherits(x, "mcmc.list")) return(x)
    if (inherits(x, "mcmc"))      return(coda::mcmc.list(x))
    if (is.null(x)) return(NULL)
    coda::mcmc.list(coda::as.mcmc(as.data.frame(x)))
  }

  s1 <- to_mlist(res$samples)
  s2 <- to_mlist(res$samples2)

  if (is.null(s1) && is.null(s2)) stop("No samples found (samples/samples2).")
  if (is.null(s2)) return(s1)
  if (is.null(s1)) return(s2)

  stopifnot(length(s1) == length(s2))
  coda::mcmc.list(lapply(seq_along(s1), function(i) {
    m1 <- as.matrix(s1[[i]])
    m2 <- as.matrix(s2[[i]])
    coda::mcmc(
      cbind(m1, m2),
      start = attr(s1[[i]], "mcpar")[1],
      end   = attr(s1[[i]], "mcpar")[2],
      thin  = attr(s1[[i]], "mcpar")[3]
    )
  }))
}


#' Family-level ESS bar plot (harmonised style)
#'
#' Produce a bar plot showing the median Effective Sample Size (ESS)
#' for each parameter family. A horizontal dashed line indicates the
#' 5 percent quantile of ESS across all families.
#'
#' Families are defined by the substring before the first "[" in
#' each parameter name. If \code{sampler_params} is provided, only
#' parameters included in that vector are used when forming families.
#'
#' @param samples An object that can be converted to an \code{mcmc.list}
#'   (e.g., \code{mcmc.list}, \code{mcmc}, \code{data.frame}, or
#'   \code{matrix}).
#' @param runtime_s Numeric scalar giving the total runtime in seconds.
#'   This value is not used directly inside the plot, but is kept for
#'   API consistency.
#' @param sampler_params Optional character vector; if supplied, only
#'   parameters listed here are used when computing family-level ESS.
#'
#' @return A \code{ggplot} object representing median ESS per family.
#'
#' @export
#' @importFrom stats quantile median
plot_family_ess_bar <- function(samples, runtime_s,
                                sampler_params = NULL,
                                model = NULL,
                                mcmc_conf = NULL,
                                samplers_df = NULL,
                                ignore_patterns = c("^lifted_","^logProb_"),
                                strict_sampler_only = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_family_ess_bar: requires ggplot2.")
  ap <- assess_performance(samples, runtime_s)
  pp <- as.data.frame(ap$per_param)

  if (is.null(sampler_params) || !length(sampler_params)) {
    sampler_params <- derive_sampler_params(samples, model, mcmc_conf, samplers_df,
                                            include_data = FALSE,
                                            ignore_patterns = ignore_patterns)
  }
  if (strict_sampler_only) {
    if (!length(sampler_params)) {
      stop("plot_family_ess_bar: strict_sampler_only=TRUE but no sampler-attached parameters could be derived. ",
           "Pass sampler_params or provide model/mcmc_conf/samplers_df.")
    }
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  } else if (length(sampler_params)) {
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  }

  family_of <- function(x) sub("\\[.*$", "", x)
  pp$Family <- family_of(pp$parameter)

  ess_family <- stats::aggregate(list(median_ESS = pp$ESS), by = list(Family = pp$Family),
                          FUN = stats::median, na.rm = TRUE)
  ess_quantile_5 <- stats::quantile(pp$ESS, probs = 0.05, na.rm = TRUE, names = FALSE, type = 7)

  ggplot2::ggplot(ess_family, ggplot2::aes(x = stats::reorder(Family, median_ESS), y = median_ESS)) +
    ggplot2::geom_bar(stat = "identity", fill = "skyblue", color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = round(median_ESS, 1)), vjust = -0.5, size = 3) +
    ggplot2::geom_hline(yintercept = ess_quantile_5, linetype = "dashed", color = "red", size = 1) +
    ggplot2::labs(
      title = "Quantile 5% Effective Sample Size",
      subtitle = paste("5th quantile =", round(ess_quantile_5, 1)),
      x = "Nodes (families)",
      y = "Median Effective Sample Size (ESS)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Histograms for ESS, CE, AE, and Rhat
#'
#' Produce histogram plots (via \code{ggplot2}) for a single metric such as
#' Effective Sample Size (ESS), Computational Efficiency (CE = ESS/s),
#' Algorithmic Efficiency (AE = ESS per draw), or the Gelman-Rubin statistic
#' Rhat (often shifted as \code{Rhat - 1} for visual clarity).
#'
#' The function expects a data frame containing at least one numeric column
#' whose name is supplied via \code{metric}. A histogram is then generated
#' using a harmonised style consistent with other diagnostic plots in the
#' package.
#'
#' @param df A data frame containing at least one numeric column representing
#'   a diagnostic metric.
#' @param metric Character scalar giving the column name in \code{df} to be
#'   plotted.
#' @param bins Integer number of histogram bins (default typically 30).
#' @param vline Optional numeric value at which to draw a vertical reference
#'   line; set to \code{NULL} (default) to omit the line.
#'
#' @return A \code{ggplot} histogram object.
#'
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_histogram labs geom_vline scale_x_continuous
plot_mcmc_histograms <- function(samples, runtime_s, rhat_thresh = 1.01, bins = 30, log_x = TRUE){
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_mcmc_histograms: requires ggplot2.")
  ap   <- assess_performance(samples, runtime_s, rhat_thresh = rhat_thresh)
  diag <- ap$per_param

  ggh <- function(df, xvar, title, vline = NULL, log_x = FALSE){
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = xvar)) +
      ggplot2::geom_histogram(bins = bins, alpha = 0.9) +
      ggplot2::labs(title = title, x = xvar, y = "count")
    if (!is.null(vline)) p <- p + ggplot2::geom_vline(xintercept = vline, linetype = 2)
    if (isTRUE(log_x))   p <- p + ggplot2::scale_x_continuous(trans = "log10")
    p
  }

  plots <- list()
  plots$ESS  <- ggh(diag, "ESS",  "ESS per parameter",      log_x = log_x)
  plots$CE   <- ggh(diag, "CE",   "ESS per second (CE)",    log_x = log_x)
  plots$AE   <- ggh(diag, "AE",   "ESS per post-draw (AE)", log_x = log_x)
  plots$Rhat <- if (any(is.finite(diag$Rhat))) {
    ggh(diag[is.finite(diag$Rhat), , drop = FALSE],
        "Rhat", paste0("Rhat (threshold ", rhat_thresh, ")"),
        vline = rhat_thresh, log_x = FALSE)
  } else NULL
  plots$summary <- ap$summary
  plots
}


#' Derive sampler-attached parameter names (internal)
#'
#' Internal helper used to infer the set of parameters that are effectively
#' attached to samplers, based on the MCMC samples and, optionally, a model
#' object and its MCMC configuration.
#'
#' Parameters whose names match any of the regular expressions in
#' \code{ignore_patterns} (for example technical nodes such as \code{"^lifted_"}
#' or \code{"^logProb_"}) are dropped. When \code{auto_configure} is \code{TRUE},
#' the function may inspect \code{model} and \code{mcmc_conf} to refine the
#' list of sampler-attached targets.
#'
#' @param samples An object containing MCMC samples (for example an
#'   \code{mcmc.list}, \code{mcmc}, \code{data.frame}, or \code{matrix}).
#' @param model Optional model object used to refine the mapping between
#'   samplers and parameter names (for example a nimble model).
#' @param mcmc_conf Optional MCMC configuration object associated with
#'   \code{model}.
#' @param ignore_patterns Character vector of regular expressions; any parameter
#'   whose name matches at least one pattern is excluded from the result.
#' @param auto_configure Logical flag; if \code{TRUE}, additional information
#'   from \code{model} and \code{mcmc_conf} can be used to derive sampler
#'   parameters automatically.
#' @param include_data Logical flag; if \code{FALSE} (default), observed/data
#'   nodes are excluded from the derived parameter set.
#'
#' @return A character vector of parameter names considered as
#'   sampler-attached targets, possibly with attributes used internally.
#'
#' @keywords internal
#' @noRd
#' @importFrom stats median
.derive_sampler_params_auto <- function(samples,
                                        model = NULL,
                                        mcmc_conf = NULL,
                                        ignore_patterns = c("^lifted_","^logProb_"),
                                        auto_configure = TRUE,
                                        include_data = FALSE) {
  # columns actually present in samples
  ml <- as_mcmc_list(samples)
  cols <- colnames(as.data.frame(ml[[1]]))

  # 1) use provided mcmc_conf
  if (!is.null(mcmc_conf)) {
    sl <- try(mcmc_conf$getSamplers(), silent = TRUE)
    if (!inherits(sl, "try-error") && length(sl)) {
      targets <- unique(unlist(lapply(sl, function(s) s$target %||% character(0)), use.names = FALSE))
      targets <- targets[!is.na(targets)]
      targets <- targets[!.is_ignored(targets, ignore_patterns)]
      keep <- intersect(cols, targets)
      if (length(keep)) return(keep)
    }
  }

  # 2) if model provided and auto_configure, try to nimble::configureMCMC(model)
  if (is.null(mcmc_conf) && !is.null(model) && isTRUE(auto_configure)) {
    conf_try <- try(nimble::configureMCMC(model), silent = TRUE)
    if (!inherits(conf_try, "try-error")) {
      sl <- try(conf_try$getSamplers(), silent = TRUE)
      if (!inherits(sl, "try-error") && length(sl)) {
        targets <- unique(unlist(lapply(sl, function(s) s$target %||% character(0)), use.names = FALSE))
        targets <- targets[!is.na(targets)]
        targets <- targets[!.is_ignored(targets, ignore_patterns)]
        keep <- intersect(cols, targets)
        if (length(keep)) return(keep)
      }
    }
  }

  # 3) fallback: stochastic nodes from model
  if (!is.null(model)) {
    stoch <- try(model$getNodeNames(stochOnly = TRUE, includeData = include_data), silent = TRUE)
    if (!inherits(stoch, "try-error") && length(stoch)) {
      stoch <- stoch[!.is_ignored(stoch, ignore_patterns)]
      keep <- intersect(cols, stoch)
      if (length(keep)) return(keep)
    }
  }

  # 4) nothing found
  character(0)
}


  .derive_sampler_params_from_conf <- function(mcmc_conf, samples, ignore_patterns, model, auto_configure) {
    cfg <- mcmc_conf
    if (is.null(cfg) && isTRUE(auto_configure) && !is.null(model)) {
      cfg <- try(nimble::configureMCMC(model), silent = TRUE)
      if (inherits(cfg, "try-error")) cfg <- NULL
    }
    tgts <- character(0)
    if (!is.null(cfg) && is.function(cfg$getSamplers)) {
      sams <- try(cfg$getSamplers(), silent = TRUE)
      if (!inherits(sams, "try-error") && length(sams)) {
        tgts <- unique(unlist(lapply(sams, function(s) s$target), use.names = FALSE))
      }
    }
    # nettoie et restreint aux colonnes présentes
    cols <- .sample_cols(samples)
    if (length(ignore_patterns)) {
      re <- paste(ignore_patterns, collapse = "|")
      tgts <- tgts[!grepl(re, tgts, perl = TRUE)]
    }
    if (length(cols)) tgts <- intersect(tgts, cols)
    unique(tgts)
  }


.family_of <- function(x) sub("\\[.*$", "", x)

#' Derive sampler-attached parameter names (internal)
#'
#' Internal helper used to infer the parameter names that are effectively
#' associated with samplers. The function operates on MCMC samples and may
#' optionally use a model object and its MCMC configuration for refinement.
#'
#' Parameter names matching any pattern in \code{ignore_patterns} (for example
#' \code{"^lifted_"} or \code{"^logProb_"}) are removed. If
#' \code{auto_configure = TRUE}, additional inspection of \code{model} or
#' \code{mcmc_conf} may be performed to determine sampler targets. When
#' \code{include_data = FALSE}, observed/data nodes are excluded.
#'
#' @param samples An object containing MCMC draws (for example an
#'   \code{mcmc.list}, \code{mcmc}, \code{data.frame}, or \code{matrix}).
#' @param model Optional model object used to refine sampler-target detection.
#' @param mcmc_conf Optional MCMC configuration associated with \code{model}.
#' @param ignore_patterns Character vector of regular expressions specifying
#'   parameter-name patterns to exclude.
#' @param auto_configure Logical; if TRUE, the function attempts to derive
#'   sampler-attached parameters using additional model information.
#' @param include_data Logical; if FALSE (default), observed/data nodes are
#'   excluded from the result.
#'
#' @return A character vector of inferred sampler-attached parameter names.
#'
#' @keywords internal
#' @noRd
#' @importFrom stats median
.derive_sampler_params_auto <- function(samples,
                                        model = NULL,
                                        mcmc_conf = NULL,
                                        ignore_patterns = c("^lifted_","^logProb_"),
                                        auto_configure = TRUE,
                                        include_data = FALSE) {
  # columns actually present in samples
  ml <- as_mcmc_list(samples)
  cols <- colnames(as.data.frame(ml[[1]]))

  # 1) use provided mcmc_conf
  if (!is.null(mcmc_conf)) {
    sl <- try(mcmc_conf$getSamplers(), silent = TRUE)
    if (!inherits(sl, "try-error") && length(sl)) {
      targets <- unique(unlist(lapply(sl, function(s) s$target %||% character(0)), use.names = FALSE))
      targets <- targets[!is.na(targets)]
      targets <- targets[!.is_ignored(targets, ignore_patterns)]
      keep <- intersect(cols, targets)
      if (length(keep)) return(keep)
    }
  }

  # 2) if model provided and auto_configure, try to nimble::configureMCMC(model)
  if (is.null(mcmc_conf) && !is.null(model) && isTRUE(auto_configure)) {
    conf_try <- try(nimble::configureMCMC(model), silent = TRUE)
    if (!inherits(conf_try, "try-error")) {
      sl <- try(conf_try$getSamplers(), silent = TRUE)
      if (!inherits(sl, "try-error") && length(sl)) {
        targets <- unique(unlist(lapply(sl, function(s) s$target %||% character(0)), use.names = FALSE))
        targets <- targets[!is.na(targets)]
        targets <- targets[!.is_ignored(targets, ignore_patterns)]
        keep <- intersect(cols, targets)
        if (length(keep)) return(keep)
      }
    }
  }

  # 3) fallback: stochastic nodes from model
  if (!is.null(model)) {
    stoch <- try(model$getNodeNames(stochOnly = TRUE, includeData = include_data), silent = TRUE)
    if (!inherits(stoch, "try-error") && length(stoch)) {
      stoch <- stoch[!.is_ignored(stoch, ignore_patterns)]
      keep <- intersect(cols, stoch)
      if (length(keep)) return(keep)
    }
  }

  # 4) nothing found
  character(0)
}


.is_ignored <- function(x, patterns) {
  if (length(patterns) == 0) return(rep(FALSE, length(x)))
  Reduce("|", lapply(patterns, function(p) grepl(p, x)))
}


  .reorder_cols_family <- function(df) {
    cols <- c("family","n_members","ESS_med","CE_med","AE_med","Rhat_med",
              "slow_node_time","meet_target","below_ess_per_s",
              "CE_rank","AE_rank","TIME_rank")
    df[, intersect(cols, names(df)), drop = FALSE]
  }


  .reorder_cols_param <- function(df) {
    cols <- c("parameter","ESS","CE","AE","Rhat","slow_node_time",
              "meet_target","below_ess_per_s","CE_rank","AE_rank","TIME_rank")
    df[, intersect(cols, names(df)), drop = FALSE]
  }


  .sample_cols <- function(smp) {
    if (is.null(smp)) return(character(0))
    x <- try(as.data.frame(smp[[1]]), silent = TRUE)
    if (inherits(x, "try-error")) {
      x <- try(as.data.frame(smp), silent = TRUE)
      if (inherits(x, "try-error")) return(character(0))
    }
    colnames(x)
  }


  agg <- function(df) {
    data.frame(
      ESS_med   = stats::median(df$ESS, na.rm = TRUE),
      AE_med    = stats::median(df$AE,  na.rm = TRUE),
      CE_med    = stats::median(df$CE,  na.rm = TRUE),
      Rhat_med  = if (all(is.na(df$Rhat))) NA_real_ else stats::median(df$Rhat, na.rm = TRUE),
      n_members = nrow(df),
      stringsAsFactors = FALSE
    )
  }

  family_of <- function(x) sub("\\[.*$", "", x)

  #' Family-level R-hat bar plot (harmonised style)
  #'
  #' Builds a bar plot of median R-hat by family, using the transformation
  #' \code{median_Rhat - 1} on the y-axis. With this transformation, the usual
  #' convergence threshold \code{rhat_threshold = 1.05} appears at \code{0.05}
  #' on the vertical axis.
  #'
  #' Families are defined by the prefix before the first \code{"["} in the
  #' parameter name (for example \code{"beta[1]"} and \code{"beta[2]"} both
  #' map to the family \code{"beta"}).
  #'
  #' @param samples An object containing MCMC draws (for example an
  #'   \code{mcmc.list}, \code{mcmc}, \code{data.frame}, or \code{matrix}).
  #'   It must be compatible with \code{assess_performance()} via
  #'   \code{as_mcmc_list()}.
  #' @param runtime_s Numeric scalar; total runtime in seconds for the MCMC
  #'   run. It is passed to \code{assess_performance()} for consistency.
  #' @param sampler_params Optional character vector of parameter names to keep.
  #'   If provided, only those parameters are used to compute family-level
  #'   median R-hat. If \code{NULL} (default), all parameters in
  #'   \code{assess_performance(samples, runtime_s)$per_param} are used.
  #' @param rhat_threshold Numeric convergence threshold for R-hat (default
  #'   \code{1.05}). A dashed horizontal line is drawn at
  #'   \code{rhat_threshold - 1} on the transformed scale.
  #'
  #' @return
  #' A \pkg{ggplot2} object representing the bar plot, or \code{NULL} if no
  #' finite R-hat values are available.
  #'
  #' @export
  #' @importFrom stats median aggregate reorder
  plot_family_rhat_bar <- function(samples, runtime_s,
                                   sampler_params = NULL,
                                   rhat_threshold = 1.05) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("plot_family_rhat_bar: requires the 'ggplot2' package.")
    }

    ap <- assess_performance(samples, runtime_s)
    pp <- as.data.frame(ap$per_param)

    if (!("Rhat" %in% names(pp))) {
      stop("plot_family_rhat_bar: 'per_param' must contain a 'Rhat' column.")
    }
    if (!any(is.finite(pp$Rhat))) {
      return(NULL)
    }

    # Optional restriction to a subset of parameters
    if (!is.null(sampler_params) && length(sampler_params)) {
      pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
      if (!nrow(pp)) {
        stop("plot_family_rhat_bar: no parameters left after applying sampler_params filter.")
      }
    }

    # Family = prefix before first '['
    family_of <- function(x) sub("\\[.*$", "", x)
    pp$Family <- family_of(pp$parameter)

    # Median R-hat by family
    rhat_family <- stats::aggregate(
      list(median_Rhat = pp$Rhat),
      by   = list(Family = pp$Family),
      FUN  = function(z) stats::median(z, na.rm = TRUE)
    )

    # Transform to median_Rhat - 1 for plotting
    rhat_family$y <- pmax(0, rhat_family$median_Rhat - 1)

    ggplot2::ggplot(
      rhat_family,
      ggplot2::aes(x = stats::reorder(Family, median_Rhat), y = y)
    ) +
      ggplot2::geom_bar(stat = "identity", fill = "orange", color = "black") +
      ggplot2::geom_hline(
        yintercept = rhat_threshold - 1,
        linetype   = "dashed",
        color      = "red",
        linewidth  = 1
      ) +
      ggplot2::xlab("Nodes (families)") +
      ggplot2::ylab("Median Gelman-Rubin R-hat") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 90, hjust = 1, vjust = 0.5
        )
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, rhat_threshold - 1),
        labels = function(y) y + 1
      )
  }

