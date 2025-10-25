#' Convert an object to `mcmc.list`
#'
#' Converts an `mcmc.list` / `mcmc` / `data.frame` / `matrix` into a `coda`-compatible `mcmc.list`.
#'
#' @param x An object of class `mcmc.list`, `mcmc`, `data.frame`, or `matrix`.
#' @return An object of class `mcmc.list`.
#' @export
#' @keywords internal
#' @importFrom coda mcmc mcmc.list
as_mcmc_list <- function(x){
  if (inherits(x, 'mcmc.list')) return(x)
  if (inherits(x, 'mcmc'))      return(coda::mcmc.list(x))
  x_df <- try(as.data.frame(x), silent = TRUE)
  if (inherits(x_df, "try-error"))
    stop("as_mcmc_list: unable to convert `samples` to data.frame.")
  coda::mcmc.list(coda::mcmc(x_df))
}

#' ESS, Rhat, AE, CE and run-level summary
#'
#' Computes, per parameter: ESS, \eqn{\hat{R}}, **AE = ESS / (iters * chains)**,
#' **CE = ESS / seconds**, and returns a global summary (including total ESS and ESS/s).
#'
#' @param samples `mcmc.list` / `mcmc` / `data.frame` compatible with `coda`.
#' @param runtime_s Total wall time in seconds (numeric).
#' @param rhat_thresh \eqn{\hat{R}} threshold used for the "OK" proportion (default 1.01).
#' @return A list with `summary` (tibble) and `per_param` (tibble).
#' @export
#' @keywords internal
#' @importFrom stats median
assess_performance <- function(samples, runtime_s, rhat_thresh = 1.01){
  ml <- as_mcmc_list(samples)
  n_chains <- length(ml)
  n_iter   <- nrow(ml[[1]])
  stopifnot(is.finite(runtime_s), runtime_s >= 0)

  # ESS per parameter (combined chains)
  ess <- coda::effectiveSize(ml)
  ess <- ess[is.finite(ess)]
  if (length(ess) == 0)
    stop("assess_performance: ESS not found (no finite values).")

  # Rhat per parameter (if >= 2 chains)
  rhat_vec <- rep(NA_real_, length(ess)); names(rhat_vec) <- names(ess)
  if (n_chains >= 2) {
    gd <- coda::gelman.diag(ml, autoburnin = FALSE, transform = FALSE, multivariate = FALSE)
    rhat_in <- try(as.numeric(gd$psrf[, 1]), silent = TRUE)
    if (!inherits(rhat_in, "try-error")) {
      names(rhat_in) <- rownames(gd$psrf)
      nm <- intersect(names(ess), names(rhat_in))
      rhat_vec[nm] <- rhat_in[nm]
    }
  }

  # AE = ESS / (retained iterations * chains) ; CE = ESS / seconds
  N_post <- as.numeric(n_iter) * as.numeric(n_chains)
  AE <- ess / max(1, N_post)
  CE <- ess / max(1e-9, runtime_s)

  per_param <- tibble::tibble(
    parameter = names(ess),
    ESS       = as.numeric(ess),
    Rhat      = as.numeric(rhat_vec),
    AE        = as.numeric(AE),
    CE        = as.numeric(CE)
  )

  ESS_total  <- sum(per_param$ESS, na.rm = TRUE)
  ESS_per_s  <- ESS_total / max(1e-9, runtime_s)
  AE_mean    <- mean(per_param$AE, na.rm = TRUE)
  AE_median  <- stats::median(per_param$AE, na.rm = TRUE)
  CE_mean    <- mean(per_param$CE, na.rm = TRUE)
  CE_median  <- stats::median(per_param$CE, na.rm = TRUE)
  n_params   <- nrow(per_param)
  rhat_ok    <- mean(per_param$Rhat < rhat_thresh, na.rm = TRUE)

  summary <- tibble::tibble(
    runtime_s    = runtime_s,
    n_chains     = n_chains,
    n_iter       = n_iter,
    n_params     = n_params,
    ESS_total    = ESS_total,
    ESS_per_s    = ESS_per_s,
    AE_mean      = AE_mean,
    AE_median    = AE_median,
    CE_mean      = CE_mean,
    CE_median    = CE_median,
    prop_rhat_ok = if (n_chains >= 2) rhat_ok else NA_real_
  )

  list(summary = summary, per_param = per_param)
}

#' Identify per-parameter bottlenecks (low ESS/s, low ESS/post-draw, long time-to-target)
#'
#' Ranks parameters by:
#'  - **CE** (ESS/s): low is worse,
#'  - **AE** (ESS/post-draw): low is worse,
#'  - **slow_node_time** (seconds to reach `ess_threshold` at current CE): high is worse.
#'
#' Degenerate parameters (non-finite or non-positive ESS, AE, or CE) are listed
#' in `degenerate` and excluded from rankings.
#'
#' If you pass `sampler_params`, results are restricted to those parameters
#' (useful to exclude deterministic monitors and keep stochastic nodes with samplers only).
#'
#' @param samples `mcmc.list`/`mcmc`/`matrix`/`data.frame`.
#' @param runtime_s numeric(1) wall time in seconds.
#' @param ess_threshold numeric(1) target ESS per parameter (default 1000).
#' @param sampler_params character() optional vector of parameter names to keep (stochastic nodes with samplers).
#' @param rhat_threshold numeric(1) kept for API symmetry (not used for ranking).
#' @param ess_per_s_min numeric(1) optional CE threshold; flags params below it (0 = inactive).
#' @return list(type="ok" or "degenerate_only",
#'              details=list(ce=..., ae=..., time=..., degenerate=...),
#'              per_param=..., summary=..., top3=data.frame)
#' @export
#' @keywords internal
#' @importFrom stats median
`%||%` <- function(a,b) if (is.null(a)) b else a
.is_ignored <- function(x, patterns) {
  if (length(patterns) == 0) return(rep(FALSE, length(x)))
  Reduce("|", lapply(patterns, function(p) grepl(p, x)))
}
.family_of <- function(x) sub("\\[.*$", "", x)

# ---- core metrics (unchanged) -----------------------------------------------
# assumes you already have as_mcmc_list() and assess_performance() defined as before

# ---- internal: derive sampler-attached parameters ---------------------------
#'@importFrom stats median
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


# ===============================
# Per-parameter bottlenecks (strict sampler-only)
# ===============================
#' @export
#' @keywords internal
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

  ap <- assess_performance(samples, runtime_s, rhat_thresh = rhat_threshold)
  pp <- as.data.frame(ap$per_param)

  if (is.null(sampler_params) || !length(sampler_params)) {
    sampler_params <- .derive_sampler_params_auto(
      samples       = samples,
      model         = model,
      mcmc_conf     = mcmc_conf,
      ignore_patterns = ignore_patterns,
      auto_configure  = auto_configure
    )
  }
  if (isTRUE(strict_sampler_only)) {
    if (!length(sampler_params)) {
      stop("identify_bottlenecks: strict_sampler_only=TRUE but no sampler-attached parameters could be derived. ",
           "Provide `model` or `mcmc_conf`, or set strict_sampler_only=FALSE (not recommended).")
    }
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  } else if (length(sampler_params)) {
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  }

  if (!nrow(pp)) {
    stop("identify_bottlenecks: no parameters remain after sampler-only filtering. ",
         "Check that your samples contain the monitored stochastic nodes.")
  }

  pp$slow_node_time <- ifelse(pp$CE > 0, ess_threshold / pp$CE, Inf)
  pp$meet_target    <- is.finite(pp$slow_node_time) & (pp$slow_node_time <= runtime_s)
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

  top_k    <- min(20L, nrow(valid))
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

  # column order CE -> AE -> slow_node_time
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
  if (!is.null(top3)) rownames(top3) <- NULL

  list(
    type      = "ok",
    details   = list(ce = worst_ce, ae = worst_ae, time = worst_time, degenerate = degenerate),
    per_param = valid,
    summary   = summary,
    top3      = top3
  )
}
#' Histograms for ESS / CE / AE / Rhat (ggplot2)
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

#' Merge `res$samples` and `res$samples2` into a single `mcmc.list`
#'
#' @export
#' @importFrom coda as.mcmc mcmc mcmc.list
#' @keywords internal
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

# ============================
# Harmonised family-level plots
# ============================

#' Family-level ESS bar plot (harmonised style)
#'
#' Builds a bar plot of median ESS by family with a dashed horizontal line at the 5% quantile.
#' The style mirrors the example you provided (fill, labels, themes, rotations).
#'
#' @param samples mcmc.list/mcmc/data.frame/matrix
#' @param runtime_s numeric(1) total runtime in seconds (used to compute CE/AE if needed elsewhere)
#' @param sampler_params optional character() to keep only sampler-attached parameters
#' @return a `ggplot` object
#' @export
#' @importFrom stats quantile median
# ===============================
# Harmonised family plots (sampler-only)
# ===============================
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

#' Family-level Rhat bar plot (harmonised style)
#'
#' Plots median(Rhat) by family on the y-axis transformed as median_Rhat - 1,
#' so that the standard convergence threshold 1.05 appears at 0.05 on the axis.
#'
#' @param samples mcmc.list/mcmc/data.frame/matrix
#' @param runtime_s numeric(1) total runtime in seconds
#' @param sampler_params optional character() to keep only sampler-attached parameters
#' @param rhat_threshold numeric(1) convergence threshold (default 1.05)
#' @return a `ggplot` object (or `NULL` if no finite Rhat values)
#' @export
#' @importFrom stats median
plot_family_rhat_bar <- function(samples, runtime_s,
                                 sampler_params = NULL,
                                 model = NULL,
                                 mcmc_conf = NULL,
                                 samplers_df = NULL,
                                 ignore_patterns = c("^lifted_","^logProb_"),
                                 strict_sampler_only = TRUE,
                                 rhat_threshold = 1.05) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_family_rhat_bar: requires ggplot2.")
  ap <- assess_performance(samples, runtime_s)
  pp <- as.data.frame(ap$per_param)

  if (all(!is.finite(pp$Rhat))) return(NULL)

  if (is.null(sampler_params) || !length(sampler_params)) {
    sampler_params <- derive_sampler_params(samples, model, mcmc_conf, samplers_df,
                                            include_data = FALSE,
                                            ignore_patterns = ignore_patterns)
  }
  if (strict_sampler_only) {
    if (!length(sampler_params)) {
      stop("plot_family_rhat_bar: strict_sampler_only=TRUE but no sampler-attached parameters could be derived. ",
           "Pass sampler_params or provide model/mcmc_conf/samplers_df.")
    }
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  } else if (length(sampler_params)) {
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  }

  family_of <- function(x) sub("\\[.*$", "", x)
  pp$Family <- family_of(pp$parameter)

  rhat_family <- stats::aggregate(list(median_Rhat = pp$Rhat), by = list(Family = pp$Family),
                           FUN = function(z) stats::median(z, na.rm = TRUE))
  rhat_family$y <- pmax(0, rhat_family$median_Rhat - 1)

  ggplot2::ggplot(rhat_family, ggplot2::aes(x = stats::reorder(Family, median_Rhat), y = y)) +
    ggplot2::geom_bar(stat = "identity", fill = "orange", color = "black") +
    ggplot2::geom_hline(yintercept = rhat_threshold - 1, linetype = "dashed", color = "red", size = 1) +
    ggplot2::xlab("Nodes (families)") +
    ggplot2::ylab("Median Gelman-Rubin Rhat") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::scale_y_continuous(limits = c(0, rhat_threshold - 1), labels = function(y) y + 1)
}
# ============================
# Family-level bottlenecks (medians)
# ============================

#' Identify bottlenecks by parameter stats::family(medians within families)
#'
#' Parameters are grouped by stats::family(prefix before the first `[`), then family-level
#' **medians** are computed:
#'   - AE_med  = median(AE)            (low = worse),
#'   - CE_med  = median(CE)            (low = worse, CE = ESS/s),
#'   - ESS_med = median(ESS),
#'   - Rhat_med= median(Rhat, na.rm=TRUE).
#'
#' Derived:
#'   - slow_node_time = ess_threshold / CE_med  (seconds to target; high = worse),
#'   - meet_target    = slow_node_time <= runtime_s.
#'
#' Families with degenerate metrics (non-finite or non-positive) are listed in `degenerate`
#' and excluded from rankings.
#'
#' If you pass `sampler_params`, only parameters belonging to those names are used to form families.
#'
#' @param samples mcmc.list/mcmc/matrix/data.frame.
#' @param runtime_s numeric(1) wall time in seconds.
#' @param ess_threshold numeric(1) target ESS per stats::family(default 1000).
#' @param sampler_params character() optional vector of parameter names to keep.
#' @param rhat_threshold numeric(1) kept for API symmetry (not used in ranks).
#' @param ess_per_s_min numeric(1) optional CE threshold to flag (0 = inactive).
#'
#' @return list(type="ok" or "degenerate_only",
#'              details=list(ce=..., ae=..., time=..., degenerate=...),
#'              per_family=..., summary=..., top3=data.frame)
#'identify_bottlenecks_family
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

  # 0) base metrics
  ap <- assess_performance(samples, runtime_s, rhat_thresh = rhat_threshold)
  pp <- as.data.frame(ap$per_param)

  # A) auto-derive sampler params if not provided
  if (is.null(sampler_params) || !length(sampler_params)) {
    sampler_params <- .derive_sampler_params_auto(
      samples       = samples,
      model         = model,
      mcmc_conf     = mcmc_conf,
      ignore_patterns = ignore_patterns,
      auto_configure  = auto_configure
    )
  }

  # Enforce sampler-only if requested
  if (isTRUE(strict_sampler_only)) {
    if (!length(sampler_params)) {
      stop("identify_bottlenecks_family: strict_sampler_only=TRUE but no sampler-attached parameters could be derived. ",
           "Provide `model` or `mcmc_conf`, or set strict_sampler_only=FALSE (not recommended).")
    }
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  } else if (length(sampler_params)) {
    pp <- pp[pp$parameter %in% sampler_params, , drop = FALSE]
  }

  if (!nrow(pp)) {
    stop("identify_bottlenecks_family: no parameters remain after sampler-only filtering. ",
         "Check that your samples contain the monitored stochastic nodes.")
  }

  # 1) add family
  pp$family <- .family_of(pp$parameter)

  # 2) family medians
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

  # 3) time-to-target
  per_family$slow_node_time <- ifelse(per_family$CE_med > 0, ess_threshold / per_family$CE_med, Inf)
  per_family$meet_target     <- is.finite(per_family$slow_node_time) & (per_family$slow_node_time <= runtime_s)
  per_family$below_ess_per_s <- if (isTRUE(ess_per_s_min > 0)) per_family$CE_med < ess_per_s_min else FALSE

  # 4) degenerates
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

  # 5) ranks (worse: CE_med low, AE_med low, slow_node_time high)
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

  # column order (CE -> AE -> slow_node_time)
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
  if (!is.null(top3)) rownames(top3) <- NULL

  list(
    type       = "ok",
    details    = list(ce = worst_ce, ae = worst_ae, time = worst_time, degenerate = degenerate),
    per_family = valid,
    summary    = summary,
    top3       = top3
  )
}
