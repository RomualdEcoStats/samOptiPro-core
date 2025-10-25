# --------------------------------------------------------------
# PLOTS BOTTLENECKS
# --------------------------------------------------------------

# ==================================================================
# samOptiPro -- Diagnostics plotting utilities (English-only)
# ==================================================================

# --------------------------------------------------------------
# PLOT BOTTLENECKS (targets & families): AE/CE, time, R-hat
# --------------------------------------------------------------

#' Plot bottlenecks (targets & families) with AE/CE, total time, and R-hat
#'
#' Creates only the requested figures (toggles). If a toggle is FALSE, the
#' corresponding figure is neither computed nor saved. ECDF and joint
#' scatter plots are intentionally removed.
#'
#' @param diag_tbl Data frame with columns (flexible names supported):
#'   \itemize{
#'     \item \code{target} (if absent, rownames are used),
#'     \item \code{AE} or \code{AE_ESS_per_it} (ESS/iter; small = worse),
#'     \item \code{CE} or \code{ESS_per_sec} or \code{ess_per_s} (ESS/s; small = worse),
#'     \item \code{ESS} or \code{ess} (optional, absolute ESS),
#'     \item \code{time_s} (optional; total seconds per target),
#'     \item \code{Rhat} (optional),
#'     \item \code{Family} (otherwise derived from target prefix before '[').
#'   }
#' @param out_dir Output directory (created if missing).
#' @param top_k Number of bars to show (worst/top).
#' @param rhat_ref Reference R-hat threshold line (default 1.01).
#'
#' @param make_time_targets  Barplot of top targets by total time.
#' @param make_esss_targets  Barplot of worst targets by CE (ESS/s).
#' @param make_esss_families Barplot of worst families by median CE.
#' @param make_time_families Barplot of top families by total time.
#' @param make_rhat_hist_targets Histogram of R-hat (targets).
#' @param make_rhat_worst_targets Barplot of worst targets by R-hat.
#' @param make_rhat_median_families Barplot of median R-hat by family.
#' @param make_hist_ae_families Histogram of family medians (AE).
#' @param make_hist_ce_families Histogram of family medians (CE).
#'
#' @return (Invisibly) a named list of ggplot objects actually created.
#' @export
plot_bottlenecks <- function(diag_tbl,
                             out_dir = "outputs/diagnostics",
                             top_k   = 20L,
                             rhat_ref = 1.05,  # ton template montre 1.05 (-> 0.05 apres -1)
                             # ---- toggles (only TRUE are produced) ----
                             make_time_targets            = TRUE,
                             make_esss_targets            = TRUE,
                             make_esss_families           = TRUE,  # -> AE bar (steelblue)
                             make_time_families           = TRUE,  # -> CE bar (green)
                             make_rhat_hist_targets       = TRUE,  # -> Rhat-1 by Family (bar)
                             make_rhat_worst_targets      = TRUE,  # -> Worst targets by Rhat (bar)
                             make_rhat_median_families    = TRUE,  # (alias of the same family bar)
                             make_hist_ae_families        = FALSE, # (desactives: on suit tes barplots)
                             make_hist_ce_families        = FALSE) {

  if (is.null(diag_tbl) || !NROW(diag_tbl)) return(invisible(NULL))
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---------- utils ----------
  as_num <- function(x) suppressWarnings(as.numeric(x))
  save_fig <- function(p, base, w = 8, h = 5.2, dpi = 180) {
    ggplot2::ggsave(file.path(out_dir, paste0(base, ".pdf")), p, width = w, height = h)
    ggplot2::ggsave(file.path(out_dir, paste0(base, ".png")), p, width = w, height = h, dpi = dpi)
  }
  agg_safe <- function(formula, data, FUN, ...) {
    if (!is.data.frame(data) || !NROW(data)) return(data.frame())
    out <- try(stats::aggregate(formula, data = data, FUN = FUN, ...), silent = TRUE)
    if (inherits(out, "try-error") || !NROW(out)) data.frame() else out
  }
  headn <- function(x, n) utils::head(x, n)

  # ---------- normalize columns ----------
  d <- as.data.frame(diag_tbl, stringsAsFactors = FALSE)
  if (!("target" %in% names(d))) {
    rn <- rownames(d); if (is.null(rn)) stop("plot_bottlenecks: missing 'target' column.")
    d$target <- rn
  }
  if (!("Family" %in% names(d))) d$Family <- sub("\\[.*", "", d$target)

  # AE = ESS/iter
  d$AE <- if ("AE_ESS_per_it" %in% names(d)) as_num(d$AE_ESS_per_it) else if ("AE" %in% names(d)) as_num(d$AE) else NA_real_
  # CE = ESS/s
  d$CE <- if ("ESS_per_sec" %in% names(d)) as_num(d$ESS_per_sec) else if ("ess_per_s" %in% names(d)) as_num(d$ess_per_s) else if ("CE" %in% names(d)) as_num(d$CE) else NA_real_
  # ESS absolu
  d$ESS_abs <- if ("ESS" %in% names(d)) as_num(d$ESS) else if ("ess" %in% names(d)) as_num(d$ess) else NA_real_

  # temps
  d$time_s <- as_num(if ("time_s" %in% names(d)) d$time_s else NA_real_)
  infer_ok <- is.finite(d$ESS_abs) & is.finite(d$CE) & d$CE > 0 & !is.finite(d$time_s)
  d$time_s[infer_ok] <- d$ESS_abs[infer_ok] / d$CE[infer_ok]

  # Rhat
  d$Rhat <- if ("Rhat" %in% names(d)) as_num(d$Rhat) else NA_real_

  # ---------- family-level aggregates ----------
  fam_ce   <- agg_safe(CE ~ Family,     d[is.finite(d$CE), ], median, na.rm = TRUE); names(fam_ce)[2]   <- "CE_median"
  fam_ae   <- agg_safe(AE ~ Family,     d[is.finite(d$AE), ], median, na.rm = TRUE); names(fam_ae)[2]   <- "AE_median"
  fam_time <- agg_safe(time_s ~ Family, d[is.finite(d$time_s), ], sum,    na.rm = TRUE); names(fam_time)[2] <- "time_sum"

  # ---------- metrics barplots per your templates ----------
  res <- list()
  grouped_data <- merge(fam_ae, fam_ce, by = "Family", all = TRUE)
  names(grouped_data)[names(grouped_data) == "AE_median"] <- "MedianAlgorithmicEfficiency"
  names(grouped_data)[names(grouped_data) == "CE_median"] <- "MedianComputationalEfficiency_tot"

  if (isTRUE(make_esss_families) && NROW(grouped_data)) {
    algorithmic_plot_raw <- ggplot2::ggplot(grouped_data,
                                            ggplot2::aes(x = stats::reorder(Family, MedianAlgorithmicEfficiency),
                                                         y = MedianAlgorithmicEfficiency)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::labs(title = "Median Algorithmic Efficiency by Node Family",
                    x = "Node Family", y = "Median Algorithmic Efficiency") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(algorithmic_plot_raw, "bar_family_algorithmic_eff")
    res$bar_family_algorithmic_eff <- algorithmic_plot_raw
  }

  if (isTRUE(make_time_families) && NROW(grouped_data)) {
    computational_plot_raw_tot <- ggplot2::ggplot(grouped_data,
                                                  ggplot2::aes(x = stats::reorder(Family, MedianComputationalEfficiency_tot),
                                                               y = MedianComputationalEfficiency_tot)) +
      ggplot2::geom_bar(stat = "identity", fill = "green") +
      ggplot2::labs(title = "Median Computational Efficiency (Total Time) by Node Family",
                    x = "Node Family", y = "Median Computational Efficiency (Total Time)") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(computational_plot_raw_tot, "bar_family_computational_eff_tot")
    res$bar_family_computational_eff_tot <- computational_plot_raw_tot
  }

  # ---------- targets: CE/time barplots (on garde le tri utile) ----------
  d_ce <- d[is.finite(d$CE), , drop = FALSE]
  d_tm <- d[is.finite(d$time_s), , drop = FALSE]

  if (isTRUE(make_esss_targets) && NROW(d_ce)) {
    d_worst <- headn(d_ce[order(d_ce$CE, -d_ce$time_s), ], min(top_k, NROW(d_ce)))
    p <- ggplot2::ggplot(d_worst,
                         ggplot2::aes(x = stats::reorder(target, CE), y = CE)) +
      ggplot2::geom_bar(stat = "identity", fill = "grey60") +
      ggplot2::labs(title = "Worst Targets by CE", x = "Targets", y = "ESS/s") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "bar_target_CE_worst")
    res$bar_target_CE <- p
  }

  if (isTRUE(make_time_targets) && NROW(d_tm)) {
    d_top <- headn(d_tm[order(-d_tm$time_s, d_tm$CE), ], min(top_k, NROW(d_tm)))
    p <- ggplot2::ggplot(d_top,
                         ggplot2::aes(x = stats::reorder(target, time_s), y = time_s)) +
      ggplot2::geom_bar(stat = "identity", fill = "grey40") +
      ggplot2::labs(title = "Top Targets by Time", x = "Targets", y = "Time (s)") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "bar_target_time_top")
    res$bar_target_time <- p
  }

  # ---------- RHAT (all "hist" now use your family bar template) ----------
  d_r <- d[is.finite(d$Rhat), , drop = FALSE]
  if (NROW(d_r)) {
    fam_rhat <- agg_safe(Rhat ~ Family, data = d_r, FUN = stats::median, na.rm = TRUE)
    names(fam_rhat)[2] <- "median_Rhat"

    if (isTRUE(make_rhat_hist_targets) && NROW(fam_rhat)) {
      p_rhat <- ggplot2::ggplot(fam_rhat,
                                ggplot2::aes(x = stats::reorder(Family, median_Rhat), y = median_Rhat - 1)) +
        ggplot2::geom_bar(stat = "identity", fill = "orange", color = "black") +
        ggplot2::geom_hline(yintercept = (rhat_ref - 1), linetype = "dashed", color = "red", size = 1) +
        ggplot2::xlab("Nodes") +
        ggplot2::ylab("Median of Gelman-Rubin Rhat") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ggplot2::scale_y_continuous(limits = c(0, rhat_ref - 1), labels = function(y) y + 1)
      save_fig(p_rhat, "rhat_family_bars_template")
      res$rhat_family_template <- p_rhat
    }

    if (isTRUE(make_rhat_median_families) && NROW(fam_rhat)) {
      # alias (meme figure, nom different pour compat)
      p_rhat2 <- ggplot2::ggplot(fam_rhat,
                                 ggplot2::aes(x = stats::reorder(Family, median_Rhat), y = median_Rhat - 1)) +
        ggplot2::geom_bar(stat = "identity", fill = "orange", color = "black") +
        ggplot2::geom_hline(yintercept = (rhat_ref - 1), linetype = "dashed", color = "red", size = 1) +
        ggplot2::xlab("Nodes") +
        ggplot2::ylab("Median of Gelman-Rubin Rhat") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ggplot2::scale_y_continuous(limits = c(0, rhat_ref - 1), labels = function(y) y + 1)
      save_fig(p_rhat2, "rhat_median_families_template")
      res$rhat_family_median <- p_rhat2
    }

    if (isTRUE(make_rhat_worst_targets)) {
      ord <- order(-d_r$Rhat, d_r$CE)
      rhat_worst <- headn(d_r[ord, , drop = FALSE], min(top_k, NROW(d_r)))
      p <- ggplot2::ggplot(rhat_worst,
                           ggplot2::aes(x = stats::reorder(target, Rhat), y = Rhat - 1)) +
        ggplot2::geom_bar(stat = "identity", fill = "orange", color = "black") +
        ggplot2::geom_hline(yintercept = (rhat_ref - 1), linetype = "dashed", color = "red", size = 1) +
        ggplot2::xlab("Nodes") +
        ggplot2::ylab("Median of Gelman-Rubin Rhat") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ggplot2::scale_y_continuous(limits = c(0, rhat_ref - 1), labels = function(y) y + 1)
      save_fig(p, "rhat_worst_targets_template")
      res$rhat_worst_targets <- p
    }
  }

  invisible(res)
}
# --------------------------------------------------------------
# PLOT CONVERGENCE CHECKS (R-hat + traces)
# --------------------------------------------------------------

#' Plot convergence diagnostics (R-hat & traces) with family-level R-hat bars
#'
#' ECDF plots are removed. Only requested items are produced and saved.
#'
#' @param samples \code{coda::mcmc.list}.
#' @param out_dir Output directory (PDF & PNG).
#' @param top_k_rhat Number of targets for "worst R-hat" traces.
#' @param top_k_aelow Number of targets for "worst AE (low ESS/s)" traces.
#' @param runtime_s Total runtime in seconds (optional; used to compute CE if needed).
#' @param rhat_ref Reference R-hat threshold (default 1.01).
#' @param make_rhat_hist Logical -- global R-hat histogram.
#' @param make_traces_rhat Logical -- traces/densities for worst R-hat.
#' @param make_traces_ae Logical -- traces/densities for worst AE (low ESS/s).
#' @param make_rhat_family_bars Logical -- bars of median R-hat by family.
#'
#' @return (Invisibly) a list with created ggplots/tables.
#' @export
plot_convergence_checks <- function(samples,
                                    out_dir = "outputs/diagnostics",
                                    top_k_rhat = 12L,
                                    top_k_aelow = 12L,
                                    runtime_s = NULL,
                                    rhat_ref = 1.05,  # coherent avec ton trait 0.05
                                    make_rhat_hist   = TRUE,   # -> utilisera le BAR "template"
                                    make_traces_rhat = TRUE,
                                    make_traces_ae   = TRUE,
                                    make_rhat_family_bars = TRUE) {

  stopifnot(inherits(samples, "mcmc.list"))
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  suppressWarnings(ok_grid <- requireNamespace("gridExtra", quietly = TRUE))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  save_fig <- function(p, base, w = 7.5, h = 5.0, dpi = 180) {
    ggplot2::ggsave(file.path(out_dir, paste0(base, ".pdf")), p, width = w, height = h)
    ggplot2::ggsave(file.path(out_dir, paste0(base, ".png")), p, width = w, height = h, dpi = dpi)
  }

  mcpar <- attr(samples[[1]], "mcpar")
  start <- if (!is.null(mcpar)) mcpar[1] else 1
  thin  <- if (!is.null(mcpar)) mcpar[3] else 1
  n_row <- nrow(samples[[1]])
  iter_seq <- seq(from = start, by = thin, length.out = n_row)

  if (is.null(runtime_s)) runtime_s <- NA_real_
  diag_tbl <- compute_diag_from_mcmc(samples, runtime_s = ifelse(is.finite(runtime_s), runtime_s, 1))
  if (!("ess_per_s" %in% names(diag_tbl)) && ("AE" %in% names(diag_tbl))) diag_tbl$ess_per_s <- diag_tbl$AE

  d_r <- diag_tbl[is.finite(diag_tbl$Rhat), , drop = FALSE]

  # ----- Global "hist" -> your FAMILY BAR TEMPLATE -----
  p_rhat_hist <- NULL
  if (isTRUE(make_rhat_hist) && NROW(d_r) > 0) {
    if (!("Family" %in% colnames(d_r))) d_r$Family <- sub("\\[.*$", "", as.character(d_r$target))
    rhat_family <- stats::aggregate(Rhat ~ Family, data = d_r, FUN = function(x) stats::median(x, na.rm = TRUE))
    names(rhat_family)[names(rhat_family) == "Rhat"] <- "median_Rhat"

    p_rhat_hist <- ggplot2::ggplot(
      rhat_family,
      ggplot2::aes(x = stats::reorder(Family, median_Rhat), y = median_Rhat - 1)
    ) +
      ggplot2::geom_bar(stat = "identity", fill = "orange", color = "black") +
      ggplot2::geom_hline(yintercept = (rhat_ref - 1), linetype = "dashed", color = "red", size = 1) +
      ggplot2::xlab("Nodes") +
      ggplot2::ylab("Median of Gelman-Rubin Rhat") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggplot2::scale_y_continuous(limits = c(0, rhat_ref - 1), labels = function(y) y + 1)
    save_fig(p_rhat_hist, "rhat_family_bars_overall")
  }

  # ----- Worst targets by R-hat (kept, styled like template) -----
  worst_rhat <- character(0)
  if (NROW(d_r) > 0) {
    has_ae <- "ess_per_s" %in% names(d_r)
    if (has_ae) {
      d_r$ess_per_s <- suppressWarnings(as.numeric(d_r$ess_per_s))
      if (!any(is.finite(d_r$ess_per_s))) has_ae <- FALSE
    }
    ord <- if (isTRUE(has_ae)) order(-d_r$Rhat, d_r$ess_per_s) else order(-d_r$Rhat)
    worst_rhat <- utils::head(d_r[ord, "target"], min(top_k_rhat, NROW(d_r)))
  }

  p_rhat_family_bars <- NULL
  rhat_family <- NULL
  if (isTRUE(make_rhat_family_bars) && NROW(d_r) > 0) {
    if (!("Family" %in% colnames(d_r))) d_r$Family <- sub("\\[.*$", "", as.character(d_r$target))
    rhat_family <- stats::aggregate(Rhat ~ Family, data = d_r, FUN = function(x) stats::median(x, na.rm = TRUE))
    names(rhat_family)[names(rhat_family) == "Rhat"] <- "median_Rhat"

    p_rhat_family_bars <- ggplot2::ggplot(
      rhat_family,
      ggplot2::aes(x = stats::reorder(Family, median_Rhat), y = median_Rhat - 1)
    ) +
      ggplot2::geom_bar(stat = "identity", fill = "orange", color = "black") +
      ggplot2::geom_hline(yintercept = (rhat_ref - 1), linetype = "dashed", color = "red", size = 1) +
      ggplot2::xlab("Nodes") +
      ggplot2::ylab("Median of Gelman-Rubin Rhat") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggplot2::scale_y_continuous(limits = c(0, rhat_ref - 1), labels = function(y) y + 1)
    ggplot2::ggsave(file.path(out_dir, "rhat_family_bars.pdf"), p_rhat_family_bars, width = 8, height = 5.2)
    ggplot2::ggsave(file.path(out_dir, "rhat_family_bars.png"), p_rhat_family_bars, width = 8, height = 5.2, dpi = 180)
  }

  # ---- Worst AE (low ESS/s) traces selection ----
  worst_ae <- character(0)
  has_ae_all <- ("ess_per_s" %in% names(diag_tbl)) && any(is.finite(diag_tbl$ess_per_s))
  if (has_ae_all) {
    d_ae <- diag_tbl[is.finite(diag_tbl$ess_per_s), , drop = FALSE]
    ord_ae <- if ("time_s" %in% names(d_ae)) order(d_ae$ess_per_s, -d_ae$time_s) else order(d_ae$ess_per_s)
    worst_ae <- utils::head(d_ae[ord_ae, "target"], min(top_k_aelow, NROW(d_ae)))
  }

  # ---- Helpers for trace/density ----
  get_param_df <- function(param) {
    dfs <- lapply(seq_along(samples), function(k) {
      mat <- as.matrix(samples[[k]])
      if (!(param %in% colnames(mat))) return(NULL)
      data.frame(iter = iter_seq, value = mat[, param], chain = paste0("chain", k),
                 stringsAsFactors = FALSE)
    })
    dfs <- Filter(Negate(is.null), dfs)
    if (!length(dfs)) return(NULL)
    do.call(rbind, dfs)
  }

  mk_trace_density <- function(param, out_base) {
    df <- get_param_df(param); if (is.null(df)) return(NULL)
    p_trace <- ggplot2::ggplot(df, ggplot2::aes(x = iter, y = value, color = chain, group = chain)) +
      ggplot2::geom_line(alpha = 0.85) +
      ggplot2::labs(title = paste("Trace --", param), x = "Iteration", y = param) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "none")
    p_dens <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = chain)) +
      ggplot2::geom_density(alpha = 0.35) +
      ggplot2::labs(title = paste("Density --", param), x = param, y = "Density") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "none")

    if (ok_grid) {
      gp <- gridExtra::arrangeGrob(p_trace, p_dens, ncol = 1)
      ggplot2::ggsave(file.path(out_dir, paste0(out_base, ".pdf")), gp, width = 7.2, height = 7.8)
      ggplot2::ggsave(file.path(out_dir, paste0(out_base, ".png")), gp, width = 7.2, height = 7.8, dpi = 180)
    } else {
      ggplot2::ggsave(file.path(out_dir, paste0(out_base, "_trace.pdf")), p_trace, width = 7.2, height = 4.2)
      ggplot2::ggsave(file.path(out_dir, paste0(out_base, "_trace.png")), p_trace, width = 7.2, height = 4.2, dpi = 180)
      ggplot2::ggsave(file.path(out_dir, paste0(out_base, "_dens.pdf")),  p_dens,  width = 7.2, height = 3.6)
      ggplot2::ggsave(file.path(out_dir, paste0(out_base, "_dens.png")),  p_dens,  width = 7.2, height = 3.6, dpi = 180)
    }
    list(trace = p_trace, density = p_dens)
  }

  traces_rhat <- list()
  if (isTRUE(make_traces_rhat) && length(worst_rhat)) {
    traces_rhat <- stats::setNames(vector("list", length(worst_rhat)), worst_rhat)
    for (p in worst_rhat) {
      traces_rhat[[p]] <- mk_trace_density(
        p, out_base = paste0("trace_rhat_worst_", gsub("[^A-Za-z0-9_\\-]", "_", p))
      )
    }
  }

  traces_ae <- list()
  if (isTRUE(make_traces_ae) && length(worst_ae)) {
    traces_ae <- stats::setNames(vector("list", length(worst_ae)), worst_ae)
    for (p in worst_ae) {
      traces_ae[[p]] <- mk_trace_density(
        p, out_base = paste0("trace_ae_worst_", gsub("[^A-Za-z0-9_\\-]", "_", p))
      )
    }
  }

  invisible(list(
    rhat_hist          = p_rhat_hist,       # now the family bar template
    rhat_family_bars   = p_rhat_family_bars,
    rhat_family_table  = rhat_family,
    traces_rhat        = traces_rhat,
    traces_ae          = traces_ae
  ))
}
#' Plot diagnostics for stochastic sampler targets (strict) + all-nodes panels
#'
#' @description
#' - **Sampler-only core**: strictly restrict AE/CE (and related "sampled" panels)
#'   to targets having an MCMC sampler, discovered *inside* via
#'   `conf.mcmc$getSamplers()` (or `attr(diag_tbl,"mcmc_conf")`, or explicit
#'   `attr(diag_tbl,"sampled_targets")`). No "Family" logic anywhere.
#' - **All-nodes replacements (this change)**: replace the former sampled-only
#'   panels you listed with **all-nodes** counterparts for the **same metrics**:
#'     * AE median by Target — **all nodes** (steelblue)
#'     * CE median by Target — **all nodes** (green)
#'     * Worst CE — **all nodes** (grey60 + labels)
#'     * Top Time — **all nodes** (grey40 + labels)
#'
#' @param diag_tbl data.frame with per-target diagnostics; must contain `target` or rownames.
#'        Expected: AE (`AE_ESS_per_it` or `AE`), CE (`ESS_per_sec` or `ess_per_s` or `CE`);
#'        optional `ESS`/`ess`, `time_s`, `Rhat`.
#' @param out_dir Output directory for figures (PDF + PNG). Default: "outputs/diagnostics".
#' @param top_k   Bars shown in the “worst/top” all-nodes panels. Default: 20L.
#' @param n_worst Rows returned in worst tables. Default: 20L.
#' @param bins    Kept for API compatibility (AE/CE “hist” are barplots). Default: 30L.
#' @param rhat_ref Rhat reference (Rhat-1 axis). Default: 1.05.
#' @param make_bar_ae_median_all, make_bar_ce_median_all,
#'        make_bar_ce_worst_all, make_bar_time_top_all
#'        Toggles for the **all-nodes** panels (default TRUE).
#' @param make_bar_rhat_worst, make_bar_ess_q5,
#'        make_hist_ce, make_hist_ae Toggles for the existing **sampled-only** panels (default TRUE).
#'
#' @return Invisibly, a list with ggplot objects and worst-node tables.
#' @export
plot_bottlenecks_index <- function(
    diag_tbl,
    out_dir      = "outputs/diagnostics",
    top_k        = 20L,
    n_worst      = 20L,
    bins         = 30L,
    rhat_ref     = 1.05,
    # ---- ALL-NODES (new/replacement) ----
    make_bar_ae_median_all = TRUE,
    make_bar_ce_median_all = TRUE,
    make_bar_ce_worst_all  = TRUE,
    make_bar_time_top_all  = TRUE,
    # ---- SAMPLED-ONLY (unchanged) ----
    make_bar_rhat_worst = TRUE,
    make_bar_ess_q5     = TRUE,
    make_hist_ce        = TRUE,  # full sampled set, CE median bars (green)
    make_hist_ae        = TRUE   # full sampled set, AE median bars (steelblue)
) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ----- utils -----
  as_num <- function(x) suppressWarnings(as.numeric(x))
  save_fig <- function(p, base, w = 8, h = 5.2, dpi = 180) {
    ggplot2::ggsave(file.path(out_dir, paste0(base, ".pdf")), p, width = w, height = h)
    ggplot2::ggsave(file.path(out_dir, paste0(base, ".png")), p, width = w, height = h, dpi = dpi)
  }
  flatten_chr <- function(x) if (is.list(x)) unique(unlist(lapply(x, flatten_chr), use.names = FALSE)) else as.character(x)

  # ----- normalize diag_tbl -----
  if (is.null(diag_tbl) || !NROW(diag_tbl)) stop("plot_bottlenecks_index: 'diag_tbl' is empty.")
  d <- as.data.frame(diag_tbl, stringsAsFactors = FALSE)
  if (!("target" %in% names(d))) {
    rn <- rownames(d); if (is.null(rn)) stop("plot_bottlenecks_index: missing 'target' and rownames.")
    d$target <- rn
  }

  # Raw metrics (will be used for both “all-nodes” and “sampled-only”)
  d$AE      <- if ("AE_ESS_per_it" %in% names(d)) as_num(d$AE_ESS_per_it)
  else if ("AE" %in% names(d))       as_num(d$AE) else NA_real_
  d$CE      <- if ("ESS_per_sec" %in% names(d))   as_num(d$ESS_per_sec)
  else if ("ess_per_s" %in% names(d))as_num(d$ess_per_s)
  else if ("CE" %in% names(d))       as_num(d$CE) else NA_real_
  d$ESS     <- if ("ESS" %in% names(d))           as_num(d$ESS)
  else if ("ess" %in% names(d))      as_num(d$ess) else NA_real_
  d$time_s  <- as_num(if ("time_s" %in% names(d)) d$time_s else NA_real_)
  d$Rhat    <- if ("Rhat" %in% names(d))          as_num(d$Rhat) else NA_real_

  # If needed, infer time = ESS/CE
  infer_ok <- is.finite(d$ESS) & is.finite(d$CE) & d$CE > 0 & !is.finite(d$time_s)
  d$time_s[infer_ok] <- d$ESS[infer_ok] / d$CE[infer_ok]

  # =========================
  # ALL-NODES AGGREGATIONS
  # =========================
  agg_all <- stats::aggregate(
    cbind(AE = d$AE, CE = d$CE, ESS = d$ESS, time_s = d$time_s, Rhat = d$Rhat),
    by = list(target = d$target),
    FUN = function(z) suppressWarnings(stats::median(z, na.rm = TRUE))
  )
  # Keep targets with any finite metric (per-plot filters still apply)
  keep_all <- is.finite(agg_all$AE) | is.finite(agg_all$CE) | is.finite(agg_all$ESS) |
    is.finite(agg_all$time_s) | is.finite(agg_all$Rhat)
  agg_all <- agg_all[keep_all, , drop = FALSE]

  # Sorted frames for "all-nodes" panels
  ae_all_ord <- agg_all[order(agg_all$AE, agg_all$target), c("target","AE"), drop = FALSE]
  ce_all_ord <- agg_all[order(agg_all$CE, agg_all$target), c("target","CE"), drop = FALSE]
  tm_all_ord <- agg_all[order(-agg_all$time_s, agg_all$target), c("target","time_s","CE"), drop = FALSE]
  ce_all_worst <- utils::head(ce_all_ord, min(top_k, nrow(ce_all_ord)))
  tm_all_top   <- utils::head(tm_all_ord, min(top_k, nrow(tm_all_ord)))

  # ============================
  # SAMPLED-ONLY AGGREGATIONS
  # ============================
  # INTERNAL routine: extract sampler targets
  get_sampler_targets <- function(d_tbl) {
    pf  <- parent.frame()
    cfg <- get0("conf.mcmc", envir = pf, inherits = TRUE)
    if (is.null(cfg)) cfg <- get0("conf.mcmc.model", envir = pf, inherits = TRUE)
    if (is.null(cfg)) cfg <- attr(d_tbl, "mcmc_conf", exact = TRUE)

    if (!is.null(cfg)) {
      sams <- try(cfg$getSamplers(), silent = TRUE)
      if (!inherits(sams, "try-error") && length(sams)) {
        tgts <- unique(unlist(lapply(sams, function(s) s$target), use.names = FALSE))
        tgts <- tgts[!grepl("^lifted_|^logProb_", tgts)]
        return(intersect(tgts, d_tbl$target))
      }
    }
    tgts_attr <- attr(d_tbl, "sampled_targets", exact = TRUE)
    if (!is.null(tgts_attr)) {
      tgts <- flatten_chr(tgts_attr)
      tgts <- tgts[!grepl("^lifted_|^logProb_", tgts)]
      return(intersect(tgts, d_tbl$target))
    }
    character(0)
  }

  sampled <- get_sampler_targets(d)
  if (!length(sampled)) {
    stop("plot_bottlenecks_index: no sampler targets detected. Ensure 'conf.mcmc' (or 'conf.mcmc.model') is visible, or attach it via attr(diag_tbl,'mcmc_conf').")
  }

  # HARD FILTER to samplers only (for sampled-only panels)
  ds <- d[d$target %in% sampled, , drop = FALSE]
  if (!NROW(ds)) stop("plot_bottlenecks_index: after filtering, no sampled targets remain.")

  agg_s <- stats::aggregate(
    cbind(AE = ds$AE, CE = ds$CE, ESS = ds$ESS, time_s = ds$time_s, Rhat = ds$Rhat),
    by = list(target = ds$target),
    FUN = function(z) suppressWarnings(stats::median(z, na.rm = TRUE))
  )
  keep_s <- is.finite(agg_s$AE) | is.finite(agg_s$CE) | is.finite(agg_s$ESS) |
    is.finite(agg_s$time_s) | is.finite(agg_s$Rhat)
  agg_s <- agg_s[keep_s, , drop = FALSE]
  if (!NROW(agg_s)) stop("plot_bottlenecks_index: no finite metrics on sampled targets.")

  # Sorted frames for sampled-only panels
  ae_s_ord  <- agg_s[order(agg_s$AE,  agg_s$target), c("target","AE"),  drop = FALSE]
  ce_s_ord  <- agg_s[order(agg_s$CE,  agg_s$target), c("target","CE"),  drop = FALSE]
  tm_s_ord  <- agg_s[order(-agg_s$time_s, agg_s$target), c("target","time_s","CE"), drop = FALSE]
  rh_s_ord  <- agg_s[order(-agg_s$Rhat,    agg_s$target), c("target","Rhat","CE"), drop = FALSE]
  ess_s_ord <- agg_s[order(agg_s$ESS,      agg_s$target), c("target","ESS","CE"),  drop = FALSE]

  res <- list()

  # =========================
  # NEW / REPLACEMENT PANELS
  # =========================

  # (ALL NODES) AE median bar — full set (steelblue)
  if (make_bar_ae_median_all && nrow(ae_all_ord)) {
    p <- ggplot2::ggplot(ae_all_ord,
                         ggplot2::aes(x = stats::reorder(target, AE), y = AE)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::labs(
        title = "Median ESS/N by Target — All Nodes",
        subtitle = paste0("All nodes with finite AE (n = ", nrow(ae_all_ord), ")"),
        x = "Targets", y = "Median ESS/N"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "AE_median_by_target__all_nodes")
    res$AE_median_by_target__all_nodes <- p
  }

  # (ALL NODES) CE median bar — full set (green)
  if (make_bar_ce_median_all && nrow(ce_all_ord)) {
    p <- ggplot2::ggplot(ce_all_ord,
                         ggplot2::aes(x = stats::reorder(target, CE), y = CE)) +
      ggplot2::geom_bar(stat = "identity", fill = "green") +
      ggplot2::labs(
        title = "Median CE ESS/s by Target — All Nodes",
        subtitle = paste0("All nodes with finite CE (n = ", nrow(ce_all_ord), ")"),
        x = "Targets", y = "Median ESS/s"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "CE_median_by_target__all_nodes")
    res$CE_median_by_target__all_nodes <- p
  }

  # (ALL NODES) Worst CE — top_k (grey60 + labels)
  if (make_bar_ce_worst_all && nrow(ce_all_worst)) {
    p <- ggplot2::ggplot(ce_all_worst,
                         ggplot2::aes(x = stats::reorder(target, CE), y = CE)) +
      ggplot2::geom_bar(stat = "identity", fill = "orange") +
      ggplot2::geom_text(ggplot2::aes(label = round(CE, 3)), vjust = -0.5, size = 3) +
      ggplot2::labs(
        title = "Targets by CE ESS/s) — All Nodes",
        subtitle = paste0("Showing up to ", nrow(ce_all_worst), " worst targets"),
        x = "Targets", y = "Median ESS/s"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "CE_worst_targets__all_nodes")
    res$CE_worst_targets__all_nodes <- p
  }

  # (ALL NODES) Top Time — top_k (grey40 + labels)
  if (make_bar_time_top_all && nrow(tm_all_top)) {
    p <- ggplot2::ggplot(tm_all_top,
                         ggplot2::aes(x = stats::reorder(target, time_s), y = time_s)) +
      ggplot2::geom_bar(stat = "identity", fill = "red") +
      ggplot2::geom_text(ggplot2::aes(label = round(time_s, 2)), vjust = -0.5, size = 3) +
      ggplot2::labs(
        title = "Top Targets by Total Compute Time — All Nodes",
        subtitle = paste0("Showing up to ", nrow(tm_all_top), " highest times"),
        x = "Targets", y = "Total Time (s)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "Time_top_targets__all_nodes")
    res$Time_top_targets__all_nodes <- p
  }

  # =================
  # SAMPLED-ONLY PANELS (unchanged, still strict)
  # =================

  # Worst Rhat (Rhat-1) — sampled only (orange/black + dashed ref)
  if (make_bar_rhat_worst && nrow(rh_s_ord)) {
    ymax <- max(rh_s_ord$Rhat - 1, na.rm = TRUE)
    p <- ggplot2::ggplot(rh_s_ord,
                         ggplot2::aes(x = stats::reorder(target, Rhat), y = Rhat - 1)) +
      ggplot2::geom_bar(stat = "identity", fill = "orange", color = "black") +
      ggplot2::geom_hline(yintercept = (rhat_ref - 1), linetype = "dashed", color = "red", size = 1) +
      ggplot2::labs(
        title = "Rhat (Rhat − 1) — Sampled Only",
        subtitle = paste0("Dashed line at Rhat = ", rhat_ref),
        x = "Targets", y = "Median Rhat − 1"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggplot2::scale_y_continuous(limits = c(0, max(rhat_ref - 1, ymax)),
                                  labels = function(y) y + 1)
    save_fig(p, "Rhat_sampled_only")
    res$Rhat_worst_targets__sampled_only <- p
  }

  # ESS bar with 5% quantile — sampled only (skyblue/black + labels)
  if (make_bar_ess_q5 && nrow(ess_s_ord)) {
    ess_q5 <- tryCatch(stats::quantile(agg_s$ESS, probs = 0.05, na.rm = TRUE), error = function(e) NA_real_)
    p <- ggplot2::ggplot(utils::head(ess_s_ord, min(top_k, nrow(ess_s_ord))),
                         ggplot2::aes(x = stats::reorder(target, ESS), y = ESS)) +
      ggplot2::geom_bar(stat = "identity", fill = "skyblue", color = "black") +
      ggplot2::geom_text(ggplot2::aes(label = round(ESS, 1)), vjust = -0.5, size = 3) +
      ggplot2::geom_hline(yintercept = ess_q5, linetype = "dashed", color = "red", size = 1) +
      ggplot2::labs(
        title = "ESS Targets with 5th Quantile Line (Sampled Only)",
        subtitle = if (is.finite(ess_q5)) paste("5th quantile =", round(ess_q5, 1)) else NULL,
        x = "Targets", y = "Effective Sample Size (ESS)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "ESS_worst_targets_q5_line__sampled_only")
    res$ESS_worst_targets_q5_line__sampled_only <- p
  }

  # “Hist” toggles (sampled-only, full set) = AE/CE barplots for the sampled set
  if (make_hist_ce && nrow(ce_s_ord)) {
    p <- ggplot2::ggplot(ce_s_ord,
                         ggplot2::aes(x = stats::reorder(target, CE), y = CE)) +
      ggplot2::geom_bar(stat = "identity", fill = "green") +
      ggplot2::labs(
        title = "Median ESS/s by Target — Sampled Only",
        subtitle = paste0("All sampled targets (n = ", nrow(ce_s_ord), ")"),
        x = "Targets", y = "Median ESS/s"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "CE_median_by_target__full_sampled_set")
    res$CE_median_by_target__full_sampled_set <- p
  }
  if (make_hist_ae && nrow(ae_s_ord)) {
    p <- ggplot2::ggplot(ae_s_ord,
                         ggplot2::aes(x = stats::reorder(target, AE), y = AE)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::labs(
        title = "Median ESS/N by Target — Sampled Only ",
        subtitle = paste0("All sampled targets (n = ", nrow(ae_s_ord), ")"),
        x = "Targets", y = "Median ESS/N"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "AE_median_by_target__full_sampled_set")
    res$AE_median_by_target__full_sampled_set <- p
  }

  # ---- worst-node tables (you still have both universes available if needed) ----
  res$worst_nodes_all_computational      <- ce_all_worst
  res$worst_nodes_all_time_top           <- tm_all_top
  res$worst_nodes_sampled_rhat           <- rh_s_ord
  res$worst_nodes_sampled_ess            <- ess_s_ord

  invisible(res)
}


