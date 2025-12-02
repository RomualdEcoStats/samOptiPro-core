
## ---- Exported functions (6) ----

#' @export
#' @keywords internal
enrich_hmc_diag_tbl_for_plots <- function(res) {
  if (is.null(res$diag_tbl)) {
    stop("L'objet 'res' ne contient pas de diag_tbl.")
  }
  diag_tbl <- res$diag_tbl

  # 1) AE_ESS_per_it : on le mappe sur AE_worst (ou AE_total)
  if (!"AE_ESS_per_it" %in% names(diag_tbl)) {
    if ("AE_worst" %in% names(diag_tbl)) {
      diag_tbl$AE_ESS_per_it <- diag_tbl$AE_worst
    } else if ("AE_total" %in% names(diag_tbl)) {
      diag_tbl$AE_ESS_per_it <- diag_tbl$AE_total
    } else {
      diag_tbl$AE_ESS_per_it <- NA_real_
    }
  }

  # 2) CE_ESS_per_s : on le mappe sur ESS_per_sec_worst / total / ESS_per_sec / ESS_per_s
  if (!"CE_ESS_per_s" %in% names(diag_tbl)) {
    if ("ESS_per_sec_worst" %in% names(diag_tbl)) {
      diag_tbl$CE_ESS_per_s <- diag_tbl$ESS_per_sec_worst
    } else if ("ESS_per_sec_total" %in% names(diag_tbl)) {
      diag_tbl$CE_ESS_per_s <- diag_tbl$ESS_per_sec_total
    } else if ("ESS_per_sec" %in% names(diag_tbl)) {
      diag_tbl$CE_ESS_per_s <- diag_tbl$ESS_per_sec
    } else if ("ESS_per_s" %in% names(diag_tbl)) {
      diag_tbl$CE_ESS_per_s <- diag_tbl$ESS_per_s
    } else {
      diag_tbl$CE_ESS_per_s <- NA_real_
    }
  }

  # 3) Rhat générique (à partir de classic / split si dispo)
  if (!"Rhat" %in% names(diag_tbl)) {
    rhat <- rep(NA_real_, nrow(diag_tbl))

    if ("Rhat_classic" %in% names(diag_tbl)) {
      rhat <- diag_tbl$Rhat_classic
    }
    if ("Rhat_split" %in% names(diag_tbl)) {
      idx <- is.na(rhat) & !is.na(diag_tbl$Rhat_split)
      rhat[idx] <- diag_tbl$Rhat_split[idx]
    }

    diag_tbl$Rhat <- rhat
  }

  res$diag_tbl <- diag_tbl
  res
}


#' Plot MCMC Bottlenecks by Node or Family
#'
#' @description
#' Generates a comprehensive diagnostic panel of MCMC bottlenecks using
#' efficiency and convergence metrics computed per node or node family.
#' The function can optionally restrict analyses to nodes that are
#' effectively sampled (i.e. have associated samplers in the NIMBLE
#' configuration), identified automatically via `conf.mcmc$getSamplers()`.
#'
#' It produces publication-ready figures of:
#' \itemize{
#'   \item Median *Algorithmic Efficiency* (AE = ESS/iter) by node family;
#'   \item Median *Computational Efficiency* (CE = ESS/s) by node family;
#'   \item Worst targets by CE (lowest ESS/s);
#'   \item Median or worst \eqn{\hat{R}} (Gelman–Rubin) by node or family.
#' }
#'
#' The function saves each plot both as PDF and PNG in the specified output directory.
#' Bar widths and spacing are optimized for compact presentation.
#'
#' @param diag_tbl `data.frame` or `tibble` containing diagnostics per target node.
#'   Must include columns such as `target`, `Family`, `ESS`, `ESS_per_sec`, and optionally `Rhat`.
#' @param out_dir Character string; path to output directory for saving figures
#'   (default: `"outputs/diagnostics"`). Will be created recursively if missing.
#' @param top_k Integer; number of worst or best nodes to display (default: `20L`).
#' @param rhat_ref Numeric; reference threshold for Gelman–Rubin \eqn{\hat{R}} (default: `1.05`).
#' @param sampled_only Logical; if `TRUE`, restricts plots to nodes that have an
#'   explicit sampler in `conf.mcmc` and are present in `samples_ml`. Default: `FALSE`.
#' @param conf.mcmc NIMBLE MCMC configuration object, typically produced by
#'   `configureMCMC(model, ...)` or stored as `build_fn()$conf`. Used to extract
#'   sampler-attached target nodes when `sampled_only = TRUE`.
#' @param samples_ml Optional MCMC list (as returned by `runMCMC(..., nchains > 1)`),
#'   used to match sampler targets with actual sample column names.
#'
#' @param make_esss_targets Logical; if `TRUE`, produces barplot of worst targets
#'   by computational efficiency (default: `TRUE`).
#' @param make_esss_families Logical; if `TRUE`, produces barplot of median
#'   algorithmic efficiency (AE) by node family (default: `TRUE`).
#' @param make_time_families Logical; if `TRUE`, produces barplot of median
#'   computational efficiency (CE) by node family (default: `TRUE`).
#' @param make_rhat_hist_targets Logical; if `TRUE`, produces barplot of median
#'   Rhat per family (default: `TRUE`).
#' @param make_rhat_worst_targets Logical; if `TRUE`, produces barplot of the
#'   worst Rhat targets (default: `TRUE`).
#' @param make_rhat_median_families Logical; if `TRUE`, produces an alias
#'   of the median Rhat-by-family plot (default: `TRUE`).
#'
#' @details
#' This function is a core visualization tool for diagnosing performance
#' bottlenecks in large hierarchical Bayesian models (e.g., SAM-like or
#' GEREM-type stock assessment models). It integrates runtime, efficiency,
#' and convergence metrics in a standardized panel of plots, suitable for
#' benchmarking, model comparison, or publication figures.
#'
#' When `sampled_only = TRUE`, it automatically extracts the list of stochastic
#' nodes (targets) from `conf.mcmc$getSamplers()` and intersects them with the
#' variable names present in `samples_ml`. This ensures only stochastically
#' sampled nodes are visualized, excluding lifted or deterministic intermediates.
#'
#' @return
#' Invisibly returns a named `list` of ggplot objects:
#' \itemize{
#'   \item `bar_family_algorithmic_eff` — Median AE by family;
#'   \item `bar_family_computational_eff` — Median CE by family;
#'   \item `bar_target_CE` — Worst targets by CE;
#'   \item `rhat_family_template` — Median Rhat by family;
#'   \item `rhat_worst_targets` — Worst targets by Rhat.
#' }
#' Each plot is also saved in `out_dir` as both `.pdf` and `.png`.
#'
#' @seealso
#' \code{\link{identify_bottlenecks_family}},
#' \code{\link{profile_sampler_times}},
#' \code{\link{run_structure_and_hmc_test}}
#'
#' @examples
#' \dontrun{
#' # Example assuming an existing NIMBLE configuration and MCMC results:
#' res <- plot_bottlenecks(
#'   diag_tbl     = diag_tbl,
#'   conf.mcmc    = conf.mcmc,
#'   samples_ml   = samples_ml,
#'   sampled_only = TRUE,
#'   out_dir      = "outputs/diagnostics"
#' )
#' }
#'
#'
#' @export
plot_bottlenecks <- function(diag_tbl,
                             out_dir = "outputs/diagnostics",
                             top_k   = 20L,
                             rhat_ref = 1.05,
                             sampled_only = FALSE,
                             conf.mcmc    = NULL,
                             samples_ml   = NULL,
                             make_esss_targets            = TRUE,
                             make_esss_families           = TRUE,
                             make_time_families           = TRUE,
                             make_rhat_hist_targets       = TRUE,
                             make_rhat_worst_targets      = TRUE,
                             make_rhat_median_families    = TRUE,
                             make_hist_ae_families        = FALSE,
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

  # ---------- normalisation colonnes ----------
  d <- as.data.frame(diag_tbl, stringsAsFactors = FALSE)
  if (!("target" %in% names(d))) {
    rn <- rownames(d); if (is.null(rn)) stop("plot_bottlenecks: missing 'target' column.")
    d$target <- rn
  }
  if (!("Family" %in% names(d))) d$Family <- sub("\\[.*", "", d$target)

  # ---------- filtre samplers only ----------
  if (isTRUE(sampled_only)) {
    sampler_targets <- character(0)
    if (!is.null(conf.mcmc) && is.function(conf.mcmc$getSamplers)) {
      smp <- conf.mcmc$getSamplers()
      sampler_targets <- unique(unlist(lapply(smp, function(s) s$target)))
    }
    if (length(sampler_targets)) {
      sampler_targets <- gsub("^lifted_|^logProb_", "", sampler_targets, perl = TRUE)
    }
    if (!is.null(samples_ml) && length(samples_ml) >= 1L) {
      cn <- try(colnames(as.data.frame(samples_ml[[1]])), silent = TRUE)
      if (!inherits(cn, "try-error") && length(cn)) {
        sampler_targets <- intersect(sampler_targets, cn)
      }
    }
    if (length(sampler_targets)) {
      d <- d[d$target %in% sampler_targets, , drop = FALSE]
    } else {
      warning("sampled_only=TRUE mais aucune cible échantillonnée trouvée.")
    }
  }

  # ---------- métriques ----------
  d$AE <- if ("AE_ESS_per_it" %in% names(d)) as_num(d$AE_ESS_per_it) else if ("AE" %in% names(d)) as_num(d$AE) else NA_real_
  d$CE <- if ("ESS_per_sec" %in% names(d)) as_num(d$ESS_per_sec) else if ("ess_per_s" %in% names(d)) as_num(d$ess_per_s) else if ("CE" %in% names(d)) as_num(d$CE) else NA_real_
  d$ESS_abs <- if ("ESS" %in% names(d)) as_num(d$ESS) else if ("ess" %in% names(d)) as_num(d$ess) else NA_real_
  d$time_s <- as_num(if ("time_s" %in% names(d)) d$time_s else NA_real_)
  infer_ok <- is.finite(d$ESS_abs) & is.finite(d$CE) & d$CE > 0 & !is.finite(d$time_s)
  d$time_s[infer_ok] <- d$ESS_abs[infer_ok] / d$CE[infer_ok]
  d$Rhat <- if ("Rhat" %in% names(d)) as_num(d$Rhat) else NA_real_

  # ---------- agrégations par famille ----------
  fam_ce   <- agg_safe(CE ~ Family,     d[is.finite(d$CE), ], median, na.rm = TRUE); names(fam_ce)[2]   <- "CE_median"
  fam_ae   <- agg_safe(AE ~ Family,     d[is.finite(d$AE), ], median, na.rm = TRUE); names(fam_ae)[2]   <- "AE_median"
  fam_time <- agg_safe(time_s ~ Family, d[is.finite(d$time_s), ], sum,    na.rm = TRUE); names(fam_time)[2] <- "time_sum"

  res <- list()
  grouped_data <- merge(fam_ae, fam_ce, by = "Family", all = TRUE)
  names(grouped_data)[names(grouped_data) == "AE_median"] <- "MedianAlgorithmicEfficiency"
  names(grouped_data)[names(grouped_data) == "CE_median"] <- "MedianComputationalEfficiency_tot"

  # ---------- AE et CE par famille (espacement réduit) ----------
  if (isTRUE(make_esss_families) && NROW(grouped_data)) {
    p1 <- ggplot2::ggplot(grouped_data,
                          ggplot2::aes(x = stats::reorder(Family, MedianAlgorithmicEfficiency),
                                       y = MedianAlgorithmicEfficiency)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +  # width réduit
      ggplot2::labs(title = "Median Algorithmic Efficiency by Node Family",
                    x = "Node Family", y = "Median Algorithmic Efficiency") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p1, "bar_family_algorithmic_eff")
    res$bar_family_algorithmic_eff <- p1
  }

  if (isTRUE(make_time_families) && NROW(grouped_data)) {
    p2 <- ggplot2::ggplot(grouped_data,
                          ggplot2::aes(x = stats::reorder(Family, MedianComputationalEfficiency_tot),
                                       y = MedianComputationalEfficiency_tot)) +
      ggplot2::geom_bar(stat = "identity", fill = "darkgreen", width = 0.8) +
      ggplot2::labs(title = "Median Computational Efficiency by Node Family",
                    x = "Node Family", y = "Median Computational Efficiency (ESS/s)") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p2, "bar_family_computational_eff")
    res$bar_family_computational_eff <- p2
  }

  # ---------- CE par target (espacement réduit) ----------
  d_ce <- d[is.finite(d$CE), , drop = FALSE]
  if (isTRUE(make_esss_targets) && NROW(d_ce)) {
    d_worst <- headn(d_ce[order(d_ce$CE, -d_ce$time_s), ], min(top_k, NROW(d_ce)))
    p3 <- ggplot2::ggplot(d_worst,
                          ggplot2::aes(x = stats::reorder(target, CE), y = CE)) +
      ggplot2::geom_bar(stat = "identity", fill = "grey50", width = 0.8) +
      ggplot2::labs(title = "Worst Targets by Computational Efficiency (ESS/s)",
                    x = "Targets", y = "ESS/s") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p3, "bar_target_CE_worst")
    res$bar_target_CE <- p3
  }

  # ---------- RHAT ----------
  d_r <- d[is.finite(d$Rhat), , drop = FALSE]
  if (NROW(d_r)) {
    fam_rhat <- agg_safe(Rhat ~ Family, data = d_r, FUN = stats::median, na.rm = TRUE)
    names(fam_rhat)[2] <- "median_Rhat"

    if (isTRUE(make_rhat_hist_targets) && NROW(fam_rhat)) {
      p4 <- ggplot2::ggplot(fam_rhat,
                            ggplot2::aes(x = stats::reorder(Family, median_Rhat), y = median_Rhat - 1)) +
        ggplot2::geom_bar(stat = "identity", fill = "orange", color = "black", width = 0.8) +
        ggplot2::geom_hline(yintercept = (rhat_ref - 1),
                            linetype = "dashed", color = "red", size = 1) +
        ggplot2::xlab("Node Family") + ggplot2::ylab("Median Rhat - 1") +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
      save_fig(p4, "rhat_family_bars_template")
      res$rhat_family_template <- p4
    }

    if (isTRUE(make_rhat_worst_targets)) {
      ord <- order(-d_r$Rhat, d_r$CE)
      rhat_worst <- headn(d_r[ord, , drop = FALSE], min(top_k, NROW(d_r)))
      p5 <- ggplot2::ggplot(rhat_worst,
                            ggplot2::aes(x = stats::reorder(target, Rhat), y = Rhat - 1)) +
        ggplot2::geom_bar(stat = "identity", fill = "orange", color = "black", width = 0.8) +
        ggplot2::geom_hline(yintercept = (rhat_ref - 1),
                            linetype = "dashed", color = "red", size = 1) +
        ggplot2::xlab("Targets") + ggplot2::ylab("Rhat - 1") +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
      save_fig(p5, "rhat_worst_targets_template")
      res$rhat_worst_targets <- p5
    }
  }

  invisible(res)
}

#' Fast bottleneck plots (samplers-only) for large models
#'
#' Produces **colored bar charts** for median **AE (ESS/N)** and **CE (ESS/s)**
#' by node family, **restricted to nodes that are actually sampled** by the MCMC
#' configuration (strict "samplers only"). Designed for large models where
#' memory/time constraints require fast post-processing on a reduced subset.
#'
#' @param diag_tbl A `data.frame` with at least column `target`. If present and finite,
#'   metrics are taken from `AE_worst` (fallback to `AE_total`) and
#'   `ESS_per_sec_worst` (fallback to `ESS_per_sec_total`).
#' @param sampled_only Logical (default `TRUE`). When `TRUE`, keep **only** nodes
#'   that are targets of configured samplers and also present in `samples_ml` columns.
#'   If no node remains, the function raises an error (no permissive fallback).
#' @param conf.mcmc A Nimble MCMC configuration (result of `nimble::configureMCMC()`),
#'   required when `sampled_only=TRUE`.
#' @param samples_ml A `coda::mcmc.list` (or compatible) holding posterior draws,
#'   required when `sampled_only=TRUE`.
#' @param out_dir Optional output directory. When provided, saves:
#'   `AE_median_by_family__samplers_only.png` and
#'   `CE_median_by_family__samplers_only.png`.
#'
#' @return (invisible) A list with ggplot objects and summary tables:
#'   `list(p_ae, p_ce, table_ae, table_ce)`.
#'
#' @examples
#' \dontrun{
#' conf <- nimble::configureMCMC(m)
#' res <- plot_bottlenecks_fast(
#'   diag_tbl   = diag_tbl,
#'   sampled_only = TRUE,
#'   conf.mcmc  = conf,
#'   samples_ml = samples_mla,
#'   out_dir    = "outputs/diagnostics_fast"
#' )
#' }
#' @export
plot_bottlenecks_fast <- function(
    diag_tbl,
    sampled_only = TRUE,
    conf.mcmc = NULL,
    samples_ml = NULL,
    out_dir = NULL
) {
  stopifnot(is.data.frame(diag_tbl), "target" %in% names(diag_tbl))

  # --- Pick AE/CE with finiteness checks
  pick_vec <- function(primary, fallback) {
    if (primary %in% names(diag_tbl) && any(is.finite(diag_tbl[[primary]]))) {
      as.numeric(diag_tbl[[primary]])
    } else if (fallback %in% names(diag_tbl) && any(is.finite(diag_tbl[[fallback]]))) {
      as.numeric(diag_tbl[[fallback]])
    } else {
      rep(NA_real_, nrow(diag_tbl))
    }
  }
  AE <- pick_vec("AE_worst", "AE_total")
  CE <- pick_vec("ESS_per_sec_worst", "ESS_per_sec_total")

  df <- data.frame(
    target = diag_tbl$target,
    AE = AE,
    CE = CE,
    stringsAsFactors = FALSE
  )
  df$Family <- .family_from_fast(df$target)

  # --- Strict samplers-only filter
  if (isTRUE(sampled_only)) {
    if (is.null(conf.mcmc) || is.null(samples_ml)) {
      stop("plot_bottlenecks_fast: sampled_only=TRUE requires 'conf.mcmc' and 'samples_ml'.")
    }
    sampler_targets <- .derive_sampler_targets_fast(conf.mcmc, samples_ml)
    df <- df[df$target %in% sampler_targets, , drop = FALSE]
    if (nrow(df) == 0L) {
      stop("plot_bottlenecks_fast: no sampled nodes found in 'diag_tbl' after strict filtering.")
    }
  }

  # --- Keep finite values only for aggregation
  dff <- df[is.finite(df$AE) | is.finite(df$CE), , drop = FALSE]
  if (nrow(dff) == 0L) stop("plot_bottlenecks_fast: no finite AE/CE metrics after filtering.")

  # --- Medians by family
  fam_ae <- stats::aggregate(AE ~ Family, data = dff, FUN = function(x) stats::median(x, na.rm = TRUE))
  names(fam_ae)[2] <- "AE_median"
  fam_ce <- stats::aggregate(CE ~ Family, data = dff, FUN = function(x) stats::median(x, na.rm = TRUE))
  names(fam_ce)[2] <- "CE_median"

  # --- Colored bar charts
  p_ae <- ggplot2::ggplot(fam_ae, ggplot2::aes(x = reorder(Family, AE_median), y = AE_median, fill = Family)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(title = "Median AE by Family (samplers only)", x = "Family", y = "AE (ESS / N)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")

  p_ce <- ggplot2::ggplot(fam_ce, ggplot2::aes(x = reorder(Family, CE_median), y = CE_median, fill = Family)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(title = "Median CE by Family (samplers only)", x = "Family", y = "CE (ESS / s)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")

  # --- Save if requested
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(out_dir, "AE_median_by_family__samplers_only.png"), p_ae, width = 12, height = 8, dpi = 300)
    ggplot2::ggsave(file.path(out_dir, "CE_median_by_family__samplers_only.png"), p_ce, width = 12, height = 8, dpi = 300)
  }

  invisible(list(p_ae = p_ae, p_ce = p_ce, table_ae = fam_ae, table_ce = fam_ce))
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
  ce_all_worst <- utils::head(ce_all_ord, min(top_k, nrow(ce_all_ord)))

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
  rh_s_ord  <- agg_s[order(-agg_s$Rhat,    agg_s$target), c("target","Rhat","CE"), drop = FALSE]
  ess_s_ord <- agg_s[order(agg_s$ESS,      agg_s$target), c("target","ESS","CE"),  drop = FALSE]

  res <- list()

  # =========================
  # NEW / REPLACEMENT PANELS
  # =========================

  # (ALL NODES) AE  bar — full set (steelblue)
  if (make_bar_ae_median_all && nrow(ae_all_ord)) {
    p <- ggplot2::ggplot(ae_all_ord,
                         ggplot2::aes(x = stats::reorder(target, AE), y = AE)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::labs(
        title = "ESS/N by Target — All Nodes",
        subtitle = paste0("All nodes with finite AE (n = ", nrow(ae_all_ord), ")"),
        x = "Targets", y = "ESS/N"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "AE_by_target__all_nodes")
    res$AE_median_by_target__all_nodes <- p
  }

  # (ALL NODES) CE bar — full set (green)
  if (make_bar_ce_median_all && nrow(ce_all_ord)) {
    p <- ggplot2::ggplot(ce_all_ord,
                         ggplot2::aes(x = stats::reorder(target, CE), y = CE)) +
      ggplot2::geom_bar(stat = "identity", fill = "green") +
      ggplot2::labs(
        title = "CE ESS/s by Target — All Nodes",
        subtitle = paste0("All nodes with finite CE (n = ", nrow(ce_all_ord), ")"),
        x = "Targets", y = "ESS/s"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "CE_by_target__all_nodes")
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
        x = "Targets", y = "ESS/s"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "CE_worst_targets__all_nodes")
    res$CE_worst_targets__all_nodes <- p
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
        x = "Targets", y = "Rhat − 1"
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
        title = "ESS/s by Target — Sampled Only",
        subtitle = paste0("All sampled targets (n = ", nrow(ce_s_ord), ")"),
        x = "Targets", y = "ESS/s"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "CE_by_target__full_sampled_set")
    res$CE_median_by_target__full_sampled_set <- p
  }
  if (make_hist_ae && nrow(ae_s_ord)) {
    p <- ggplot2::ggplot(ae_s_ord,
                         ggplot2::aes(x = stats::reorder(target, AE), y = AE)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::labs(
        title = "ESS/N by Target — Sampled Only ",
        subtitle = paste0("All sampled targets (n = ", nrow(ae_s_ord), ")"),
        x = "Targets", y = "ESS/N"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(p, "AE_by_target__full_sampled_set")
    res$AE_median_by_target__full_sampled_set <- p
  }

  # ---- worst-node tables (you still have both universes available if needed) ----
  res$worst_nodes_all_computational      <- ce_all_worst
  res$worst_nodes_sampled_rhat           <- rh_s_ord
  res$worst_nodes_sampled_ess            <- ess_s_ord

  invisible(res)
}
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

#' Fast convergence plots for R-hat (samplers-only), thresholded
#'
#' Builds **samplers-only** convergence diagnostics based on **R̂** with a
#' configurable threshold (default 1.01). Prefers `Rhat_split` if available,
#' falls back to `Rhat_classic` otherwise. Produces:
#' 1) median R̂ by family (bar chart with threshold line), and
#' 2) top-k worst nodes by R̂ (bar chart with threshold line).
#'
#' @param diag_tbl A `data.frame` with at least column `target` and one of
#'   `Rhat_split` or `Rhat_classic` having finite values.
#' @param threshold Numeric; R̂ reference line (default `1.01`).
#' @param sampled_only Logical (default `TRUE`). When `TRUE`, keep **only** nodes
#'   that are targets of configured samplers and also present in `samples_ml`.
#'   Raises an error if none remain (strict behavior).
#' @param conf.mcmc Nimble MCMC configuration (required when `sampled_only=TRUE`).
#' @param samples_ml A `coda::mcmc.list` providing column names to intersect with sampler targets.
#' @param out_dir Optional output directory. Saves:
#'   `Rhat_median_by_family__samplers_only_thr_*.png`,
#'   `TopK_worst_Rhat__samplers_only_thr_*.png`,
#'   and two CSV summaries.
#' @param top_k Integer; number of worst nodes by R̂ to display (default `20L`).
#' @param prefer_split Logical; if `TRUE` (default), use `Rhat_split` when available/finite,
#'   otherwise fallback to `Rhat_classic`.
#'
#' @return (invisible) A list with ggplot objects and tables:
#'   `list(p_family, p_top, table_family, table_top, threshold)`.
#'
#' @examples
#' \dontrun{
#' conf <- nimble::configureMCMC(m)
#' rhat_res <- plot_convergence_rhat_sampled_fast(
#'   diag_tbl     = diag_tbl,
#'   threshold    = 1.01,
#'   sampled_only = TRUE,
#'   conf.mcmc    = conf,
#'   samples_ml   = samples_mla,
#'   out_dir      = "outputs/diagnostics_fast_rhat",
#'   top_k        = 20L,
#'   prefer_split = TRUE
#' )
#' }
#' @export
plot_convergence_rhat_sampled_fast <- function(
    diag_tbl,
    threshold = 1.01,
    sampled_only = TRUE,
    conf.mcmc = NULL,
    samples_ml = NULL,
    out_dir = NULL,
    top_k = 20L,
    prefer_split = TRUE
) {
  stopifnot(is.data.frame(diag_tbl), "target" %in% names(diag_tbl))

  # --- choose R-hat column
  pick_rhat <- function(dt, prefer_split = TRUE) {
    r_split_ok   <- "Rhat_split"   %in% names(dt) && any(is.finite(dt$Rhat_split))
    r_classic_ok <- "Rhat_classic" %in% names(dt) && any(is.finite(dt$Rhat_classic))
    if (prefer_split && r_split_ok) return(as.numeric(dt$Rhat_split))
    if (r_classic_ok)               return(as.numeric(dt$Rhat_classic))
    if (r_split_ok)                 return(as.numeric(dt$Rhat_split))
    stop("plot_convergence_rhat_sampled_fast: no finite Rhat (classic/split) found in 'diag_tbl'.")
  }
  Rhat <- pick_rhat(diag_tbl, prefer_split = prefer_split)

  df <- data.frame(
    target = diag_tbl$target,
    Rhat   = Rhat,
    stringsAsFactors = FALSE
  )
  df$Family <- .family_from_fast(df$target)

  # --- Strict samplers-only filter
  if (isTRUE(sampled_only)) {
    if (is.null(conf.mcmc) || is.null(samples_ml)) {
      stop("plot_convergence_rhat_sampled_fast: sampled_only=TRUE requires 'conf.mcmc' and 'samples_ml'.")
    }
    sampler_targets <- .derive_sampler_targets_fast(conf.mcmc, samples_ml)
    df <- df[df$target %in% sampler_targets, , drop = FALSE]
    if (nrow(df) == 0L) {
      stop("plot_convergence_rhat_sampled_fast: no sampled nodes found in 'diag_tbl' after strict filtering.")
    }
  }

  # --- Keep finite R-hat only
  dff <- df[is.finite(df$Rhat), , drop = FALSE]
  if (nrow(dff) == 0L) stop("plot_convergence_rhat_sampled_fast: no finite R-hat after filtering.")

  # --- Family summary and top-k worst nodes
  fam_rhat <- stats::aggregate(Rhat ~ Family, data = dff, FUN = function(x) stats::median(x, na.rm = TRUE))
  names(fam_rhat)[2] <- "Rhat_median"

  dff_ord <- dff[order(dff$Rhat, decreasing = TRUE), , drop = FALSE]
  k <- max(1L, min(as.integer(top_k), nrow(dff_ord)))
  top_tbl <- utils::head(dff_ord, k)
  top_tbl$target <- factor(top_tbl$target, levels = top_tbl$target)

  # --- Plots
  p_fam <- ggplot2::ggplot(fam_rhat, ggplot2::aes(x = reorder(Family, Rhat_median), y = Rhat_median, fill = Family)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    ggplot2::labs(title = sprintf("Median R\u0302 by Family (samplers only) — threshold %.3f", threshold),
                  x = "Family", y = "Median R\u0302") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")

  p_top <- ggplot2::ggplot(top_tbl, ggplot2::aes(x = target, y = Rhat, fill = Family)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    ggplot2::labs(title = sprintf("Top-%d worst nodes by R\u0302 (samplers only)", k),
                  x = "Target", y = "R\u0302") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # --- Save if requested
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    thr_tag <- gsub("\\.", "_", as.character(threshold))
    ggplot2::ggsave(file.path(out_dir, sprintf("Rhat_median_by_family__samplers_only_thr_%s.png", thr_tag)),
                    p_fam, width = 12, height = 8, dpi = 300)
    ggplot2::ggsave(file.path(out_dir, sprintf("TopK_worst_Rhat__samplers_only_thr_%s.png", thr_tag)),
                    p_top, width = 14, height = 8, dpi = 300)
    utils::write.table(top_tbl,  file = file.path(out_dir, "Rhat_top_worst__samplers_only.csv"),
                       sep = ",", row.names = FALSE)
    utils::write.table(fam_rhat, file = file.path(out_dir, "Rhat_median_by_family__samplers_only.csv"),
                       sep = ",", row.names = FALSE)
  }

  invisible(list(
    p_family = p_fam,
    p_top    = p_top,
    table_family = fam_rhat,
    table_top    = top_tbl,
    threshold    = threshold
  ))
}


## ---- Internal functions (20) ----


#' @keywords internal
.derive_sampler_targets_fast <- function(conf.mcmc, samples_ml) {
  stopifnot(!is.null(conf.mcmc), !is.null(samples_ml))
  samps <- conf.mcmc$getSamplers()
  raw_targets <- unique(unlist(lapply(samps, function(s) s$target), use.names = FALSE))
  cn <- colnames(as.data.frame(samples_ml[[1]]))
  intersect(raw_targets, cn)  # keep only targets that exist in the MCMC samples
}

#' @keywords internal
.family_from_fast <- function(x) sub("\\[.*\\]", "", x)

  mk_plot_all <- function(metric, ylab, file) {
    if (!nrow(Agg_top)) return(invisible(NULL))
    data.table::setorder(Agg_top, -get(metric))
    Agg_top[, key_f := factor(key, levels = rev(unique(key)))]
    p <- ggplot2::ggplot(Agg_top, ggplot2::aes(x = key_f, y = .data[[metric]], fill = .data$`__sampler__`)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.85)) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~ `.__step__`, scales = "free_y") +
      ggplot2::labs(
        x = NULL, y = ylab, fill = "Sampler",
        title = sprintf("%s – strategy comparison (per = %s, top_by = %s, k = %d)",
                        metric, per, top_by, as.integer(top_k))
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold"),
                     legend.position = "bottom")
    .ggsave(p, file.path(out_dir, file))
    invisible(p)
  }


      mk_plot_step <- function(metric, ylab, file) {
        if (!nrow(sub)) return(invisible(NULL))
        data.table::setorder(sub, -get(metric))
        sub[, key_f := factor(key, levels = rev(unique(key)))]
        p <- ggplot2::ggplot(sub, ggplot2::aes(x = key_f, y = .data[[metric]], fill = .data$`__sampler__`)) +
          ggplot2::geom_col() +
          ggplot2::coord_flip() +
          ggplot2::labs(x = NULL, y = ylab, fill = "Sampler",
                        title = sprintf("%s – %s", metric, st)) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                         plot.title = ggplot2::element_text(face = "bold"),
                         legend.position = "bottom")
        .ggsave(p, file.path(sdir, file))
        invisible(p)
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


  pick_rhat <- function(dt, prefer_split = TRUE) {
    r_split_ok   <- "Rhat_split"   %in% names(dt) && any(is.finite(dt$Rhat_split))
    r_classic_ok <- "Rhat_classic" %in% names(dt) && any(is.finite(dt$Rhat_classic))
    if (prefer_split && r_split_ok) return(as.numeric(dt$Rhat_split))
    if (r_classic_ok)               return(as.numeric(dt$Rhat_classic))
    if (r_split_ok)                 return(as.numeric(dt$Rhat_split))
    stop("plot_convergence_rhat_sampled_fast: no finite Rhat (classic/split) found in 'diag_tbl'.")
  }


  pick_vec <- function(primary, fallback) {
    if (primary %in% names(diag_tbl) && any(is.finite(diag_tbl[[primary]]))) {
      as.numeric(diag_tbl[[primary]])
    } else if (fallback %in% names(diag_tbl) && any(is.finite(diag_tbl[[fallback]]))) {
      as.numeric(diag_tbl[[fallback]])
    } else {
      rep(NA_real_, nrow(diag_tbl))
    }
  }


#' Plot strategy comparisons from test_strategy_family_fast results (fast path)
#'
#' @description
#' Fast plotting utility to compare sampler strategies (scalar and block plans,
#' including optional full-model HMC) using median efficiency metrics computed
#' per target or per family. Designed to work with objects returned by
#' [test_strategy_family_fast()] and [configure_hmc_safely()].
#'
#' @details
#' The function consumes `res` objects that contain:
#' - `baseline$diag_tbl` **or** (`baseline$samples`, `baseline$runtime_s`)
#' - a list of `steps` with per-step diagnostics under `steps[[i]]$res$dg`
#' - optionally, `mode == "HMC_full"` with `hmc$diag_tbl`.
#'
#' When diagnostics are missing but samples/runtime are provided, the function
#' prefers `compute_diag_from_mcmc_vect()` (vectorized) and falls back to
#' `compute_diag_from_mcmc()` if needed (both assumed available in your package).
#'
#' Filtering removes likelihood-like targets (`logLik`, `log_lik`, `logdens`,
#' `lpdf`) and internal/generated nodes (`lifted_`, `logProb_`). If the `Family`
#' column is absent, it is derived as the root prefix of `target`.
#'
#' Plots include:
#' - Global bar charts across steps for AE (ESS/iteration), CE (ESS/second),
#'   and \eqn{\hat{R}} (max per key);
#' - Per-step bar charts for the same metrics.
#'
#' @param res A result object from [test_strategy_family_fast()] or
#'   [configure_hmc_safely()] containing baseline/steps/HMC diagnostics.
#' @param out_dir Output directory for plots (created if missing).
#' @param per Aggregation level: `"target"` or `"family"`.
#' @param top_k Integer; number of keys (targets or families) to keep for plotting
#'   based on `top_by`. Use a large number to keep all.
#' @param top_by Selection metric for top-k: `"CE"` or `"AE"`.
#' @param export_csv Logical; if `TRUE`, writes a CSV summary at `out_dir`.
#' @param width,height,dpi Plot export parameters passed to `ggsave()`.
#'
#' @return (Invisibly) a list with:
#' \describe{
#'   \item{data}{`data.table` of aggregated metrics used for plotting.}
#'   \item{out_dir}{Output directory path.}
#'   \item{plots_all}{A list of ggplot objects for AE/CE/Rhat global charts.}
#' }
#'
#' @section Metrics:
#' - AE = median(ESS per iteration) over the aggregation key;
#' - CE = median(ESS per second) over the aggregation key;
#' - Rhat = max(\eqn{\hat{R}}) over the aggregation key.
#'
#' @seealso [test_strategy_family_fast()], [configure_hmc_safely()]
#'
#' @importFrom data.table as.data.table rbindlist setnames setorder fwrite
#' @importFrom ggplot2 ggplot aes geom_col coord_flip facet_wrap labs theme_minimal
#'   theme element_blank element_text ggsave
#'
#'
#' @examples
#' \dontrun{
#' res <- test_strategy_family_fast(build_fn = build_M, nbot = 2, try_hmc = TRUE)
#' plot_strategies_from_test_result_fast(res, per = "family", top_k = 30, top_by = "CE")
#' }
#'  @export
plot_strategies_from_test_result_fast <- function(
    res,
    out_dir = "outputs/strategies_plots",
    per = c("target", "family"),
    top_k = 40L,
    top_by = c("CE","AE"),
    export_rds = TRUE,
    width = 12, height = 7, dpi = 300
) {
  ## ------------------------------------------------------------
  ## 0) Préliminaires & vérifs
  ## ------------------------------------------------------------
  per    <- match.arg(per)
  top_by <- match.arg(top_by)

  pkgs <- c("data.table","ggplot2")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing required packages: ", paste(miss, collapse = ", "))

  `%||%` <- function(x, y) if (is.null(x)) y else x
  is_ll  <- function(x) grepl("^logLik(\\[.*\\])?$|log_?lik|logdens|lpdf", x, perl = TRUE, ignore.case = TRUE)
  root_of <- function(x) sub("\\[.*", "", x)
  .mk_dir <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE) }
  .ggsave <- function(p, f) try(ggplot2::ggsave(filename = f, plot = p, width = width, height = height, dpi = dpi), silent = TRUE)

  .mk_dir(out_dir)

  ## ------------------------------------------------------------
  ## 1) Extraction / harmonisation des diagnostics
  ## ------------------------------------------------------------
  pieces <- list()

  # Cas : sortie configure_hmc_safely_bis()
  if (!is.null(res$samples) && !is.null(res$runtime_s)) {
    if (!is.null(res$diag_tbl)) {
      di <- res$diag_tbl
    } else {
      diag_fun <- if (exists("compute_diag_from_mcmc_vect", mode = "function")) {
        compute_diag_from_mcmc_vect
      } else if (exists("compute_diag_from_mcmc", mode = "function")) {
        compute_diag_from_mcmc
      } else {
        stop("No diagnostic function found: compute_diag_from_mcmc_vect or compute_diag_from_mcmc.")
      }
      di <- diag_fun(res$samples, runtime_s = res$runtime_s, compute_rhat = "both", ess_for = "both")
    }
    di$.__step__    <- "HMC_full"
    di$.__sampler__ <- "NUTS_full"
    pieces[[length(pieces)+1]] <- di
  }

  # Cas : présence d'autres structures (baseline, steps)
  if (!is.null(res$baseline) && !is.null(res$baseline$diag_tbl)) {
    di <- res$baseline$diag_tbl
    di$.__step__    <- "baseline"
    di$.__sampler__ <- "baseline"
    pieces[[length(pieces)+1]] <- di
  }

  if (!is.null(res$steps) && length(res$steps)) {
    for (st in res$steps) {
      di <- try(st$res$dg, silent = TRUE)
      if (!inherits(di, "try-error") && !is.null(di)) {
        di$.__step__    <- as.character(st$level %||% "step")
        di$.__sampler__ <- as.character(st$sampler %||% "NA")
        pieces[[length(pieces)+1]] <- di
      }
    }
  }

  if (!length(pieces))
    stop("No diagnostic tables (`diag_tbl`) detected in input object `res`.")

  ## ------------------------------------------------------------
  ## 2) Fusion et nettoyage
  ## ------------------------------------------------------------
  needed_primary <- c("target","AE_ESS_per_it","ESS_per_sec","Rhat",".__step__","Family")
  needed_sampler <- c("__sampler__", ".__sampler__")

  DT <- data.table::rbindlist(
    lapply(pieces, function(d) {
      cols <- intersect(colnames(d), c(needed_primary, needed_sampler))
      data.table::as.data.table(d[, cols, drop = FALSE])
    }),
    use.names = TRUE, fill = TRUE
  )

  if (!"__sampler__" %in% names(DT) && ".__sampler__" %in% names(DT))
    data.table::setnames(DT, ".__sampler__", "__sampler__")
  if (!"__sampler__" %in% names(DT))
    DT[, `__sampler__` := "unknown"]

  DT <- DT[!is_ll(target)]
  DT <- DT[!grepl("^lifted_|^logProb_", target)]
  if (!"Family" %in% names(DT)) DT[, Family := root_of(target)]

  DT[, AE := AE_ESS_per_it]
  DT[, CE := ESS_per_sec]
  if (!"AE" %in% names(DT)) DT[, AE := NA_real_]
  if (!"CE" %in% names(DT)) DT[, CE := NA_real_]
  if (!"Rhat" %in% names(DT)) DT[, Rhat := NA_real_]

  ## ------------------------------------------------------------
  ## 3) Agrégation (family/target)
  ## ------------------------------------------------------------
  key_col <- if (per == "family") "Family" else "target"

  Agg <- DT[
    , .(
      AE   = suppressWarnings(stats::median(AE, na.rm = TRUE)),
      CE   = suppressWarnings(stats::median(CE, na.rm = TRUE)),
      Rhat = suppressWarnings(max(Rhat, na.rm = TRUE))
    ),
    by = c("__sampler__", ".__step__", key_col)
  ]

  Agg[is.infinite(Rhat), Rhat := NA_real_]
  data.table::setnames(Agg, key_col, "key")

  ## ------------------------------------------------------------
  ## 4) Sélection top-k
  ## ------------------------------------------------------------
  if (!is.finite(top_k) || top_k <= 0L) top_k <- nrow(Agg)
  if (all(is.na(Agg[[top_by]]))) {
    keys_keep <- unique(Agg$key)[1:min(length(unique(Agg$key)), as.integer(top_k))]
  } else {
    data.table::setorder(Agg, -get(top_by))
    keys_keep <- unique(Agg[1:min(.N, as.integer(top_k)), key])
  }
  Agg_top <- Agg[key %in% keys_keep]

  ## ------------------------------------------------------------
  ## 5) Graphiques globaux
  ## ------------------------------------------------------------
  mk_plot_all <- function(metric, ylab, file) {
    if (!nrow(Agg_top)) return(invisible(NULL))
    data.table::setorder(Agg_top, -get(metric))
    Agg_top[, key_f := factor(key, levels = rev(unique(key)))]
    p <- ggplot2::ggplot(Agg_top, ggplot2::aes(x = key_f, y = .data[[metric]], fill = .data$`__sampler__`)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.85)) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~ `.__step__`, scales = "free_y") +
      ggplot2::labs(
        x = NULL, y = ylab, fill = "Sampler",
        title = sprintf("%s – strategy comparison (per = %s, top_by = %s, k = %d)",
                        metric, per, top_by, as.integer(top_k))
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold"),
                     legend.position = "bottom")
    .ggsave(p, file.path(out_dir, file))
    invisible(p)
  }

  pAE   <- mk_plot_all("AE",   "AE (ESS / iteration)", "strategy_bar_AE.png")
  pCE   <- mk_plot_all("CE",   "CE (ESS / second)",    "strategy_bar_CE.png")
  pRhat <- mk_plot_all("Rhat", "R\u0302 (max)",        "strategy_bar_Rhat.png")

  ## ------------------------------------------------------------
  ## 6) Graphiques par step
  ## ------------------------------------------------------------
  if (nrow(Agg_top)) {
    for (st in unique(Agg_top$`.__step__`)) {
      sub <- Agg_top[`.__step__` == st]
      sdir <- file.path(out_dir, paste0("step_", gsub("[^A-Za-z0-9_]+","_", st)))
      .mk_dir(sdir)
      mk_plot_step <- function(metric, ylab, file) {
        if (!nrow(sub)) return(invisible(NULL))
        data.table::setorder(sub, -get(metric))
        sub[, key_f := factor(key, levels = rev(unique(key)))]
        p <- ggplot2::ggplot(sub, ggplot2::aes(x = key_f, y = .data[[metric]], fill = .data$`__sampler__`)) +
          ggplot2::geom_col() +
          ggplot2::coord_flip() +
          ggplot2::labs(x = NULL, y = ylab, fill = "Sampler",
                        title = sprintf("%s – %s", metric, st)) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                         plot.title = ggplot2::element_text(face = "bold"),
                         legend.position = "bottom")
        .ggsave(p, file.path(sdir, file))
        invisible(p)
      }
      mk_plot_step("AE",   "AE (ESS / iteration)", file = "bar_AE.png")
      mk_plot_step("CE",   "CE (ESS / second)",    file = "bar_CE.png")
      mk_plot_step("Rhat", "R\u0302 (max)",        file = "bar_Rhat.png")
    }
  }

  ## ------------------------------------------------------------
  ## 7) Export RDS
  ## ------------------------------------------------------------
  if (isTRUE(export_rds) && nrow(Agg_top)) {
    f_rds <- file.path(out_dir, "strategy_summary_top.rds")
    saveRDS(Agg_top, f_rds)
    message("[IO] Saved summary RDS: ", f_rds)
  }

  ## ------------------------------------------------------------
  ## 8) Retour
  ## ------------------------------------------------------------
  invisible(list(
    data      = Agg_top,
    out_dir   = out_dir,
    plots_all = list(AE = pAE, CE = pCE, Rhat = pRhat)
  ))
}


#'  @export
plot_strategies_from_test_result_hmc_fast <- function(
    res,
    out_dir        = NULL,
    per            = c("target", "family"),
    top_k          = 40,
    top_by         = c("CE", "AE", "Rhat"),
    show_rhat_label = TRUE
) {
  if (is.null(res$diag_tbl)) {
    stop("L'objet 'res' ne contient pas de diag_tbl.")
  }
  diag_tbl <- res$diag_tbl

  per    <- match.arg(per)
  top_by <- match.arg(top_by)

  # --- Vérif colonnes essentielles (on suppose enrich_hmc_diag_tbl_for_plots déjà appelé) ---
  needed <- c("AE_ESS_per_it", "CE_ESS_per_s", "Rhat")
  missing <- setdiff(needed, names(diag_tbl))
  if (length(missing) > 0L) {
    stop("Colonnes manquantes dans diag_tbl : ", paste(missing, collapse = ", "),
         ". As-tu bien appelé enrich_hmc_diag_tbl_for_plots() ?")
  }

  # --- choix de la metric de ranking ---
  metric_col <- switch(
    top_by,
    "CE"   = "CE_ESS_per_s",
    "AE"   = "AE_ESS_per_it",
    "Rhat" = "Rhat"
  )

  # Nettoyage NA/Inf/NaN
  metric <- diag_tbl[[metric_col]]
  metric[!is.finite(metric)] <- NA_real_
  diag_tbl[[metric_col]] <- metric

  # --- catégorisation de Rhat pour la couleur ---
  rhat <- diag_tbl$Rhat
  rhat_cat <- rep(NA_character_, length(rhat))

  rhat_cat[is.na(rhat)]                  <- "Rhat NA"
  rhat_cat[!is.na(rhat) & rhat < 1.05]  <- "< 1.05"
  rhat_cat[!is.na(rhat) & rhat >= 1.05 & rhat < 1.1] <- "1.05–1.10"
  rhat_cat[!is.na(rhat) & rhat >= 1.1]  <- ">= 1.10"

  diag_tbl$rhat_cat <- factor(
    rhat_cat,
    levels = c("< 1.05", "1.05–1.10", ">= 1.10", "Rhat NA")
  )

  # --- Construction du tableau pour le plot ---
  if (per == "target") {

    if (!"target" %in% names(diag_tbl)) {
      stop("diag_tbl ne contient pas de colonne 'target'.")
    }

    plot_tbl <- diag_tbl

    # filtre NA sur la metric choisie
    ok <- !is.na(plot_tbl[[metric_col]])
    plot_tbl <- plot_tbl[ok, , drop = FALSE]

    if (nrow(plot_tbl) == 0L) {
      stop("Aucune ligne valide (métrique NA partout) pour per = 'target'.")
    }

    # tri décroissant
    o <- order(-plot_tbl[[metric_col]])
    plot_tbl <- plot_tbl[o, , drop = FALSE]

    # top_k
    if (nrow(plot_tbl) > top_k) {
      plot_tbl <- plot_tbl[seq_len(top_k), , drop = FALSE]
    }

    # ordre des cibles sur le graphique
    plot_tbl$target <- factor(
      plot_tbl$target,
      levels = rev(plot_tbl$target[order(plot_tbl[[metric_col]])])
    )

    title_txt <- sprintf("Top %d targets par %s (HMC/NUTS)", top_k, metric_col)

    p <- ggplot2::ggplot(
      plot_tbl,
      ggplot2::aes_string(x = "target", y = metric_col, fill = "rhat_cat")
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = "Target",
        y = metric_col,
        fill = "Rhat",
        title = title_txt,
        subtitle = "Couleur = classe de Rhat"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 6)
      )

    if (show_rhat_label) {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes_string(
            x = "target",
            y = metric_col,
            label = "sprintf('Rhat=%.2f', Rhat)"
          ),
          hjust = -0.05,
          size  = 2.5
        ) +
        ggplot2::expand_limits(y = max(plot_tbl[[metric_col]], na.rm = TRUE) * 1.1)
    }

  } else if (per == "family") {

    if (!"Family" %in% names(diag_tbl)) {
      stop("per='family' demandé mais diag_tbl ne contient pas 'Family'.")
    }

    # agrégation par famille avec base R
    ae_by_family <- aggregate(
      diag_tbl$AE_ESS_per_it,
      by = list(Family = diag_tbl$Family),
      FUN = function(z) mean(z, na.rm = TRUE)
    )
    names(ae_by_family)[names(ae_by_family) == "x"] <- "AE_ESS_per_it"

    ce_by_family <- aggregate(
      diag_tbl$CE_ESS_per_s,
      by = list(Family = diag_tbl$Family),
      FUN = function(z) mean(z, na.rm = TRUE)
    )
    names(ce_by_family)[names(ce_by_family) == "x"] <- "CE_ESS_per_s"

    rhat_by_family <- aggregate(
      diag_tbl$Rhat,
      by = list(Family = diag_tbl$Family),
      FUN = function(z) mean(z, na.rm = TRUE)
    )
    names(rhat_by_family)[names(rhat_by_family) == "x"] <- "Rhat"

    plot_tbl <- Reduce(function(x, y) merge(x, y, by = "Family", all = TRUE),
                       list(ae_by_family, ce_by_family, rhat_by_family))

    # recalc rhat_cat par famille
    rhat_f <- plot_tbl$Rhat
    rhat_cat_f <- rep(NA_character_, length(rhat_f))
    rhat_cat_f[is.na(rhat_f)]                        <- "Rhat NA"
    rhat_cat_f[!is.na(rhat_f) & rhat_f < 1.05]       <- "< 1.05"
    rhat_cat_f[!is.na(rhat_f) & rhat_f >= 1.05 & rhat_f < 1.1] <- "1.05–1.10"
    rhat_cat_f[!is.na(rhat_f) & rhat_f >= 1.1]       <- ">= 1.10"

    plot_tbl$rhat_cat <- factor(
      rhat_cat_f,
      levels = c("< 1.05", "1.05–1.10", ">= 1.10", "Rhat NA")
    )

    # filtre NA sur la metric choisie
    ok <- !is.na(plot_tbl[[metric_col]])
    plot_tbl <- plot_tbl[ok, , drop = FALSE]

    if (nrow(plot_tbl) == 0L) {
      stop("Aucune famille valide (métrique NA partout) pour per = 'family'.")
    }

    # tri décroissant
    o <- order(-plot_tbl[[metric_col]])
    plot_tbl <- plot_tbl[o, , drop = FALSE]

    if (nrow(plot_tbl) > top_k) {
      plot_tbl <- plot_tbl[seq_len(top_k), , drop = FALSE]
    }

    plot_tbl$Family <- factor(
      plot_tbl$Family,
      levels = rev(plot_tbl$Family[order(plot_tbl[[metric_col]])])
    )

    title_txt <- sprintf("Top %d families par %s (HMC/NUTS)", top_k, metric_col)

    p <- ggplot2::ggplot(
      plot_tbl,
      ggplot2::aes_string(x = "Family", y = metric_col, fill = "rhat_cat")
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = "Family",
        y = metric_col,
        fill = "Rhat",
        title = title_txt,
        subtitle = "Couleur = classe de Rhat (moyenne par famille)"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 7)
      )

    if (show_rhat_label) {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes_string(
            x = "Family",
            y = metric_col,
            label = "sprintf('Rhat=%.2f', Rhat)"
          ),
          hjust = -0.05,
          size  = 2.5
        ) +
        ggplot2::expand_limits(y = max(plot_tbl[[metric_col]], na.rm = TRUE) * 1.1)
    }
  } else {
    stop("Valeur de 'per' non gérée : ", per)
  }

  # --- export éventuel ---
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    fname <- file.path(
      out_dir,
      sprintf("diag_hmc_%s_top%d_by_%s.png", per, top_k, top_by)
    )

    ggplot2::ggsave(
      filename = fname,
      plot     = p,
      width    = 8,
      height   = 6,
      dpi      = 300
    )
  }

  invisible(list(plot = p, data = plot_tbl))
}


  root_of <- function(x) sub("\\[.*", "", x)
  .mk_dir <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE) }


  save_fig <- function(p, base, w = 8, h = 5.2, dpi = 180) {
    ggplot2::ggsave(file.path(out_dir, paste0(base, ".pdf")), p, width = w, height = h)
    ggplot2::ggsave(file.path(out_dir, paste0(base, ".png")), p, width = w, height = h, dpi = dpi)
  }


