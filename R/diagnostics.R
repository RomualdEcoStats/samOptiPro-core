# diagnostics.R -- target-level diagnostics (time + step proxy)
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

# ----------------------------------------------------------------------
# Step proxy: sd(diff(chain)) on columns corresponding to a pattern
# ----------------------------------------------------------------------

#' Step proxy: sd(diff(chain)) on columns matching a pattern
#'
#' Computes, for each parameter whose name matches `pattern`, the standard
#' deviation of successive differences (a simple step-size proxy).
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

# --------------------------------------------------------------
# Sampler-level time profiling (compiled) aligned with conf$...
# --------------------------------------------------------------
#' Profile sampler times per node/family
#'
#' @param cmodel A compiled \code{nimbleModel}.
#' @param mcmc_conf An MCMC configuration (e.g., \code{nimble::configureMCMC(model)}). If \code{NULL}, an automatic attempt is made.
#' @param include_data Logical; include data nodes (default \code{FALSE}).
#' @param removed_nodes Character vector; nodes to exclude (e.g., with ESS = 0).
#' @param ignore_patterns Character vector of regex patterns to ignore.
#' @param plot_dependencies Logical; if \code{TRUE}, makes a ggplot barplot of medians per family.
#' @param output_dir Character or \code{NULL}; directory where CSV/PNG are exported. If \code{NULL}, nothing is written.
#' @param save_csv Logical; if \code{TRUE}, writes \code{dependencies_per_node.csv} to \code{output_dir}.
#' @param node_of_interest Character or \code{NULL}; optionally filter downstream dependencies for this node.
#' @param np Numeric in (0,1]; kept for API compatibility (default \code{0.10}).
#' @param ... Additional arguments passed to \code{nimble::configureMCMC()}.
#'
#' @return A list with (at least) data frames per node/family and (optionally) ggplot objects.
#' @examples
#' \dontrun{
#'   # p <- profile_sampler_times(cmodel, mcmc_conf)
#' }
#' @export
#' @keywords internal
profile_sampler_times <- function(conf, cmodel, niter = 5e4) {
  m  <- nimble::buildMCMC(conf)
  cm <- nimble::compileNimble(m, project = cmodel, resetFunctions = TRUE)
  cm$run(niter, time = TRUE)
  tm <- cm$getTimes()
  ns <- length(conf$getSamplers())
  if (length(tm) > ns) tm <- tm[seq_len(ns)]
  if (length(tm) < ns) tm <- c(tm, rep(NA_real_, ns - length(tm)))
  tm
}

# --------------------------------------------------------------
# Target-level diagnostics (time + step proxy when samples provided)
# --------------------------------------------------------------

#' Target-level diagnostics (time + optional step-size proxy)
#'
#' Builds a robust MCMC configuration, profiles per-sampler time with
#' \code{\link{profile_sampler_times}}, maps samplers to their targets, and
#' optionally computes a step-size proxy per target via \code{\link{proxy_step_sd}}
#' if \code{samples} are provided.
#'
#' Note: uses a tolerant default monitor set followed by a robust configuration
#' (no \code{monitors=} passed directly to \code{configureMCMC}).
#'
#' @param build_fn A builder function returning \code{list(model=, cmodel=, monitors=?)}.
#' @param opts     A list of options (e.g., from \code{samOptiPro_options()}).
#' @param niter    Integer; iterations for time profiling (defaults to \code{opts$time_profile_n}).
#' @param samples  Optional \code{coda::mcmc.list} used to compute the step proxy.
#'
#' @return A \code{data.frame} with columns \code{target}, \code{type}, \code{time_s}, \code{step_sd}.
#' @export
#' @keywords internal
diagnostics_by_target <- function(build_fn, opts = samOptiPro_options(),
                                  niter = opts$time_profile_n, samples = NULL) {
  built <- build_fn()

  # Monitors (tolerant) -> then ROBUST configuration
  mons <- default_monitors(built$model, opts)
  mons <- ensure_monitors_exist(built$model, mons)
  conf <- .configure_with_monitors(built$model, monitors = mons)

  time_s <- profile_sampler_times(conf, built$cmodel, niter = niter)

  df_map <- sampler_df_from_conf(conf)  # name, type, target, scale (fallback)
  df_map$time_s <- time_s

  if (!is.null(samples)) {
    pats <- unique(gsub("\\[.*$", "[", df_map$target))
    step_df <- do.call(rbind, lapply(pats, function(p) {
      vals <- proxy_step_sd(samples, p)
      if (!length(vals)) return(NULL)
      data.frame(target = names(vals), step_sd = as.numeric(vals), stringsAsFactors = FALSE)
    }))
    if (!is.null(step_df)) df_map <- merge(df_map, step_df, by = "target", all.x = TRUE, sort = FALSE)
  } else {
    df_map$step_sd <- NA_real_
  }

  df_map[, c("target","type","time_s","step_sd")]
}

# ----------------------------------------------------------------------
# Full execution: structure diagnostics & HMC/NUTS smoke test
# ----------------------------------------------------------------------

#' Run structural diagnostics and (optional) HMC/NUTS smoke test
#'
#' Inspects a NIMBLE model produced by \code{build_fn()} to:
#' \itemize{
#'   \item count stochastic vs. deterministic nodes (optionally including data),
#'   \item parse the model code to detect non-differentiable functions and BUGS-style truncation,
#'   \item rebuild a per-node table with distribution, support, and bounds,
#'   \item tag \strong{HMC showstoppers} (e.g., discrete latents, simplex constraints, truncation,
#'         non-differentiable deterministic ops feeding latents),
#'   \item optionally attempt a short HMC/NUTS run (via \pkg{nimbleHMC}) to check practical feasibility.
#' }
#'
#' Compared to \code{diagnose_model_structure()}, this function does not filter nodes by
#' ignore patterns or explicit removals and does not compute dependency fan-out.
#' Its emphasis is \emph{HMC suitability} rather than structural load.
#'
#' @param build_fn     A zero-arg function returning a list with at least
#'   \code{model} (\code{nimbleModel}); optional elements: \code{cmodel}, \code{monitors},
#'   \code{code_text} (or \code{code}) used to strengthen non-diff detection.
#' @param include_data Logical; include data nodes in counts and scans. Default \code{FALSE}.
#' @param try_hmc      Logical; if \code{TRUE}, attempt a brief HMC/NUTS run. Default \code{TRUE}.
#' @param niter        Integer; total iterations for the HMC test. Default \code{50L}.
#' @param nburnin      Integer; burn-in for the HMC test. Default \code{10L}.
#' @param seed         Integer; RNG seed for the HMC run. Default \code{1L}.
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
#' run_structure_and_hmc_test
#' @export
run_structure_and_hmc_test <- function(build_fn,
                                       include_data = FALSE,
                                       try_hmc     = TRUE,
                                       niter       = 50L,
                                       nburnin     = 10L,
                                       seed        = 1L) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
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

  # --- Normalize distribution names (alias -> canonical)
  normalize_dist <- function(d) {
    d <- tolower(trimws(d %||% ""))
    if (d == "") return(NA_character_)
    # Bernoulli / binomial / multinomial aliases
    d <- sub("^dbernoulli$", "dbern", d)
    d <- sub("^dbin$",       "dbinom", d)
    d <- sub("^dmulti$",     "dmultinom", d)
    # NegBin (JAGS alias)
    d <- sub("^dnbinom$",    "dnegbin", d)
    # Dirichlet variants
    d <- sub("^ddirich(let)?$", "ddirichlet", d)
    d <- sub("^ddirch(let)?$",  "ddirichlet", d)
    # Hypergeometric variants
    d <- sub("^dhyperg(eometric)?$", "dhypergeom", d)
    d <- sub("^dhyper$",             "dhypergeom", d)
    # Poisson alias
    d <- sub("^dpois_rate$", "dpois", d)
    # Uniform / triangle / half
    d <- sub("^dtri(angle)?$",   "dtriangle",  d)
    d <- sub("^dhalfnorm(al)?$", "dhalfnorm",  d)
    d <- sub("^dhalfcauchy$",    "dhalfcauchy", d)
    d
  }

  # Families (canonical names after normalize_dist)
  discrete_dists <- c("dbern","dbinom","dcat","dmultinom","dgeom","dnegbin","dpois","dhypergeom","dinterval")
  simplex_dists  <- c("ddirichlet")
  bounded_support_dists <- c("dbeta","dunif","dtriangle","dhalfnorm","dhalfcauchy")

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

  cat("\n[STRUCTURE]\n")
  cat(sprintf("- # stochastic nodes   : %d\n", length(m$getNodeNames(stochOnly = TRUE, includeData = include_data))))
  cat(sprintf("- # deterministic nodes: %d\n",
              length(m$getNodeNames(stochOnly = FALSE, includeData = include_data)) -
                length(m$getNodeNames(stochOnly = TRUE, includeData = include_data))))

  # --- Non-differentiable operators (code scan + introspection)
  .sop_detect_nondiff_functions <- function(model,
                                            nondiff_candidates = c(
                                              "round","floor","ceiling","trunc",
                                              "abs","max","min","step","ifElse","ifelse","equals"
                                            ),
                                            code_text_override = NULL) {
    out <- character(0)
    code_txt <- ""
    if (is.character(code_text_override) && length(code_text_override)) {
      code_txt <- paste(code_text_override, collapse = "\n")
    } else {
      code_try <- try({
        code_txt <- paste0(paste(deparse(model$modelDef$code, width.cutoff = 500), collapse = "\n"), "\n")
      }, silent = TRUE)
      if (inherits(code_try, "try-error") || is.null(code_txt) || nchar(code_txt) == 0) {
        code_txt <- tryCatch(paste(utils::capture.output(print(model)), collapse = "\n"), error = function(e) "")
      }
    }
    if (nchar(code_txt)) {
      for (fn in nondiff_candidates) {
        pat1 <- paste0("\\b", fn, "\\s*\\(")
        pat2 <- paste0("(?i)\\b", fn, "\\s*\\(")
        if (grepl(pat1, code_txt, perl = TRUE) || grepl(pat2, code_txt, perl = TRUE)) out <- c(out, fn)
      }
    }
    # Introspection maps
    pull_char <- function(x) {
      x <- try(x, silent = TRUE)
      if (inherits(x, "try-error") || is.null(x)) return(character(0))
      as.character(unlist(x, use.names = FALSE))
    }
    maps <- try(model$modelDef$maps, silent = TRUE)
    if (!inherits(maps, "try-error") && !is.null(maps)) {
      out <- c(out, pull_char(maps$nodeFxnNames))
      out <- c(out, pull_char(maps$allFxnNames))
      out <- c(out, pull_char(maps$fxnNamesByNode))
      out <- c(out, pull_char(maps$deterministicFxnNames))
    }
    out <- unique(tolower(out))
    hits <- intersect(out, tolower(nondiff_candidates))
    list(hits = hits, code_txt = code_txt)
  }

  nd <- .sop_detect_nondiff_functions(m, code_text_override = code_override)
  nondiff_hits     <- nd$hits
  code_txt_scanned <- nd$code_txt

  # Distribution scan + BUGS truncation "T(a,b)"
  truncation_detected_bugst <- grepl("\\)\\s*T\\s*\\(", code_txt_scanned, perl = TRUE)
  dists_in_code_raw <- unique(tolower(unlist(regmatches(
    code_txt_scanned, gregexpr("\\bd[a-zA-Z_]+\\b", code_txt_scanned, perl = TRUE)
  ))))
  dists_in_code <- unique(vapply(dists_in_code_raw, normalize_dist, character(1)))

  diag_raw <- try(diagnose_model_structure(m, include_data = include_data), silent = TRUE)

  rebuild_nodes_df <- function(model, include_data) {
    all_nodes <- model$getNodeNames(includeData = include_data)
    stoch     <- model$getNodeNames(stochOnly = TRUE, includeData = include_data)

    is_data_vec <- stats::setNames(rep(FALSE, length(all_nodes)), all_nodes)
    for (n in all_nodes) is_data_vec[n] <- tryCatch(model$isData(n), error = function(e) FALSE)

    get_dist_raw  <- function(node) tolower(tryCatch(model$getDistribution(node), error = function(e) NA_character_))
    get_bound <- function(node, side) suppressWarnings(tryCatch(model$getBound(node, side), error = function(e) NA_real_))

    nodes_df <- do.call(rbind, lapply(all_nodes, function(node) {
      is_stoch <- node %in% stoch
      dist_raw <- if (is_stoch) get_dist_raw(node) else NA_character_
      dist     <- if (is_stoch) normalize_dist(dist_raw) else NA_character_
      lb       <- if (is_stoch) get_bound(node, "lower") else NA_real_
      ub       <- if (is_stoch) get_bound(node, "upper") else NA_real_
      sup      <- if (is_stoch) support_of(dist, lb, ub) else NA_character_
      data.frame(
        node      = node,
        var       = strip_index(node),
        is_stoch  = is_stoch,
        is_data   = is_data_vec[[node]],
        dist_raw  = dist_raw,
        dist      = dist,
        support   = sup,
        lower     = lb,
        upper     = ub,
        stringsAsFactors = FALSE
      )
    }))

    nodes_df$var      <- as.character(nodes_df$var)
    nodes_df$is_stoch <- to_logical(nodes_df$is_stoch)
    nodes_df$is_data  <- to_logical(nodes_df$is_data)

    fin_lb <- is.finite(as.numeric(nodes_df$lower))
    fin_ub <- is.finite(as.numeric(nodes_df$upper))

    # Mark truncation via bounds only when the natural support is unbounded continuous
    trunc_by_bounds <- nodes_df$is_stoch &
      (nodes_df$support == "continuous") &
      fin_lb & fin_ub

    nodes_df$is_truncated <- trunc_by_bounds | (nodes_df$is_stoch & isTRUE(truncation_detected_bugst))

    nodes_df$hmc_showstopper_reason <- NA_character_
    nodes_df$hmc_showstopper_reason[nodes_df$is_stoch & !nodes_df$is_data & nodes_df$support == "discrete"] <- "discrete-latent"
    nodes_df$hmc_showstopper_reason[nodes_df$is_stoch & !nodes_df$is_data & nodes_df$support == "simplex" &
                                      is.na(nodes_df$hmc_showstopper_reason)] <- "simplex-constraint"
    nodes_df$hmc_showstopper_reason[nodes_df$is_stoch & !nodes_df$is_data & nodes_df$is_truncated &
                                      is.na(nodes_df$hmc_showstopper_reason)] <- "truncation"
    if (length(nondiff_hits)) {
      nodes_df$hmc_showstopper_reason[!nodes_df$is_stoch & is.na(nodes_df$hmc_showstopper_reason)] <- "non-diff-deterministic-op"
    }

    # dims
    vars <- unique(nodes_df$var)
    dims <- lapply(vars, function(v) {
      vi <- tryCatch(model$getVarInfo(v), error = function(e) NULL)
      if (is.null(vi) || is.null(vi$nDim)) NA_integer_ else as.integer(vi$nDim)
    })
    names(dims) <- vars
    nodes_df$dims <- unname(vapply(nodes_df$var, function(v) {
      if (v %in% names(dims)) dims[[v]] else NA_integer_
    }, integer(1)))

    nodes_df
  }

  nodes <- if (!inherits(diag_raw, "try-error") && is.list(diag_raw) && is.data.frame(diag_raw$nodes)) {
    nds <- diag_raw$nodes
    if ("dist" %in% names(nds))     nds$dist     <- vapply(nds$dist, normalize_dist, character(1))
    if ("dist_raw" %in% names(nds)) nds$dist_raw <- tolower(nds$dist_raw)

    nds$var      <- as.character(nds$var)
    nds$is_stoch <- to_logical(nds$is_stoch)
    nds$is_data  <- to_logical(nds$is_data %||% FALSE)

    fin_lb <- is.finite(as.numeric(nds$lower %||% NA_real_))
    fin_ub <- is.finite(as.numeric(nds$upper %||% NA_real_))

    trunc_by_bounds <- nds$is_stoch &
      (mapply(support_of, nds$dist, nds$lower, nds$upper, SIMPLIFY = TRUE) == "continuous") &
      fin_lb & fin_ub

    nds$is_truncated <- trunc_by_bounds | (nds$is_stoch & isTRUE(truncation_detected_bugst))

    # Recompute support defensively
    nds$support <- mapply(support_of, nds$dist, nds$lower, nds$upper, SIMPLIFY = TRUE, USE.NAMES = FALSE)

    if (!"hmc_showstopper_reason" %in% names(nds)) nds$hmc_showstopper_reason <- NA_character_
    nds$hmc_showstopper_reason[nds$is_stoch & !nds$is_data & nds$support == "discrete"] <- "discrete-latent"
    nds$hmc_showstopper_reason[nds$is_stoch & !nds$is_data & nds$support == "simplex" &
                                 is.na(nds$hmc_showstopper_reason)] <- "simplex-constraint"
    nds$hmc_showstopper_reason[nds$is_stoch & !nds$is_data & nds$is_truncated &
                                 is.na(nds$hmc_showstopper_reason)] <- "truncation"
    if (length(nondiff_hits)) {
      nds$hmc_showstopper_reason[!nds$is_stoch & is.na(nds$hmc_showstopper_reason)] <- "non-diff-deterministic-op"
    }

    # Fill dims if missing
    if (!"dims" %in% names(nds)) {
      vars <- unique(nds$var)
      dims <- lapply(vars, function(v) {
        vi <- tryCatch(m$getVarInfo(v), error = function(e) NULL)
        if (is.null(vi) || is.null(vi$nDim)) NA_integer_ else as.integer(vi$nDim)
      })
      names(dims) <- vars
      nds$dims <- unname(vapply(nds$var, function(v) {
        if (v %in% names(dims)) dims[[v]] else NA_integer_
      }, integer(1)))
    }

    nds
  } else {
    rebuild_nodes_df(m, include_data)
  }

  # Scan summaries
  cat("\n[NON-DIFF INDICATORS]\n")
  cat(sprintf("- Non-diff functions detected: %s\n",
              if (length(nondiff_hits)) paste(nondiff_hits, collapse = ",") else "None"))
  cat(sprintf("- Distributions found in code : %s\n",
              if (length(dists_in_code)) paste(sort(unique(dists_in_code)), collapse = ", ") else "None"))

  bounded_latent_trunc <- any(nodes$is_stoch & !nodes$is_data &
                                (nodes$support == "continuous") &
                                is.finite(nodes$lower) & is.finite(nodes$upper), na.rm = TRUE)

  cat(sprintf("- BUGS-style truncation 'T(a,b)' spotted: %s\n",
              if (isTRUE(truncation_detected_bugst)) "Yes" else "No"))
  cat(sprintf("- Bounded latent nodes (implicit truncation via finite bounds on continuous support): %s\n",
              if (bounded_latent_trunc) "Yes" else "No"))

  # ---- Targeted tagging of non-diff deterministic ancestors feeding latents (for 'round')
  if (length(nondiff_hits)) {
    code_lines <- unlist(strsplit(code_txt_scanned, "\n", fixed = TRUE))
    round_lines <- grep("\\bround\\s*\\(", code_lines, value = TRUE, perl = TRUE)

    vars_in_lines <- unique(gsub("\\[.*?\\]$", "", unlist(regmatches(
      round_lines, gregexpr("\\b[A-Za-z_][A-Za-z0-9_]*(\\[[^\\]]+\\])?", round_lines, perl = TRUE)
    ))))

    vars_in_model <- unique(nodes$var)
    seed_vars <- intersect(vars_in_lines, vars_in_model)

    latents <- nodes$node[nodes$is_stoch & !nodes$is_data]
    det_anc_all <- unique(unlist(lapply(latents, function(nd)
      tryCatch(m$getDependencies(target = nd, upstream = TRUE, downstream = FALSE,
                                 includeData = FALSE, stochOnly = FALSE),
               error = function(e) character(0)))))

    det_anc_vars <- gsub("\\[.*?\\]$", "", det_anc_all)
    to_tag <- det_anc_all[det_anc_vars %in% seed_vars]

    is_det <- !nodes$is_stoch
    nodes$hmc_showstopper_reason[is_det] <- ifelse(
      nodes$node[is_det] %in% to_tag, "non-diff-deterministic-op", NA_character_
    )
  }

  cat("\n[DEBUG before summarise_showstoppers]\n")
  cat("has hmc_showstopper_reason? ", "hmc_showstopper_reason" %in% names(nodes), "\n")
  if ("hmc_showstopper_reason" %in% names(nodes)) {
    cat("non-NA reasons (non-data): ",
        sum(!nodes$is_data & !is.na(nodes$hmc_showstopper_reason)), "\n")
    print(utils::head(nodes[!nodes$is_data & !is.na(nodes$hmc_showstopper_reason),
                            c("node","is_stoch","hmc_showstopper_reason")], 10), row.names=FALSE)
  }

  # ---- Summarize HMC showstoppers
  summarise_showstoppers <- function(nodes_df) {
    ss <- nodes_df
    if (!"hmc_showstopper_reason" %in% names(ss)) ss$hmc_showstopper_reason <- NA_character_
    if (!"is_data" %in% names(ss))               ss$is_data               <- FALSE
    ss$is_data <- to_logical(ss$is_data)
    ss <- ss[!is.na(ss$hmc_showstopper_reason) & !ss$is_data, , drop = FALSE]

    if (!nrow(ss)) {
      return(list(summary = data.frame(reason = character(0), n_nodes = integer(0)),
                  examples = data.frame(reason = character(0), examples = character(0))))
    }

    tt <- as.data.frame(table(ss$hmc_showstopper_reason), stringsAsFactors = FALSE)
    names(tt) <- c("reason","n_nodes")
    tt$n_nodes <- as.integer(tt$n_nodes)
    tt <- tt[order(-tt$n_nodes, tt$reason), , drop = FALSE]

    ex <- do.call(rbind, lapply(split(ss, ss$hmc_showstopper_reason), function(d) {
      data.frame(
        reason   = as.character(d$hmc_showstopper_reason[1]),
        examples = paste(unique(utils::head(as.character(d$node), 5)), collapse = ","),
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
  cat(sprintf("\n- HMC globally feasible? %s\n", if (isTRUE(hmc_globally_ok)) "Yes" else "No"))

  # ------ Optional HMC/NUTS trial ------
  hmc_res <- list(ok = FALSE, error = "HMC test not requested.", details = NULL)
  if (isTRUE(try_hmc)) {
    cat("\n[HMC/NUTS SMOKE TEST]\n")
    has_hmc <- suppressWarnings(requireNamespace("nimbleHMC", quietly = TRUE))
    if (!has_hmc) {
      msg <- "nimbleHMC not available (not installed or not loaded)."
      cat(sprintf("- HMC feasible? No\n  reason: %s\n", msg))
      hmc_res <- list(ok = FALSE, error = msg, details = NULL)
    } else if (!isTRUE(hmc_globally_ok)) {
      reasons <- unique(stats::na.omit(nodes$hmc_showstopper_reason[!nodes$is_data]))
      msg <- sprintf("Model not suitable for HMC/NUTS (structural): %s", paste(reasons, collapse = ","))
      cat(sprintf("- HMC feasible? No\n  reason: %s\n", msg))
      hmc_res <- list(ok = FALSE, error = msg, details = NULL)
    } else {
      model_has_derivs <- .sop_supports_derivs(m)
      if (!isTRUE(model_has_derivs)) {
        msg <- "NUTS not attempted: the model was not built with buildDerivs=TRUE."
        cat(sprintf("- HMC feasible? No\n  reason: %s\n", msg))
        hmc_res <- list(ok = FALSE, error = msg, details = NULL)
      } else {
        base_conf <- try(nimble::configureMCMC(m), silent = TRUE)
        if (inherits(base_conf, "try-error")) {
          msg <- paste("configureMCMC failed:", as.character(base_conf))
          cat(sprintf("- HMC feasible? No\n  reason: %s\n", msg))
          hmc_res <- list(ok = FALSE, error = msg, details = NULL)
        } else {
          if (length(mons)) try(base_conf$addMonitors(mons), silent = TRUE)
          okH <- TRUE
          tryCatch({
            nimbleHMC::configureHMC(base_conf, model = m)
          }, error = function(e) {
            okH <<- FALSE
            cat("configureHMC failed: ", e$message, "\n", sep = "")
          })
          if (!okH) {
            msg <- "configureHMC failed."
            cat(sprintf("- HMC feasible? No\n  reason: %s\n", msg))
            hmc_res <- list(ok = FALSE, error = msg, details = NULL)
          } else {
            mcmc  <- try(nimble::buildMCMC(base_conf), silent = TRUE)
            if (inherits(mcmc, "try-error")) {
              msg <- paste("buildMCMC failed:", as.character(mcmc))
              cat(sprintf("- HMC feasible? No\n  reason: %s\n", msg))
              hmc_res <- list(ok = FALSE, error = msg, details = NULL)
            } else {
              cm_obj <- if (!is.null(cm)) cm else try(nimble::compileNimble(m), silent = TRUE)
              cmcmc  <- try(nimble::compileNimble(mcmc, project = cm_obj, resetFunctions = TRUE), silent = TRUE)
              if (inherits(cmcmc, "try-error")) {
                msg <- paste("compileNimble failed:", as.character(cmcmc))
                cat(sprintf("- HMC feasible? No\n  reason: %s\n", msg))
                hmc_res <- list(ok = FALSE, error = msg, details = NULL)
              } else {
                set.seed(seed)
                run <- try(nimble::runMCMC(cmcmc,
                                           niter   = as.integer(niter),
                                           nburnin = as.integer(nburnin),
                                           thin    = 1L, nchains = 1L,
                                           samplesAsCodaMCMC = FALSE,
                                           summary = FALSE, WAIC = FALSE),
                           silent = TRUE)
                if (inherits(run, "try-error")) {
                  msg <- paste("runMCMC failed:", as.character(run))
                  cat(sprintf("- HMC feasible? No\n  reason: %s\n", msg))
                  hmc_res <- list(ok = FALSE, error = msg, details = NULL)
                } else {
                  cat("- HMC feasible? Yes\n")
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
        functions_found           = sort(unique(nondiff_hits)),
        distributions_found       = sort(unique(dists_in_code)),
        distributions_found_raw   = sort(unique(dists_in_code_raw)),
        truncation_bugst          = truncation_detected_bugst,
        bounded_latent_trunc      = bounded_latent_trunc
      ),
      code_scan = code_txt_scanned,
      hmc_globally_ok = hmc_globally_ok
    ),
    hmc  = hmc_res
  ))
}


# --------------------------------------------------------------
# DIAGNOSTIC STRUCTURE (robuste + dependances + samplers, conf interne)
# --------------------------------------------------------------

#' Diagnose model structure (stochastic vs deterministic), dimensions, dependencies, and samplers
#'
#' @description
#' Analyze a NIMBLE model: stochastic/deterministic nodes, dimensions,
#' \strong{downstream dependencies} (ignoring artificial nodes), and \strong{samplers}
#' (type & targets) from an MCMC configuration \strong{built internally} via
#' \code{nimble::configureMCMC(model, ...)}.
#'
#' Artificial nodes \code{lifted_*}, \code{logProb_*}, and \code{removed_nodes} are
#' \strong{ignored} everywhere. Optionally: CSV export and ggplot2 barplot of
#' \strong{medians} of dependencies by stats::family(with \strong{labels = count}).
#'
#' @param model A \code{nimbleModel} (uncompiled).
#' @param include_data Logical; include data nodes (default \code{FALSE}).
#' @param removed_nodes Character vector; nodes to exclude (e.g., ESS = 0).
#' @param ignore_patterns Character vector; regex to ignore
#'   (default \code{c("^lifted_", "^logProb_")}).
#' @param plot_dependencies Logical; if \code{TRUE}, barplot (ggplot2) of medians per family.
#' @param output_dir Character or \code{NULL}; directory for CSV/PNG exports
#'   (\code{NULL} = no export).
#' @param save_csv Logical; if \code{TRUE}, write \code{dependencies_per_node.csv} to \code{output_dir}.
#' @param node_of_interest Character or \code{NULL}; node whose filtered downstream deps are requested.
#' @param np Numeric in (0,1]; proportion for \code{prop_worst = ceil(n_samplers * np)} (default \code{0.10}).
#' @param ... Additional arguments passed to \code{nimble::configureMCMC()} (e.g., \code{print = FALSE}, \code{useConjugacy = FALSE}).
#'
#' @return A list:
#' \preformatted{
#'   (stochastic_nodes, deterministic_nodes, dims, nodes = NULL,
#'    dependencies_df, dep_counts, dep_summary, plot_regular,
#'    node_of_interest_deps,
#'    samplers_df, n_samplers, prop_worst,
#'    mcmc_conf)  # configuration created internally
#' }
#'
# --------------------------------------------------------------
# DIAGNOSTIC STRUCTURE (dependencies + samplers + dual plots)
# --------------------------------------------------------------
#' Profile per-sampler elapsed times (uncompiled MCMC, R-level granularity)
#'
#' @description
#' Build an uncompiled MCMC (\code{nimble::buildMCMC(conf)}) and measure wall-clock
#' elapsed time spent in each sampler by instrumenting calls to \code{sampler$run()}
#' over a short run. This preserves \strong{per-sampler granularity}, which is lost in
#' compiled MCMC.
#'
#' @param model A \code{nimbleModel} (uncompiled), already initialized.
#' @param conf  An MCMC configuration created by \code{nimble::configureMCMC(model, ...)}.
#' @param niter Integer; total iterations to run (including burn-in).
#' @param burnin Integer; number of warm-up iterations to discard at the beginning.
#' @param thin Integer; thinning interval (affects only internal advancement here).
#' @param set_seed Integer or \code{NULL}; if not \code{NULL}, uses \code{set.seed(set_seed)} for reproducibility.
#' @param progress Logical; if \code{TRUE}, prints a simple progress indicator every 10\%.
#'
#' @return A numeric vector of length \code{length(conf$getSamplers())}, aligned with the
#'         sampler order returned by \code{conf$getSamplers()}. Units are seconds.
#'
#' Profile per-sampler elapsed times (uncompiled MCMC, R-level granularity)
#'
#' Builds an uncompiled MCMC and measures wall-clock time spent in each sampler
#' by instrumenting calls to each samplerFunction$run() over a short run.
#'
#' @param model \code{nimbleModel} (uncompiled), initialized.
#' @param conf  MCMC configuration from \code{nimble::configureMCMC(model, ...)}.
#' @param niter Total iterations including burn-in.
#' @param burnin Number of warm-up iterations before timing.
#' @param thin Thinning interval (\eqn{\ge} 1). If \code{> 1}, extra sweeps are advanced without timing.
#' @param set_seed Integer or \code{NULL}; if not \code{NULL}, \code{set.seed()} for reproducibility.
#' @param progress Logical; print simple 10\% progress.
#'
#' @return Numeric vector of length \code{length(conf$getSamplers())}, seconds per sampler (aligned by index).
#' profile_sampler_times
#' @export
#' @keywords internal
profile_sampler_times <- function(model,
                                  conf,
                                  niter    = 2000L,
                                  burnin   = 500L,
                                  thin     = 1L,
                                  set_seed = NULL,
                                  progress = FALSE) {
  stopifnot(niter > 0, burnin >= 0, thin >= 1)
  if (!requireNamespace("nimble", quietly = TRUE)) stop("Package 'nimble' is required.")
  if (!is.null(set_seed)) set.seed(as.integer(set_seed))

  # Build uncompiled MCMC to keep per-sampler visibility
  mcmcR <- nimble::buildMCMC(conf)

  # IMPORTANT: sampler functions live here (NOT in conf)
  samplerFns <- mcmcR$samplerFunctions
  ns <- length(samplerFns)
  if (ns == 0L) return(numeric(0))

  times <- numeric(ns)

  # Warm-up (burn-in) using whole sweeps
  if (burnin > 0L) mcmcR$run(burnin)

  n_eff <- as.integer(niter - burnin)
  if (n_eff <= 0L) return(times)

  pb_steps <- if (progress) unique(pmax(1L, floor(seq(0.1, 1.0, by = 0.1) * n_eff))) else integer(0)

  for (it in seq_len(n_eff)) {
    # one iteration = one sweep: run each sampler exactly once, timing per sampler
    for (i in seq_len(ns)) {
      t0 <- proc.time()[["elapsed"]]
      samplerFns[[i]]$run()
      times[i] <- times[i] + (proc.time()[["elapsed"]] - t0)
    }
    # thinning: advance (thin-1) extra sweeps without timing
    if (thin > 1L) {
      for (k in seq_len(thin - 1L)) {
        for (i in seq_len(ns)) samplerFns[[i]]$run()
      }
    }
    if (progress && it %in% pb_steps) cat(sprintf("... %d%%\n", round(100 * it / n_eff)))
  }

  times
}

# --------------------------------------------------------------
# DIAGNOSTIC STRUCTURE (dependencies + samplers + auto timing + dual plots)
# --------------------------------------------------------------

#' Diagnose model structure, dependencies, and per-sampler time (parameter- and family-level)
#'
#' @description
#' Inspect a NIMBLE model to extract the universe of nodes, classify stochastic vs.
#' deterministic nodes, compute downstream dependencies per node, map configured samplers
#' to their target nodes, and (optionally) profile per-sampler run time. Produces tidy
#' tables and publication-quality figures at both the \emph{parameter} level and the
#' \emph{family} level (where "family" = base variable name before any "i,j" index,
#' replacing indices by \code{[]}).
#'
#' Key features:
#' \itemize{
#'   \item Robust filtering of nodes (ignored patterns, removal list).
#'   \item Downstream dependency counts per node (+ per family via summary stats).
#'   \item Per-sampler times aggregated to parameters (+ per family with optional normalization).
#'   \item Optional auto-profiling of samplers via a short MCMC run.
#'   \item Non-regressive: original per-parameter plots are preserved by default.
#' }
#'
#' @param model                A compiled or uncompiled \code{nimbleModel} (required).
#' @param include_data         Logical; include data nodes when enumerating nodes (default \code{FALSE}).
#' @param removed_nodes        Character vector of nodes to exclude explicitly (default \code{NULL}).
#' @param ignore_patterns      Character vector of regex patterns to exclude (e.g., \code{"^lifted_"}, \code{"^logProb_"}).
#' @param make_plots           Logical; if \code{TRUE}, generate ggplot objects and optionally save them (default \code{TRUE}).
#' @param output_dir           Directory where figures/CSVs will be saved; if \code{NULL}, nothing is written (default \code{NULL}).
#' @param save_csv             Logical; if \code{TRUE}, write CSV exports (dependencies per node, family tables) (default \code{FALSE}).
#' @param node_of_interest     Optional character vector of nodes to highlight or subset (reserved for user logic) (default \code{NULL}).
#' @param sampler_times        Optional numeric vector of per-sampler times aligned with \code{nimble::configureMCMC(model)$getSamplers()}.
#' @param sampler_times_unit   Character label for time axis (e.g., \code{"seconds"}, \code{"ms"}) (default \code{"seconds"}).
#' @param auto_profile         Logical; if \code{TRUE} and \code{sampler_times} is \code{NULL}, profile sampler times automatically (default \code{TRUE}).
#' @param profile_niter        Integer; iterations used by the auto-profiler (default \code{2000L}).
#' @param profile_burnin       Integer; burn-in iterations for auto-profiler (default \code{500L}).
#' @param profile_thin         Integer; thinning for auto-profiler (default \code{1L}).
#' @param profile_seed         Optional integer seed for reproducibility (default \code{NULL}).
#' @param np                   Proportion in (0,1] used elsewhere for "worst sampler" selection (kept for API compatibility) (default \code{0.10}).
#' @param by_family            Logical; if \code{TRUE}, compute and plot family-level summaries in addition to parameter-level (default \code{TRUE}).
#' @param family_stat          One of \code{c("median","mean","sum")}; summary statistic for family-level aggregation (default \code{"median"}).
#' @param time_normalize       One of \code{c("none","per_node")}; if \code{"per_node"}, divide family time by the number of distinct nodes in the stats::family(default \code{"none"}).
#' @param only_family_plots    Logical; if \code{TRUE}, only family-level figures are exported (parameter plots not written).
#' @param ...                  Additional arguments forwarded to \code{nimble::configureMCMC()}.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{dependencies_df}        : data.frame of (node, dependency) pairs,
#'   \item \code{dep_counts}             : data.frame with per-parameter downstream dependency counts,
#'   \item \code{samplers_df}            : data.frame listing samplers and their target nodes (list-column),
#'   \item \code{per_param_times}        : data.frame with per-parameter aggregated sampler time,
#'   \item \code{deps_df}                : tidy data.frame used for parameter-level dependency plotting,
#'   \item \code{sampler_df}             : tidy data.frame used for parameter-level time plotting,
#'   \item \code{fam_deps_df}            : family-level dependency summary (statistic per family),
#'   \item \code{fam_time_df}            : family-level time summary (optionally normalized),
#'   \item \code{plots}                  : list of ggplot objects (some may be \code{NULL} if not created):
#'         \code{plot_dependencies}, \code{plot_sampler_time}, \code{plot_combined},
#'         \code{plot_dependencies_family}, \code{plot_sampler_time_family}, \code{plot_combined_family}.
#' }
#'
#' @examples
#' \dontrun{
#' res <- diagnose_model_structure(
#'   model = my_nimble_model,
#'   make_plots = TRUE,
#'   output_dir = "outputs/diagnostics",
#'   save_csv = TRUE,
#'   by_family = TRUE,
#'   family_stat = "median",
#'   time_normalize = "per_node",
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
                                     profile_niter       = 2000L,
                                     profile_burnin      = 500L,
                                     profile_thin        = 1L,
                                     profile_seed        = NULL,
                                     np                  = 0.10,
                                     by_family           = TRUE,
                                     family_stat         = c("median","mean","sum"),
                                     time_normalize      = c("none","per_node"),
                                     only_family_plots   = FALSE,
                                     ...) {
  stopifnot(!missing(model))

  # -- utilities ---------------------------------------------------------------
  `%||%` <- function(a, b) if (is.null(a)) b else a
  .is_ignored <- function(x, patterns) {
    if (length(patterns) == 0) return(rep(FALSE, length(x)))
    Reduce("|", lapply(patterns, function(p) grepl(p, x)))
  }
  .to_family <- function(x) sub("\\[.*\\]", "", x)

  family_stat    <- match.arg(family_stat)
  time_normalize <- match.arg(time_normalize)

  # (3) garantir des tracés au niveau des nœuds si by_family = FALSE
  if (!isTRUE(by_family)) {
    only_family_plots <- FALSE
  }

  # -- 1) Universe of nodes ----------------------------------------------------
  all_nodes_raw <- model$getNodeNames(includeData = include_data)
  removed_nodes <- unique(c(removed_nodes %||% character(0)))
  base_nodes    <- all_nodes_raw[!.is_ignored(all_nodes_raw, ignore_patterns)]
  base_nodes    <- setdiff(base_nodes, removed_nodes)

  # -- 2) Stochastic vs deterministic -----------------------------------------
  stoch_all           <- model$getNodeNames(stochOnly = TRUE, includeData = include_data)
  stochastic_nodes    <- intersect(stoch_all, base_nodes)
  deterministic_nodes <- setdiff(base_nodes, stochastic_nodes)

  # -- 3) Base-variable dimensions --------------------------------------------
  base_vars <- unique(gsub("\\[.*\\]", "", base_nodes))
  dims <- lapply(base_vars, function(v) {
    info <- try(model$getVarInfo(v), silent = TRUE)
    if (inherits(info, "try-error") || is.null(info)) NA_integer_ else info$nDim
  })
  names(dims) <- base_vars

  # -- 4) Downstream dependencies per parameter (filtered) ---------------------
  if (length(base_nodes) > 0) {
    num_dependencies <- sapply(base_nodes, function(node) {
      deps <- model$getDependencies(nodes = node, self = FALSE, downstream = TRUE)
      if (!is.null(deps) && length(deps)) {
        deps_clean <- deps[!.is_ignored(deps, ignore_patterns)]
        deps_clean <- intersect(deps_clean, base_nodes)
        length(deps_clean)
      } else 0L
    })

    filtered_node_names <- names(num_dependencies)[num_dependencies > 0]

    node_dependencies <- stats::setNames(vector("list", length(filtered_node_names)), filtered_node_names)
    for (i in seq_along(filtered_node_names)) {
      nd <- filtered_node_names[i]
      deps <- model$getDependencies(nodes = nd, self = FALSE, downstream = TRUE)
      deps_clean <- character(0)
      if (!is.null(deps) && length(deps)) {
        deps_clean <- deps[!.is_ignored(deps, ignore_patterns)]
        deps_clean <- intersect(deps_clean, base_nodes)
      }
      node_dependencies[[i]] <- deps_clean
    }

    dependencies_df <- if (length(node_dependencies)) {
      do.call(
        rbind,
        lapply(names(node_dependencies), function(nd) {
          deps <- node_dependencies[[nd]]
          if (length(deps) == 0L) {
            data.frame(node = nd, dependency = NA_character_, stringsAsFactors = FALSE)
          } else {
            data.frame(node = nd, dependency = deps, stringsAsFactors = FALSE)
          }
        })
      )
    } else {
      data.frame(node = character(0), dependency = character(0), stringsAsFactors = FALSE)
    }

    dep_counts <- data.frame(
      parameter = names(num_dependencies),
      total_dependencies = as.numeric(num_dependencies),
      stringsAsFactors = FALSE
    )
  } else {
    node_dependencies <- list()
    dependencies_df   <- data.frame(node = character(0), dependency = character(0), stringsAsFactors = FALSE)
    dep_counts        <- data.frame(parameter = character(0), total_dependencies = numeric(0))
  }

  # -- 5) Optional CSV export (node-level) ------------------------------------
  if (isTRUE(save_csv)) {
    if (is.null(output_dir)) {
      warning("save_csv=TRUE but output_dir is NULL: CSV not written.")
    } else {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(dependencies_df,
                       file = file.path(output_dir, "dependencies_per_node.csv"),
                       row.names = FALSE)
    }
  }

  # -- 6) Build internal MCMC config + extract samplers ------------------------
  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("Package 'nimble' is required for nimble::configureMCMC().")
  }
  mcmc_conf <- try(nimble::configureMCMC(model, ...), silent = TRUE)
  if (inherits(mcmc_conf, "try-error")) {
    stop("Failed to run nimble::configureMCMC(model, ...). Check your model and options passed via '...'.")
  }

  smp_list <- mcmc_conf$getSamplers()
  samplers_df <- if (length(smp_list) > 0) {
    do.call(
      rbind,
      lapply(smp_list, function(s) {
        tgt <- s$target %||% character(0)
        if (length(tgt)) {
          tgt <- tgt[!.is_ignored(tgt, ignore_patterns)]
          tgt <- intersect(tgt, base_nodes)
          tgt <- intersect(tgt, stochastic_nodes)   # **seuls nœuds stochastiques**
          if (length(tgt) == 0) tgt <- NA_character_
        } else {
          tgt <- NA_character_
        }
        data.frame(
          Type        = paste(class(s), collapse = "/"),
          TargetNodes = I(list(tgt)),
          stringsAsFactors = FALSE
        )
      })
    )
  } else {
    data.frame(Type = character(0), TargetNodes = I(list()), stringsAsFactors = FALSE)
  }

  n_samplers <- nrow(samplers_df)
  if (length(np) != 1 || !is.finite(np) || np <= 0 || np > 1) {
    warning("Invalid 'np'; using np = 0.10")
    np <- 0.10
  }
  prop_worst <- if (n_samplers > 0) max(1L, min(n_samplers, ceiling(n_samplers * np))) else 0L

  # Nœuds ayant effectivement un sampler (donc stochastiques)
  has_sampler_nodes <- unique(unlist(samplers_df$TargetNodes))
  has_sampler_nodes <- has_sampler_nodes[!is.na(has_sampler_nodes)]

  # -- 7) Per-sampler timing: provided or auto-profiled ------------------------
  if (is.null(sampler_times) && isTRUE(auto_profile)) {
    if (!exists("profile_sampler_times", mode = "function")) {
      warning("auto_profile=TRUE but 'profile_sampler_times()' is not available; setting sampler_times to NULL.")
    } else {
      sampler_times <- try(
        profile_sampler_times(
          model    = model,
          conf     = mcmc_conf,
          niter    = profile_niter,
          burnin   = profile_burnin,
          thin     = profile_thin,
          set_seed = profile_seed,
          progress = FALSE
        ),
        silent = TRUE
      )
      if (inherits(sampler_times, "try-error")) {
        warning("Profiling failed; proceeding without sampler times.")
        sampler_times <- NULL
      }
    }
  }

  # Agrégation des temps (uniquement nœuds stochastiques avec sampler)
  per_param_times <- data.frame(parameter = character(0), sampler_time = numeric(0))
  if (!is.null(sampler_times)) {
    if (!is.numeric(sampler_times)) {
      warning("sampler_times is not numeric; ignoring sampler times.")
      sampler_times <- NULL
    } else if (length(sampler_times) != n_samplers) {
      warning(sprintf("Length of sampler_times (%d) != number of samplers (%d); ignoring sampler times.",
                      length(sampler_times), n_samplers))
      sampler_times <- NULL
    }
  }
  if (!is.null(sampler_times) && n_samplers > 0 && length(has_sampler_nodes) > 0) {
    rows <- lapply(seq_len(n_samplers), function(i) {
      tgt <- samplers_df$TargetNodes[[i]]
      t_i <- sampler_times[i]
      if (is.null(tgt) || (length(tgt) == 1 && is.na(tgt))) return(NULL)
      tgt <- intersect(tgt, has_sampler_nodes)
      if (!length(tgt)) return(NULL)
      data.frame(parameter = tgt, sampler_time = t_i, stringsAsFactors = FALSE)
    })
    long_tbl <- do.call(rbind, rows)
    if (!is.null(long_tbl) && nrow(long_tbl) > 0) {
      per_param_times <- stats::aggregate(sampler_time ~ parameter, data = long_tbl, FUN = sum)
    }
  }

  # -- 8) Tidy frames for parameter-level plots --------------------------------
  ord_nodes <- if (nrow(dep_counts) > 0) {
    dep_counts$parameter[order(-dep_counts$total_dependencies)]
  } else base_nodes

  deps_df <- data.frame(
    parameter          = factor(ord_nodes, levels = ord_nodes),
    total_dependencies = as.numeric(stats::setNames(dep_counts$total_dependencies, dep_counts$parameter)[ord_nodes]),
    stringsAsFactors   = FALSE
  )

  # Temps : restreints aux nœuds avec sampler
  st_named    <- stats::setNames(per_param_times$sampler_time, per_param_times$parameter)
  time_nodes  <- intersect(ord_nodes, has_sampler_nodes)
  sampler_vec <- st_named[time_nodes]
  ok          <- !is.na(sampler_vec)
  time_nodes  <- time_nodes[ok]
  sampler_vec <- sampler_vec[ok]

  sampler_df <- data.frame(
    parameter    = factor(time_nodes, levels = time_nodes),
    sampler_time = as.numeric(sampler_vec),
    stringsAsFactors = FALSE
  )

  # -- 9) Node of interest ------------------------------------------------------
  node_of_interest_deps <- NULL
  if (!is.null(node_of_interest)) {
    if (node_of_interest %in% base_nodes) {
      deps <- model$getDependencies(nodes = node_of_interest, self = FALSE, downstream = TRUE)
      deps_clean <- character(0)
      if (!is.null(deps) && length(deps)) {
        deps_clean <- deps[!.is_ignored(deps, ignore_patterns)]
        deps_clean <- intersect(deps_clean, base_nodes)
      }
      node_of_interest_deps <- deps_clean
    } else {
      warning(sprintf("node_of_interest = '%s' not present after filtering; skipping.", node_of_interest))
    }
  }

  # -- 10) Family summaries -----------------------------------------------------
  # (1) Dépendances par famille = toutes les familles présentes dans le modèle
  if (nrow(deps_df) > 0) {
    deps_df$fam <- .to_family(as.character(deps_df$parameter))
    dep_by_fam <- stats::aggregate(
      total_dependencies ~ fam, data = deps_df,
      FUN = switch(family_stat,
                   median = stats::median,
                   mean   = base::mean,
                   sum    = base::sum)
    )
    names(dep_by_fam) <- c("family","deps_stat")
  } else {
    dep_by_fam <- data.frame(family = character(0), deps_stat = numeric(0))
  }

  # (2) Temps par famille = uniquement nœuds stochastiques avec sampler
  if (nrow(sampler_df) > 0) {
    sampler_df$fam <- .to_family(as.character(sampler_df$parameter))
    time_by_fam <- stats::aggregate(
      sampler_time ~ fam, data = sampler_df,
      FUN = switch(family_stat,
                   median = stats::median,
                   mean   = base::mean,
                   sum    = base::sum)
    )
    names(time_by_fam) <- c("family","time_stat")

    if (time_normalize == "per_node") {
      n_by_fam <- stats::aggregate(parameter ~ fam, data = sampler_df, FUN = function(z) length(unique(z)))
      names(n_by_fam) <- c("family","n_nodes")
      time_by_fam <- merge(time_by_fam, n_by_fam, by = "family", all.x = TRUE)
      time_by_fam$time_stat <- time_by_fam$time_stat / pmax(time_by_fam$n_nodes, 1L)
      time_by_fam$n_nodes <- NULL
    }
  } else {
    time_by_fam <- data.frame(family = character(0), time_stat = numeric(0))
  }

  # Ordres pour facettes/axes
  if (nrow(dep_by_fam) > 0) {
    ord_fam_all <- dep_by_fam$family[order(-dep_by_fam$deps_stat)]
  } else {
    ord_fam_all <- character(0)
  }

  fam_deps_df <- dep_by_fam
  if (nrow(fam_deps_df) > 0) {
    fam_deps_df$family <- factor(fam_deps_df$family, levels = ord_fam_all)
  }

  fam_time_df <- time_by_fam
  if (nrow(fam_time_df) > 0) {
    ord_fam_time <- if (length(ord_fam_all)) intersect(ord_fam_all, fam_time_df$family) else fam_time_df$family
    fam_time_df$family <- factor(fam_time_df$family, levels = ord_fam_time)
  }

  # -- 11) Plotting -------------------------------------------------------------
  plot_dependencies <- plot_sampler_time <- plot_combined <- NULL
  plot_dependencies_family <- plot_sampler_time_family <- plot_combined_family <- NULL

  if (isTRUE(make_plots)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("ggplot2 not installed: plots will not be produced.")
    } else {

      base_theme <- ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          plot.title   = ggplot2::element_text(face = "bold", hjust = 0.5),
          axis.title.x = ggplot2::element_text(face = "bold"),
          axis.title.y = ggplot2::element_text(face = "bold"),
          panel.grid.minor = ggplot2::element_blank(),
          legend.position = "top",
          legend.title = ggplot2::element_text(face = "bold"),
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
        )

      # -- parameter-level plots ----------------------------------------------
      if (!isTRUE(only_family_plots)) {
        if (nrow(deps_df) > 0) {
          plot_dependencies <- ggplot2::ggplot(deps_df,
                                               ggplot2::aes(x = parameter, y = total_dependencies, fill = "Total dependencies")) +
            ggplot2::geom_col(width = 1) +
            ggplot2::scale_x_discrete(expand = c(0, 0)) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
            ggplot2::scale_fill_manual(name = "", values = c("Total dependencies" = "orange")) +
            ggplot2::labs(title = "Number of  dependencies per node",
                          x = "Node", y = "Number of dependencies") +
            base_theme
        }

        if (nrow(sampler_df) > 0) {
          time_label <- sprintf("Time in samplers (%s)", sampler_times_unit)
          plot_sampler_time <- ggplot2::ggplot(sampler_df,
                                               ggplot2::aes(x = parameter, y = sampler_time, fill = time_label)) +
            ggplot2::geom_col(width = 1) +
            ggplot2::scale_x_discrete(expand = c(0, 0)) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
            ggplot2::scale_fill_manual(name = "", values = stats::setNames("blue", time_label)) +
            ggplot2::labs(title = "Time spent in samplers per node",
                          x = "Node", y = time_label) +
            base_theme
        }

        if (requireNamespace("patchwork", quietly = TRUE) && !is.null(plot_dependencies) && !is.null(plot_sampler_time)) {
          plot_combined <- plot_dependencies + plot_sampler_time + patchwork::plot_layout(ncol = 2)
        } else if (requireNamespace("gridExtra", quietly = TRUE) && !is.null(plot_dependencies) && !is.null(plot_sampler_time)) {
          plot_combined <- gridExtra::arrangeGrob(plot_dependencies, plot_sampler_time, ncol = 2)
        }
      }

      # -- family-level plots ---------------------------------------------------
      if (isTRUE(by_family)) {
        if (nrow(fam_deps_df) > 0) {
          lab_dep <- sprintf("Number dependencies (%s by family)", family_stat)
          plot_dependencies_family <- ggplot2::ggplot(
            fam_deps_df, ggplot2::aes(x = family, y = deps_stat, fill = "Deps by family")
          ) +
            ggplot2::geom_col(width = 1) +
            ggplot2::scale_x_discrete(expand = c(0,0)) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
            ggplot2::scale_fill_manual(name = "", values = c("Deps by family" = "orange")) +
            ggplot2::labs(title = lab_dep, x = "Family", y = family_stat) +
            base_theme +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        }

        if (nrow(fam_time_df) > 0) {
          lab_tim <- sprintf("Time in samplers (%s by family, %s)",
                             family_stat,
                             if (time_normalize == "per_node") "per node" else "raw")
          plot_sampler_time_family <- ggplot2::ggplot(
            fam_time_df, ggplot2::aes(x = family, y = time_stat, fill = "Time by family")
          ) +
            ggplot2::geom_col(width = 1) +
            ggplot2::scale_x_discrete(expand = c(0,0)) +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
            ggplot2::scale_fill_manual(name = "", values = c("Time by family" = "blue")) +
            ggplot2::labs(title = lab_tim, x = "Family",
                          y = sprintf("Time (%s)", sampler_times_unit)) +
            base_theme +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        }

        if (requireNamespace("patchwork", quietly = TRUE) &&
            !is.null(plot_dependencies_family) && !is.null(plot_sampler_time_family)) {
          plot_combined_family <- plot_dependencies_family + plot_sampler_time_family + patchwork::plot_layout(ncol = 2)
        } else if (requireNamespace("gridExtra", quietly = TRUE) &&
                   !is.null(plot_dependencies_family) && !is.null(plot_sampler_time_family)) {
          plot_combined_family <- gridExtra::arrangeGrob(plot_dependencies_family, plot_sampler_time_family, ncol = 2)
        }
      }

      # -- optional export ------------------------------------------------------
      if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

        # CSVs (family level)
        if (isTRUE(save_csv)) {
          if (nrow(fam_deps_df) > 0) {
            utils::write.csv(fam_deps_df,
                             file = file.path(output_dir, sprintf("deps_by_family_%s.csv", family_stat)),
                             row.names = FALSE)
          }
          if (nrow(fam_time_df) > 0) {
            utils::write.csv(fam_time_df,
                             file = file.path(output_dir, sprintf("sampler_time_by_family_%s_%s.csv",
                                                                  family_stat, time_normalize)),
                             row.names = FALSE)
          }
        }

        # Parameter-level plots (toujours actifs si by_family = FALSE)
        if (!isTRUE(only_family_plots)) {
          if (!is.null(plot_dependencies)) {
            ggplot2::ggsave(file.path(output_dir, "dependencies_per_parameter.png"),
                            plot = plot_dependencies, width = 12, height = 6, dpi = 300)
          }
          if (!is.null(plot_sampler_time)) {
            ggplot2::ggsave(file.path(output_dir, "sampler_time_per_parameter.png"),
                            plot = plot_sampler_time, width = 12, height = 6, dpi = 300)
          }
          if (!is.null(plot_combined)) {
            if (inherits(plot_combined, "ggplot")) {
              ggplot2::ggsave(file.path(output_dir, "deps_and_sampler_time_side_by_side.png"),
                              plot = plot_combined, width = 14, height = 6, dpi = 300)
            } else {
              grDevices::png(file.path(output_dir, "deps_and_sampler_time_side_by_side.png"),
                             width = 14, height = 6, units = "in", res = 300)
              grid::grid.draw(plot_combined)
              grDevices::dev.off()
            }
          }
        }

        # Family-level plots
        if (nrow(fam_deps_df) > 0 && !is.null(plot_dependencies_family)) {
          ggplot2::ggsave(file.path(output_dir, sprintf("deps_by_family_%s.png", family_stat)),
                          plot = plot_dependencies_family, width = 12, height = 6, dpi = 300)
        }
        if (nrow(fam_time_df) > 0 && !is.null(plot_sampler_time_family)) {
          ggplot2::ggsave(file.path(output_dir, sprintf("sampler_time_by_family_%s_%s.png",
                                                        family_stat, time_normalize)),
                          plot = plot_sampler_time_family, width = 12, height = 6, dpi = 300)
        }
        if (!is.null(plot_combined_family)) {
          if (inherits(plot_combined_family, "ggplot")) {
            ggplot2::ggsave(file.path(output_dir, sprintf("family_deps_and_time_side_by_side_%s_%s.png",
                                                          family_stat, time_normalize)),
                            plot = plot_combined_family, width = 14, height = 6, dpi = 300)
          } else {
            grDevices::png(file.path(output_dir, sprintf("family_deps_and_time_side_by_side_%s_%s.png",
                                                         family_stat, time_normalize)),
                           width = 14, height = 6, units = "in", res = 300)
            grid::grid.draw(plot_combined_family)
            grDevices::dev.off()
          }
        }
      }
    }
  }

  # -- 12) Return structured result -------------------------------------------
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
    per_param_times       = per_param_times,   # only stochastic nodes with a sampler
    deps_df               = deps_df,
    sampler_df            = sampler_df,        # only stochastic nodes with a sampler
    fam_deps_df           = fam_deps_df,       # deps: toutes familles du modèle
    fam_time_df           = fam_time_df,       # time: familles filtrées aux nœuds avec samplers
    plots = list(
      plot_dependencies          = plot_dependencies,
      plot_sampler_time          = plot_sampler_time,
      plot_combined              = plot_combined,
      plot_dependencies_family   = plot_dependencies_family,
      plot_sampler_time_family   = plot_sampler_time_family,
      plot_combined_family       = plot_combined_family
    )
  ))
}


#' Test and compare MCMC strategies on selected bottleneck nodes
#'
#' Runs a small, reproducible workflow to (i) build a model via `build_fn`,
#' (ii) compute a baseline with the default MCMC configuration, then (iii)
#' try alternative samplers in a *strict*, user-defined order on one or two
#' “bottleneck” targets (singleton and then optional block on the union).
#' Optionally attempts full-model HMC/NUTS first (if `nimbleHMC` is available
#' and the model supports derivatives). Results and diagnostic plots
#' (R-hat bars and trace/density) are written under `out_dir`.
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
#' When \code{ask = TRUE}, interactive yes/no prompts allow you to stop early.
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
#'test_strategy
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
  if (is.list(build_fn) && !is.null(build_fn$model)) { obj <- build_fn; build_fn <- function() obj }
  stopifnot(is.function(build_fn))
  if (!requireNamespace("nimble", quietly = TRUE)) stop("samOptiPro: 'nimble' is required.")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- utils
  say     <- function(...) { msg <- try(sprintf(...), silent = TRUE); if (!inherits(msg,"try-error")) cat(msg,"\n") }
  is_ll   <- function(x) grepl("^logLik(\\[.*\\])?$|log_?lik|logdens|lpdf", x, perl=TRUE, ignore.case=TRUE)

  # (A) sanitize niter/burnin once, and reuse
  .sanitize_iters <- local({
    warned <- FALSE
    function(niter, nburnin) {
      niter   <- as.integer(niter)
      nburnin <- as.integer(nburnin)
      if (is.na(niter) || niter < 1L) niter <- 1L
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

  # (B) strict yes/no with flush + loop
  .ask_yes_no_strict <- function(prompt) {
    if (!isTRUE(ask) || !interactive()) return(TRUE)
    while (TRUE) {
      cat(paste0(prompt, " (yes/no): "))
      utils::flush.console()
      ans <- try(readline(), silent = TRUE)
      if (inherits(ans, "try-error")) return(TRUE) # non-interactive fallback = proceed
      ans <- tolower(trimws(ans))
      if (ans %in% c("y","yes","o","oui")) return(TRUE)
      if (ans %in% c("n","no","non"))      return(FALSE)
      cat("Please answer 'yes' or 'no'.\n")
    }
  }

  # barrier & restore for smooth scalar->block
  .barrier_before_block <- function(tag = "scalar_to_block") {
    try(nimble::clearCompiled(), silent = TRUE); gc()
    new_bd <- file.path(tempdir(), paste0("samOptiPro_build_", gsub("[^A-Za-z0-9_]+","_", tag), "_", as.integer(stats::runif(1,1e9,9e9))))
    dir.create(new_bd, recursive = TRUE, showWarnings = FALSE)
    assign(".sop_old_builddir", nimble::nimbleOptions("buildDir"), inherits = FALSE)
    nimble::nimbleOptions(buildDir = new_bd)
    if (.Platform$OS.type == "windows") Sys.sleep(0.3) else Sys.sleep(0.05)
    invisible(TRUE)
  }
  .restore_builddir <- function() {
    old <- try(get(".sop_old_builddir", inherits = FALSE), silent = TRUE)
    if (!inherits(old, "try-error") && !is.null(old)) try(nimble::nimbleOptions(buildDir = old), silent = TRUE)
    invisible(TRUE)
  }

  .ensure_unsampled <- function(conf) {
    uns <- try(conf$getUnsampledNodes(), silent = TRUE)
    if (!inherits(uns, "try-error") && length(uns)) {
      uns <- uns[!is_ll(uns)]
      for (u in uns) conf$addSampler(u, type = "slice")
    }
  }

  .compile_and_run <- function(step_tag, build_obj, conf, niter, nburnin) {
    it <- .sanitize_iters(niter, nburnin); niter <- it$niter; nburnin <- it$nburnin
    attempts <- 0L; last_err <- NULL
    repeat {
      attempts <- attempts + 1L
      try(nimble::clearCompiled(), silent = TRUE); gc()
      res <- try({
        .ensure_unsampled(conf)
        cmcmc <- .compile_mcmc_with_build(conf, build_obj, reset = TRUE, show = FALSE)
        out   <- .run_and_collect(cmcmc, niter = niter, nburnin = nburnin, thin = thin, nchains = nchains)
        ml    <- as_mcmc_list_sop(out$samples, out$samples2, drop_loglik = FALSE, thin = thin)
        dg    <- compute_diag_from_mcmc(ml, runtime_s = out$runtime_s)
        if (!is.null(cmcmc) && is.list(cmcmc) && "unloadDLL" %in% names(cmcmc)) try(cmcmc$unloadDLL(), silent = TRUE)
        list(out=out, ml=ml, dg=dg)
      }, silent = TRUE)
      if (!inherits(res, "try-error")) return(res)
      last_err <- tryCatch(conditionMessage(attr(res, "condition")), error=function(e) as.character(res))
      if (attempts == 1L) message(sprintf("[info] First attempt failed at %s; retrying cleanly...", step_tag))
      else message(sprintf("[retry %d @ %s] %s", attempts, step_tag, last_err))
      build_obj <- .fresh_build(build_fn, monitors = monitors, thin = thin)
      conf <- build_obj$conf
      if (attempts >= 3L) stop(sprintf("Compilation failed after %d attempts at step '%s': %s", attempts, step_tag, last_err))
    }
  }

  # assigners
  .add_scalar <- function(conf, nodes, type, has_hmc) {
    if (identical(type, "NUTS")) {
      if (isTRUE(has_hmc)) { for (n in nodes) conf$addSampler(target = n, type = "NUTS") }
      else { cat("[Info] nimbleHMC not available -> falling back to slice (scalar).\n"); for (n in nodes) conf$addSampler(n, "slice", slice_control) }
      return(invisible())
    }
    if (identical(type, "slice")) { for (n in nodes) conf$addSampler(n, "slice", slice_control); return(invisible()) }
    if (identical(type, "RW"))    { for (n in nodes) conf$addSampler(n, "RW",    rw_control);    return(invisible()) }
    for (n in nodes) conf$addSampler(n, "slice", slice_control)
  }
  .add_block <- function(conf, nodes, type, has_hmc, Sc = NULL) {
    nodes <- unique(nodes)
    if (length(nodes) < 2L) {
      if (type == "NUTS_block") .add_scalar(conf, nodes, "NUTS", has_hmc) else
        if (type == "AF_slice")   for (n in nodes) conf$addSampler(n, "AF_slice", af_slice_control) else
          if (type == "RW_block")   for (n in nodes) conf$addSampler(n, "RW", rw_control) else
            for (n in nodes) conf$addSampler(n, "slice", slice_control)
      return(invisible())
    }
    if (identical(type, "NUTS_block")) {
      if (isTRUE(has_hmc)) conf$addSampler(target = nodes, type = "NUTS")
      else { cat("[Info] nimbleHMC not available -> falling back to AF_slice (block).\n"); conf$addSampler(target = nodes, type = "AF_slice", control = af_slice_control) }
      return(invisible())
    }
    if (identical(type, "AF_slice")) {
      ok <- TRUE
      tryCatch({ conf$addSampler(target = nodes, type = "AF_slice", control = af_slice_control) }, error = function(e) ok <<- FALSE)
      if (!ok) for (n in nodes) conf$addSampler(n, "slice", slice_control)
      return(invisible())
    }
    if (identical(type, "RW_block")) {
      ctrl <- rwblock_control; if (!is.null(Sc)) ctrl$propCov <- Sc
      conf$addSampler(nodes, "RW_block", ctrl); return(invisible())
    }
    for (n in nodes) conf$addSampler(n, "slice", slice_control)
  }

  .metrics_for <- function(dg, nodes) {
    nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
    keep <- dg$target %in% nodes & !is_ll(dg$target)
    list(AE = stats::median(dg$AE_ESS_per_it[keep], na.rm = TRUE),
         CE = stats::median(dg$ESS_per_sec[keep],   na.rm = TRUE),
         Rhat = if (any(keep)) suppressWarnings(max(dg$Rhat[keep], na.rm = TRUE)) else NA_real_)
  }
  .print_step <- function(title, label, nodes, sampler, runtime, met, base_dx, base_rt) {
    nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
    say("--- %s ---", title)
    if (!is.null(label)) say("Target(s): %s", label)
    say("Nodes: %s", paste(nodes, collapse=", "))
    say("Sampler: %s", sampler)
    say("Runtime_s: %.2f (baseline: %.2f)", runtime %||% NA_real_, base_rt %||% NA_real_)
    say("AE median (ESS/iter): %.3g (baseline: %.3g)", met$AE %||% NA_real_,
        stats::median(base_dx$AE_ESS_per_it, na.rm=TRUE))
    say("CE median (ESS/s):    %.3g (baseline: %.3g)", met$CE %||% NA_real_,
        stats::median(base_dx$ESS_per_sec,   na.rm=TRUE))
    say("Rhat max:             %.3g (baseline max: %.3g)",
        met$Rhat %||% NA_real_, suppressWarnings(max(base_dx$Rhat, na.rm=TRUE)))
  }

  # ---------- 0) Build + baseline ----------
  b0  <- .fresh_build(build_fn, monitors = monitors, thin = thin)
  mdl <- b0$model
  it  <- .sanitize_iters(pilot_niter, pilot_burnin)
  rb  <- run_baseline_config(build_fn, it$niter, it$nburnin, thin, monitors, nchains)
  base_ml <- as_mcmc_list_sop(rb$samples, rb$samples2, drop_loglik = FALSE, thin = thin)
  base_dg <- compute_diag_from_mcmc(base_ml, runtime_s = rb$runtime_s)
  base_dx <- base_dg[!is_ll(base_dg$target), , drop = FALSE]
  say("Baseline runtime_s: %.2f s", rb$runtime_s %||% NA_real_)
  say("Baseline median AE(ESS/iter): %.3g", stats::median(base_dx$AE_ESS_per_it, na.rm = TRUE))
  say("Baseline median CE(ESS/s):    %.3g", stats::median(base_dx$ESS_per_sec,   na.rm = TRUE))

  # ---------- 1) Optional full-model HMC ----------
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
    if (.ask_yes_no_strict(sprintf("Baseline ready. Runtime=%.2fs; median AE=%.3g; median CE=%.3g\nProceed with full-model HMC/NUTS?",
                                   rb$runtime_s %||% NA_real_,
                                   median(base_dx$AE_ESS_per_it, na.rm = TRUE),
                                   median(base_dx$ESS_per_sec,   na.rm = TRUE)))) {
      hmc_try <- configure_hmc_safely(build_fn = build_fn, niter = it$niter, nburnin = it$nburnin,
                                      thin = thin, monitors = monitors, nchains = nchains,
                                      out_dir = file.path(out_dir, "HMC_full"))
      if (isTRUE(hmc_try$ok)) {
        dg <- hmc_try$diag_tbl
        cat(sprintf("HMC runtime_s: %.3f\n", hmc_try$res$runtime_s %||% NA_real_))
        cat(sprintf("HMC median AE: %.3g ; CE: %.3g ; max Rhat: %.3g\n",
                    median(dg$AE_ESS_per_it, na.rm=TRUE), median(dg$ESS_per_sec, na.rm=TRUE),
                    if (all(is.na(dg$Rhat))) NA_real_ else max(dg$Rhat, na.rm=TRUE)))
        return(list(mode="HMC_full",
                    baseline=list(runtime_s=rb$runtime_s, samples=base_ml, diag_tbl=base_dg),
                    hmc=hmc_try, messages="Full HMC completed."))
      } else {
        cat("[Warn] Full-model HMC failed -> continuing with surgical singleton plan.\n")
      }
    } else {
      cat("User declined full-model HMC/NUTS. Switching to surgical singleton plan.\n")
    }
  }

  # ---------- 2) Choose singleton bottlenecks (robust + forcing) ----------
  stoch_nodes <- mdl$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  stoch_nodes <- stoch_nodes[!is_ll(stoch_nodes)]

  if (is.null(force_union_nodes) && !is.null(force_union)) force_union_nodes <- force_union

  pick_nodes <- character(0)
  # From forcing
  if (!is.null(force_singletons) && length(force_singletons)) {
    pn <- unique(force_singletons)
    present <- intersect(pn, intersect(stoch_nodes, unique(colnames(as.matrix(base_ml[[1]])))))
    pick_nodes <- utils::head(present, nbot)
  }
  # From diagnostics at node-level
  if (!length(pick_nodes)) {
    node_like <- grepl("\\[", base_dx$target) # heuristic
    if (any(node_like)) {
      ord <- order(base_dx$ESS_per_sec[node_like], decreasing = FALSE)
      cand <- base_dx$target[node_like][ord]
      pick_nodes <- intersect(utils::head(cand, nbot), stoch_nodes)
    }
  }
  # From samples ESS/s
  if (!length(pick_nodes)) {
    if (!requireNamespace("coda", quietly = TRUE)) stop("samOptiPro: 'coda' package required for node-level ESS.")
    common_cols <- Reduce(intersect, lapply(base_ml, function(m) colnames(as.matrix(m))))
    node_cols   <- intersect(common_cols, stoch_nodes)
    if (length(node_cols)) {
      ess_per_node <- sapply(node_cols, function(nm) {
        ess_ch <- sapply(base_ml, function(m) { X <- as.matrix(m); if (nm %in% colnames(X)) coda::effectiveSize(X[, nm]) else NA_real_ })
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
  if (!length(pick_nodes)) stop("No stochastic singleton targets to operate on.")
  say("Selected singleton bottlenecks (nbot=%d): %s", length(pick_nodes), paste(pick_nodes, collapse = ", "))

  steps <- list()

  # ======================= CASE nbot = 1 =======================
  if (nbot == 1L) {
    node1 <- pick_nodes[1]
    for (samp in unique(as.character(strict_scalar_seq))) {
      bS <- .fresh_build(build_fn, monitors = monitors, thin = thin)
      confS <- bS$conf
      try(confS$removeSamplers(node1), silent = TRUE)
      .add_scalar(confS, node1, samp, has_hmc)

      itS  <- .sanitize_iters(pilot_niter, pilot_burnin)
      resS <- .compile_and_run(paste0("nbot1_", gsub("[^A-Za-z0-9_]+","_", node1), "_", samp), bS, confS, itS$niter, itS$nburnin)
      metS <- .metrics_for(resS$dg, node1)
      .print_step("Scalar plan on singleton", node1, node1, samp, resS$out$runtime_s, metS, base_dx, rb$runtime_s)

      pdir <- file.path(out_dir, sprintf("singleton_%s_%s", gsub("[^A-Za-z0-9_]", "_", node1), samp))
      if (!dir.exists(pdir)) dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
      .plot_rhat_bar(resS$dg, nodes = node1, out_file = file.path(pdir, "rhat_bar.png"))
      .plot_traces(resS$ml, nodes = node1, out_file_prefix = file.path(pdir, "trace_"))

      steps <- c(steps, list(list(level="singleton-scalar", nodes=node1, sampler=samp, res=resS, dir=pdir)))
      if (!.ask_yes_no_strict(sprintf("Proceed to the next sampler for '%s'?", node1))) break
    }
    return(list(status="completed", mode="surgical_nbot1_singleton", baseline=rb, targets=pick_nodes, steps=steps))
  }

  # ======================= CASE nbot >= 2 =======================
  # Step 1 -- singleton 1 scalar strict
  node1 <- pick_nodes[1]
  for (samp in unique(as.character(strict_scalar_seq))) {
    bS <- .fresh_build(build_fn, monitors = monitors, thin = thin)
    confS <- bS$conf
    try(confS$removeSamplers(node1), silent = TRUE)
    .add_scalar(confS, node1, samp, has_hmc)

    itS  <- .sanitize_iters(pilot_niter, pilot_burnin)
    resS <- .compile_and_run(paste0("nbot2_scalar_", gsub("[^A-Za-z0-9_]+","_", node1), "_", samp), bS, confS, itS$niter, itS$nburnin)
    metS <- .metrics_for(resS$dg, node1)
    .print_step("Scalar plan on singleton", node1, node1, samp, resS$out$runtime_s, metS, base_dx, rb$runtime_s)

    pdir <- file.path(out_dir, sprintf("singleton_%s_%s", gsub("[^A-Za-z0-9_]", "_", node1), samp))
    if (!dir.exists(pdir)) dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    .plot_rhat_bar(resS$dg, nodes = node1, out_file = file.path(pdir, "rhat_bar.png"))
    .plot_traces(resS$ml, nodes = node1, out_file_prefix = file.path(pdir, "trace_"))

    steps <- c(steps, list(list(level="singleton-scalar", nodes=node1, sampler=samp, res=resS, dir=pdir)))
    if (!.ask_yes_no_strict(sprintf("Proceed to the next sampler for '%s'?", node1))) break
  }

  # Step 2 -- Block strict on union {node1 ? node2} (or forced union)
  node2 <- pick_nodes[2] %||% NA_character_
  # build union safely (no NA)
  union_nodes <- if (!is.null(force_union_nodes) && length(force_union_nodes) >= 2L) {
    present <- intersect(unique(force_union_nodes),
                         intersect(mdl$getNodeNames(stochOnly=TRUE, includeData=FALSE),
                                   unique(colnames(as.matrix(base_ml[[1]])))))
    unique(present)
  } else unique(c(node1, node2))
  union_nodes <- union_nodes[!is.na(union_nodes) & nzchar(union_nodes)]
  union_nodes <- unique(union_nodes)
  if (length(union_nodes) > block_max) union_nodes <- utils::head(union_nodes, block_max)

  if (length(union_nodes) < 2L) {
    message(sprintf("[info] Skipping block phase: need at least 2 valid nodes, got %d (%s).",
                    length(union_nodes),
                    paste(union_nodes, collapse = ", ")))
    return(list(status="completed", mode="surgical_nbot2_singleton_no_block",
                baseline=rb, targets=pick_nodes, steps=steps))
  }

  .barrier_before_block("scalar_to_block_singletons"); on.exit(.restore_builddir(), add = TRUE)

  # optional PD propCov from baseline
  propCov <- NULL
  M_full <- do.call(rbind, lapply(base_ml, function(m) {
    X <- as.matrix(m); keep <- intersect(colnames(X), union_nodes); X[, keep, drop = FALSE]
  }))
  if (!is.null(M_full) && is.matrix(M_full) && ncol(M_full) >= 2L) {
    S <- try(stats::cov(M_full, use = "pairwise.complete.obs"), silent = TRUE)
    if (!inherits(S,"try-error")) {
      pc <- try(.sop_make_propCov_PD(S), silent = TRUE)
      if (!inherits(pc,"try-error")) propCov <- pc
    }
  }

  for (samp in unique(as.character(strict_block_seq))) {
    bB <- .fresh_build(build_fn, monitors = monitors, thin = thin)
    confB <- bB$conf
    for (n in union_nodes) try(confB$removeSamplers(n), silent = TRUE)
    .add_block(confB, union_nodes, samp, has_hmc, Sc = propCov)

    itB  <- .sanitize_iters(pilot_niter, pilot_burnin)
    resB <- .compile_and_run(paste0("nbot2_block_", paste(gsub("[^A-Za-z0-9_]+","_", union_nodes), collapse = "_"), "_", samp),
                             bB, confB, itB$niter, itB$nburnin)
    metB <- .metrics_for(resB$dg, union_nodes)
    .print_step("Block plan on union (singletons)", paste(union_nodes, collapse = " + "),
                union_nodes, samp, resB$out$runtime_s, metB, base_dx, rb$runtime_s)

    pdir <- file.path(out_dir, sprintf("block_union_singletons_%s",
                                       paste(gsub("[^A-Za-z0-9_]+","_", union_nodes), collapse = "_")))
    if (!dir.exists(pdir)) dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    .plot_rhat_bar(resB$dg, nodes = union_nodes, out_file = file.path(pdir, "rhat_bar.png"))
    .plot_traces(resB$ml, nodes = union_nodes, out_file_prefix = file.path(pdir, "trace_"))

    steps <- c(steps, list(list(level="singletons-block", nodes=union_nodes, sampler=samp, res=resB, dir=pdir)))
    if (!.ask_yes_no_strict("Proceed to the next BLOCK sampler for the union?")) break
  }

  return(list(status="completed", mode="surgical_nbot2_singleton",
              baseline=rb, targets=pick_nodes, steps=steps))
}

# ======================================================================
# test_strategy_family(): HMC if possible, otherwise surgical per-bottleneck
# ======================================================================

#' Family-based sampler strategy: full HMC if allowed, else surgical on bottlenecks
#'
#' @description
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
#' test_strategy_family
#' @export
test_strategy_family <- function(build_fn,
                                 monitors            = NULL,   # optional, just passed through
                                 try_hmc             = TRUE,   # only used for full-model path; surgical ignores
                                 nchains             = 3L,
                                 pilot_niter         = 4000L,
                                 pilot_burnin        = 1000L,
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
  if (is.list(build_fn) && !is.null(build_fn$model)) { obj <- build_fn; build_fn <- function() obj }
  stopifnot(is.function(build_fn))
  if (!requireNamespace("nimble", quietly = TRUE)) stop("samOptiPro: 'nimble' is required.")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  say     <- function(...) { msg <- try(sprintf(...), silent = TRUE); if(!inherits(msg,"try-error")) cat(msg,"\n") }
  is_ll   <- function(x) grepl("^logLik(\\[.*\\])?$|log_?lik|logdens|lpdf", x, perl=TRUE, ignore.case=TRUE)
  root_of <- function(x) sub("\\[.*", "", x)

  .prompt_info <- function(txt) {
    if (!isTRUE(ask) || !interactive()) return(invisible(NULL))
    try(readline(paste0(txt, " (yes/no): ")), silent = TRUE); invisible(NULL)
  }

  # --------- compile/run with robust unload & retries ----------
  .ensure_unsampled <- function(conf) {
    uns <- try(conf$getUnsampledNodes(), silent = TRUE)
    if (!inherits(uns, "try-error") && length(uns)) {
      uns <- uns[!is_ll(uns)]
      for (u in uns) conf$addSampler(u, type = "slice")
    }
  }

  .compile_and_run <- function(step_tag, build_obj, conf) {
    attempts <- 0L
    last_err <- NULL
    repeat {
      attempts <- attempts + 1L
      # clean state between attempts
      try(nimble::clearCompiled(), silent = TRUE)
      gc()
      res <- try({
        .ensure_unsampled(conf)
        cmcmc <- .compile_mcmc_with_build(conf, build_obj, reset = TRUE, show = FALSE)
        out   <- .run_and_collect(cmcmc, niter = pilot_niter, nburnin = pilot_burnin,
                                  thin = thin, nchains = nchains)
        ml    <- as_mcmc_list_sop(out$samples, out$samples2, drop_loglik = FALSE, thin = thin)
        dg    <- compute_diag_from_mcmc(ml, runtime_s = out$runtime_s)
        # unload DLL to avoid lock for next step
        if (!is.null(cmcmc) && is.list(cmcmc) && "unloadDLL" %in% names(cmcmc)) {
          try(cmcmc$unloadDLL(), silent = TRUE)
        }
        list(out=out, ml=ml, dg=dg)
      }, silent = TRUE)

      if (!inherits(res, "try-error")) return(res)

      last_err <- tryCatch(conditionMessage(attr(res, "condition")), error=function(e) as.character(res))
      cat(sprintf("[Retry %d @ %s] %s\n", attempts, step_tag, last_err))

      # fresh rebuild for next attempt
      build_obj <- .fresh_build(build_fn, monitors = monitors, thin = thin)
      conf <- build_obj$conf

      if (attempts >= 3L) stop(sprintf("Compilation failed after %d attempts at step '%s': %s",
                                       attempts, step_tag, last_err))
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
    if (identical(type, "slice")) { for (n in nodes) conf$addSampler(n, "slice", slice_control); return(invisible()) }
    if (identical(type, "RW"))    { for (n in nodes) conf$addSampler(n, "RW",    rw_control);    return(invisible()) }
    for (n in nodes) conf$addSampler(n, "slice", slice_control) # fallback
  }

  .add_block_family <- function(conf, nodes, type, has_hmc, Sc = NULL) {
    nodes <- unique(nodes)
    if (length(nodes) < 2L) {
      # graceful fallback if block is degenerate
      if (type == "NUTS_block") .add_scalar_family(conf, nodes, "NUTS", has_hmc) else
        if (type == "AF_slice")   for (n in nodes) conf$addSampler(n, "AF_slice", af_slice_control) else
          if (type == "RW_block")   for (n in nodes) conf$addSampler(n, "RW", rw_control) else
            for (n in nodes) conf$addSampler(n, "slice", slice_control)
      return(invisible())
    }
    if (identical(type, "NUTS_block")) {
      if (isTRUE(has_hmc)) conf$addSampler(target = nodes, type = "NUTS") else {
        cat("[Info] nimbleHMC not available -> falling back to AF_slice (block).\n")
        conf$addSampler(target = nodes, type = "AF_slice", control = af_slice_control)
      }
      return(invisible())
    }
    if (identical(type, "AF_slice")) {
      ok <- TRUE
      tryCatch({ conf$addSampler(target = nodes, type = "AF_slice", control = af_slice_control) },
               error = function(e) ok <<- FALSE)
      if (!ok) for (n in nodes) conf$addSampler(n, "slice", slice_control)
      return(invisible())
    }
    if (identical(type, "RW_block")) {
      ctrl <- rwblock_control; if (!is.null(Sc)) ctrl$propCov <- Sc
      conf$addSampler(nodes, "RW_block", ctrl); return(invisible())
    }
    for (n in nodes) conf$addSampler(n, "slice", slice_control)
  }

  .fam_metrics <- function(dg, nodes) {
    keep <- dg$target %in% nodes & !is_ll(dg$target)
    list(AE = stats::median(dg$AE_ESS_per_it[keep], na.rm = TRUE),
         CE = stats::median(dg$ESS_per_sec[keep],   na.rm = TRUE),
         Rhat = if (any(keep)) suppressWarnings(max(dg$Rhat[keep], na.rm = TRUE)) else NA_real_)
  }

  .print_step <- function(title, fam_label, nodes, sampler, runtime, met, base_dx, base_rt) {
    say("--- %s ---", title)
    if (!is.null(fam_label)) say("Family: %s", fam_label)
    say("Nodes: %s", paste(nodes, collapse=", "))
    say("Sampler: %s", sampler)
    say("Runtime_s: %.2f (baseline: %.2f)", runtime %||% NA_real_, base_rt %||% NA_real_)
    say("AE median (ESS/iter): %.3g (baseline: %.3g)", met$AE %||% NA_real_,
        stats::median(base_dx$AE_ESS_per_it, na.rm=TRUE))
    say("CE median (ESS/s):    %.3g (baseline: %.3g)", met$CE %||% NA_real_,
        stats::median(base_dx$ESS_per_sec,   na.rm=TRUE))
    say("Rhat max:             %.3g (baseline max: %.3g)",
        met$Rhat %||% NA_real_, suppressWarnings(max(base_dx$Rhat, na.rm=TRUE)))
  }

  # ---------- 0) Build + baseline ----------
  bld <- .fresh_build(build_fn, monitors = monitors, thin = thin)
  mdl <- bld$model
  rb  <- run_baseline_config(build_fn, pilot_niter, pilot_burnin, thin, monitors, nchains)
  base_ml <- as_mcmc_list_sop(rb$samples, rb$samples2, drop_loglik = FALSE, thin = thin)
  base_dg <- compute_diag_from_mcmc(base_ml, runtime_s = rb$runtime_s)
  base_dx <- base_dg[!is_ll(base_dg$target), , drop = FALSE]
  say("Baseline runtime_s: %.2f s", rb$runtime_s %||% NA_real_)
  say("Baseline median AE(ESS/iter): %.3g", stats::median(base_dx$AE_ESS_per_it, na.rm = TRUE))
  say("Baseline median CE(ESS/s):    %.3g", stats::median(base_dx$ESS_per_sec,   na.rm = TRUE))

  # ---------- 1) Full-model HMC (optionnel) ----------
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
    if (isTRUE(ask) && isTRUE(ask_before_hmc) && isTRUE(interactive())) {
      cat(sprintf("Baseline ready. Runtime=%.2fs; median AE=%.3g; median CE=%.3g\n",
                  rb$runtime_s %||% NA_real_,
                  median(base_dx$AE_ESS_per_it, na.rm = TRUE),
                  median(base_dx$ESS_per_sec,   na.rm = TRUE)))
      ans <- try(readline("Proceed with full-model HMC/NUTS? (yes/no): "), silent = TRUE)
      if (!inherits(ans, "try-error") && tolower(trimws(ans)) %in% c("no","n")) {
        cat("User declined full-model HMC/NUTS. Switching to surgical family plan.\n")
      } else {
        hmc_try <- configure_hmc_safely(
          build_fn = build_fn, niter = pilot_niter, nburnin = pilot_burnin,
          thin = thin, monitors = monitors, nchains = nchains,
          out_dir = file.path(out_dir, "HMC_full"))
        if (isTRUE(hmc_try$ok)) {
          dg <- hmc_try$diag_tbl
          cat(sprintf("HMC runtime_s: %.3f\n", hmc_try$res$runtime_s %||% NA_real_))
          cat(sprintf("HMC median AE: %.3g ; CE: %.3g ; max Rhat: %.3g\n",
                      median(dg$AE_ESS_per_it, na.rm=TRUE),
                      median(dg$ESS_per_sec,   na.rm=TRUE),
                      if (all(is.na(dg$Rhat))) NA_real_ else max(dg$Rhat, na.rm=TRUE)))
          return(list(mode="HMC_full",
                      baseline=list(runtime_s=rb$runtime_s, samples=base_ml, diag_tbl=base_dg),
                      hmc=hmc_try, messages="Full HMC completed."))
        } else cat("[Warn] Full-model HMC failed -> continuing with surgical family plan.\n")
      }
    }
  }

  # ---------- 2) Families selection (respecting forcing) ----------
  stoch_nodes <- mdl$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  stoch_nodes <- stoch_nodes[!is_ll(stoch_nodes)]
  if (!"Family" %in% names(base_dx)) base_dx$Family <- root_of(base_dx$target)

  fam_med <- stats::aggregate(ESS_per_sec ~ Family, data = base_dx, median)
  fam_med <- fam_med[order(fam_med$ESS_per_sec, decreasing = FALSE), , drop = FALSE]

  pick_fams <- character(0)
  if (!is.null(force_families) && length(force_families)) {
    pf <- unique(force_families)
    present <- pf[pf %in% unique(base_dx$Family)]
    pick_fams <- utils::head(present, nbot)
    if (!length(pick_fams)) warning("force_families provided but none matched diagnostics; using automatic selection.")
  }
  if (!length(pick_fams)) {
    for (fam in fam_med$Family) {
      fam_nodes_all <- intersect(stoch_nodes[root_of(stoch_nodes) == fam], unique(base_dx$target))
      if (length(fam_nodes_all)) pick_fams <- unique(c(pick_fams, fam))
      if (length(pick_fams) >= nbot) break
    }
    pick_fams <- utils::head(pick_fams, nbot)
  }
  if (!length(pick_fams)) stop("No stochastic bottleneck families found.")
  say("Selected bottleneck families (nbot=%d): %s", nbot, paste(pick_fams, collapse = ", "))

  get_family_nodes <- function(fam) {
    allf <- intersect(stoch_nodes[root_of(stoch_nodes) == fam], unique(base_dx$target))
    forced <- force_nodes[[fam]] %||% NULL
    if (!is.null(forced)) {
      forced <- intersect(forced, allf)
      if (!length(forced)) warning(sprintf("force_nodes for '%s' not found in diagnostics; using detected nodes.", fam))
      return(if (length(forced)) forced else allf)
    }
    allf
  }

  steps <- list()

  # ========================== CASE nbot = 1 ==========================
  if (nbot == 1L) {
    fam1 <- pick_fams[1]
    fam1_nodes <- get_family_nodes(fam1)
    for (samp in unique(as.character(strict_scalar_seq))) {
      bS <- .fresh_build(build_fn, monitors = monitors, thin = thin)
      confS <- bS$conf
      for (n in fam1_nodes) try(confS$removeSamplers(n), silent = TRUE)
      .add_scalar_family(confS, fam1_nodes, samp, has_hmc)

      resS <- .compile_and_run(paste0("nbot1_", fam1, "_", samp), bS, confS)
      metS <- .fam_metrics(resS$dg, fam1_nodes)
      .print_step("Scalar plan on family", fam1, fam1_nodes, samp,
                  resS$out$runtime_s, metS, base_dx, rb$runtime_s)

      pdir <- file.path(out_dir, sprintf("scalar_family_%s_%s", gsub("[^A-Za-z0-9_]", "_", fam1), samp))
      if (!dir.exists(pdir)) dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
      .plot_rhat_bar(resS$dg, nodes = fam1_nodes, out_file = file.path(pdir, "rhat_bar.png"))
      .plot_traces(resS$ml, nodes = fam1_nodes, out_file_prefix = file.path(pdir, "trace_"))

      steps <- c(steps, list(list(level="family-scalar", family=fam1, nodes=fam1_nodes,
                                  sampler=samp, res=resS, dir=pdir)))

      .prompt_info(sprintf("Proceed to the next sampler for family '%s'?", fam1))
    }
    return(list(status="completed", mode="surgical_nbot1",
                baseline=rb, families=pick_fams, steps=steps))
  }

  # ========================== CASE nbot >= 2 ==========================
  # Step 1 -- Family 1 scalar strict
  fam1 <- pick_fams[1]
  fam1_nodes <- get_family_nodes(fam1)
  for (samp in unique(as.character(strict_scalar_seq))) {
    bS <- .fresh_build(build_fn, monitors = monitors, thin = thin)
    confS <- bS$conf
    for (n in fam1_nodes) try(confS$removeSamplers(n), silent = TRUE)
    .add_scalar_family(confS, fam1_nodes, samp, has_hmc)

    resS <- .compile_and_run(paste0("nbot2_scalar_", fam1, "_", samp), bS, confS)
    metS <- .fam_metrics(resS$dg, fam1_nodes)
    .print_step("Scalar plan on family", fam1, fam1_nodes, samp,
                resS$out$runtime_s, metS, base_dx, rb$runtime_s)

    pdir <- file.path(out_dir, sprintf("scalar_family_%s_%s", gsub("[^A-Za-z0-9_]", "_", fam1), samp))
    if (!dir.exists(pdir)) dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    .plot_rhat_bar(resS$dg, nodes = fam1_nodes, out_file = file.path(pdir, "rhat_bar.png"))
    .plot_traces(resS$ml, nodes = fam1_nodes, out_file_prefix = file.path(pdir, "trace_"))

    steps <- c(steps, list(list(level="family-scalar", family=fam1, nodes=fam1_nodes,
                                sampler=samp, res=resS, dir=pdir)))
    .prompt_info(sprintf("Proceed to the next sampler for family '%s'?", fam1))
  }

  # Step 2 -- Block strict on union
  fam2 <- pick_fams[2]
  fam2_nodes <- get_family_nodes(fam2)
  if (!is.null(force_union) && length(force_union) >= 2L) {
    union_fams  <- unique(force_union)
    union_nodes <- unique(unlist(lapply(union_fams, get_family_nodes)))
    fam_label   <- paste(union_fams, collapse = " + ")
  } else {
    union_nodes <- unique(c(fam1_nodes, fam2_nodes))
    fam_label   <- paste(fam1, "+", fam2)
  }
  if (length(union_nodes) > block_max) union_nodes <- utils::head(union_nodes, block_max)

  # optional PD propCov from baseline
  propCov <- NULL
  if (length(union_nodes) >= 2L) {
    M_full <- do.call(rbind, lapply(base_ml, function(m) {
      X <- as.matrix(m); keep <- intersect(colnames(X), union_nodes); X[, keep, drop = FALSE]
    }))
    if (!is.null(M_full) && is.matrix(M_full) && ncol(M_full) >= 2L) {
      S <- try(stats::cov(M_full, use = "pairwise.complete.obs"), silent = TRUE)
      if (!inherits(S,"try-error")) {
        pc <- try(.sop_make_propCov_PD(S), silent = TRUE)
        if (!inherits(pc,"try-error")) propCov <- pc
      }
    }
  }

  for (samp in unique(as.character(strict_block_seq))) {
    bB <- .fresh_build(build_fn, monitors = monitors, thin = thin)
    confB <- bB$conf
    for (n in union_nodes) try(confB$removeSamplers(n), silent = TRUE)
    .add_block_family(confB, union_nodes, samp, has_hmc, Sc = propCov)

    resB <- .compile_and_run(paste0("nbot2_block_", gsub("[^A-Za-z0-9_]+","_", fam_label), "_", samp), bB, confB)
    metB <- .fam_metrics(resB$dg, union_nodes)
    .print_step("Block plan on families union", fam_label, union_nodes, samp,
                resB$out$runtime_s, metB, base_dx, rb$runtime_s)

    pdir <- file.path(out_dir, sprintf("block_union_%s_%s", gsub("[^A-Za-z0-9_]", "_", fam_label), samp))
    if (!dir.exists(pdir)) dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    .plot_rhat_bar(resB$dg, nodes = union_nodes, out_file = file.path(pdir, "rhat_bar.png"))
    .plot_traces(resB$ml, nodes = union_nodes, out_file_prefix = file.path(pdir, "trace_"))

    steps <- c(steps, list(list(level="families-block", families=fam_label, nodes=union_nodes,
                                sampler=samp, res=resB, dir=pdir)))
    .prompt_info(sprintf("Proceed to the next BLOCK sampler for union '%s'?", fam_label))
  }

  return(list(status="completed", mode=if (nbot==1L) "surgical_nbot1" else "surgical_nbot2",
              baseline=rb, families=pick_fams, steps=steps))
}



# --------------------------------------------------------------
# DIAGNOSTICS
# --------------------------------------------------------------

#' Build a standard diagnostics table from a \code{coda::mcmc.list}
#'
#' Columns: \code{target}, \code{ESS}, \code{AE_ESS_per_it} (ESS/iterations),
#' \code{ESS_per_sec}, \code{time_s_per_ESS}, \code{Rhat}, \code{Family}.
#'
#' @param samples An \code{mcmc.list} (or single \code{mcmc} / \code{matrix} accepted).
#' @param runtime_s Numeric(1); wall-clock runtime in seconds.
#'
#' @return A data.frame \code{diag_tbl}.
#' @export
compute_diag_from_mcmc <- function(samples, runtime_s) {
  stopifnot(is.numeric(runtime_s), length(runtime_s) == 1L, is.finite(runtime_s))

  # Normalisation d'entree
  if (inherits(samples, "mcmc"))  samples <- coda::mcmc.list(samples)
  if (is.matrix(samples))         samples <- coda::mcmc.list(coda::mcmc(samples))
  if (!inherits(samples, "mcmc.list")) stop("compute_diag_from_mcmc: 'samples' must be mcmc.list/mcmc/matrix.")

  # mcpar / iterations gardees par chaine (1ere chaine comme reference)
  mcpar <- attr(samples[[1]], "mcpar")
  n_keep <- if (!is.null(mcpar) && length(mcpar) >= 3L) {
    as.integer((mcpar[2] - mcpar[1]) / mcpar[3] + 1L)
  } else nrow(samples[[1]])

  mats <- lapply(samples, function(m) as.matrix(m))
  params <- Reduce(intersect, lapply(mats, colnames))
  if (is.null(params) || !length(params)) stop("compute_diag_from_mcmc: no common parameters between channels.")
  params <- as.character(params)

  ess_min <- numeric(length(params)); names(ess_min) <- params
  rhat    <- rep(NA_real_, length(params)); names(rhat) <- params

  # ESS min et Rhat robuste par parametre
  for (p in params) {
    ess_chain <- vapply(mats, function(M) {
      v <- M[, p]
      as.numeric(tryCatch(coda::effectiveSize(v), error = function(e) NA_real_))
    }, numeric(1))
    ess_min[p] <- suppressWarnings(min(ess_chain, na.rm = TRUE))

    ml_param <- coda::mcmc.list(lapply(seq_along(samples), function(k) {
      m <- as.matrix(samples[[k]])
      coda::mcmc(m[, p, drop = FALSE],
                 start = if (!is.null(mcpar)) mcpar[1] else 1,
                 end   = if (!is.null(mcpar)) mcpar[2] else (nrow(m) + (mcpar[1] %||% 1) - 1),
                 thin  = if (!is.null(mcpar)) mcpar[3] else 1)
    }))
    rhat[p] <- tryCatch({
      gd <- coda::gelman.diag(ml_param, autoburnin = FALSE)$psrf
      as.numeric(gd[1])  # plus robuste que nommer la colonne
    }, error = function(e) NA_real_)
  }

  ess_min[!is.finite(ess_min)] <- NA_real_

  # AE = ESS / iterations gardees ; ESS/s ; cout temps par ESS
  AE           <- ess_min / max(n_keep, 1L)
  ESS_per_sec  <- if (is.finite(runtime_s) && runtime_s > 0) ess_min / runtime_s else NA_real_
  time_s_per_ESS <- runtime_s / pmax(ess_min, .Machine$double.eps)

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
  out <- out[is.finite(out$ESS) & out$ESS > 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}

