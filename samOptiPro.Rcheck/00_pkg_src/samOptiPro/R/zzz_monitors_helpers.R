#' importFrom(utils, sessionInfo)
#' importFrom(utils, capture.output)
#' importFrom(stats, reorder)
#' importFrom(stats, aggregate)
#' Internal helpers for samOptiPro
#'
#' Central helper functions used internally by samOptiPro for model
#' interrogation, monitor configuration, plotting helpers and various
#' low-level utilities. Not intended to be called directly by end-users.
#'
#' @keywords internal
#' @importFrom utils globalVariables sessionInfo capture.output
#' @importFrom stats runif aggregate reorder
NULL

## Global variables for NSE (dplyr/data.table/ggplot2)
utils::globalVariables(c(
  "target", "Rhat", "Family", "MedianAlgorithmicEfficiency",
  "MedianComputationalEfficiency_tot", "CE", "time_s", "median_Rhat",
  "iter", "value", "chain", "nodes", "parameter", "total_dependencies",
  "sampler_time", "family", "deps_stat", "time_stat", "y",
  ".", ".N", ".__step__", ".data", ".ggsave",
  "AE", "AE_ESS_per_it", "AE_median", "AE_rank",
  "Agg_top", "CE_median", "CE_rank", "ESS", "ESS_per_sec",
  "Rhat_median", "TIME_rank", "__sampler__", "all_nodes",
  "base", "bins", "built", "diag_tbl", "key", "key_f", "m",
  "ok_grid", "out_dir", "per", "sdir", "st", "top_by", "top_k",
  "valid"
))

## ---- Exported functions (1) ----

#' Compile a nimble MCMC for a given build object
#'
#' Always compiles the MCMC against the same R-level model used during the
#' build step (\code{build_obj$model}). This avoids subtle mismatches.
#'
#' @param conf A \code{configureMCMC} object.
#' @param build_obj A list returned by your build function (must contain
#'   \code{$model}).
#' @param reset Logical; passed to \code{nimble::compileNimble()}.
#' @param show Logical; show compiler output.
#'
#' @return A compiled MCMC object as returned by \code{nimble::compileNimble()}.
#' @export
#' @keywords internal
.compile_mcmc_with_build <- function(conf, build_obj, reset = TRUE, show = FALSE) {
  mcmc <- nimble::buildMCMC(conf)
  nimble::compileNimble(
    mcmc,
    project            = build_obj$model,
    resetFunctions     = reset,
    showCompilerOutput = show
  )
}


#' Available variable roots present in a model
#' @keywords internal
.avail_vars <- function(model) {
  stopifnot(.sop_is_model(model))
  m  <- .sop_get_uncompiled_model(model)
  nn <- try(m$getNodeNames(includeData = TRUE), silent = TRUE)
  if (inherits(nn, "try-error") || is.null(nn)) return(character(0))
  unique(sub("\\[.*\\]$", "", nn))
}

#' Temporarily switch nimble buildDir (clear compiled & GC)
#'
#' @param tag Character tag used in the temporary directory name.
#' @return TRUE (invisibly).
#' @keywords internal
.barrier_before_block <- function(tag = "scalar_to_block") {
  try(nimble::clearCompiled(), silent = TRUE)
  gc()
  new_bd <- file.path(
    tempdir(),
    paste0(
      "samOptiPro_build_",
      gsub("[^A-Za-z0-9_]+", "_", tag),
      "_",
      as.integer(stats::runif(1, 1e9, 9e9))
    )
  )
  dir.create(new_bd, recursive = TRUE, showWarnings = FALSE)
  assign(".sop_old_builddir", nimble::nimbleOptions("buildDir"), inherits = FALSE)
  nimble::nimbleOptions(buildDir = new_bd)
  if (.Platform$OS.type == "windows") Sys.sleep(0.3) else Sys.sleep(0.05)
  invisible(TRUE)
}

#' Build a configureMCMC with sanitized/expanded monitors
#'
#' @param model Nimble model (R-level or compiled).
#' @param monitors Monitors roots (optional).
#' @param thin,thin2 Thinning for \code{monitors} / \code{monitors2}.
#' @param opts Optional package options.
#'
#' @return A \code{configureMCMC} object with thin/thin2 harmonised and the
#'   attributes \code{._sop_monitors_roots}, \code{._sop_monitors2_roots},
#'   \code{._sop_monitors_nodes}, \code{._sop_monitors2_nodes}.
#' @keywords internal
.configure_with_monitors <- function(model,
                                     monitors = NULL,
                                     thin = 1L,
                                     thin2 = NULL,
                                     opts = samOptiPro_options()) {
  if (!.sop_is_model(model)) {
    stop(
      "samOptiPro: '.configure_with_monitors' did not receive a nimble model. Classes: ",
      paste(class(model), collapse = "/")
    )
  }
  m <- .sop_get_uncompiled_model(model)
  `%||%` <- function(x, y) if (is.null(x)) y else x

  thin  <- as.integer(thin  %||% 1L); if (thin  < 1L) thin  <- 1L
  thin2 <- as.integer(thin2 %||% thin); if (thin2 < 1L) thin2 <- 1L

  mons_in <- monitors %||% try(default_monitors(m, opts), silent = TRUE)
  if (inherits(mons_in, "try-error") || is.null(mons_in)) mons_in <- character(0)
  mons_in <- .sop_sanitize_roots(mons_in)

  is_ll      <- grepl("^logLik(\\[.*\\])?$|log_?lik|logdens|lpdf",
                      mons_in, ignore.case = TRUE, perl = TRUE)
  roots_main <- mons_in[!is_ll]
  roots_ll   <- mons_in[is_ll]

  nodes_main <- .sop_expand_roots_to_nodes(m, roots_main)
  nodes_ll   <- .sop_expand_roots_to_nodes(m, roots_ll)

  if (length(nodes_main) == 0 && length(nodes_ll) == 0) {
    conf <- nimble::configureMCMC(m, thin = thin)
  } else {
    conf <- nimble::configureMCMC(
      m,
      monitors  = nodes_main,
      monitors2 = nodes_ll,
      thin      = thin
    )
  }

  if (is.function(try(conf$setThin,  silent = TRUE))) {
    try(conf$setThin(as.integer(thin)),  silent = TRUE)
  }
  if (is.function(try(conf$setThin2, silent = TRUE))) {
    try(conf$setThin2(as.integer(thin2)), silent = TRUE)
  }

  attr(conf, "._sop_monitors_roots")  <- roots_main
  attr(conf, "._sop_monitors2_roots") <- roots_ll
  attr(conf, "._sop_monitors_nodes")  <- nodes_main
  attr(conf, "._sop_monitors2_nodes") <- nodes_ll
  conf
}

#' Save a bar chart of R-hat for selected nodes
#' @keywords internal
.plot_rhat_bar <- function(diag_tbl, nodes, out_file) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  df <- diag_tbl[diag_tbl$target %in% nodes, c("target", "Rhat"), drop = FALSE]
  if (!nrow(df)) return(invisible(NULL))
  df$target <- factor(df$target, levels = unique(df$target))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = target, y = Rhat, fill = "Rhat")) +
    ggplot2::geom_col(width = 0.9) +
    ggplot2::geom_hline(yintercept = 1.01, linetype = 2) +
    ggplot2::scale_fill_manual(name = "Legend", values = c(Rhat = "grey40")) +
    ggplot2::labs(
      title = "Gelman–Rubin diagnostic (R-hat) for selected node(s)",
      x     = "Node",
      y     = "R-hat"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "top",
      legend.title    = ggplot2::element_text(face = "bold")
    )
  ggplot2::ggsave(out_file, plot = p, width = 8, height = 4.5, dpi = 300)
  invisible(p)
}

#' Save trace plots for a set of nodes
#' @keywords internal
.plot_traces <- function(mcmc_list, nodes, out_file_prefix) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  for (node in nodes) {
    mats <- lapply(mcmc_list, function(m) {
      mm <- as.matrix(m)
      if (!node %in% colnames(mm)) return(NULL)
      mm[, node, drop = FALSE]
    })
    if (any(vapply(mats, is.null, logical(1)))) next
    df <- do.call(rbind, lapply(seq_along(mats), function(i) {
      data.frame(
        iter  = seq_len(nrow(mats[[i]])),
        value = as.numeric(mats[[i]][, 1]),
        chain = paste0("chain", i),
        stringsAsFactors = FALSE
      )
    }))
    if (is.null(df) || !nrow(df)) next
    p <- ggplot2::ggplot(df, ggplot2::aes(x = iter, y = value, color = chain)) +
      ggplot2::geom_line(linewidth = 0.4, alpha = 0.9) +
      ggplot2::labs(
        title = sprintf("Trace plot for %s", node),
        x     = "Iteration",
        y     = "Value",
        color = "Chain"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(face = "bold", hjust = 0.5),
        legend.position = "top",
        legend.title    = ggplot2::element_text(face = "bold")
      )
    ggplot2::ggsave(
      paste0(out_file_prefix, gsub("[^A-Za-z0-9_]", "_", node), ".png"),
      plot = p, width = 8, height = 4.5, dpi = 300
    )
  }
  invisible(TRUE)
}

#' Pretty-print selected monitors (optional)
#' @keywords internal
.print_monitors <- function(title, nodes, thin = 1L) {
  cat(title, "\n", sep = "")
  if (!length(nodes)) {
    cat("(None)\n")
  } else {
    cat(sprintf("thin = %d: %s\n", as.integer(thin), paste(nodes, collapse = ",")))
  }
}

#' Restore previous nimble buildDir if it was changed
#' @keywords internal
.restore_builddir <- function() {
  old <- try(get(".sop_old_builddir", inherits = FALSE), silent = TRUE)
  if (!inherits(old, "try-error") && !is.null(old)) {
    try(nimble::nimbleOptions(buildDir = old), silent = TRUE)
  }
  invisible(TRUE)
}

#' Bind two matrices by row count (truncate to common minimum)
#' @keywords internal
.sop_cbind_align <- function(m1, m2) {
  if (is.null(m2)) return(m1)
  n1 <- NROW(m1); n2 <- NROW(m2)
  if (isTRUE(n1 == n2)) return(cbind(m1, m2))
  n <- min(n1, n2)
  if (n <= 0) stop("as_mcmc_list_sop: empty or incompatible matrices (0 rows).")
  message(sprintf(
    "as_mcmc_list_sop: align %d vs %d rows -> truncated at %d.",
    n1, n2, n
  ))
  cbind(
    m1[seq_len(n), , drop = FALSE],
    m2[seq_len(n), , drop = FALSE]
  )
}

#' Detect chain prefix "chainK." in column names
#' @keywords internal
.sop_detect_chain_prefix <- function(nms) {
  m <- regexpr("^chain(\\d+)\\.", nms, perl = TRUE)
  as.integer(ifelse(m > 0, sub("^chain(\\d+)\\..*$", "\\1", nms), NA_character_))
}

#' Expand many roots to nodes
#' @keywords internal
.sop_expand_roots_to_nodes <- function(mdl, roots) {
  roots <- unique(trimws(as.character(roots)))
  roots <- roots[nzchar(roots)]
  unique(unlist(
    lapply(roots, function(v) .sop_expand_var_nodes(mdl, v)),
    use.names = FALSE
  ))
}

#' Expand one root to nodes via getVarInfo()
#' @keywords internal
.sop_expand_var_nodes <- function(mdl, v) {
  stopifnot(.sop_is_model(mdl))
  m  <- .sop_get_uncompiled_model(mdl)
  vi <- try(m$getVarInfo(v), silent = TRUE)
  if (inherits(vi, "try-error") || is.null(vi)) return(character(0))
  if (isTRUE(vi$nDim == 0L)) return(v)
  mins <- as.integer(vi$mins)
  maxs <- as.integer(vi$maxs)
  idx  <- lapply(seq_len(vi$nDim), function(d) seq.int(mins[d], maxs[d]))
  grid <- do.call(
    expand.grid,
    c(idx, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  )
  apply(
    grid,
    1L,
    function(row) paste0(v, "[", paste(row, collapse = ","), "]")
  )
}

#' Return the R-level (uncompiled) nimble model when given either R or compiled model
#' @keywords internal
.sop_get_uncompiled_model <- function(x) {
  if (is.environment(x)) {
    getM <- try(x$getModel, silent = TRUE)
    if (is.function(getM)) {
      m <- try(getM(), silent = TRUE)
      if (!inherits(m, "try-error") && is.environment(m)) return(m)
    }
  }
  x
}

#' Heuristic check that an object exposes nimble model API
#' @keywords internal
.sop_has_model_api <- function(x) {
  m <- .sop_get_uncompiled_model(x)
  if (!is.environment(m)) return(FALSE)
  has_fun <- function(obj, nm) {
    f <- tryCatch(obj[[nm]], error = function(e) NULL)
    is.function(f)
  }
  all(c(
    has_fun(m, "getNodeNames"),
    has_fun(m, "getVarInfo"),
    has_fun(m, "getVarNames"),
    has_fun(m, "getDependencies"),
    has_fun(m, "calculate"),
    has_fun(m, "getLogProb")
  ))
}

#' Cholesky check
#' @keywords internal
.sop_is_chol_ok <- function(S) {
  ok <- TRUE
  tryCatch({ chol(S) }, error = function(e) ok <<- FALSE)
  ok
}

#' Is this a nimble model?
#' @keywords internal
.sop_is_model <- function(x) {
  is.environment(x) && .sop_has_model_api(x)
}


#' Try to make a symmetric matrix positive definite
#'
#' Attempts to repair a symmetric matrix by:
#' \enumerate{
#'   \item applying \code{Matrix::nearPD()} if available;
#'   \item shrinking towards the diagonal with a grid of weights;
#'   \item adding diagonal jitter.
#' }
#'
#' @param S Symmetric numeric matrix.
#' @param jitter_seq Numeric vector of diagonal jitter values.
#' @param shrink_grid Numeric vector of shrinkage weights toward the diagonal.
#'
#' @return A positive definite numeric matrix, or \code{NULL} on failure.
#' @keywords internal
.sop_make_propCov_PD <- function(S,
                                 jitter_seq  = c(0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2),
                                 shrink_grid = c(0, 0.05, 0.1, 0.2, 0.4)) {
  S <- as.matrix(S)
  p <- ncol(S)
  if (p == 1L) {
    return(matrix(max(S[1, 1], .Machine$double.eps), 1L, 1L))
  }

  # 1) nearPD if Matrix is available
  if (requireNamespace("Matrix", quietly = TRUE)) {
    npd <- try(Matrix::nearPD(S, corr = FALSE), silent = TRUE)
    if (!inherits(npd, "try-error")) {
      Snpd <- as.matrix(npd$mat)
      for (j in jitter_seq) {
        Sc <- Snpd + diag(j, p)
        if (.sop_is_chol_ok(Sc)) return(Sc)
      }
    }
  }

  # 2) shrink towards diag(S)
  D <- diag(diag(S), p)
  for (alpha in shrink_grid) {
    Ssh <- (1 - alpha) * S + alpha * D
    for (j in jitter_seq) {
      Sc <- Ssh + diag(j, p)
      if (.sop_is_chol_ok(Sc)) return(Sc)
    }
  }

  # 3) fall back to diagonal with jitter
  for (j in jitter_seq) {
    Sc <- diag(diag(S), p) + diag(max(j, 1e-8), p)
    if (.sop_is_chol_ok(Sc)) return(Sc)
  }

  NULL
}

#' Sanitize monitor roots (drop numerics and control tokens)
#'
#' Removes numeric-only tokens and simple control tokens such as \code{"thin"},
#' \code{"="}, \code{":"}, and \code{","}.
#'
#' @param x Character vector.
#'
#' @return Cleaned character vector.
#' @keywords internal
.sop_sanitize_roots <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[nzchar(x)]
  bad <- grepl("^\\d+(?:\\.\\d*)?$", x) | x %in% c("thin", "=", ":", ",")
  unique(x[!bad])
}

#' Correlation-based clustering into compact blocks
#'
#' Partitions variables into blocks using hierarchical clustering on a
#' correlation-based distance. Intended to build compact blocks for samplers.
#'
#' @param M Numeric matrix with variables in columns.
#' @param max_block Integer; maximum block size.
#' @param min_corr Reserved numeric argument for a minimal correlation target
#'   (currently not enforced).
#'
#' @return A list of character vectors (column names per block).
#' @keywords internal
.sop_split_into_blocks <- function(M, max_block = 40L, min_corr = 0.2) { # min_corr kept for API
  M <- as.matrix(M)
  if (is.null(colnames(M))) {
    colnames(M) <- paste0("V", seq_len(ncol(M)))
  }
  p <- ncol(M)
  if (p <= 1L) return(list(colnames(M)))
  if (p <= max_block) return(list(colnames(M)))

  keep <- which(apply(M, 2L, function(x) sum(is.finite(x)) > 1L))
  if (!length(keep)) {
    idx <- split(seq_len(p), ceiling(seq_along(seq_len(p)) / max_block))
    return(lapply(idx, function(ii) colnames(M)[ii]))
  }

  M2 <- M[, keep, drop = FALSE]
  R  <- try(stats::cor(M2, use = "pairwise.complete.obs"), silent = TRUE)
  if (inherits(R, "try-error") || !is.matrix(R)) {
    idx <- split(seq_len(p), ceiling(seq_along(seq_len(p)) / max_block))
    return(lapply(idx, function(ii) colnames(M)[ii]))
  }

  R <- pmax(pmin(R, 1), -1)
  diag(R) <- 1

  dist_mat <- sqrt(pmax(0, 0.5 * (1 - R)))
  d  <- try(stats::as.dist(dist_mat), silent = TRUE)
  hc <- try(stats::hclust(d, method = "average"), silent = TRUE)
  if (inherits(hc, "try-error")) {
    idx <- split(seq_len(p), ceiling(seq_along(seq_len(p)) / max_block))
    return(lapply(idx, function(ii) colnames(M)[ii]))
  }

  K   <- ceiling(ncol(M2) / max_block)
  grp <- stats::cutree(hc, k = K)
  blocks <- split(colnames(M2), grp)

  rest <- setdiff(colnames(M), colnames(M2))
  if (length(rest)) {
    idx <- split(seq_along(rest), ceiling(seq_along(rest) / max_block))
    blocks <- c(blocks, lapply(idx, function(ii) rest[ii]))
  }

  Filter(function(v) length(v) > 0L, blocks)
}

#' Strip "chainK." prefix from column names
#'
#' @param nms Character vector of column names.
#'
#' @return Character vector with chain prefixes removed.
#' @keywords internal
.sop_strip_chain_prefix <- function(nms) {
  sub("^chain\\d+\\.", "", nms, perl = TRUE)
}

# --- Derivatives capability probe -------------------------------------------

#' Check if derivatives/HMC are supported for a nimble model
#'
#' Uses model flags if available, otherwise tries a small \code{configureMCMC}
#' + \code{nimbleHMC::configureHMC()} call to detect derivative support.
#'
#' @param model Nimble model.
#'
#' @return Logical.
#' @keywords internal
.sop_supports_derivs <- function(model) {
  bd1 <- try(model$modelDef$buildDerivs, silent = TRUE)
  if (!inherits(bd1, "try-error") && !is.null(bd1)) {
    return(isTRUE(bd1))
  }
  bd2 <- try(model$buildDerivs, silent = TRUE)
  if (!inherits(bd2, "try-error") && !is.null(bd2)) {
    return(isTRUE(bd2))
  }

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

#' Expand a compact node specification into explicit indices (numeric only)
#'
#' For example, \code{"p[1, 1:3]"} is expanded to
#' \code{c("p[1,1]", "p[1,2]", "p[1,3]")}, provided all indices are numeric.
#'
#' @param node Single node string.
#'
#' @return Character vector of expanded nodes, or \code{node} if not expandable.
#' @keywords internal
.us_expand_compact_node <- function(node) {
  lb <- regexpr("\\[", node, fixed = TRUE)
  rb <- regexpr("\\]", node, fixed = TRUE)
  if (lb < 0 || rb < 0 || rb < lb) return(node)

  base   <- substr(node, 1L, lb - 1L)
  inside <- substr(node, lb + 1L, rb - 1L)
  dims   <- strsplit(inside, ",", fixed = TRUE)[[1]]
  dims   <- trimws(dims)

  dim_idx <- lapply(dims, function(d) {
    if (grepl(":", d, fixed = TRUE)) {
      ab <- strsplit(d, ":", fixed = TRUE)[[1]]
      ab <- trimws(ab)
      if (all(grepl("^-?\\d+$", ab))) {
        a <- as.integer(ab[1]); b <- as.integer(ab[2])
        if (is.na(a) || is.na(b)) return(d)
        if (a <= b) seq.int(a, b) else seq.int(a, b, by = -1L)
      } else {
        d
      }
    } else if (grepl("^-?\\d+$", d)) {
      as.integer(d)
    } else {
      d
    }
  })

  if (any(!vapply(dim_idx, is.integer, logical(1)))) return(node)

  grid <- do.call(
    expand.grid,
    c(dim_idx, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  )
  apply(
    grid,
    1L,
    function(row) paste0(base, "[", paste(row, collapse = ","), "]")
  )
}

#' Expand variable roots to explicit node names using model metadata
#'
#' Uses \code{getVarInfo()} and model node names to expand roots such as
#' \code{"beta"} into all corresponding node names present in the model.
#'
#' @param model Nimble model (R-level or compiled).
#' @param vars Character vector of variable roots.
#'
#' @return Character vector of node names present in the model.
#' @keywords internal
.us_expand_vars_to_nodes <- function(model, vars) {
  if (!length(vars)) return(character(0))
  stopifnot(.sop_is_model(model))

  m <- .sop_get_uncompiled_model(model)
  vars <- unique(trimws(as.character(vars)))
  vars <- vars[nzchar(vars)]

  all_nodes <- try(
    m$getNodeNames(stochOnly = FALSE, includeRHSonly = FALSE),
    silent = TRUE
  )
  if (inherits(all_nodes, "try-error") || is.null(all_nodes)) {
    all_nodes <- character(0)
  }

  var_to_nodes <- function(v) {
    vi <- try(m$getVarInfo(v), silent = TRUE)
    if (!inherits(vi, "try-error") && length(vi$nDim)) {
      if (vi$nDim == 0L) {
        cand <- v
      } else {
        mins <- as.integer(vi$mins)
        maxs <- as.integer(vi$maxs)
        if (length(mins) == vi$nDim && length(maxs) == vi$nDim) {
          idx  <- lapply(seq_len(vi$nDim), function(d) seq.int(mins[d], maxs[d]))
          grid <- do.call(
            expand.grid,
            c(idx, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
          )
          cand <- apply(
            grid,
            1L,
            function(row) paste0(v, "[", paste(row, collapse = ","), "]")
          )
        } else {
          cand <- character(0)
        }
      }
    } else {
      cand <- all_nodes[
        all_nodes == v | startsWith(all_nodes, paste0(v, "["))
      ]
      cand <- unlist(
        lapply(cand, .us_expand_compact_node),
        use.names = FALSE
      )
    }
    cand[cand %in% all_nodes]
  }

  unique(unlist(lapply(vars, var_to_nodes), use.names = FALSE))
}

#' Drop numeric junk and control tokens from monitors
#'
#' @param x Character vector.
#'
#' @return Cleaned character vector.
#' @keywords internal
.us_sanitize_monitors <- function(x) {
  x <- .us_tokenize_monitors(x)
  bad <- grepl("^\\d+(?:\\.\\d*)?$", x) | x %in% c("thin", "=", ":", ",")
  unique(x[!bad])
}

#' Tokenize a monitors specification into clean symbols
#'
#' @param x Character vector or scalar.
#'
#' @return Character vector of tokens, no empty strings.
#' @keywords internal
.us_tokenize_monitors <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  x <- gsub("[=:\\|;]", " ", x, perl = TRUE)
  parts <- unlist(strsplit(x, "[,\\s]+", perl = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

#' Convert various sample formats to \code{coda::mcmc.list}
#'
#' Optionally merges \code{samples2} into \code{samples}, aligning by row count
#' when needed and stripping log-likelihood columns if requested.
#'
#' @param samples An object containing MCMC draws (list, matrix, \code{mcmc},
#'   or \code{mcmc.list}).
#' @param samples2 Optional object of the same structure to be merged.
#' @param drop_loglik Logical; drop log-likelihood columns.
#' @param thin Integer; thinning factor for returned \code{mcmc} objects.
#'
#' @return A \code{coda::mcmc.list} object.
#' @keywords internal
as_mcmc_list_sop <- function(samples,
                             samples2    = NULL,
                             drop_loglik = FALSE,
                             thin        = 1L) {

  strip_loglik_cols <- function(M) {
    if (is.null(M)) return(NULL)
    keep <- !grepl(
      "^logLik(\\[.*\\])?$|log_?lik|logdens|lpdf",
      colnames(M),
      perl        = TRUE,
      ignore.case = TRUE
    )
    M[, keep, drop = FALSE]
  }

  # Case 1: mcmc.list
  if (inherits(samples, "mcmc.list")) {
    if (!is.null(samples2)) {
      if (!inherits(samples2, "mcmc.list")) {
        stop("samples2 must be an mcmc.list if samples is one.")
      }
      stopifnot(length(samples) == length(samples2))
      out <- coda::mcmc.list(lapply(seq_along(samples), function(i) {
        m1 <- as.matrix(samples[[i]])
        m2 <- as.matrix(samples2[[i]])
        if (drop_loglik) {
          m1 <- strip_loglik_cols(m1)
          m2 <- strip_loglik_cols(m2)
        }
        colnames(m1) <- .sop_strip_chain_prefix(colnames(m1))
        colnames(m2) <- .sop_strip_chain_prefix(colnames(m2))
        coda::mcmc(.sop_cbind_align(m1, m2), thin = thin)
      }))
      return(out)
    }
    return(samples)
  }

  # Case 2: single mcmc
  if (inherits(samples, "mcmc")) {
    return(coda::mcmc.list(samples))
  }

  # Case 3: list of matrices (chains)
  if (is.list(samples) &&
      length(samples) &&
      all(vapply(samples, is.matrix, TRUE))) {

    if (!is.null(samples2)) {
      if (!(is.list(samples2) &&
            length(samples2) == length(samples) &&
            all(vapply(samples2, is.matrix, TRUE)))) {
        stop("samples2 must be a list of matrices, with the same length as samples.")
      }
    }

    out <- coda::mcmc.list(lapply(seq_along(samples), function(i) {
      m1 <- samples[[i]]
      if (drop_loglik) m1 <- strip_loglik_cols(m1)
      if (!is.null(samples2)) {
        m2 <- samples2[[i]]
        if (drop_loglik) m2 <- strip_loglik_cols(m2)
        m  <- .sop_cbind_align(m1, m2)
      } else {
        m <- m1
      }
      colnames(m) <- .sop_strip_chain_prefix(colnames(m))
      coda::mcmc(m, thin = thin)
    }))
    return(out)
  }

  # Case 4: big matrix, possibly with chainK. prefixes
  if (is.matrix(samples)) {
    nms <- colnames(samples)
    cid <- .sop_detect_chain_prefix(nms)

    if (all(is.na(cid))) {
      m <- if (drop_loglik) strip_loglik_cols(samples) else samples
      colnames(m) <- .sop_strip_chain_prefix(colnames(m))
      return(coda::mcmc.list(coda::mcmc(m, thin = thin)))
    }

    K    <- max(cid, na.rm = TRUE)
    cid2 <- if (is.matrix(samples2)) .sop_detect_chain_prefix(colnames(samples2)) else NULL

    out_list <- vector("list", K)
    for (k in seq_len(K)) {
      keep1 <- which(cid == k)
      if (!length(keep1)) next
      mk1 <- samples[, keep1, drop = FALSE]
      if (drop_loglik) mk1 <- strip_loglik_cols(mk1)

      if (is.matrix(samples2)) {
        keep2 <- which(cid2 == k)
        mk2 <- if (length(keep2)) samples2[, keep2, drop = FALSE] else NULL
        if (!is.null(mk2) && drop_loglik) mk2 <- strip_loglik_cols(mk2)
        mk  <- if (!is.null(mk2)) .sop_cbind_align(mk1, mk2) else mk1
      } else {
        mk <- mk1
      }

      colnames(mk) <- .sop_strip_chain_prefix(colnames(mk))
      out_list[[k]] <- coda::mcmc(mk, thin = thin)
    }

    out_list <- Filter(Negate(is.null), out_list)
    if (!length(out_list)) {
      stop("Unable to extract chains from the sample matrix.")
    }
    return(coda::mcmc.list(out_list))
  }

  stop("Sample format not supported by as_mcmc_list_sop().")
}
#' Internal ggplot2 histogram helper
#'
#' Convenience wrapper to build a 1D histogram with an optional vertical line
#' and optional log10 transformation on the x-axis.
#'
#' Note: this function relies on a global object \code{bins} (number of bins),
#' which is typically set in the plotting context and declared via
#' \code{utils::globalVariables()}.
#'
#' @param df Data frame containing the variable to plot.
#' @param xvar Character scalar; name of the column in \code{df} to use on x.
#' @param title Plot title.
#' @param vline Optional numeric; if non-NULL, draws a vertical reference line.
#' @param log_x Logical; if TRUE, apply a \code{log10} transform to the x-axis.
#'
#' @return A \code{ggplot} object.
#' @keywords internal
ggh <- function(df, xvar, title, vline = NULL, log_x = FALSE) {
  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = xvar)) +
    ggplot2::geom_histogram(bins = bins, alpha = 0.9) +
    ggplot2::labs(title = title, x = xvar, y = "count")

  if (!is.null(vline)) {
    p <- p + ggplot2::geom_vline(xintercept = vline, linetype = 2)
  }
  if (isTRUE(log_x)) {
    p <- p + ggplot2::scale_x_continuous(trans = "log10")
  }
  p
}

#' Build a compact "top-3" summary table
#'
#' Extracts the top 3 rows of a data frame and keeps only a subset of
#' columns relevant for bottleneck summaries.
#'
#' @param df Data frame with at least columns \code{parameter}, \code{ESS},
#'   \code{CE}, \code{AE}, \code{slow_node_time}, \code{CE_rank},
#'   \code{AE_rank}, \code{TIME_rank} (only the intersection is kept).
#' @param axis_label Character scalar added as a new \code{axis} column.
#'
#' @return A data frame with up to 3 rows, or \code{NULL} if \code{df} is empty.
#' @keywords internal
mk_top3 <- function(df, axis_label) {
  if (is.null(df) || !nrow(df)) return(NULL)
  k   <- min(3L, nrow(df))
  out <- df[seq_len(k), , drop = FALSE]
  out$axis <- axis_label
  out[, intersect(
    c("axis", "parameter", "ESS", "CE", "AE", "slow_node_time",
      "CE_rank", "AE_rank", "TIME_rank"),
    names(out)
  ), drop = FALSE]
}

#' Select top-k rows given an ordering index
#'
#' Utility used in bottleneck summaries to pick the top-k rows from
#' a "valid" data frame according to a precomputed order.
#'
#' Note: this helper relies on the existence of global objects
#' \code{valid}, \code{CE_rank}, \code{AE_rank}, and \code{TIME_rank}
#' (typically created in the surrounding analysis function and declared
#' via \code{utils::globalVariables()}).
#'
#' @param ord Integer vector of row indices (order).
#' @param k   Integer scalar; desired number of rows.
#'
#' @return Subset of \code{valid} with rank columns attached.
#' @keywords internal
take <- function(ord, k) {
  idx <- ord[seq_len(min(k, length(ord)))]
  if (length(idx) == 0L) return(valid[integer(0), ])
  out <- valid[idx, , drop = FALSE]
  out$CE_rank   <- CE_rank[idx]
  out$AE_rank   <- AE_rank[idx]
  out$TIME_rank <- TIME_rank[idx]
  out
}

#' Coerce objects to \code{coda::mcmc.list}
#'
#' Lightweight helper to normalise various sample formats to a single
#' \code{mcmc.list} representation.
#'
#' @param x An object of class \code{mcmc.list}, \code{mcmc}, or something
#'   coercible to a data frame.
#'
#' @return An \code{mcmc.list} object, or \code{NULL} if \code{x} is \code{NULL}.
#' @keywords internal
to_mlist <- function(x) {
  if (inherits(x, "mcmc.list")) return(x)
  if (inherits(x, "mcmc"))      return(coda::mcmc.list(x))
  if (is.null(x)) return(NULL)
  coda::mcmc.list(coda::as.mcmc(as.data.frame(x)))
}
