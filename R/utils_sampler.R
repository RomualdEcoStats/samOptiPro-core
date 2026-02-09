# utils_sampler.R -- helpers for NIMBLE (ASCII only)
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

# --- helpers internes robustes pour les monitors ------------------------

#' @keywords internal
# eclate "p[1, 1:3]" -> c("p[1,1]","p[1,2]","p[1,3]") (si bornes numeriques)

## ---- Exported functions (9) ----

#' Return a data.frame (name, type, target, scale) from a sampler container.
#' Accepts: samplerConf list (conf$getSamplers()) or other sampler-like lists.
#' @keywords internal
#' @noRd
sampler_df <- function(samplers) {
  L <- .as_list_safe(samplers)
  n <- length(L)
  nm <- paste0('sampler_', seq_len(n))
  types <- vapply(seq_len(n), function(i) .get_type(L[[i]]), character(1), USE.NAMES = FALSE)
  targets <- vapply(seq_len(n), function(i) .get_targets(L[[i]]), character(1), USE.NAMES = FALSE)
  scales  <- vapply(seq_len(n), function(i) sampler_scale(L[[i]]), numeric(1), USE.NAMES = FALSE)
  data.frame(name = nm, type = types, target = targets, scale = scales, stringsAsFactors = FALSE)
}


#' Convenience: directly return sampler_df from an existing conf.
#' @export
#' @keywords internal
sampler_df_from_conf <- function(conf) {
  sampler_df(conf$getSamplers())
}


#' Long data.frame of sampler env numeric fields
#' @param cm compiled MCMC
#' @return data.frame(sampler_index, field, value)
#' @export
#' @keywords internal
sampler_env_dump <- function(cm) {
  lst <- sampler_env_numeric_fields(cm)
  idx <- rep.int(seq_along(lst), vapply(lst, length, integer(1)))
  fld <- unlist(lapply(lst, names), use.names = FALSE)
  val <- unlist(lapply(lst, unlist), use.names = FALSE)
  if (!length(idx)) {
    return(data.frame(sampler_index = integer(0), field = character(0), value = numeric(0)))
  }
  data.frame(sampler_index = idx, field = fld, value = val, row.names = NULL)
}


#' Numeric scalar fields in the environment of each compiled sampler function
#' @param cm compiled MCMC (from compileNimble(buildMCMC(...)))
#' @return list of named lists (one per sampler)
#' @export
#' @keywords internal
sampler_env_numeric_fields <- function(cm) {
  out <- list()
  sfun <- tryCatch(as.list(cm$samplerFunctions), error = function(e) cm$samplerFunctions)
  for (i in seq_along(sfun)) {
    s <- sfun[[i]]
    fields <- list()
    e <- try(environment(s), silent = TRUE)
    if (!inherits(e, "try-error") && !is.null(e)) {
      nms <- ls(e, all.names = TRUE)
      for (nm in nms) {
        val <- try(get(nm, envir = e), silent = TRUE)
        if (!inherits(val, "try-error") && is.numeric(val) && length(val) == 1L && is.finite(val)) {
          fields[[nm]] <- val
        }
      }
    }
    out[[i]] <- fields
  }
  out
}


#' Build and compile an MCMC from a conf to access samplerFunctions.
#' @export
#' @keywords internal
sampler_functions_from_conf <- function(conf, cmodel) {
  m  <- nimble::buildMCMC(conf)
  cm <- nimble::compileNimble(m, project = cmodel, resetFunctions = TRUE)
  cm$samplerFunctions
}


#' Extract a numeric scale proxy from a sampler object.
#' For RW: control$scale ; For slice: control$width (returned as 'scale').
#' Fallbacks: $scale (numeric or function) or environment(s)$scale.
#' @export
#' @keywords internal
sampler_scale <- function(s) {
  val <- NA_real_
  try({
    ctrl <- s$control
    if (!is.null(ctrl) && !is.null(ctrl$scale)) val <- ctrl$scale
  }, silent = TRUE)
  if (is.na(val)) try({
    ctrl <- s$control
    tp   <- try(s$type, silent = TRUE)
    if (!inherits(tp, 'try-error') && !is.null(tp) && identical(tolower(tp), 'slice')) {
      if (!is.null(ctrl) && !is.null(ctrl$width)) val <- ctrl$width
    }
  }, silent = TRUE)
  if (is.na(val)) try({ sc <- s$scale; if (!is.null(sc)) val <- if (is.function(sc)) sc() else sc }, silent = TRUE)
  if (is.na(val)) try({ e <- try(environment(s), silent = TRUE); if (!inherits(e, 'try-error') && !is.null(e) && !is.null(e$scale)) val <- e$scale }, silent = TRUE)
  val
}


#' Extract sampler scales after a short adaptive run.
#' Runs a short MCMC to let samplers adapt, then inspects samplerFunctions envs.
#' @param build_fn function returning list(model, cmodel, monitors)
#' @param niter total iterations (default 2000)
#' @param nburnin burn-in (default 500)
#' @param thin thinning (kept for API symmetry; not used here)
#' @return data.frame(name, type, target, scale)
#' @export
#' @keywords internal
sampler_scales_after_run <- function(build_fn, niter = 2000, nburnin = 500, thin = 1) {
  built <- build_fn()
  conf  <- .us_configure_with_monitors(built$model, built$monitors %||% NULL)
  df_conf <- sampler_df_from_conf(conf)
  m  <- nimble::buildMCMC(conf)
  cm <- nimble::compileNimble(m, project = built$cmodel, resetFunctions = TRUE)
  cm$run(niter)
  n <- nrow(df_conf)
  env_scale <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    v <- NA_real_
    s <- try(cm$samplerFunctions[[i]], silent = TRUE)
    if (!inherits(s, 'try-error') && !is.null(s)) {
      e <- try(environment(s), silent = TRUE)
      if (!inherits(e, 'try-error') && !is.null(e)) {
        if (!is.null(e$scale))                 v <- e$scale
        if (is.na(v) && !is.null(e$currentScale)) v <- e$currentScale
        if (is.na(v) && !is.null(e$propScale))    v <- e$propScale
      }
    }
    env_scale[i] <- v
  }
  df_conf$scale <- ifelse(is.na(env_scale), df_conf$scale, env_scale)
  df_conf
}


#' Preferred: read scales from a fresh samplerConf list (pre-run).
#' @export
#' @keywords internal
sampler_scales_from_run <- function(run_result, build_fn) {
  built <- build_fn()
  # ? NE PAS passer monitors=... directement a configureMCMC : on pose explicitement
  conf  <- .us_configure_with_monitors(built$model, built$monitors %||% NULL)
  sampler_df(conf$getSamplers())
}


#' Targets (node names) for each sampler in a MCMC configuration
#' @param conf a nimble MCMC configuration (from configureMCMC)
#' @return character vector (length = length(conf$getSamplers()))
#' @export
#' @keywords internal
sampler_targets <- function(conf) {
  s <- conf$getSamplers()
  vapply(s, function(si) {
    t <- try(si$target, silent = TRUE)
    if (inherits(t, "try-error") || is.null(t)) "" else as.character(t)
  }, character(1))
}


## ---- Internal functions (10) ----

#' @keywords internal
.as_list_safe <- function(x) {
  out <- tryCatch(as.list(x), error = function(e) x)
  if (!is.list(out)) out <- list(out)
  out
}


#' @keywords internal
.get_targets <- function(s) {
  tgt <- NA_character_
  try({
    ctrl <- s$control
    if (!is.null(ctrl)) {
      if (!is.null(ctrl$target))      tgt <- paste(ctrl$target, collapse = ',')
      if (!is.null(ctrl$targetNode))  tgt <- paste(ctrl$targetNode, collapse = ',')
      if (!is.null(ctrl$targetNodes)) tgt <- paste(ctrl$targetNodes, collapse = ',')
    }
  }, silent = TRUE)
  if (is.na(tgt)) try({ tgt <- paste(s$target, collapse = ',') }, silent = TRUE)
  tgt
}


#' @keywords internal
.get_type <- function(s) {
  t <- NA_character_
  try({ t <- s$type }, silent = TRUE)
  if (is.null(t) || is.na(t)) try({ t <- class(s)[1] }, silent = TRUE)
  t
}


#' @keywords internal
.us_configure_with_monitors <- function(model, monitors = NULL) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  sanitize_roots <- function(x) {
    if (is.null(x)) return(character(0))
    x <- as.character(x); x <- trimws(x); x <- x[nzchar(x)]
    bad <- grepl("^\\d+(?:\\.\\d*)?$", x) | x %in% c("thin","=",":",",")
    unique(x[!bad])
  }

  mons_in <- sanitize_roots(monitors %||% character(0))

  is_loglike <- grepl("^logLik(\\[.*\\])?$|log_?lik|logdens|lpdf",
                      mons_in, ignore.case = TRUE, perl = TRUE)
  roots_main <- mons_in[!is_loglike]
  roots_ll   <- mons_in[ is_loglike]

  all_vars <- unique(sub("\\[.*\\]$", "", model$getNodeNames(stochOnly = FALSE, includeData = FALSE)))
  roots_main <- roots_main[roots_main %in% all_vars]
  roots_ll   <- roots_ll  [roots_ll   %in% all_vars]

  conf <- nimble::configureMCMC(model)
  try(silent = TRUE, conf$clearMonitors()); try(silent = TRUE, conf$clearMonitors2())

  # DEBUG anti "1"
  cat("DEBUG(us) roots_main =", paste(roots_main, collapse="|"), "\n")
  stopifnot(!any(roots_main %in% c("1","thin","=",":",",")))
  cat("DEBUG(us) roots_ll   =", paste(roots_ll, collapse="|"), "\n")

  nimble::configureMCMC(model, monitors = nodes, monitors2 = nodes)


  attr(conf, "._us_monitors")  <- roots_main
  attr(conf, "._us_monitors2") <- roots_ll
  conf
}


.us_expand_compact_node <- function(node) {
  lb <- regexpr("\\[", node, fixed = TRUE)
  rb <- regexpr("\\]", node, fixed = TRUE)
  if (lb < 0 || rb < 0 || rb < lb) return(node)
  base <- substr(node, 1L, lb - 1L)
  inside <- substr(node, lb + 1L, rb - 1L)
  dims <- strsplit(inside, ",", fixed = TRUE)[[1]]
  dims <- trimws(dims)
  dim_idx <- lapply(dims, function(d) {
    if (grepl(":", d, fixed = TRUE)) {
      ab <- strsplit(d, ":", fixed = TRUE)[[1]]; ab <- trimws(ab)
      if (all(grepl("^-?\\d+$", ab))) {
        a <- as.integer(ab[1]); b <- as.integer(ab[2])
        if (is.na(a) || is.na(b)) return(d)
        if (a <= b) seq.int(a, b) else seq.int(a, b, by = -1L)
      } else d
    } else if (grepl("^-?\\d+$", d)) {
      as.integer(d)
    } else d
  })
  if (any(!vapply(dim_idx, is.integer, logical(1)))) return(node)
  grid <- do.call(expand.grid, c(dim_idx, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  apply(grid, 1L, function(row) paste0(base, "[", paste(row, collapse=","), "]"))
}


.us_expand_vars_to_nodes <- function(model, vars) {
  if (!length(vars)) return(character(0))
  vars <- unique(trimws(as.character(vars))); vars <- vars[nzchar(vars)]

  # all known nodes, for final filtering
  all_nodes <- try(model$getNodeNames(stochOnly = FALSE, includeRHSonly = FALSE), silent = TRUE)
  if (inherits(all_nodes, "try-error") || is.null(all_nodes)) all_nodes <- character(0)

  # try to use shape metadata (min/max) when available
  var_to_nodes <- function(v) {
    vi <- try(model$getVarInfo(v), silent = TRUE)
    if (!inherits(vi, "try-error") && !is.null(vi) && length(vi$nDim)) {
      ndim <- vi$nDim
      if (ndim == 0L) {
        cand <- v
      } else {
        mins <- as.integer(vi$mins); maxs <- as.integer(vi$maxs)
        if (length(mins) == ndim && length(maxs) == ndim) {
          idx <- lapply(seq_len(ndim), function(d) seq.int(mins[d], maxs[d]))
          grid <- do.call(expand.grid, c(idx, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
          cand <- apply(grid, 1L, function(row) paste0(v, "[", paste(row, collapse=","), "]"))
        } else cand <- character(0)
      }
    } else {
      # fallback: takes everything that starts with “v[” or is equal to “v”
      cand <- all_nodes[ all_nodes == v | startsWith(all_nodes, paste0(v, "[")) ]
      # break up compact ratings if necessary
      cand <- unlist(lapply(cand, .us_expand_compact_node), use.names = FALSE)
    }
    # keep only what the model knows
    cand[cand %in% all_nodes]
  }

  out <- unique(unlist(lapply(vars, var_to_nodes), use.names = FALSE))
  out
}


#' @keywords internal
.us_sanitize_monitors <- function(x) {
  x <- .us_tokenize_monitors(x)
  # remove bare numbers (e.g., “1,” “2.”, “3.0”) AND log tokens
  bad <- grepl("^\\d+(?:\\.\\d*)?$", x) | x %in% c("thin", "=", ":", ",")
  unique(x[!bad])
}


#' @keywords internal
.us_tokenize_monitors <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  # replace common separators with spaces
  x <- gsub("[=:\\|;]", " ", x, perl = TRUE)
  # split on comma OR contiguous spaces
  parts <- unlist(strsplit(x, "[,\\s]+", perl = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts[nzchar(parts)]
}


  sanitize_roots <- function(x) {
    if (is.null(x)) return(character(0))
    x <- as.character(x); x <- trimws(x); x <- x[nzchar(x)]
    bad <- grepl("^\\d+(?:\\.\\d*)?$", x) | x %in% c("thin","=",":",",")
    unique(x[!bad])
  }


  var_to_nodes <- function(v) {
    vi <- try(model$getVarInfo(v), silent = TRUE)
    if (!inherits(vi, "try-error") && !is.null(vi) && length(vi$nDim)) {
      ndim <- vi$nDim
      if (ndim == 0L) {
        cand <- v
      } else {
        mins <- as.integer(vi$mins); maxs <- as.integer(vi$maxs)
        if (length(mins) == ndim && length(maxs) == ndim) {
          idx <- lapply(seq_len(ndim), function(d) seq.int(mins[d], maxs[d]))
          grid <- do.call(expand.grid, c(idx, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
          cand <- apply(grid, 1L, function(row) paste0(v, "[", paste(row, collapse=","), "]"))
        } else cand <- character(0)
      }
    } else {
      # fallback: takes everything that starts with “v[” or is equal to “v”
      cand <- all_nodes[ all_nodes == v | startsWith(all_nodes, paste0(v, "[")) ]
      # break up compact ratings if necessary
      cand <- unlist(lapply(cand, .us_expand_compact_node), use.names = FALSE)
    }
    # keep only what the model knows
    cand[cand %in% all_nodes]
  }

  #' Extract sampled-only monitors from a configureMCMC object
  #'
  #' @param conf nimble::MCMCconf returned by nimble::configureMCMC()
  #' @param monitors character vector of base node names (your monitor vector)
  #' @param strict logical; if TRUE, error if none are sampled
  #' @param keep_order logical; if TRUE keep original order
  #' @return character vector monitors_sampled with attributes for diagnostics
  #' @export
  monitors_sampled_from_conf <- function(conf,
                                         monitors,
                                         strict = FALSE,
                                         keep_order = TRUE) {

    stopifnot(is.character(monitors))
    if (!length(monitors)) stop("`monitors` is empty.")

    # NIMBLE conf is often a reference object, not a list
    if (is.null(conf) || is.null(conf$getSamplers) || !is.function(conf$getSamplers)) {
      stop("`conf` must be a NIMBLE configureMCMC object exposing `$getSamplers()`.")
    }

    samplers <- conf$getSamplers()
    if (!length(samplers)) stop("No samplers found in the MCMC configuration.")

    # Robust extraction of 'target' across possible sampler object types
    get_target <- function(s) {
      if (is.null(s)) return(NULL)
      # Most common: list-like with $target
      if (!is.null(s$target)) return(s$target)
      # Sometimes: list-like with [["target"]]
      if (!is.null(s[["target"]])) return(s[["target"]])
      # Rare: accessor method
      if (!is.null(s$getTarget) && is.function(s$getTarget)) return(s$getTarget())
      NULL
    }

    targets_list <- lapply(samplers, get_target)
    sampler_targets <- unique(unlist(targets_list, use.names = FALSE))

    if (!length(sampler_targets) || all(is.na(sampler_targets))) {
      stop("Sampler targets could not be extracted from `conf$getSamplers()` output.")
    }

    # Base-node reduction: "x[1,2]" -> "x"
    sampler_nodes <- unique(sub("\\[.*$", "", sampler_targets))

    out <- monitors[monitors %in% sampler_nodes]
    if (!keep_order) out <- sort(out)

    if (isTRUE(strict) && !length(out)) {
      stop("None of the requested monitors are sampled according to this MCMC configuration.")
    }

    attr(out, "sampler_targets") <- sampler_targets
    attr(out, "sampler_nodes")   <- sampler_nodes
    attr(out, "dropped")         <- setdiff(monitors, out)
    attr(out, "n_kept")          <- length(out)
    attr(out, "n_dropped")       <- length(setdiff(monitors, out))

    out
  }


  #' One-pass build (uncompiled) + configureMCMC + monitors_sampled
  #'
  #' @param monitors character vector of requested monitors
  #' @param chain_id integer(1)
  #' @param strict logical
  #' @param keep_order logical
  #' @param buildDerivs logical passed to nimbleModel
  #' @return list(model=m, conf=conf.mcmc, monitors_sampled=..., chain_id=...)
  #' @export
  monitors_sampled_safe <- function(monitors,
                                    chain_id = 1L,
                                    strict = FALSE,
                                    keep_order = TRUE,
                                    buildDerivs = TRUE) {

    stopifnot(is.character(monitors))
    chain_id <- as.integer(chain_id)

    if (!exists("model.nimble", inherits = TRUE)) stop("Object `model.nimble` not found.")
    if (!exists("Const_nimble", inherits = TRUE)) stop("Object `Const_nimble` not found.")
    if (!exists("Data_nimble", inherits = TRUE))  stop("Object `Data_nimble` not found.")
    if (!exists("inits_nimble", inherits = TRUE)) stop("Object `inits_nimble` not found.")

    inits <- get("inits_nimble", inherits = TRUE)
    stopifnot(is.list(inits), chain_id >= 1L, chain_id <= length(inits))

    m <- nimble::nimbleModel(
      code        = get("model.nimble",  inherits = TRUE),
      name        = sprintf("Model (chain %d)", chain_id),
      constants   = get("Const_nimble",  inherits = TRUE),
      data        = get("Data_nimble",   inherits = TRUE),
      inits       = inits[[chain_id]],
      buildDerivs = isTRUE(buildDerivs)
    )
    m$initializeInfo()

    conf.mcmc <- nimble::configureMCMC(m)

    monitors_sampled <- monitors_sampled_from_conf(
      conf       = conf.mcmc,
      monitors   = monitors,
      strict     = strict,
      keep_order = keep_order
    )

    list(
      model            = m,
      conf             = conf.mcmc,
      monitors_sampled = monitors_sampled,
      chain_id         = chain_id
    )
  }

