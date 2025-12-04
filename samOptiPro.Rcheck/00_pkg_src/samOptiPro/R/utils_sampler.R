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

  # tous les noeuds connus, pour filtrage final
  all_nodes <- try(model$getNodeNames(stochOnly = FALSE, includeRHSonly = FALSE), silent = TRUE)
  if (inherits(all_nodes, "try-error") || is.null(all_nodes)) all_nodes <- character(0)

  # essaie d'utiliser la metadonnee de forme (mins/maxs) quand dispo
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
      # fallback: prend tout ce qui commence par "v[" ou est egal a "v"
      cand <- all_nodes[ all_nodes == v | startsWith(all_nodes, paste0(v, "[")) ]
      # eclate les notations compactes si necessaire
      cand <- unlist(lapply(cand, .us_expand_compact_node), use.names = FALSE)
    }
    # garde seulement ce que le modele connait
    cand[cand %in% all_nodes]
  }

  out <- unique(unlist(lapply(vars, var_to_nodes), use.names = FALSE))
  out
}


#' @keywords internal
.us_sanitize_monitors <- function(x) {
  x <- .us_tokenize_monitors(x)
  # vire nombres nus (ex "1", "2.", "3.0") ET tokens de log
  bad <- grepl("^\\d+(?:\\.\\d*)?$", x) | x %in% c("thin", "=", ":", ",")
  unique(x[!bad])
}


#' @keywords internal
.us_tokenize_monitors <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  # remplace separateurs courants par des espaces
  x <- gsub("[=:\\|;]", " ", x, perl = TRUE)
  # split sur virgule OU espaces contigus
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
      # fallback: prend tout ce qui commence par "v[" ou est egal a "v"
      cand <- all_nodes[ all_nodes == v | startsWith(all_nodes, paste0(v, "[")) ]
      # eclate les notations compactes si necessaire
      cand <- unlist(lapply(cand, .us_expand_compact_node), use.names = FALSE)
    }
    # garde seulement ce que le modele connait
    cand[cand %in% all_nodes]
  }


