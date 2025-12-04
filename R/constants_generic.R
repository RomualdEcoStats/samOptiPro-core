# constants_generic.R -- model-agnostic constants normalization
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Normalize constants in a model-agnostic way
#' @param Const list of constants
#' @param n_rule "min"|"nyear"|"max"
#' @param pad_rule "repeat_last"|"zero"|"NA"
#' @param trim_rule logical (garder TRUE)
#' @return constants list normalized (with integer 'n')
#' @export
#' @keywords internal
normalize_constants_generic <- function(Const,
                                        n_rule = "min",
                                        pad_rule = "repeat_last",
                                        trim_rule = TRUE) {
  C <- Const
  lens <- .collect_lengths(C)
  if (!length(lens)) stop("normalize_constants_generic: cannot infer series lengths from constants.")

  if (!is.null(C$nyear) && n_rule == "nyear") {
    n <- as.integer(C$nyear)
  } else if (n_rule == "max") {
    n <- as.integer(max(lens))
  } else {
    n <- as.integer(min(lens))
  }
  C$n <- n

  # clean unfinished scalars
  for (f in names(C)) if (is.numeric(C[[f]]) && length(C[[f]]) == 1) {
    if (is.na(C[[f]]) || !is.finite(C[[f]])) C[[f]] <- 0
  }

  # crop vector and matrix
  for (nm in names(C)) {
    x <- C[[nm]]
    if (is.numeric(x) || is.integer(x)) {
      if (length(x) > 1) C[[nm]] <- .pad_vec(x, n, pad_rule)
    } else if (is.matrix(x)) {
      C[[nm]] <- .pad_mat_cols(x, n, pad_rule)
    }
  }

  C
}


## ---- Internal functions (4) ----

#' @keywords internal
.collect_lengths <- function(Const) {
  lens <- integer()
  push <- function(v) { if (length(v)) return(length(v)); NA_integer_ }
  for (nm in names(Const)) {
    x <- Const[[nm]]
    if (is.numeric(x) || is.integer(x)) {
      lens <- c(lens, push(x))
    } else if (is.matrix(x)) {
      lens <- c(lens, ncol(x))
    }
  }
  lens[is.finite(lens) & lens > 0]
}


#' @keywords internal
.pad_mat_cols <- function(M, n, pad_rule) {
  if (!is.matrix(M)) M <- as.matrix(M)
  p <- ncol(M)
  if (p == n) return(M)
  if (p >  n) return(M[, seq_len(n), drop = FALSE])
  last <- M[, p, drop = FALSE]
  add  <- if (pad_rule == "NA") matrix(NA_real_, nrow(M), n - p)
  else if (pad_rule == "zero") matrix(0, nrow(M), n - p)
  else last[, rep(1, n - p), drop = FALSE]
  cbind(M, add)
}


#' @keywords internal
.pad_vec <- function(v, n, pad_rule) {
  if (length(v) == n) return(v)
  if (length(v) >  n) return(v[seq_len(n)])
  if (length(v) == 0) {
    if (pad_rule == "NA")   return(rep(NA_real_, n))
    if (pad_rule == "zero") return(rep(0, n))
    return(rep(0, n))
  }
  last <- utils::tail(v, 1)
  if (pad_rule == "NA")   return(c(v, rep(NA_real_, n - length(v))))
  if (pad_rule == "zero") return(c(v, rep(0,        n - length(v))))
  c(v, rep(last, n - length(v))) # repeat_last
}


  push <- function(v) { if (length(v)) return(length(v)); NA_integer_ }


