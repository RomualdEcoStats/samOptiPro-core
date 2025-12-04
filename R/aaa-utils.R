# utils.R (par ex.)

## ---- Internal helpers: max / timing ----

#' Internal numeric max helper
#'
#' Small wrapper around [base::max()] that is NA-safe and returns `NA_real_`
#' for empty or all-NA inputs. Used internally in diagnostic summaries.
#'
#' @param x Numeric vector.
#'
#' @return A length-one numeric value. Returns `NA_real_` if `x` is empty
#'   or if all entries are `NA`.
#'
#' @keywords internal
.sop_max <- function(x) {
  if (length(x) == 0L) return(NA_real_)
  if (all(is.na(x)))   return(NA_real_)
  max(x, na.rm = TRUE)
}

#' Internal helper to time an expression
#'
#' Evaluate an expression in the parent frame and return both the result
#' and the elapsed time in seconds. Intended for internal profiling of
#' MCMC runs and diagnostic steps.
#'
#' @param expr Expression to be evaluated.
#'
#' @return A list with two components:
#' \describe{
#'   \item{result}{The evaluated value of \code{expr}.}
#'   \item{runtime_s}{Elapsed time in seconds (numeric scalar).}
#' }
#'
#' @keywords internal
time_it <- function(expr) {
  s <- proc.time()[["elapsed"]]
  res <- eval.parent(substitute(expr))
  e <- proc.time()[["elapsed"]]
  list(result = res, runtime_s = as.numeric(e - s))
}



