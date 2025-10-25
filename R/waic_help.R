#' Compute WAIC from coda samples that include logLik[] columns
#' @param samples coda::mcmc.list
#' @export
#' @keywords internal
compute_WAIC <- function(samples) {
  stopifnot(requireNamespace("coda", quietly = TRUE))
  M <- as.matrix(samples[[1]])
  cols <- grep("^logLik\\[", colnames(M))
  if (!length(cols)) stop("No logLik[...] columns found.")
  ll <- M[, cols, drop = FALSE]
  lppd <- sum(log(colMeans(exp(ll), na.rm = TRUE)))
  p_waic <- 2 * sum(log(colMeans(exp(ll), na.rm = TRUE)) - colMeans(ll, na.rm = TRUE))
  waic <- -2 * (lppd - p_waic)
  list(summary = data.frame(WAIC = waic, lppd = lppd, p_waic = p_waic),
       pointwise = data.frame(loglik_col = colnames(ll),
                              lppd = log(colMeans(exp(ll), na.rm = TRUE)),
                              p_waic = 2 * (log(colMeans(exp(ll), na.rm = TRUE)) - colMeans(ll, na.rm = TRUE))))
}
