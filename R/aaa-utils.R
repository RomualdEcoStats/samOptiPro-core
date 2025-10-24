time_it <- function(expr){s<-proc.time()[['elapsed']]; res<-eval.parent(substitute(expr)); e<-proc.time()[['elapsed']]; list(result=res, runtime_s=as.numeric(e-s))}
# utils.R (par ex.)
#' @export
.sop_max <- function(x) {
  if (length(x) == 0L) return(NA_real_)
  if (all(is.na(x)))   return(NA_real_)
  max(x, na.rm = TRUE)
}
