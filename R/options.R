# options.R -- package-wide options helper for samOptiPro

#' samOptiPro options constructor
#'
#' Build a consistent options list used internally by samOptiPro.
#'
#' This helper gathers three sources of information, in increasing order
#' of priority:
#' \enumerate{
#'   \item Package defaults (hard-coded in this function),
#'   \item Global R option \code{"samOptiPro.options"} (a named list),
#'   \item Explicit arguments passed to \code{samOptiPro_options()}.
#' }
#'
#' The returned object is a named list. Unknown fields are kept as-is so
#' that new components can be added in the future without breaking the API.
#'
#' Typical consumers only need a few fields:
#' \itemize{
#'   \item \code{include_data}      : logical; whether to include data nodes
#'         when discovering monitors (rarely needed; default \code{FALSE}).
#'   \item \code{include_logLik}    : logical; whether to include a
#'         \code{"logLik"} variable in default monitors (default \code{FALSE}).
#'   \item \code{extra_monitors}    : character vector of additional variable
#'         roots to monitor in the main samples.
#'   \item \code{extra_monitors2}   : character vector of additional variable
#'         roots to monitor in the secondary samples (if any).
#' }
#'
#' Internal functions such as \code{default_monitors()},
#' \code{.configure_with_monitors()}, or the various runner helpers are
#' expected to accept a generic list \code{opts} and use only the fields
#' they understand (ignoring unknown components).
#'
#' @param include_data logical; default \code{FALSE}. If \code{TRUE}, data
#'   nodes may be included when discovering monitors.
#' @param include_logLik logical; default \code{FALSE}. If \code{TRUE}, a
#'   variable called \code{"logLik"} is added to default monitors when present.
#' @param extra_monitors character vector of additional variable roots to
#'   monitor in the main samples. Default \code{NULL}.
#' @param extra_monitors2 character vector of additional variable roots to
#'   monitor in the secondary samples. Default \code{NULL}.
#' @param ... Named components that either extend or override entries in the
#'   options list. This is mainly intended for future extensions and for
#'   advanced users.
#'
#' @return A named list of options, to be passed around as \code{opts}.
#' @examples
#' # Default options
#' samOptiPro_options()
#'
#' # Override a few fields explicitly
#' samOptiPro_options(include_logLik = TRUE,
#'                    extra_monitors = c("logLik", "sigma_proc"))
#'
#' # Use a global R option and then refine locally
#' old <- options(samOptiPro.options = list(include_logLik = TRUE))
#' samOptiPro_options()  # will have include_logLik = TRUE
#' samOptiPro_options(include_logLik = FALSE)  # local override
#' options(old)
#'
#' @export
samOptiPro_options <- function(include_data      = FALSE,
                               include_logLik    = FALSE,
                               extra_monitors    = NULL,
                               extra_monitors2   = NULL,
                               ...) {
  # 1) Hard-coded defaults used by the package
  defaults <- list(
    include_data    = include_data,
    include_logLik  = include_logLik,
    extra_monitors  = extra_monitors,
    extra_monitors2 = extra_monitors2
  )

  # 2) Global options, if any: options(samOptiPro.options = list(...))
  global_opts <- getOption("samOptiPro.options", NULL)
  if (is.list(global_opts) && length(global_opts)) {
    for (nm in names(global_opts)) {
      defaults[[nm]] <- global_opts[[nm]]
    }
  }

  # 3) Explicit arguments in ... have highest priority
  dots <- list(...)
  if (length(dots)) {
    # If names are missing, we ignore them silently
    has_name <- nzchar(names(dots))
    if (any(has_name)) {
      for (nm in names(dots)[has_name]) {
        defaults[[nm]] <- dots[[nm]]
      }
    }
  }

  defaults
}
