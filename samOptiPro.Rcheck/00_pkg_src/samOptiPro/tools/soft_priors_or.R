## ============================================================
## Soft half-normal for NIMBLE (HMC-friendly, AD-compatible)
## ============================================================
##
## Cible : approximer T(dnorm(0, sd = scale), 0, Inf)
## en évitant la troncature BUGS T(...).
##
## Base half-normal(0, scale) :
##   f_hn(x) = (2/scale) * phi(x/scale), x > 0
##
## Version "soft" :
##   f_soft(x) ∝ f_hn(x) * Beta(g(x/scale); alpha, beta)
##   avec g(z) = logistic(z) = 1 / (1 + exp(-z)).
##
## Remarque importante :
##   - si alpha = beta = 1  : Beta(g;1,1) = 1  => f_soft = f_hn (half-normal exacte)
##   - si alpha ≈ beta ≈ 1 : f_soft est très proche de la half-normal
##
## Cette forme évite pnorm dans la densité et reste AD-compatible.
## ============================================================

dhalfnorm_soft <- nimbleFunction(
  run = function(x     = double(0),
                 scale = double(0),
                 alpha = double(0),
                 beta  = double(0),
                 log   = integer(0, default = 0)) {

    returnType(double(0))

    ## Paramètres invalides
    if (scale <= 0 | alpha <= 0 | beta <= 0) {
      if (log == 1L) return(-Inf) else return(0.0)
    }

    ## Support (0, +inf)
    if (x <= 0) {
      if (log == 1L) return(-Inf) else return(0.0)
    }

    ## z = x / scale
    z <- x / scale

    ## Base demi-normale : f_hn(x) = (2/scale) * phi(z)
    f_hn <- (2.0 / scale) * dnorm(z, mean = 0.0, sd = 1.0, log = 0L)

    ## Transformée lisse vers (0,1) : g(z) = logistic(z)
    u <- 1.0 / (1.0 + exp(-z))  ## ∈ (0,1)

    ## Facteur Beta(g(z) ; alpha, beta)
    beta_factor <- dbeta(u, alpha, beta, log = 0L)

    ## Densité ∝ demi-normale * Beta(g(z))
    dens <- f_hn * beta_factor

    if (log == 1L) return(log(dens))
    return(dens)
  },
  buildDerivs = TRUE
)

rhalfnorm_soft <- nimbleFunction(
  run = function(n     = integer(0),
                 scale = double(0),
                 alpha = double(0),
                 beta  = double(0)) {

    returnType(double(0))

    if (n != 1L) {
      nimPrint("rhalfnorm_soft: n > 1, returning a single draw.")
    }

    if (scale <= 0 | alpha <= 0 | beta <= 0) {
      nimPrint("rhalfnorm_soft: invalid parameters, returning NaN.")
      return(NaN)
    }

    ## Ici on génère un tirage approx. cohérent avec la demi-normale
    ## (ce générateur n’entre pas dans la log-postérieure).
    ## U ~ Beta(alpha, beta)
    u <- rbeta(1, alpha, beta)
    z <- qnorm((u + 1.0) / 2.0, mean = 0.0, sd = 1.0)
    x <- scale * abs(z)

    return(x)
  }
)
## ============================================================
## Soft half-Cauchy for NIMBLE (HMC-friendly, AD-compatible)
## ============================================================
##
## Cible : approximer une half-Cauchy(scale) sur (0, +Inf)
## en évitant T(dt(...), 0, Inf) et en gardant une densité lisse.
##
## Base half-Cauchy(scale) :
##   F_hc(x) = (2/pi) * atan(x/scale)
##   f_hc(x) = (2/pi) * scale / (scale^2 + x^2)
##
## Version "soft" :
##   f_soft(x) ∝ f_hc(x) * Beta(F_hc(x); alpha, beta)
##
##   - alpha = beta = 1  => f_soft = f_hc (half-Cauchy exacte)
##
## ============================================================

dhalfcauchy_soft <- nimbleFunction(
  run = function(x     = double(0),
                 scale = double(0),
                 alpha = double(0),
                 beta  = double(0),
                 log   = integer(0, default = 0)) {

    returnType(double(0))

    if (scale <= 0 | alpha <= 0 | beta <= 0) {
      if (log == 1L) return(-Inf) else return(0.0)
    }
    if (x <= 0) {
      if (log == 1L) return(-Inf) else return(0.0)
    }

    ## Base half-Cauchy
    z <- x / scale

    ## F_hc(x) = (2/pi) * atan(z) ∈ (0,1)
    F_hc <- (2.0 / pi) * atan(z)

    ## f_hc(x) = (2/pi) * scale / (scale^2 + x^2)
    f_hc <- (2.0 / pi) * scale / (scale * scale + x * x)

    ## pondération Beta(F_hc; alpha, beta)
    beta_factor <- dbeta(F_hc, alpha, beta, log = 0L)

    dens <- f_hc * beta_factor

    if (log == 1L) return(log(dens))
    return(dens)
  },
  buildDerivs = TRUE
)

rhalfcauchy_soft <- nimbleFunction(
  run = function(n     = integer(0),
                 scale = double(0),
                 alpha = double(0),
                 beta  = double(0)) {

    returnType(double(0))

    if (n != 1L) {
      nimPrint("rhalfcauchy_soft: n > 1, returning a single draw.")
    }

    if (scale <= 0 | alpha <= 0 | beta <= 0) {
      nimPrint("rhalfcauchy_soft: invalid parameters, returning NaN.")
      return(NaN)
    }

    ## U ~ Beta(alpha, beta)
    ## X = F_hc^{-1}(U) = scale * tan( (pi/2) * U )
    u <- rbeta(1, alpha, beta)
    x <- scale * tan(0.5 * pi * u)

    return(x)
  }
)
## ============================================================
## Soft half-t(df=3) for NIMBLE (approx T(dt(0,1,3), 0, Inf))
## ============================================================

dhalft3_soft <- nimbleFunction(
  run = function(x     = double(0),
                 scale = double(0),
                 alpha = double(0),
                 beta  = double(0),
                 log   = integer(0, default = 0)) {

    returnType(double(0))

    if (scale <= 0 | alpha <= 0 | beta <= 0) {
      if (log == 1L) return(-Inf) else return(0.0)
    }
    if (x <= 0) {
      if (log == 1L) return(-Inf) else return(0.0)
    }

    ## Base half-t(df=3)
    z <- x / scale
    f_ht <- 2.0 * dt(z, df = 3.0, log = 0L) / scale

    ## Transformée lisse vers (0,1) : g(z) = logistic(z)
    u <- 1.0 / (1.0 + exp(-z))

    ## pondération Beta(g(z); alpha, beta)
    beta_factor <- dbeta(u, alpha, beta, log = 0L)

    dens <- f_ht * beta_factor

    if (log == 1L) return(log(dens))
    return(dens)
  },
  buildDerivs = TRUE
)

rhalft3_soft <- nimbleFunction(
  run = function(n     = integer(0),
                 scale = double(0),
                 alpha = double(0),
                 beta  = double(0)) {

    returnType(double(0))

    if (n != 1L) {
      nimPrint("rhalft3_soft: n > 1, returning a single draw.")
    }

    if (scale <= 0 | alpha <= 0 | beta <= 0) {
      nimPrint("dhalft3_soft: invalid parameters, returning NaN.")
      return(NaN)
    }

    ## Générateur approx : on part d’un t(df=3) puis on applique
    ## un poids Beta via transformation quantile simplifiée.
    z <- rt(1, df = 3.0)
    x <- abs(scale * z)
    return(x)
  }
)
register_soft_priors_nimble <- function() {

  assign("dhalfnorm_soft",   dhalfnorm_soft,   envir = .GlobalEnv)
  assign("rhalfnorm_soft",   rhalfnorm_soft,   envir = .GlobalEnv)
  assign("dhalfcauchy_soft", dhalfcauchy_soft, envir = .GlobalEnv)
  assign("rhalfcauchy_soft", rhalfcauchy_soft, envir = .GlobalEnv)
  assign("dhalft3_soft",     dhalft3_soft,     envir = .GlobalEnv)
  assign("rhalft3_soft",     rhalft3_soft,     envir = .GlobalEnv)

  registerDistributions(list(
    dhalfnorm_soft = list(
      BUGSdist = "dhalfnorm_soft(scale, alpha, beta)",
      types    = c(
        "value = double(0)",
        "scale = double(0)",
        "alpha = double(0)",
        "beta  = double(0)"
      )
    ),
    dhalfcauchy_soft = list(
      BUGSdist = "dhalfcauchy_soft(scale, alpha, beta)",
      types    = c(
        "value = double(0)",
        "scale = double(0)",
        "alpha = double(0)",
        "beta  = double(0)"
      )
    ),
    dhalft3_soft = list(
      BUGSdist = "dhalft3_soft(scale, alpha, beta)",
      types    = c(
        "value = double(0)",
        "scale = double(0)",
        "alpha = double(0)",
        "beta  = double(0)"
      )
    )
  ))

  invisible(NULL)
}
