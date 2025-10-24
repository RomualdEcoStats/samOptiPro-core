## Rscript/0_generate_inits_base.R
## Inits robustes: tolérant aux constantes manquantes/hors bornes, sans NA

if (!requireNamespace("nimble", quietly = TRUE)) {
  stop("Le package 'nimble' est requis pour générer des inits.")
}

.base_name <- function(x) sub("\\[.*\\]$", "", x)

.get_consts <- function(constants = NULL) {
  if (!is.null(constants)) return(constants)
  if (exists("constants", inherits = TRUE)) return(get("constants", inherits = TRUE))
  if (exists("Const_nimble", inherits = TRUE)) return(get("Const_nimble", inherits = TRUE))
  stop("Aucune liste de constantes trouvée (ni `constants` ni `Const_nimble`).")
}

# borne et remplace NA/inf
.clamp01 <- function(x, eps = 1e-6) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return(0.5)
  x <- max(eps, min(1 - eps, x))
  x
}
.pos_num <- function(x, default, eps = 1e-6) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x) || x <= 0) return(max(default, eps))
  x
}
.num <- function(x, default = 0) {
  x <- suppressWarnings(as.numeric(x))
  if (!is.finite(x)) return(default)
  x
}

# Transforme un vecteur nommé (values(model, nodes)) -> liste par familles
.values_to_inits_list <- function(named_vec) {
  if (length(named_vec) == 0) return(list())
  fams <- .base_name(names(named_vec))
  spl  <- split(named_vec, fams)
  lapply(spl, function(v) if (length(v) > 1) as.numeric(v) else as.numeric(v[[1]]))
}

make_inits <- function(constants = NULL,
                       sizes = list(N3_tot = 26L, N2 = 24L, N6 = 26L, N9 = 27L,
                                    theta3 = 25L, theta4f = 26L, theta4m = 26L, theta4 = 26L,
                                    theta5 = 26L, theta8 = 26L, prop3f = 25L, prop6f = 26L, prop9f = 27L,
                                    p_smolt = 24L)) {
  Const <- .get_consts(constants)
  
  # --- assainissement des hyperparamètres
  theta1_max <- .clamp01( Const$theta1_max %||% Const$theta1.prior.max %||% 0.8 )
  nsample    <- .pos_num( Const$nsample_theta1 %||% Const$theta1.nsample %||% 50, default = 50 )
  a_shape    <- max(1e-6, theta1_max * nsample)
  b_shape    <- max(1e-6, (1 - theta1_max) * nsample)
  
  logk_pr    <- .num(Const$logk_pr %||% Const$k.log.prior %||% 0, default = 0)
  sdlog_k    <- .pos_num(Const$sdlog_k %||% 1, default = 1)
  
  # petite fonction infix pour valeurs par défaut
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # logit = qlogis
  logit <- stats::qlogis
  
  # tailles
  n_N3_tot  <- as.integer(sizes$N3_tot)
  n_N2      <- as.integer(sizes$N2)
  n_N6      <- as.integer(sizes$N6)
  n_N9      <- as.integer(sizes$N9)
  n_theta3  <- as.integer(sizes$theta3)
  n_theta4f <- as.integer(sizes$theta4f)
  n_theta4m <- as.integer(sizes$theta4m)
  n_theta4  <- as.integer(sizes$theta4)
  n_theta5  <- as.integer(sizes$theta5)
  n_theta8  <- as.integer(sizes$theta8)
  n_prop3f  <- as.integer(sizes$prop3f)
  n_prop6f  <- as.integer(sizes$prop6f)
  n_prop9f  <- as.integer(sizes$prop9f)
  n_psmolt  <- as.integer(sizes$p_smolt)
  
  # Générations robustes
  fill_psmolt <- runif(n_psmolt, 0.4, 0.6)
  
  # sds positives
  sd_N <- 1
  big_mu <- log(1e6)
  
  # formes de beta > 0 (déjà bornées ci-dessus)
  alpha_draw <- suppressWarnings(rbeta(1, a_shape, b_shape))
  if (!is.finite(alpha_draw)) {
    warning("rbeta a produit NA (a=", a_shape, ", b=", b_shape, "). Utilisation de 0.8 par défaut.")
    alpha_draw <- 0.8
  }
  
  k_draw <- suppressWarnings(rlnorm(1, meanlog = logk_pr, sdlog = sdlog_k))
  if (!is.finite(k_draw)) {
    warning("rlnorm invalide (meanlog=", logk_pr, ", sdlog=", sdlog_k, "). Utilisation de exp(0)=1 par défaut.")
    k_draw <- 1
  }
  
  myinits <- list(
    N3_tot = rlnorm(n = n_N3_tot, meanlog = big_mu, sdlog = sd_N),
    N2     = rlnorm(n = n_N2,     meanlog = big_mu, sdlog = sd_N),
    N6     = rlnorm(n = n_N6,     meanlog = big_mu, sdlog = sd_N),
    N9     = rlnorm(n = n_N9,     meanlog = big_mu, sdlog = sd_N),
    
    p_smolt = cbind(fill_psmolt, 1 - fill_psmolt),  # si le modèle attend 2 colonnes
    
    logit_theta3  = runif(n_theta3,  -0.5, 0.5),
    logit_theta4f = runif(n_theta4f, -0.5, 0.5),
    logit_theta4m = runif(n_theta4m, -0.5, 0.5),
    logit_theta4  = runif(n_theta4,  -0.5, 0.5),
    logit_theta5  = rep(logit(0.8), n_theta5),
    logit_theta8  = rep(logit(0.8), n_theta8),
    
    prop3f = runif(n_prop3f, 0.4, 0.6),
    prop6f = runif(n_prop6f, 0.4, 0.6),
    prop9f = runif(n_prop9f, 0.4, 0.6),
    
    logN2_sd = .pos_num(runif(1, 0, 5), default = 1),
    
    alpha = alpha_draw,
    k     = k_draw
  )
  
  # Diagnostics rapides
  message(sprintf("Inits: theta1_max=%.3f, nsample=%g, a=%.4f, b=%.4f; logk_pr=%.3f, sdlog_k=%.3f",
                  theta1_max, nsample, a_shape, b_shape, logk_pr, sdlog_k))
  
  inits <- myinits
  inits
}

# Si on source le fichier et qu'on a des constantes dans l'env, on produit `inits`
try({
  if (exists("constants", inherits = TRUE) || exists("Const_nimble", inherits = TRUE)) {
    inits <- make_inits()
  }
}, silent = TRUE)
