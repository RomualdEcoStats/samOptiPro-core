## ------------------------------------------------------------------
#devtools::load_all("samOptiPro")
devtools::load_all("~/Scorff LCM_model1/samOptiPro")
devtools::document()
devtools::load_all()

####################################################################################"
## ======================================================================
## 0) Préambule
## ======================================================================
stopifnot(requireNamespace("nimble", quietly = TRUE))
stopifnot(requireNamespace("samOptiPro", quietly = TRUE))

library(nimble)
library(nimbleHMC)
library(samOptiPro)
ilogit <- plogis; logit <- qlogis
`%||%` <- function(x, y) if (is.null(x)) y else x

set.seed(123)

out_dir <- "outputs/tutorial_nondiff"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ======================================================================
## 1) Données synthétiques
##    - z[t] processus RW (aléatoire)
##    - p[t,·] probas de K catégories ~ Dirichlet(α)
##    - counts[t,·] ~ Multinomial(N[t], p[t,·])
##    - y[t] observation gaussienne autour de z[t]
## ======================================================================
n  <- 80         # longueur de série
K  <- 3          # catégories
N_counts <- 50   # taille d’échantillon multinomial par t
sd_proc_true <- 0.4
sd_obs_true  <- 0.3

z_true <- numeric(n)
z_true[1] <- rnorm(1, 0, sd_proc_true)
for (t in 2:n) z_true[t] <- rnorm(1, z_true[t-1], sd_proc_true)
y <- rnorm(n, z_true, sd_obs_true)

# Probabilités (douces) puis arrondies pour injecter de la non-diff plus loin
base_logits <- matrix(rnorm(n*K, sd = 0.7), n, K)
p_true <- t(apply(base_logits, 1, function(v) { u <- exp(v - max(v)); u/sum(u) }))

counts <- t(apply(p_true, 1, function(p) rmultinom(1, size = N_counts, prob = p)))

# NIMBLE attend un tableau |counts| (n x K), probas |p| (n x K)
Data_nimble <- list(
  y          = y,
  counts     = counts
)

Constants <- list(
  n              = n,
  K              = K,
  N_counts       = rep(N_counts, n),
  sd_obs_const   = sd_obs_true,
  dirich_alpha   = rep(2, K)    # concentration modérée
)

# Inits grossiers
myinits <- list(
  sd_proc = runif(1, 0.1, 1),
  z       = c(rnorm(1, 0, 1), diff(y)),
  # init p en Dirichlet (toutes lignes identiques, NIMBLE adaptera)
  p       = matrix(1/K, n, K)
)

## ======================================================================
## 2) Modèle nimble (avec non-diff)
##    - Dirichlet + Multinomial
##    - round() introduit une non-linéarité non différentiable
##    - logLik[t] pour WAIC
## ======================================================================
model.nondiff <- nimble::nimbleCode({
  sd_proc ~ dunif(0, 5)
  
  z[1] ~ dnorm(0, sd = sd_proc)
  for (t in 2:n) {
    z[t] ~ dnorm(z[t-1], sd = sd_proc)
  }
  
  # Observation continue
  for (t in 1:n) {
    y[t] ~ dnorm(z[t], sd = sd_obs_const)
  }
  
  # Couche discrète : Dirichlet -> Multinomial (non-diff pour HMC)
  for (t in 1:n) {
    p[t, 1:K] ~ ddirich(dirich_alpha[1:K])
    counts[t, 1:K] ~ dmulti(p[t, 1:K], size = N_counts[t])
    
    # Non-diff volontaire via round ; utilisé dans une eq. déterministe
    zr[t] <- round(z[t])             # rupture de différentiabilité
    # Petit “indice” corrélé à z : sert juste à rajouter une dépendance
    idx[t] <- zr[t] + 10
    # log-vraisemblance agrégée (avertissement NIMBLE attendu)
    logLik[t] <- dnorm(y[t], mean = z[t], sd = sd_obs_const, log = 1) +
      dmulti(counts[t,1:K], prob = p[t,1:K], size = N_counts[t], log = 1)
  }
})

## ======================================================================
## 3) Construction & compilation
## ======================================================================
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = FALSE)

m  <- nimbleModel(code = model.nondiff, name = "nondiff",
                  constants = Constants, data = Data_nimble, inits = myinits)
cm <- compileNimble(m)

## Fabrique build_fn attendue par samOptiPro
monitors = c("sd_proc","z","p","logLik")
build_nondiff <- function() list(
  model    = m,
  cmodel   = cm,
  monitors = monitors,
  code_text = paste(deparse(model.nondiff), collapse = "\n")
)

## ======================================================================
## 4) Diagnostic structure & test HMC
## ======================================================================
cat("\n[MODEL STRUCTURE CHECK]\n")
diag_s <- diagnose_model_structure(model =m,
                                   include_data        = FALSE,
                                   removed_nodes       = NULL,
                                   ignore_patterns     = c("^lifted_", "^logProb_"),
                                   make_plots          = TRUE,
                                   output_dir          ="outputs/tutorial_nondiff",
                                   save_csv            = TRUE,
                                   node_of_interest    = NULL,
                                   sampler_times       = NULL,
                                   sampler_times_unit  = "seconds",
                                   auto_profile        = TRUE,
                                   profile_niter       = 2000L,
                                   profile_burnin      = 500L,
                                   profile_thin        = 1L,
                                   profile_seed        = NULL,
                                   np                  = 0.10,
                                   by_family           = TRUE,
                                   family_stat         = c("median","mean","sum"),
                                   time_normalize      = c("none","per_node"),
                                   only_family_plots   = TRUE)
diag_s$plots 
cat(sprintf("- Stochastic nodes   : %d\n", length(diag_s$stochastic_nodes)))
cat(sprintf("- Deterministic nodes: %d\n", length(diag_s$deterministic_nodes)))
out <- run_structure_and_hmc_test(build_nondiff, include_data = FALSE)
#Baseline MCMC, Bottelenecks and Performance Assessment
n.iter   <- 4000
n.burnin <- 1000
n.thin   <- 2
n.chains <- 3

res_b <- run_baseline_config(
  build_nondiff,
  niter   = n.iter,
  nburnin = n.burnin,
  nchains = n.chains,
  thin    = n.thin,
  monitors = monitors
)

samples_ml <- as_mcmc_list_sop(res_b$samples, res_b$samples2,
                               drop_loglik = FALSE, thin = n.thin)

runtime_s <- res_b$runtime_s   # Time total
ap  <- assess_performance(samples_ml, runtime_s)

bot <- identify_bottlenecks_family(samples_ml,
                                   runtime_s,
                                   ess_threshold = 1000,
                                   time_threshold = "auto",
                                   rhat_threshold = 1.01,
                                   ess_per_s_min = 0)
runtime_s
ap$summary
bot$top3

# Plots écrits dans outputs/diagnostics/{baseline|strategy}
diff <- test_strategy_family(
  build_fn     = build_nondiff,
  try_hmc      = FALSE,
  nbot         = 2,
  pilot_niter  = 5000,
  pilot_burnin = 1500,
  thin         = 2,
  ask          = TRUE,
  out_dir      = "outputs/tutorial_nondiff",
  order_scalar = c("NUTS","slice","RW"),
  order_block  = c("NUTS_block","AF_slice","RW_block")
)











## ======================================================================
## 5) Paramètres MCMC
## ======================================================================
n.iter   <- 6e4
n.burnin <- 1e4
n.thin   <- 20
ESS_goal <- 1000

## ======================================================================
## 6) Baseline
## ======================================================================
cat("\n[BASELINE] run_baseline_config…\n")
res_b <- run_baseline_config(build_nondiff, niter = n.iter, nburnin = n.burnin,nchains = 3, thin = n.thin,monitors = c("sd_proc", "z", "p", "logLik"))
samples   <- res_b$samples  # <-- fusionne samples + samples2

runtime_s <- res_b$runtime_s
# Fusion propre en mcmc.list multi-chaînes (et concat samples2)
samples_ml <- as_mcmc_list_sop(res_b$samples, res_b$samples2, drop_loglik = FALSE, thin = n.thin)

ap  <- assess_performance(samples_ml, runtime_s)
#bot <- identify_bottlenecks(samples_ml, runtime_s)


bot2 <- identify_bottlenecks(samples_ml, res_b$runtime_s, top_k = 5)
bot2$details$algo   # 15 pires en AE (faibles AE)
bot2$details$comp   # 15 pires en CE (coût élevé)
bot2$details$joint  # 15 pires "simultanément" (agrégation de rangs)

diag_tbl <- compute_diag_from_mcmc(samples_ml, runtime_s = res_b$runtime_s)
plots_bn <- plot_bottlenecks(diag_tbl, out_dir = "outputs/diagnostics")
plots_cv <- plot_convergence_checks(samples_ml, out_dir = "outputs/diagnostics",
                                    top_k_rhat = 12, top_k_aelow = 12,
                                    runtime_s = res_b$runtime_s)


perf_b <- tryCatch(safe_assess_performance(res_b$samples, res_b$runtime_s), error = function(e) NULL)

cat(sprintf("- runtime baseline: %.2f s\n", res_b$runtime_s))
if (!is.null(perf_b) && is.finite(perf_b$ESS_per_s)) {
  cat(sprintf("- ESS/s baseline : %.1f  (temps pour ESS=%d ≈ %.1f s)\n",
              perf_b$ESS_per_s, ESS_goal, ESS_goal / perf_b$ESS_per_s))
}
######Mode en parallele et plusieurs chaines pour le test###############################################
## ======================================================================
## 6bis) Baseline en PARALLÈLE — 3 chaînes
## ======================================================================
## ======================================================================
## 6bis) Baseline en PARALLÈLE — 3 chaînes (robuste OS-adaptive)
## ======================================================================
stopifnot(requireNamespace("parallel", quietly = TRUE))
library(parallel)

n_chains <- 3
monitors <- c("sd_proc", "z", "p", "logLik")

## exports en formats simples (typage strict)
y_export         <- as.numeric(Data_nimble$y)
counts_export    <- as.matrix(Data_nimble$counts)
storage.mode(counts_export) <- "integer"   # multinomial attend des entiers
n_export         <- as.integer(Constants$n)
K_export         <- as.integer(Constants$K)
N_counts_export  <- as.integer(Constants$N_counts)
sd_obs_export    <- as.numeric(Constants$sd_obs_const)
dirich_export    <- as.numeric(Constants$dirich_alpha)

## worker (commun aux deux modes); zéro '$' dans le code exécuté
## -- Worker v4 : répertoire de compilation unique + logs détaillés
worker_core <- function(seed_offset) {
  step <- "init"
  tryCatch({
    set.seed(123 + seed_offset)
    
    y_loc        <- y_export
    counts_loc   <- counts_export
    n_loc        <- n_export
    K_loc        <- K_export
    N_counts_loc <- N_counts_export
    sd_obs_loc   <- sd_obs_export
    dirich_loc   <- dirich_export
    
    if (!is.integer(counts_loc)) storage.mode(counts_loc) <- "integer"
    stopifnot(is.numeric(y_loc), is.matrix(counts_loc),
              length(N_counts_loc) == n_loc,
              length(dirich_loc) == K_loc)
    
    # --- code modèle ---
    model_code_w <- nimble::nimbleCode({
      sd_proc ~ dunif(0, 5)
      
      z[1] ~ dnorm(0, sd = sd_proc)
      for (t in 2:n) {
        z[t] ~ dnorm(z[t-1], sd = sd_proc)
      }
      
      for (t in 1:n) {
        y[t] ~ dnorm(z[t], sd = sd_obs_const)
      }
      
      for (t in 1:n) {
        p[t, 1:K] ~ ddirich(dirich_alpha[1:K])
        counts[t, 1:K] ~ dmulti(p[t, 1:K], size = N_counts[t])
        
        zr[t] <- round(z[t])
        idx[t] <- zr[t] + 10
        logLik[t] <- dnorm(y[t], mean = z[t], sd = sd_obs_const, log = 1) +
          dmulti(counts[t,1:K], prob = p[t,1:K], size = N_counts[t], log = 1)
      }
    })
    
    Data_w <- list(y = y_loc, counts = counts_loc)
    Constants_w <- list(
      n = n_loc, K = K_loc, N_counts = N_counts_loc,
      sd_obs_const = sd_obs_loc, dirich_alpha = dirich_loc
    )
    
    z_init <- c(stats::rnorm(1, 0, 1), diff(y_loc))
    p_init <- matrix(1 / K_loc, n_loc, K_loc)
    myinits_w <- list(sd_proc = stats::runif(1, 0.1, 1), z = z_init, p = p_init)
    
    # === Répertoire de compilation unique par worker ===
    comp_dir <- file.path(tempdir(), sprintf("nimbleW_%s_%d", Sys.getpid(), seed_offset))
    if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
    
    # --- build/compile modèle ---
    step <- "nimbleModel"
    m_w  <- nimbleModel(model_code_w, constants = Constants_w, data = Data_w, inits = myinits_w)
    
    step <- "compileNimble(model)"
    cm_w <- compileNimble(m_w, dirName = comp_dir, showCompilerOutput = TRUE)
    
    # --- configure/build/compile MCMC ---
    step <- "configureMCMC"
    conf_w <- configureMCMC(m_w)
    if (length(monitors)) conf_w$addMonitors(monitors)
    
    step <- "buildMCMC"
    mcmc_w  <- buildMCMC(conf_w)
    
    step <- "compileNimble(mcmc)"
    cmcmc_w <- compileNimble(mcmc_w, project = cm_w, dirName = comp_dir, resetFunctions = TRUE, showCompilerOutput = TRUE)
    
    if (!inherits(cmcmc_w, "CompiledObject")) {
      # Essaie de sortir des logs utiles
      log_files <- list.files(comp_dir, pattern = "\\.(out|log|err)$", full.names = TRUE)
      msg <- "Compilation MCMC échouée sans message C++ capturé."
      if (length(log_files)) {
        tail_txt <- try(paste(utils::tail(readLines(log_files[1]), 50), collapse = "\n"), silent = TRUE)
        if (!inherits(tail_txt, "try-error")) msg <- paste(msg, "\n--- Tail log ---\n", tail_txt)
      }
      stop(msg)
    }
    
    # --- run ---
    step <- "runMCMC"
    res_w <- nimble::runMCMC(cmcmc_w,
                             niter = as.integer(n.iter),
                             nburnin = as.integer(n.burnin),
                             thin   = as.integer(n.thin),
                             nchains = 1L,
                             samplesAsCodaMCMC = TRUE,
                             summary = FALSE,
                             WAIC = FALSE)
    
    list(error = FALSE, samples = res_w$samples, runtime_s = res_w$runtime_s)
  }, error = function(e) list(error = TRUE, stage = step, message = conditionMessage(e)))
}

## --- lance (comme avant) ---
if (.Platform$OS.type != "windows") {
  res_list <- parallel::mclapply(1:n_chains, worker_core, mc.cores = n_chains)
} else {
  cl <- parallel::makeCluster(min(n_chains, parallel::detectCores()))
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  parallel::clusterEvalQ(cl, {
    library(nimble)
    nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
    nimbleOptions(MCMCsaveHistory = FALSE)
    RNGkind("L'Ecuyer-CMRG")
    NULL
  })
  parallel::clusterExport(cl,
                          varlist = c("y_export","counts_export","n_export","K_export","N_counts_export",
                                      "sd_obs_export","dirich_export","n.iter","n.burnin","n.thin","monitors",
                                      "worker_core"),
                          envir = environment()
  )
  res_list <- parallel::parLapply(cl, 1:n_chains, worker_core)
}

## ======================================================================
## 7) Évaluation des perfs + export (commune séquentiel / parallèle)
## ======================================================================

## Choix des objets d'évaluation :
## - Si tu as exécuté la 6bis (parallèle), 'samples' et 'runtime_s' ont été
##   redéfinis en fin de 6bis pour pointer vers samples_par / runtime_s_par.
## - Sinon, on garde la baseline séquentielle (res_b).
if (!exists("samples") || !exists("runtime_s")) {
  samples   <- res_b$samples
  runtime_s <- res_b$runtime_s
}

cat("\n[ASSESS] Calcul ESS, Rhat, AE, CE…\n")
ap <- assess_performance(samples, runtime_s)

cat("\n[ASSESS] Résumé global\n")
print(ap$summary)

cat("\n[ASSESS] 10 premiers paramètres (tri par ESS croissant)\n")
print(utils::head(ap$per_param[order(ap$per_param$ESS), ], 10))

## Export CSV (résumé + par paramètre)
summary_csv   <- file.path(out_dir, "performance_summary.csv")
perparam_csv  <- file.path(out_dir, "performance_per_param.csv")
try(utils::write.csv(ap$summary,    summary_csv,  row.names = FALSE), silent = TRUE)
try(utils::write.csv(ap$per_param,  perparam_csv, row.names = FALSE), silent = TRUE)
cat(sprintf("\n[EXPORT] %s\n[EXPORT] %s\n", summary_csv, perparam_csv))

## ======================================================================
## 8) Bottlenecks (convergence/efficacité)
## ======================================================================
cat("\n[BOTTLENECKS] Recherche de goulets…\n")
bot <- identify_bottlenecks(
  samples, runtime_s,
  ess_threshold  = 100,   # ESS mini par paramètre
  time_threshold = 60,    # runtime trop long (s)
  rhat_threshold = 1.01,  # seuil Gelman-Rubin
  ess_per_s_min  = 1.0    # CE mini (ESS/s)
)

cat(sprintf("- Type: %s\n", bot$type))
if (length(bot$details)) {
  cat("[BOTTLENECKS] Détails:\n")
  print(bot$details)
}
if (!is.null(bot$summary)) {
  cat("\n[BOTTLENECKS] Rappel résumé global:\n")
  print(bot$summary)
}

## ======================================================================
## 9) Plots (ESS, CE, AE, Rhat) + export PNG
## ======================================================================
have_gg <- requireNamespace("ggplot2", quietly = TRUE)
if (isTRUE(have_gg)) {
  cat("\n[PLOTS] Génération des histogrammes…\n")
  plots <- plot_mcmc_histograms(samples, runtime_s, rhat_thresh = 1.01, bins = 35)
  
  ## Sauvegarde PNG (si un plot est NULL, on saute)
  save_plot <- function(p, name) {
    if (!is.null(p)) {
      fp <- file.path(out_dir, paste0(name, ".png"))
      ggplot2::ggsave(filename = fp, plot = p, width = 7, height = 5, dpi = 150)
      cat(sprintf("[PLOTS] Exporté: %s\n", fp))
    }
  }
  save_plot(plots$ESS,  "hist_ESS")
  save_plot(plots$CE,   "hist_CE")
  save_plot(plots$AE,   "hist_AE")
  save_plot(plots$Rhat, "hist_Rhat")
  
  ## Affichage interactif (décommente si tu veux voir à l’écran)
  # print(plots$ESS); print(plots$CE); print(plots$AE)
  # if (!is.null(plots$Rhat)) print(plots$Rhat)
} else {
  cat("\n[PLOTS] ggplot2 non disponible — install.packages('ggplot2') pour les figures.\n")
}

## ======================================================================
## 10) Indicateurs rapides style “baseline”
## ======================================================================
cat("\n[KPIs] Indicateurs rapides\n")
ESS_goal <- ESS_goal %||% 1000
ESS_per_s <- ap$summary$ESS_per_s
if (is.finite(ESS_per_s)) {
  t_for_goal <- ESS_goal / ESS_per_s
  cat(sprintf("- runtime total : %.2f s\n", runtime_s))
  cat(sprintf("- ESS/s global  : %.2f\n", ESS_per_s))
  cat(sprintf("- Temps estimé pour ESS=%d : %.1f s (%.1f min)\n",
              ESS_goal, t_for_goal, t_for_goal/60))
} else {
  cat("- ESS/s non calculable (vérifier les samples / runtime)\n")
}


# 7) Goulets
bot_par <- identify_bottlenecks(
  samples_par, runtime_s_par,
  ess_threshold = 100,
  time_threshold = 60,
  rhat_threshold = 1.01,
  ess_per_s_min = 1.0
)
cat("\n[PARALLEL x3] Bottleneck type:", bot_par$type, "\n")
if (length(bot_par$details)) print(bot_par$details)

# 8) Histogrammes (si tu veux les afficher)
# plots_par <- plot_mcmc_histograms(samples_par, runtime_s_par, rhat_thresh = 1.01, bins = 35)
# print(plots_par$ESS); print(plots_par$CE); print(plots_par$AE); if (!is.null(plots_par$Rhat)) print(plots_par$Rhat)

# 9) (Option) remplacer les objets de la baseline par la version parallèle pour la suite du script
samples   <- samples_par
runtime_s <- runtime_s_par

## ======================================================================
## 7) Diagnostics par cible (temps & ESS/s), plots, sélection goulots
## ======================================================================
cat("\n[DIAGNOSTICS cibles] diagnostics_by_target…\n")
diag_tbl <- tryCatch(
  diagnostics_by_target(build_nondiff, niter = min(2e4, n.iter), samples = res_b$samples),
  error = function(e) NULL
)

# Plots robustes (recalcule ESS/s si non fourni)
diag_tbl2 <- plot_bottlenecks(
  diag_tbl,
  out_dir   = file.path(out_dir, "diagnostics"),
  samples   = res_b$samples,
  runtime_s = res_b$runtime_s
)

# Remplir les métriques si besoin
diag_tbl2 <- fill_bottleneck_metrics(diag_tbl, samples = res_b$samples, runtime_s = res_b$runtime_s)

# Option 1: pick_bottlenecks classique (si ta version ne gère que 2 args)
worst <- pick_bottlenecks(diag_tbl2, top_k = 4)

# Option 2: version robuste (accepte samples/runtime_s)
# worst <- pick_bottlenecks_robust(diag_tbl, top_k = 4, samples = res_b$samples, runtime_s = res_b$runtime_s)

cat("\nCibles retenues pour action:\n")
print(worst)

## ======================================================================
## 8) Plan d’action samplers + application + nouveau run
##    - 1ère moitié -> slice
##    - 2ème moitié -> RW, avec RW_block par paires
## ======================================================================
plan <- propose_sampler_actions(worst, block = TRUE)
cat("\nPlan de samplers proposé:\n")
print(plan)

patched <- apply_sampler_plan(build_nondiff, plan)
res_opt <- run_baseline_config(function() patched, niter = n.iter, nburnin = n.burnin, thin = n.thin)
perf_opt <- tryCatch(safe_assess_performance(res_opt$samples, res_opt$runtime_s), error = function(e) NULL)

cat("\n[APRES OPTIMISATION]\n")
cat(sprintf("- runtime: %.2f s\n", res_opt$runtime_s))
if (!is.null(perf_opt) && is.finite(perf_opt$ESS_per_s)) {
  cat(sprintf("- ESS/s  : %.1f  (temps pour ESS=%d ≈ %.1f s)\n",
              perf_opt$ESS_per_s, ESS_goal, ESS_goal / perf_opt$ESS_per_s))
}

## ======================================================================
## 9) Tentative HMC (attendue négative ici) — mais tolérante
## ======================================================================
cat("\n[HMC] run_hmc_all_nodes (tolérant)…\n")
res_h <- tryCatch(run_hmc_all_nodes(build_nondiff, niter = n.iter, nburnin = n.burnin, thin = n.thin),
                  error = function(e) NULL)
if (!is.null(res_h)) {
  perf_h <- tryCatch(safe_assess_performance(res_h$samples, res_h$runtime_s), error = function(e) NULL)
  cat(sprintf("- runtime HMC: %.2f s\n", res_h$runtime_s))
  if (!is.null(perf_h) && is.finite(perf_h$ESS_per_s)) {
    cat(sprintf("- ESS/s HMC : %.1f  (temps pour ESS=%d ≈ %.1f s)\n",
                perf_h$ESS_per_s, ESS_goal, ESS_goal / perf_h$ESS_per_s))
  }
} else {
  cat("- HMC indisponible (modèle non différentiable), ce qui est attendu ici.\n")
}

## ======================================================================
## 10) WAIC (grâce à logLik[t])
## ======================================================================
waic_b   <- tryCatch(compute_WAIC(res_b$samples), error = function(e) NULL)
waic_opt <- tryCatch(compute_WAIC(res_opt$samples), error = function(e) NULL)
if (!is.null(waic_b))   cat(sprintf("\nWAIC baseline     : %.2f\n", waic_b$summary$WAIC))
if (!is.null(waic_opt)) cat(sprintf("WAIC optimisé      : %.2f\n", waic_opt$summary$WAIC))

## ======================================================================
## 11) Extraction infos de samplers et temps (profil)
## ======================================================================
cat("\n[Profil samplers]\n")
# Table des samplers de la conf baseline
df_conf <- tryCatch(sampler_df_from_conf(patched$conf), error = function(e) NULL)
if (is.null(df_conf)) {
  # fallback: fabrique simple
  sfun <- patched$conf$getSamplers()
  nm <- names(sfun); if (is.null(nm)) nm <- paste0("sampler_", seq_along(sfun))
  df_conf <- data.frame(name = nm,
                        type = vapply(sfun, function(s) class(s)[1], character(1)),
                        target = vapply(sfun, function(s) paste0(s$target, collapse=","), character(1)),
                        stringsAsFactors = FALSE)
}

# Temps par sampler (alignés par index) via un mini run
m2  <- nimble::buildMCMC(patched$conf)
cm2 <- nimble::compileNimble(m2, project = patched$cmodel, resetFunctions = TRUE)
cm2$run(2e4)
times <- cm2$getTimes()
df_conf$time_s <- times[seq_len(nrow(df_conf))]

head(df_conf)

## ======================================================================
## 12) Exports
## ======================================================================
if (!is.null(perf_b))   readr::write_csv(perf_b,   file.path(out_dir, "perf_baseline.csv"))
if (!is.null(perf_opt)) readr::write_csv(perf_opt, file.path(out_dir, "perf_optim.csv"))
if (!is.null(diag_tbl)) readr::write_csv(diag_tbl, file.path(out_dir, "diagnostics_by_target_raw.csv"))
if (!is.null(diag_tbl2)) readr::write_csv(diag_tbl2, file.path(out_dir, "diagnostics_by_target_filled.csv"))
if (!is.null(waic_b))   readr::write_csv(waic_b$pointwise,   file.path(out_dir, "waic_pointwise_baseline.csv"))
if (!is.null(waic_opt)) readr::write_csv(waic_opt$pointwise, file.path(out_dir, "waic_pointwise_optim.csv"))
readr::write_csv(df_conf, file.path(out_dir, "samplers_profile.csv"))

cat(sprintf("\nFichiers écrits dans: %s\n", normalizePath(out_dir)))

