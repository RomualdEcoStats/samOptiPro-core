pkgname <- "samOptiPro"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "samOptiPro-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('samOptiPro')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("assess_performance")
### * assess_performance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: assess_performance
### Title: Assess MCMC performance metrics (ESS, R-hat, AE, CE)
### Aliases: assess_performance

### ** Examples

## Not run: 
##D res <- run_baseline_config(build_M, niter = 2000, nburnin = 500, thin = 2)
##D perf <- assess_performance(res$samples, runtime_s = res$runtime_s)
##D perf$summary
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("assess_performance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_diag_from_mcmc_vect")
### * compute_diag_from_mcmc_vect

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_diag_from_mcmc_vect
### Title: Scalable diagnostics from MCMC samples (block-wise)
### Aliases: compute_diag_from_mcmc_vect

### ** Examples

## Not run: 
##D res_diag <- compute_diag_from_mcmc_vect(
##D   samples   = my_mcmc,
##D   runtime_s = 5400,
##D   compute_rhat         = "both",
##D   ess_for              = "both",
##D   target_block_ram_gb  = 2
##D )
##D head(res_diag)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_diag_from_mcmc_vect", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("configure_hmc_safely")
### * configure_hmc_safely

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: configure_hmc_safely
### Title: Configure and run HMC/NUTS safely
### Aliases: configure_hmc_safely
### Keywords: internal

### ** Examples

## Not run: 
##D   out <- configure_hmc_safely(
##D     build_fn = my_build_fn, niter = 2000, nburnin = 1000, thin = 1,
##D     monitors = c("beta[1]", "sigma"), nchains = 4, out_dir = "outputs"
##D   )
##D   out$diag_tbl
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("configure_hmc_safely", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("diagnose_model_structure")
### * diagnose_model_structure

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: diagnose_model_structure
### Title: Diagnose model structure, dependencies, and sampler time
###   (parameter and family levels)
### Aliases: diagnose_model_structure

### ** Examples

## Not run: 
##D res <- diagnose_model_structure(
##D   model            = my_nimble_model,
##D   make_plots       = TRUE,
##D   output_dir       = "outputs/diagnostics",
##D   save_csv         = TRUE,
##D   by_family        = TRUE,
##D   family_stat      = "median",
##D   time_normalize   = "per_node",
##D   only_family_plots = FALSE
##D )
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("diagnose_model_structure", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_bottlenecks")
### * plot_bottlenecks

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_bottlenecks
### Title: Plot MCMC Bottlenecks by Node or Family
### Aliases: plot_bottlenecks

### ** Examples

## Not run: 
##D # Example assuming an existing NIMBLE configuration and MCMC results:
##D res <- plot_bottlenecks(
##D   diag_tbl     = diag_tbl,
##D   conf.mcmc    = conf.mcmc,
##D   samples_ml   = samples_ml,
##D   sampled_only = TRUE,
##D   out_dir      = "outputs/diagnostics"
##D )
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_bottlenecks", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_bottlenecks_fast")
### * plot_bottlenecks_fast

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_bottlenecks_fast
### Title: Fast bottleneck plots (samplers-only) for large models
### Aliases: plot_bottlenecks_fast

### ** Examples

## Not run: 
##D conf <- nimble::configureMCMC(m)
##D res <- plot_bottlenecks_fast(
##D   diag_tbl   = diag_tbl,
##D   sampled_only = TRUE,
##D   conf.mcmc  = conf,
##D   samples_ml = samples_mla,
##D   out_dir    = "outputs/diagnostics_fast"
##D )
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_bottlenecks_fast", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_convergence_rhat_sampled_fast")
### * plot_convergence_rhat_sampled_fast

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_convergence_rhat_sampled_fast
### Title: Fast convergence plots for R-hat (samplers-only), thresholded
### Aliases: plot_convergence_rhat_sampled_fast

### ** Examples

## Not run: 
##D conf <- nimble::configureMCMC(m)
##D rhat_res <- plot_convergence_rhat_sampled_fast(
##D   diag_tbl     = diag_tbl,
##D   threshold    = 1.01,
##D   sampled_only = TRUE,
##D   conf.mcmc    = conf,
##D   samples_ml   = samples_mla,
##D   out_dir      = "outputs/diagnostics_fast_rhat",
##D   top_k        = 20L,
##D   prefer_split = TRUE
##D )
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_convergence_rhat_sampled_fast", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_strategies_from_test_result_fast")
### * plot_strategies_from_test_result_fast

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_strategies_from_test_result_fast
### Title: Plot strategy comparisons from test_strategy_family_fast results
###   (fast path)
### Aliases: plot_strategies_from_test_result_fast

### ** Examples

## Not run: 
##D res <- test_strategy_family_fast(build_fn = build_M, nbot = 2, try_hmc = TRUE)
##D plot_strategies_from_test_result_fast(res, per = "family", top_k = 30, top_by = "CE")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_strategies_from_test_result_fast", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("profile_sampler_times")
### * profile_sampler_times

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: profile_sampler_times
### Title: Profile sampler times per MCMC sampler
### Aliases: profile_sampler_times

### ** Examples

## Not run: 
##D   model  <- nimbleModel(code, constants = Const, data = Data, inits = Inits)
##D   cmodel <- compileNimble(model)
##D   conf   <- configureMCMC(model)
##D   ts <- profile_sampler_times(conf, cmodel, niter = 1e4)
##D   print(ts)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("profile_sampler_times", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("run_structure_and_hmc_test")
### * run_structure_and_hmc_test

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: run_structure_and_hmc_test
### Title: Run structural diagnostics and (optional) HMC/NUTS smoke test
### Aliases: run_structure_and_hmc_test

### ** Examples

## Not run: 
##D out <- run_structure_and_hmc_test(my_builder, include_data = FALSE, try_hmc = TRUE)
##D if (!out$hmc$ok) message("HMC not feasible: ", out$hmc$error)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("run_structure_and_hmc_test", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("test_strategy")
### * test_strategy

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: test_strategy
### Title: Test and compare MCMC strategies on selected bottleneck nodes
### Aliases: test_strategy
### Keywords: internal

### ** Examples

## Not run: 
##D res <- test_strategy(
##D   build_fn = my_build_fn,
##D   monitors = c("theta","beta"),
##D   try_hmc  = TRUE,
##D   nbot     = 2,
##D   out_dir  = "outputs/diagnostics"
##D )
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("test_strategy", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
