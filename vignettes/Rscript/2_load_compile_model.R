nimbleOptions(unsupportedDerivativeHandling = "warn")
source(modelfile)
# library(cli) # for communication
cli_h2("Start building and compiling model {model_name} for chain {thischain}")
logfile <- glue("log/build_model_{model_name}_chain{thischain}.log")
logfile_connection <- file(logfile, open = "wt")

cli_alert_info("Log for chain {thischain} registered in {logfile}")
# File to save global times
time_log_file <- "log/times.log"
# Utility function to log timings
log_time <- function(step_name, duration, file = time_log_file) {
  log_entry <- paste(Sys.time(), "-", step_name, "- Duration (seconds):", round(duration, 2), "\n")
  write(log_entry, file = file, append = TRUE)
}

# ------------------------------------------------------------------------------
# Step 1 - Build Nimble model --------------------------------------------------
cli_progress_step("Build model {model_name} for chain {thischain}")

sink(file = logfile_connection)
sink(file = logfile_connection, type = "message", append = TRUE)
model.nimble <- nimbleModel(
  code = model.nimble,
  name = 'model.nimble',
  constants = Const_nimble,
  data = Data_nimble,
  inits = myinits, buildDerivs = TRUE
)
#warnings()
model.nimble$initializeInfo()
sink(type = "message")
sink()

# ------------------------------------------------------------------------------
# Step 2 - Compile the model ---------------------------------------------------
cli_progress_step("Compile model {model_name} for chain {thischain}")

timeStart_compile <- Sys.time()  # Start timer for compilation
sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection, type = "message", append = TRUE)
compiled.model <- compileNimble(model.nimble, showCompilerOutput = TRUE)
assign("compiled.model", compiled.model, envir = .GlobalEnv)
sink(type = "message")
sink()
timeEnd_compile <- Sys.time()  # End timer for compilation
time_compile <- difftime(timeEnd_compile, timeStart_compile, units = "sec")
cli_alert_success("Model compilation completed in {round(time_compile, 2)} seconds.")
# Save compilation time
log_time("Model Compilation", time_compile)
# ------------------------------------------------------------------------------
# Step 3 - Create a MCMC sampler for the model ---------------------------------
cli_progress_step("Configure MCMC sampler for model {model_name} for chain {thischain}")

sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection, type = "message", append = TRUE)
conf.mcmc.model <- nimbleHMC::configureHMC(model.nimble,
                        monitors = monitors,
                        enableWAIC = exists("enableWAIC"))
#conf.mcmc.model <- configureMCMC(
  #model.nimble,
  #monitors = monitors,
  #enableWAIC = TRUE
#)
# Apply samplers for specific nodes
# 1. Appliquer RW_block pour N2
#conf.mcmc.model$removeSamplers(c('N2', 'p_smolt_cohort'))
#conf.mcmc.model$addSampler(target = c('N2', 'p_smolt_cohort'), type = 'RW_block')#'AF_slice'

# 2. Appliquer NUTS pour logit_theta4f
#conf.mcmc.model$removeSamplers(c('logit_theta4f', 'sd_theta4f'))
#conf.mcmc.model$addSampler(target = c('logit_theta4f', 'sd_theta4f'), type = 'NUTS')
# 3. Appliquer NUTS pour logit_theta4f
# Supprimer les samplers existants pour 'alpha' et 'k'
#conf.mcmc.model$removeSamplers(c('alpha', 'k'))

# Ajouter un sampler pour 'alpha' et pour 'k'
#conf.mcmc.model$addSampler(target = c('alpha', 'k'), type = 'NUTS')

sink(type = "message")
sink()

# Build MCMC sampler with the above configuration
cli_progress_step("Build MCMC sampler for model {model_name} for chain {thischain}")

sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection, type = "message", append = TRUE)
MCMC.model <- buildMCMC(conf.mcmc.model)
sink(type = "message")
sink()

# ------------------------------------------------------------------------------
# Step 4 - Compile the MCMC sampler --------------------------------------------
cli_progress_step("Compile full model {model_name} for chain {thischain}")

timeStart_mcmc_compile <- Sys.time()  # Start timer for MCMC compilation
sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection, type = "message", append = TRUE)
compiled.MCMC.model <- compileNimble(MCMC.model, project = model.nimble, showCompilerOutput = TRUE)
sink(type = "message")
sink()
timeEnd_mcmc_compile <- Sys.time()  # End timer for MCMC compilation
time_mcmc_compile <- difftime(timeEnd_mcmc_compile, timeStart_mcmc_compile, units = "sec")
cli_alert_success("MCMC sampler compilation completed in {round(time_mcmc_compile, 2)} seconds.")
# Save MCMC compilation time
log_time("MCMC Sampler Compilation", time_mcmc_compile)
# ------------------------------------------------------------------------------
# Step 5 - Run the MCMC and measure execution time -----------------------------
cli_progress_step("Run MCMC for model {model_name} for chain {thischain}")

timeStart_mcmc <- Sys.time()  # Start timer for MCMC execution
compiled.MCMC.model$run(
  niter = n.iter,
  nburnin = n.burnin,
  thin = n.thin,
  time = TRUE,
  progressBar = TRUE
)
timeEnd_mcmc <- Sys.time()  # End timer for MCMC execution
time_mcmc <- difftime(timeEnd_mcmc, timeStart_mcmc, units = "sec")
cli_alert_success("MCMC execution completed in {round(time_mcmc, 2)} seconds.")

# Save MCMC execution time
log_time("MCMC Execution", time_mcmc)
# ------------------------------------------------------------------------------
# Step 6 - Get sampler times ---------------------------------------------------
cli_progress_step("Get sampler times for model {model_name} for chain {thischain}")

sampler_times <- compiled.MCMC.model$getTimes()
print(sampler_times)

# Save sampler times
sampler_times_file <- glue("outputs_hindcast/sampler_times_{model_name}_chain{thischain}.rds")
saveRDS(sampler_times, file = sampler_times_file)
cli_alert_success("Sampler times saved to {sampler_times_file}.")

cli_progress_done()




