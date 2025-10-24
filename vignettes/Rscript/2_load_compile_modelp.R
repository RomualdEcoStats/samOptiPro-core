# Define helper functions to manage sink
start_sink <- function(connection) {
  sink(file = connection, append = TRUE)
  sink(file = connection, type = "message", append = TRUE)
}

end_sink <- function() {
  sink(type = "message")
  sink()
}

# Main script
source(modelfile)
# library(cli) # for communication
cli_h2("Start building and compiling model {model_name} for chain {thischain}")
logfile <- glue("log/build_model_{model_name}_chain{thischain}.log")
logfile_connection <- file(logfile, open = "wt")

cli_alert_info("Log for chain {thischain} registered in {logfile}")

# ------------------------------------------------------------------------------
# Step 1 - Build Nimble model --------------------------------------------------
cli_progress_step("Build model {model_name} for chain {thischain}")

start_sink(logfile_connection)
model.nimble <- nimbleModel(
  code = model.nimble,
  name = 'model.nimble',
  constants = Const_nimble,
  data = Data_nimble,
  inits = myinits
)
model.nimble$initializeInfo()
end_sink()

# ------------------------------------------------------------------------------
# Step 2 - Compile the model ---------------------------------------------------
cli_progress_step("Compile model {model_name} for chain {thischain}")

timeStart_compile <- Sys.time()  # Start timer for compilation

start_sink(logfile_connection)
compiled.model <- compileNimble(model.nimble, showCompilerOutput = TRUE)
assign("compiled.model", compiled.model, envir = .GlobalEnv)
end_sink()

timeEnd_compile <- Sys.time()  # End timer for compilation
time_compile <- difftime(timeEnd_compile, timeStart_compile, units = "sec")
cli_alert_success("Model compilation completed in {round(time_compile, 2)} seconds.")

# ------------------------------------------------------------------------------
# Step 3 - Create a MCMC sampler for the model ---------------------------------
cli_progress_step("Configure MCMC sampler for model {model_name} for chain {thischain}")

start_sink(logfile_connection)
conf.mcmc.model <- configureMCMC(
  model.nimble,
  monitors = monitors,
  enableWAIC = TRUE
)
end_sink()

cli_progress_step("Build MCMC sampler for model {model_name} for chain {thischain}")

start_sink(logfile_connection)
MCMC.model <- buildMCMC(conf.mcmc.model)
end_sink()

# ------------------------------------------------------------------------------
# Step 4 - Compile the MCMC sampler --------------------------------------------
cli_progress_step("Compile full model {model_name} for chain {thischain}")

timeStart_mcmc_compile <- Sys.time()  # Start timer for MCMC compilation

start_sink(logfile_connection)
compiled.MCMC.model <- compileNimble(MCMC.model, project = model.nimble, showCompilerOutput = TRUE)
end_sink()

timeEnd_mcmc_compile <- Sys.time()  # End timer for MCMC compilation
time_mcmc_compile <- difftime(timeEnd_mcmc_compile, timeStart_mcmc_compile, units = "sec")
cli_alert_success("MCMC sampler compilation completed in {round(time_mcmc_compile, 2)} seconds.")

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
