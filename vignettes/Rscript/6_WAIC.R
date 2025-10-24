# R. Patin - 23 Nov. 2021
cli_h1("Calculate wAIC with posteriors for {n.chains} chains")

cli_progress_step("Load MCMC results")
if(!exists("mcmc_nimble")){
  mcmc_nimble <- readRDS(name_file_mcmc)
  time_mcmc <- readRDS(name_file_time)
}

# Load configuration ------------------------------------------------------
cli_progress_step("Load and compile model")
thischain <- 1
set.seed(thischain)
source(initfile)
source("Rscript/2_load_compile_model.R")

cli_progress_step("Initialize model")

compiled.MCMC.model$run(niter = 100,
                         nburnin = 0,
                         thin = 100,
                         progressBar = TRUE,
                         reset = TRUE)

lapply(mcmc_nimble,function(x){
  x[-seq(1,burnin_output),]
}) %>% 
  do.call('rbind',.) -> new_mvSamples

cli_progress_step("Load samples into the model")
nimble:::matrix2mv(new_mvSamples, compiled.MCMC.model$mvSamples)
cli_progress_step("Calculate wAIC")

sink(file = logfile_connection)
sink(file = logfile_connection, type = "message", append = TRUE)
this_wAIC <- calculateWAIC(compiled.MCMC.model)
sink(type = "message")
sink()


cli_process_done()
cli_alert_success("wAIC = {round(this_wAIC$WAIC, digits = 2)}")
saveRDS(this_wAIC, file =
          glue("outputs_hindcast/wAIC_{model_name}.rds"))
