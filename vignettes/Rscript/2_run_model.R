set.seed(42)

# Build and compile model -------------------------------------------------

if(!parallel_run){
  cli_h1("Build, compile and run model {model_name} for {n.chains} chain sequentially.")
  
  # Single core run ---------------------------------------------------------
  timeStart <- Sys.time()
  thischain <- 1
  cli_progress_step("Run model for chain {thischain}")
  for(thischain in seq_len(n.chains)){
    cli_progress_update()
    set.seed(thischain)
    source(initfile)
    source("Rscript/2_load_compile_model.R")
    mcmc_nimble <- list()
    message("Start running chain ",thischain)
    compiled.MCMC.model$run(niter = n.iter,
                            nburnin = n.burnin,
                            thin = n.thin,
                            progressBar = TRUE, 
                            reset = TRUE)
    mcmc_nimble[[thischain]] <- as.matrix(compiled.MCMC.model$mvSamples)
  }
  timeEnd <- Sys.time()
  time_mcmc <- difftime(timeEnd, timeStart, units='sec')
  time_mcmc_min <- round(time_mcmc/60, digits = 0)
  cli_progress_done()
  cli_alert_success("Model run sequentially for {n.chains} chains in {time_mcmc_min}mins")
} else {
  # Parallel Run ------------------------------------------------------------
  cli_h1("Build, compile and run model {model_name} for {n.chains} chain in parallel.")
  
  # 1 - Create a cluster with n nodes (here n = n.chains)
  cl <- makeCluster(n.chains, outfile = "")
  
  # 2 - Clone everything that is needed in each node
  # including the libraries needed to run the model
  
  # Clone all other arguments needed to run the function
  clusterExport(cl, c('modelfile',
                      'Data_nimble',
                      'Const_nimble',
                      'n.iter',
                      'n.burnin',
                      'n.thin',
                      'monitors',
                      'model_name',
                      'initfile'))
  
  for(thischain in seq_len(n.chains)){
    clusterExport(cl[thischain], c('thischain'))
  }
  
  timeStart <- Sys.time()
  clusterEvalQ(cl, {
    
    # load MCMC control run  
    source("Rscript/0_load_library.R")
    source("Rscript/0_load_functions.R")
    
    set.seed(thischain)
    source(initfile)
    
    # compile model and MCMC configuration
    source("Rscript/2_load_compile_model.R")
    cli_h1("Run full model for chain {thischain}")
    
    compiled.MCMC.model$run(niter = n.iter,
                            nburnin = n.burnin,
                            thin = n.thin,
                            time = TRUE,
                            progressBar = TRUE, 
                            reset = TRUE)
  })
  timeEnd <- Sys.time()
  time_mcmc <- difftime(timeEnd, timeStart, units='sec')
  time_mcmc_min <- round(time_mcmc/60, digits = 0)
  # time_mcmc <- 2*time_mcmc
  # print(time_mcmc)
  
  # A function to run the compiled nimble model
  # Will be run with 1 chain in each cluster, starting from different inits
  return_results <- function(thischain2)
  {
    mcmc_nimble <- as.matrix(compiled.MCMC.model$mvSamples)
    return(mcmc_nimble)
  }
  # Outputs will be aggregated as a list
  mcmc_nimble <- parLapply(cl, seq_len(n.chains),  fun = return_results)  
  
  # Colate MCMC samples of different chains in a single mcmc.list object
  
  # save results ------------------------------------------------------------
  stopCluster(cl) 
  gc()
  cli_alert_success("Model run in parallel for {n.chains} chains in {time_mcmc_min}mins")
  
}

names(mcmc_nimble) <- glue("chain{seq_len(n.chains)}")

saveRDS(mcmc_nimble, file = name_file_mcmc)
saveRDS(time_mcmc, file = name_file_time)
cli_alert_success("Matrix Results Saved")

dfname <- get_df_varname(colnames(mcmc_nimble[[1]]))
mcmc_df <- get_mcmc(dfname, mcmc_nimble, burnin = 0)
saveRDS(mcmc_df, file = name_file_mcmc_df)
cli_alert_success("Dataframe Results Saved")
#rm(list = ls())

#gc() # garbage collector to try to clean R memory as much as possible.
#gc() # dobble gc() is better

