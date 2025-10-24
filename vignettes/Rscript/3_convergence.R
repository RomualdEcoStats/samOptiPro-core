cli_h1("Assess model {model_name} convergence")

cli_progress_step("Load MCMC posteriors")
# load matrix results if need be
if(!exists("mcmc_nimble")){
  mcmc_nimble <- readRDS(name_file_mcmc)
  time_mcmc <- readRDS(name_file_time)
}

cli_progress_step("Define variables for which convergence will be assessed")

cli_progress_done()
source("Rscript/3_calc_convergence.R")



# remove results from workspace
rm(mcmc_nimble)
# garbage collector to clean memory
gc()

source("Rscript/3_plot_convergence.R")

