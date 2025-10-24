# allnames <- colnames(mcmc_nimble[[1]])
# 
# custome function to close sink() on error
handle_error <- function(...){
  # message("Additional Error Message")
  if(sink.number()>0){
    sink();sink(type = "message");
    message("\n An error was encountered, please check the relevant log file")
  }
}
options("error" = handle_error)


source("Rscript/functions/get_df_varname.R")
source("Rscript/functions/nimble_hyperg.R")
source("Rscript/functions/get_mcmc.R")