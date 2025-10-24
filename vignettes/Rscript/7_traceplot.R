cli_h1("Generating Traceplot")


logfile <- glue("log/traceplot_{model_name}.log")
logfile_connection <- file(logfile, open = "wt")
cli_alert_info("log registered in {logfile}")


cli_progress_step("Load Posteriors")
dfres <- readRDS(glue("outputs_hindcast/mcmc_results_df_{model_name}.rds"))

cli_progress_step("Selecting variable to look at")
current_var <- var_to_look[1]
whichvar <- 1
totvar <- length(var_to_look)
cli_progress_step("Traceplot for variable {whichvar}/{totvar}: {current_var}")
for(current_var in var_to_look){
  whichvar <- which(var_to_look == current_var)
  cli_progress_update()
  
  current_file <- 
    glue("../Figures/Traceplot/traceplot_{current_var}_{model_name}.html")
  sink(file = logfile_connection,append = TRUE)
  sink(file = logfile_connection, type = "message",append = TRUE)
  render(input= "Rscript/7_traceplot.Rmd",
         output_file = current_file, 
         quiet=TRUE)
  sink()
  sink(type = "message")
  
}
current_var <- "All done"
cli_progress_update()
cli_progress_done()
