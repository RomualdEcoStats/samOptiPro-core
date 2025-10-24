dfres <- readRDS(glue("outputs_hindcast/mcmc_results_df_{model_name}.rds"))

cli_h1("Assess model {model_name} fit to the data")

source("Rscript/5_calc_fit_assessment.R")


cli_progress_step("Plotting Fit Assessment")
sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection,type = "message", append = TRUE)
render(input= "Rscript/5_fit_assessment.Rmd",
       output_file = glue("../Figures/Model Fit/Fit_Assessment_{model_name}.html"), 
       # output_file is relative to where Rmd is located
       quiet=TRUE)
sink()
sink(type = "message")

cli_progress_done()
