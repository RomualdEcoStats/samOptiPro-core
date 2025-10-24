
file_median <- glue("outputs_hindcast/median_{model_name}.rds")
cli_h1("Posteriors median and quantiles")
dfres <- readRDS(glue("outputs_hindcast/mcmc_results_df_{model_name}.rds"))

cli_progress_step("Calculating median and quantiles")
dfres %>% 
  filter(iter > burnin_output) %>% 
  group_by(fullname, name, dim1, dim2, dim3) %>% 
  summarize(median = median(value),
            mean = mean(value),
            min = min(value),
            max = max(value),
            q2.5 = quantile(value, p = 0.025),
            q5 = quantile(value, p = 0.05),
            q10 = quantile(value, p = 0.10),
            q25 = quantile(value, p = 0.25),
            q75 = quantile(value, p = 0.75),
            q90 = quantile(value, p = 0.90),
            q95 = quantile(value, p = 0.95),
            q97.5 = quantile(value, p = 0.975),
            .groups = 'drop') -> dfmed
saveRDS(dfmed, file = file_median)

cli_progress_done()
cli_h2("Plotting posteriors")

logfile <- glue("log/posteriors_{model_name}.log")
logfile_connection <- file(logfile, open = "wt")
cli_alert_info("log registered in {logfile}")

# Smolt Sex Ratio & Freshwater age ----------------------------------

cli_progress_step("Plotting smolt sex ratio (prop3f)")
sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection, type = "message", append = TRUE)
render(input= "Rscript/4_plot_smolt_sex_ratio.Rmd",
       output_file = glue("../Figures/Parameter Estimation/smolt_sex_ratio_{model_name}.html"), 
       # output_file is relative to where Rmd is located
       quiet=TRUE)
sink()
sink(type = "message")

# Stock-Recruitment ----------------------------------

cli_progress_step("Plotting Stock-Recruitment (N1-N2)")
sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection, type = "message", append = TRUE)
render(input= "Rscript/4_plot_stock_recruitment.Rmd",
       output_file = glue("../Figures/Parameter Estimation/stock_recruitment_{model_name}.html"), 
       # output_file is relative to where Rmd is located
       quiet=TRUE)
sink()
sink(type = "message")
