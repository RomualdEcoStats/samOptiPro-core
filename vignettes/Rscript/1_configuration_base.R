
name_file_mcmc <- glue("outputs_hindcast/mcmc_results_{model_name}.rds")
name_file_mcmc_df <- glue("outputs_hindcast/mcmc_results_df_{model_name}.rds")
# Time
name_file_time <- glue("outputs_hindcast/time_mcmc_results_{model_name}.rds")

# Load Monitors -----------------------------------------------------------

monitors <- c(
  
  # Abundances
  "N1","N2","N3",
  "N3_tot","N4","N5","N6","N7","N8","N9","N10",
  "N3f_tot","N4f","N5f","N6f","N7f","N8f","N9f","N10f",
  "N3m_tot","N4m","N5m","N6m","N7m","N8m","N9m","N10m",

  # Stock recruitment
  "logN2_mean","logN2_sd",
  "alpha","k",
  
  # Smolt ages
  "p_smolt_cohort",
  "p_smolt1_migr",  "p_smolt_dirch",

  # Sex-ratio
  "prop3f","prop6f","prop9f",

  # Post-smolt survival
  "theta3","logit_theta3",
  "sd_theta3",

  # Probability of maturation as 1SW
  "theta4f","theta4m","theta4", 
  "logit_theta4f","logit_theta4m",
  "sd_theta4f","sd_theta4m","mu_logit_theta3",  "mu_logit_theta4f", "mu_logit_theta4m",

  # survival from end sum1 to return
  "theta5","theta8",
  "logit_theta5","logit_theta8",
  
  "h6_hw","h9_hw",
  "C6_hw","C9_hw"
)


# Only variable for which convergence should be looked (here all variables)
var_to_look <- monitors