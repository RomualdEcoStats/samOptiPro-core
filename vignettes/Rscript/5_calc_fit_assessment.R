
logfile <- glue("log/fit_assessment_{model_name}.log")
logfile_connection <- file(logfile, open = "wt")
cli_alert_info("log registered in {logfile}")

# Fit Smolt Abundance -----------------------------------------------------

cli_progress_step("Calculates fit assessment for Smolt Abundance N3_tot")
sink(file = logfile_connection)
sink(file = logfile_connection,type = "message", append = TRUE)
subdf <- filter(dfres,
                name == "N3_tot") %>% 
  mutate(year = dim1) %>% 
  filter(iter > burnin_output, 
         year >= 3)

subdf %>% 
  mutate(type = "output") %>% 
  select(type, year, value ) -> dfAbundance
list_year <- c()
list_SU <- c()
simul <- c()
nsample <- nrow(filter(subdf, year == 3))


for(thisyear in unique(dfAbundance$year)){
  
  this_sd <- Const_nimble$log_N3_sd[thisyear]
  list_year <- c(list_year,rep(thisyear,nsample))
  simul <- c(simul,rlnorm(nsample, 
                          meanlog = Data_nimble$log_N3_mu[thisyear],
                          sdlog = this_sd)
  )
}

tmp2 <- data.frame(type = "data", 
                   year = list_year,
                   value = simul)
dfAbundance <- rbind(dfAbundance, tmp2)
saveRDS(dfAbundance, 
        glue("outputs_figures/Fit Assessment/fit_N3_tot_{model_name}.rds"))
sink()
sink(type = "message")

# Fit 1SW Abundance -----------------------------------------------------

cli_progress_step("Calculates fit assessment for 1SW Abundance N6")

sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection,type = "message", append = TRUE)
subdf <- filter(dfres,
                name == "N6",
                iter > burnin_output) %>% 
  mutate(year = dim1)

subdf %>% 
  mutate(type = "output") %>% 
  select(type, year, value ) -> dfAbundance
list_year <- c()
list_SU <- c()
simul <- c()
nsample <- nrow(filter(subdf, year == 5))


for(thisyear in unique(dfAbundance$year)){
  
  this_sd <- Const_nimble$log_N6_sd[thisyear]
  list_year <- c(list_year,rep(thisyear,nsample))
  simul <- c(simul,rlnorm(nsample, 
                          meanlog = Data_nimble$log_N6_mu[thisyear],
                          sdlog = this_sd)
  )
}

tmp2 <- data.frame(type = "data", 
                   year = list_year,
                   value = simul)
dfAbundance <- rbind(dfAbundance, tmp2)
saveRDS(dfAbundance, 
        glue("outputs_figures/Fit Assessment/fit_N6_{model_name}.rds"))
sink()
sink(type = "message")


# Fit 2SW Abundance -----------------------------------------------------

cli_progress_step("Calculates fit assessment for 2SW Abundance N9")

sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection,type = "message", append = TRUE)
subdf <- filter(dfres,
                name == "N9",
                iter > burnin_output) %>% 
  mutate(year = dim1)

subdf %>% 
  mutate(type = "output") %>% 
  select(type, year, value ) -> dfAbundance
list_year <- c()
list_SU <- c()
simul <- c()
nsample <- nrow(filter(subdf, year == 5))


for(thisyear in unique(dfAbundance$year)){
  
  this_sd <- Const_nimble$log_N9_sd[thisyear]
  list_year <- c(list_year,rep(thisyear,nsample))
  simul <- c(simul,rlnorm(nsample, 
                          meanlog = Data_nimble$log_N9_mu[thisyear],
                          sdlog = this_sd)
  )
}

tmp2 <- data.frame(type = "data", 
                   year = list_year,
                   value = simul)
dfAbundance <- rbind(dfAbundance, tmp2)
saveRDS(dfAbundance, 
        glue("outputs_figures/Fit Assessment/fit_N9_{model_name}.rds"))
sink()
sink(type = "message")


# fit Smolt sex ratio -------------------------------------------------------

cli_progress_step("Calculates fit assessment for Smolt Sex Ratio")

sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection,type = "message", append = TRUE)
subdf <- filter(dfres,
                name == "prop3f") %>% 
  mutate(year = dim1) %>% 
  filter(iter > burnin_output, 
         year >= 3)
subdf %>% 
  mutate(type = "output") %>% 
  select(type, year, value ) -> dfAbundance

list_year <- c()
list_SU <- c()
simul <- c()
nsample <- nrow(filter(subdf, year == 5))

list_prob <- Data_nimble$sex_ratio_smolt/
  rowSums(Data_nimble$sex_ratio_smolt)
thisyear <- 3
for(thisyear in unique(dfAbundance$year)){
  list_year <- c(list_year,rep(thisyear,nsample))
  for(i in seq_len(nsample)){
    # simul <- c(simul,
    #            (rmulti(n = 1, 
    #                    size = Const_nimble$N_sample_smolt_sex[thisyear],
    #                    prob = list_prob[thisyear,])/
    #               Const_nimble$N_sample_smolt_sex[thisyear])[1]
    # )
    simul <- c(simul,         
               (rhyper(nn = 1,
                       k =  Const_nimble$N_sample_smolt_sex[thisyear],
                       m = list_prob[thisyear,1]*exp(Data_nimble$log_N3_mu[thisyear]),
                       n = list_prob[thisyear,2]*exp(Data_nimble$log_N3_mu[thisyear])
               )
               )/Const_nimble$N_sample_smolt_sex[thisyear]
    )
  }
}
tmp2 <- data.frame(type = "data", 
                   year = list_year,
                   value = simul)
dfAbundance <- rbind(dfAbundance, tmp2)
saveRDS(dfAbundance, 
        glue("outputs_figures/Fit Assessment/fit_sex_ratio_smolt_{model_name}.rds"))
sink()
sink(type = "message")

# fit 1SW sex ratio -------------------------------------------------------

cli_progress_step("Calculates fit assessment for 1SW Sex Ratio")

sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection,type = "message", append = TRUE)
subdf <- filter(dfres,
                name == "prop6f") %>% 
  mutate(year = dim1) %>% 
  filter(iter > burnin_output)

subdf %>% 
  mutate(type = "output") %>% 
  select(type, year, value ) -> dfAbundance

list_year <- c()
list_SU <- c()
simul <- c()
nsample <- nrow(filter(subdf, year == 5))

list_prob <- Data_nimble$sex_ratio_1SW/
  rowSums(Data_nimble$sex_ratio_1SW)

for(thisyear in unique(dfAbundance$year)){
  list_year <- c(list_year,rep(thisyear,nsample))
  for(i in seq_len(nsample)){
    # simul <- c(simul,
    #            (rmulti(n = 1, 
    #                    size = Const_nimble$N_sample_1SW[thisyear],
    #                    prob = list_prob[thisyear,])/
    #               Const_nimble$N_sample_1SW[thisyear])[1]
    # )
    tmphyper <- rhyper(nn = 1,
                       k =  Const_nimble$N_sample_1SW[thisyear],
                       m = list_prob[thisyear,1]*exp(Data_nimble$log_N6_mu[thisyear]),
                       n = list_prob[thisyear,2]*exp(Data_nimble$log_N6_mu[thisyear]))
    if(is.na(tmphyper)) message(thisyear)
    simul <- c(
      simul, 
      tmphyper/Const_nimble$N_sample_1SW[thisyear]
    )
    
  }
}
tmp2 <- data.frame(type = "data", 
                   year = list_year,
                   value = simul)
dfAbundance <- rbind(dfAbundance, tmp2)
saveRDS(dfAbundance, 
        glue("outputs_figures/Fit Assessment/fit_sex_ratio_1SW_{model_name}.rds"))
sink()
sink(type = "message")

# fit 2SW sex ratio -------------------------------------------------------

cli_progress_step("Calculates fit assessment for 2SW Sex Ratio")

sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection,type = "message", append = TRUE)
subdf <- filter(dfres,
                name == "prop9f",
                iter > burnin_output) %>% 
  mutate(year = dim1) 

subdf %>% 
  mutate(type = "output") %>% 
  select(type, year, value ) -> dfAbundance

list_year <- c()
list_SU <- c()
simul <- c()
nsample <- nrow(filter(subdf, year == 5))

list_prob <- Data_nimble$sex_ratio_2SW/
  rowSums(Data_nimble$sex_ratio_2SW)

for(thisyear in unique(dfAbundance$year)){
  list_year <- c(list_year,rep(thisyear,nsample))
  for(i in seq_len(nsample)){
    # simul <- c(simul,
    #            (rmulti(n = 1, 
    #                    size = Const_nimble$N_sample_2SW[thisyear],
    #                    prob = list_prob[thisyear,])/
    #               Const_nimble$N_sample_2SW[thisyear])[1]
    # )
    tmp <-  (rhyper(nn = 1,
                    k =  Const_nimble$N_sample_2SW[thisyear],
                    m = list_prob[thisyear,1]*exp(Data_nimble$log_N9_mu[thisyear]),
                    n = list_prob[thisyear,2]*exp(Data_nimble$log_N9_mu[thisyear])
    )
    )/Const_nimble$N_sample_2SW[thisyear]
    if(is.na(tmp))  message(thisyear, " = ", tmp)
    simul <- c(simul,tmp)
  }
}
tmp2 <- data.frame(type = "data", 
                   year = list_year,
                   value = simul)
dfAbundance <- rbind(dfAbundance, tmp2)
saveRDS(dfAbundance, 
        glue("outputs_figures/Fit Assessment/fit_sex_ratio_2SW_{model_name}.rds"))
sink()
sink(type = "message")

cli_process_done()



# Fit 1SW Homewater catches -----------------------------

cli_progress_step("Calculates fit assessment for 1SW Homewater catches C6_hw")

sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection,type = "message", append = TRUE)
subdf <- filter(dfres,
                name == "C6_hw",
                iter > burnin_output) %>% 
  mutate(year = dim1) 

subdf %>% 
  mutate(type = "output") %>% 
  select(type, year, value ) -> dfAbundance
list_year <- c()
list_SU <- c()
simul <- c()
nsample <- nrow(filter(subdf, year == 5))


for(thisyear in unique(dfAbundance$year)){
  
  this_sd <- Const_nimble$log_C6_sd[thisyear]
  list_year <- c(list_year,rep(thisyear,nsample))
  simul <- c(simul,rlnorm(nsample, 
                          meanlog = Data_nimble$log_C6_mu[thisyear],
                          sdlog = this_sd)
  )
}

tmp2 <- data.frame(type = "data", 
                   year = list_year,
                   value = simul)
dfAbundance <- rbind(dfAbundance, tmp2)
saveRDS(dfAbundance, 
        glue("outputs_figures/Fit Assessment/fit_C6_hw_{model_name}.rds"))
sink()
sink(type = "message")

# Fit 2SW Homewater catches -------------------------------

cli_progress_step("Calculates fit assessment for 2SW Homewater catches C9_hw")

sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection,type = "message", append = TRUE)
subdf <- filter(dfres,
                name == "C9_hw",
                iter > burnin_output) %>% 
  mutate(year = dim1) 

subdf %>% 
  mutate(type = "output") %>% 
  select(type, year, value ) -> dfAbundance
list_year <- c()
list_SU <- c()
simul <- c()
nsample <- nrow(filter(subdf, year == 5))


for(thisyear in unique(dfAbundance$year)){
  
  this_sd <- Const_nimble$log_C9_sd[thisyear]
  list_year <- c(list_year,rep(thisyear,nsample))
  simul <- c(simul,rlnorm(nsample, 
                          meanlog = Data_nimble$log_C9_mu[thisyear],
                          sdlog = this_sd)
  )
}

tmp2 <- data.frame(type = "data", 
                   year = list_year,
                   value = simul)
dfAbundance <- rbind(dfAbundance, tmp2)
saveRDS(dfAbundance, 
        glue("outputs_figures/Fit Assessment/fit_C9_hw_{model_name}.rds"))
sink()
sink(type = "message")


# fit Smolt Freshwater age -------------------------------------------------------

cli_progress_step("Calculates fit assessment for Smolt Freshwater age")

sink(file = logfile_connection, append = TRUE)
sink(file = logfile_connection,type = "message", append = TRUE)
subdf <- filter(dfres,
                name == "p_smolt_dirch") %>% 
  mutate(year = dim1,
         FW = dim2) %>% 
  filter(iter > burnin_output, 
                year >= 4,
                FW == 1)

subdf %>% 
  mutate(type = "output",
         value = value/1000) %>% 
  select(type, year, value ) -> dfAbundance

list_year <- c()
list_SU <- c()
simul <- c()
nsample <- nrow(filter(subdf, year == 5))

list_prob <- Data_nimble$p_smolt_migr_data
thisyear <- 3
for(thisyear in unique(dfAbundance$year)){
  list_year <- c(list_year,rep(thisyear,nsample))
  add_simul <- rbeta(n = nsample, 
                     shape1 = Const_nimble$N_sample_smolt_age *
                       list_prob[thisyear,1],
                     shape2 = Const_nimble$N_sample_smolt_age *
                       list_prob[thisyear,2])
  # if(is.na(add_simul)) break()
  simul <- c(simul,add_simul)
  
}
tmp2 <- data.frame(type = "data", 
                   year = list_year,
                   value = simul)
dfAbundance <- rbind(dfAbundance, tmp2)
saveRDS(dfAbundance, 
        glue("outputs_figures/Fit Assessment/fit_freshwater_age_smolt_{model_name}.rds"))
sink()
sink(type = "message")
