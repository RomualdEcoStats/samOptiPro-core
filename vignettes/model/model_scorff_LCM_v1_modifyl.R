model.nimble <- nimbleCode({
  
  # Observation equations ---------------------------------------------------
  
  ## Total smolt abundances estimates -----
  ## we estimate the total true abundance of smolts (N3_tot) using the log of the median, the mean and the standard deviation estimated from a previous CMR model
  
  # Fit of N3_tot to available data  with observation errors (Pseudolikelihood method)
  # from 1996 to 2019
  for (t in 3:(nyear-1)) {
    log_N3_mu[t] ~ dnorm(log(N3_tot[t]), sd = log_N3_sd[t])
  }
  
  
  ## Smolt Sex ratio -----
  ## use of the number of individuals for each sex  (SAMARCH data) to estimate the sex ratio (p_sex_smolt) (multinomial)
  
  # from 1996 to 2018
  for (t in 3:(nyear-2))  {
    ## Remplacement de dhyperg par une approximation binomiale pour assurer la différentiabilité
    #RH sex_ratio_smolt[t,1] ~ dhyperg(N1 = round(prop3f[t]*N3_tot[t]), 
    #RH                              N0 = round((1-prop3f[t])*N3_tot[t]),
    #RH                              nx = N_sample_smolt_sex[t])
    sex_ratio_smolt[t,1] ~ dbinom(prob = prop3f[t], size = N_sample_smolt_sex[t])  #RH
  }
  
  ## smolt age proportion ----------------------------------------------------
  # 1997 to 2019
  for(t in 4:(nyear-1)){
    p_smolt_dirch[t,1] <- p_smolt1_migr[t]*N_sample_smolt_age
    p_smolt_dirch[t,2] <- (1-p_smolt1_migr[t])*N_sample_smolt_age
    p_smolt_migr_data[t, 1:2] ~ ddirich(p_smolt_dirch[t,1:2])
  }
  
  
  
  ## Total adult abundances estimates ----
  ## we estimate the total true abundance of adults (N6/N9) using the log of the median, the mean and the standard deviation estimated from a previous CMR model
  
  # from 1994 to 2019
  for (t in 1:(nyear-1))  {
    log_N6_mu[t] ~ dnorm(log(N6[t]), sd = log_N6_sd[t])
  }
  # from 1994 to 2020
  for (t in 1:(nyear))  {
    log_N9_mu[t] ~ dnorm(log(N9[t]), sd = log_N9_sd[t])
  }
  
  
  
  ## adult Sex ratio -----
  ## use of the number of individuals for each sex  (SAMARCH data) to estimate the sex ratio for 1SW and 2SW (prop6f, prop9f) (multinomial)
  
  # from 1994 to 2019
  for (t in 1:(nyear-1)){
    ## Remplacement de dhyperg par une approximation binomiale pour assurer la différentiabilité
    #RH sex_ratio_1SW[t,1] ~ dhyperg(N1 = round(prop6f[t]*N6[t]), 
    #RH                              N0 = round((1-prop6f[t])*N6[t]),
    #RH                              nx = N_sample_1SW[t])
    sex_ratio_1SW[t,1] ~ dbinom(prob = prop6f[t], size = N_sample_1SW[t])  #RH
  }
  
  # from 1994 to 2020
  for (t in 1:(nyear)){
    ## Remplacement de dhyperg par une approximation binomiale pour assurer la différentiabilité
    #RH sex_ratio_2SW[t,1] ~ dhyperg(N1 = round(prop9f[t]*N9[t]), 
    #RH                              N0 = round((1-prop9f[t])*N9[t]),
    #RH                              nx = N_sample_2SW[t])
    sex_ratio_2SW[t,1] ~ dbinom(prob = prop9f[t], size = N_sample_2SW[t])  #RH
  }
  
  
  ## Homewater Catches -------------------------------------------------------
  # from 1994 to 2019
  for (t in 1:(nyear-1))  {
    log_C6_mu[t] ~ dnorm(log(C6_hw[t]), sd = log_C6_sd[t])
  }
  # from 1994 to 2020
  for (t in 1:(nyear))  {
    log_C9_mu[t] ~ dnorm(log(C9_hw[t]), sd = log_C9_sd[t])
  }
  
  
  
  # Priors ------------------------------------------------------------------
  
  ## Stock-recruitment parameters --------------------------------------------
  
  logN2_sd ~ dunif(0,5)
  alpha ~ dbeta(theta1_max*nsample_theta1, (1-theta1_max)*nsample_theta1)
  k ~ dlnorm(meanlog = logk_pr, sdlog = 1)
  
  ## Smolt Abundance & Sex Ratio  ---------------------
  
  # year 1996 
  N3_tot[3] ~ dlnorm(0,0.01)
  
  # year 1996-2018
  for (t in 3:(nyear-2)){
    prop3f[t] ~ dbeta(10,10)
  }
  
  
  ## Smolt freshwater age proportion ----------------------------------------------------
  # 1994 to 2017
  for(t in 1:(nyear-3)){
    p_smolt_cohort[t,1] ~ dbeta(2,2)
    p_smolt_cohort[t,2] <- 1-p_smolt_cohort[t,1]
  }
  
  ## Post-smolt survival  ----------------------
  # survival during first summer (probability for individuals to survive)
  # psurv_smolt_F[t,k] is for smolt survival, in year t and size class k. 
  # Here we consider only one size class so k=1.
  
  # year 1996-2018
  for (t in 3:(nyear-2)){
    theta3[t] <- ilogit(logit_theta3[t])
    # survival in logit scale
    logit_theta3[t] ~ dnorm(mu_logit_theta3,sd = sd_theta3)
  }
  sd_theta3 ~ dunif(0,5)
  mu_logit_theta3 ~ dnorm(0,0.1)
  
  ## Probability of maturation  -----------
  # maturation during first summer (probability for individuals to mature)
  # pmat1_F[t,k] is for maturation at end of first summer, in year t and size class k. 
  
  # year 1997-2019
  for (t in 4:(nyear-1)){
    theta4f[t] <- ilogit(logit_theta4f[t])
    theta4m[t] <- ilogit(logit_theta4m[t])
    theta4[t] <- prop3f[t-1]*theta4f[t]+ (1-prop3f[t-1])*theta4m[t]
    
    #maturation in logit scale: intercept + temporal trend
    logit_theta4m[t] ~ dnorm(mu_logit_theta4m,sd = sd_theta4m)
    logit_theta4f[t] ~ dnorm(mu_logit_theta4f,sd = sd_theta4f)
    
  }
  sd_theta4f ~ dunif(0,5)
  sd_theta4m ~ dunif(0,5)
  mu_logit_theta4f ~ dnorm(0,0.1)
  mu_logit_theta4m ~ dnorm(0,0.1)
  
  ## Post PFA mortality ------------------------------------------------------
  # sync with cohort maturing year t
  # year 1997-2019
  for (t in 4:(nyear-1)){
    ### Survival probabilities -----
    # s1 and s2 are constant rates.
    # s1 is the additional survival probability during the first year at sea for 1SW adults
    # s2 is the additional survival probability during the first year at sea and during the second year at sea for 2SW adults
    # M is the mortality rate per additional month at sea after the end of first summer (november). Delta_t is the number of additional month at sea. 
    # 1SW fish return on next July, 2SW fish stay an additional year et return on the next march
    logit_theta5[t] <- logit(exp(-E_M*(delta5)))
    logit_theta8[t] <-  logit(exp(-E_M*(delta8)))
    theta5[t] <- ilogit(logit_theta5[t])
    theta8[t] <- ilogit(logit_theta8[t])
    
  }
  
  
  ## Exploitation Rates ------------------------------------------------------
  
  # year 1994-2019
  for(t in 1:(nyear-1)){
    h6_hw[t] ~ dbeta(1,2)
  }
  # year 1994-2020
  for(t in 1:(nyear)){
    h9_hw[t] ~ dbeta(1,2)
  }
  
  ## Return Abundance & sex-ratio 1994-1997 ----------------------------
  
  # 1SW : year 1994-1996
  for(t in 1:3){
    prop6f[t] ~ dunif(0,1)
    N6[t] ~ dlnorm(meanlog = 6.5, sdlog = 1)
    N6f[t] <- N6[t] * prop6f[t]
    N6m[t] <- N6[t] * ( 1 - prop6f[t] )
  }
  
  # 2SW : year 1994-1997
  for(t in 1:4){
    prop9f[t] ~ dunif(0,1)
    N9[t] ~ dlnorm(meanlog = 4.5, sdlog = 1)
    N9f[t] <- N9[t] * prop9f[t]
    N9m[t] <- N9[t] * ( 1 - prop9f[t] )
  }
  
  # Life Cycle  -----------------------------------------------------------
  # (indices are available up to n.year)
  
  ## Egg deposition and Stock-Recruitment -----------------
  
  # 1994 to 2017
  for (t in 1:(nyear-3)){
    N1[t] <- 
      N7f[t] * eggs[1,t] + 
      N10f[t] * eggs[2,t]
    
    # Stock 
    logN2_mean[t] <-  log(N1[t] / (1/alpha + N1[t]/k))
    N2[t] ~ dlnorm(logN2_mean[t]-0.5*logN2_sd^2, sd = logN2_sd)
  }
  
  ## Parr-smolt transition ---------------
  # 1994 to 2017
  for(t in 1:(nyear-3)){
    for(fw in 1:2){
      N3[t+1+fw,fw] <- N2[t]*p_smolt_cohort[t,fw]
    }
  }
  
  # 1997 to 2019
  # N3_tot[1] has a direct prior
  for(t in 4:(nyear-1)){
    N3_tot[t] <- N3[t,1] + N3[t,2]
    p_smolt1_migr[t] <- N3[t,1]/N3_tot[t]
  }
  for(t in 3:(nyear-2)){
    N3f_tot[t] <- N3_tot[t] * prop3f[t]
    N3m_tot[t] <- N3_tot[t] * (1 - prop3f[t])
  }
  
  ## post-smolt survival -----------------
  ### N4 is PFA
  # N3_tot: 1996 to 2018
  for (t in 3:(nyear-2)){
    N4f[t+1] <- N3f_tot[t] * theta3[t]
    N4m[t+1] <- N3m_tot[t] * theta3[t]
    N4[t+1]  <- N4m[t+1]   + N4f[t+1]
  }
  
  ## maturation -----------------
  # 1997 to 2019
  for (t in 4:(nyear-1)){
    # Dynamic equation for maturing individuals
    # Total number of maturing fish
    N5f[t] <- N4f[t]  * theta4f[t]  
    N5m[t] <- N4m[t]  * theta4m[t] 
    N5[t] <- N5f[t]+N5m[t]
    # Total number of non maturing fish
    N8f[t] <-  N4f[t]  * (1-theta4f[t])
    N8m[t] <-  N4m[t]  * (1-theta4m[t]) 
    N8[t] <- N8f[t]+N8m[t]
  }
  
  
  ### total numbers N6f, N6m, N9f, N9m
  #total number of individuals that survived, grew and matured in each sea age and sex categories should be multiplied by the additional survival, differently for 1SW and 2SW
  
  # 1997 to 2019
  for (t in 4:(nyear-1))
  {	
    # 1SW F and M (maturing individuals * survival s1 between sum1 and return)
    N6f[t] <-  N5f[t] * theta5[t]
    N6m[t] <-  N5m[t] * theta5[t]
    N6[t] <- N6f[t]+N6m[t]
    prop6f[t] <- N6f[t]/N6[t] 
    
    # 2SW F and M (non maturing individuals * survival s2 between sum1 and return)
    N9f[t+1] <-  N8f[t] * theta8[t]
    N9m[t+1] <-  N8m[t] * theta8[t]
    N9[t+1] <- N9f[t+1]+N9m[t+1]
    prop9f[t+1] <- N9f[t+1]/N9[t+1]
  }
  
  
  ## Homewater Catches  ------------------------------------------------------
  # 1994 to 2019
  for(t in 1:(nyear-1)){
    C6_hw[t] <- N6[t]*h6_hw[t]
    N7f[t] <- N6f[t]*(1-h6_hw[t])
    N7m[t] <- N6m[t]*(1-h6_hw[t])
    N7[t] <- N7f[t]+N7m[t]
  }
  # 1994 to 2020
  for(t in 1:(nyear)){
    C9_hw[t] <- N9[t]*h9_hw[t]
    N10f[t] <- N9f[t]*(1-h9_hw[t])
    N10m[t] <- N9m[t]*(1-h9_hw[t])
    N10[t] <- N10f[t]+N10m[t]
  }
  
  
  ## End of population dynamic process  ----------------
  # end model
})
