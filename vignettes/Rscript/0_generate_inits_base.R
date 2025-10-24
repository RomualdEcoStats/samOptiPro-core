fill_psmolt <- runif(24, 0.4, 0.6)
myinits <- 
  list(
    N3_tot = rlnorm(n = 26, meanlog = log(1000000), sdlog = 1),
    N2 = rlnorm(n = 24, meanlog = log(1000000), sdlog = 1),
    N6 = rlnorm(n = 26, meanlog = log(1000000), sdlog = 1),
    N9 = rlnorm(n = 27, meanlog = log(1000000), sdlog = 1),
    p_smolt = cbind(fill_psmolt, 1-fill_psmolt),
    logit_theta3 = runif(n= 25, -0.5, 0.5),
    logit_theta4f = runif(n= 26, -0.5, 0.5),
    logit_theta4m = runif(n= 26, -0.5, 0.5),
    logit_theta4 = runif(n= 26, -0.5, 0.5),
    logit_theta5 = rep(logit(0.8), 26),
    logit_theta8 = rep(logit(0.8), 26),
    prop3f = runif(n = 25, 0.4,0.6),
    prop6f = runif(n = 26, 0.4,0.6),
    prop9f = runif(n = 27, 0.4,0.6),
    logN2_sd = runif(n = 1, 0,5),
    alpha = rbeta(1,
                  Const_nimble$theta1_max*Const_nimble$nsample_theta1,
                  (1-Const_nimble$theta1_max)*Const_nimble$nsample_theta1),
    k = rlnorm(1, meanlog = Const_nimble$logk_pr, sdlog = 1)
    
  )







