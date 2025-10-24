# Use data and Simplified back-calculation from model equation to find appropriate 
# fixed values to constrain Nimble simulation in order to create appropriate inits 
# for the life cycle
# 14/12/2021
# R. Patin

generate_fixed_values <- function(Data_nimble, Const_nimble, n.chains = 3,
																	logfile_connection){
	set.seed(42)
	# For loop
	
	nyear     <- Const_nimble$nyear
	nSm   <- Const_nimble$nSm
	
	Fixed_for_simul <- list()
	
	for(thischain in seq_len(n.chains)){
		set.seed(thischain)
		cli_process_start("Generate fixed values for chain {thischain}")
		sink(file = logfile_connection, append = TRUE)
		sink(file = logfile_connection,type = "message",append = TRUE)
		
		
		# Setup -------------------------------------------------------------------
		
		
		## Variable to fill --------------------------------------------------------
		
		
		# Abundances
		N1       <- rep(NA, nyear-4)
		logN2 <- rep(NA, nyear-4)
		N2 		<- rep(NA, nyear-4)
		N3     <- matrix(NA,nyear+3,nSm)
		N3_tot <- rep(NA,nyear-2)
		N5     <- rep(NA,nyear-1)
		N5_1   <- rep(NA,nyear-1)
		N7     <- rep(NA,nyear-1)
		N8     <- rep(NA,nyear-1)
		N8_1   <- rep(NA,nyear-1)
		N10    <- rep(NA,nyear)
		
		# Harvest rates
		
		h6_hw       <- rep(NA,nyear-1)
		h9_hw       <- rep(NA,nyear)
		h5       <- rep(NA,nyear-1)
		h8       <- rep(NA,nyear-1)
		
		# Captures
		
		
		## sample from data distribution -------------------------------------------
		
		C5  <- rlnorm_vec(Data_nimble$log_C5_mu, Const_nimble$log_C5_sd)
		C8  <- rlnorm_vec(Data_nimble$log_C8_mu, Const_nimble$log_C8_sd)
		sd_hw <- sqrt(log(Const_nimble$CV_hw*Const_nimble$CV_hw + 1))
		C6_hw <- rlnorm_vec(Data_nimble$log_C6_mu, Const_nimble$log_C6_sd)
		C9_hw <- rlnorm_vec(Data_nimble$log_C9_mu, Const_nimble$log_C9_sd)
		N6 <- rlnorm_vec(Data_nimble$log_N6_mu, Const_nimble$log_N6_sd)
		N9 <-rlnorm_vec(Data_nimble$log_N9_mu, Const_nimble$log_N9_sd)
		C6_hw <- pmin(0.95*N6,C6_hw)
		C9_hw <- pmin(0.95*N9,C9_hw)
		
		
		## Initialization variables ----------------------------------------------------------
		
		N1_pr <- rep(NA,(nSm +2))
		N2_pr <- rep(NA,(nSm +2))
		N3_pr <- matrix(NA,2*nSm+3,nSm)
		log_N1_pr     <- rep(NA,(nSm +2))
		
		p_smolt_pr <- Const_nimble$p_smolt_pr
		p_smolt <- Const_nimble$p_smolt
		
		## Survival & maturation -------------------------------------------------------------------
		
		M <- rep(Const_nimble$E_M, nyear)
		# Survival
		theta8   <- exp(-  M  * (Const_nimble$deltat8))
		theta5   <- exp(-  M  * (Const_nimble$deltat5))
		theta4 <- rep(NA,nyear-1)
		theta3 <- rep(NA,nyear-2)
		
		
		
		# Start life-cycle --------------------------------------------------------
		## From N9 to N3 (forward) -------------------------------------------------
		### N7 & N10 ----------------------------------------------
		h6_hw  <- C6_hw[1:(nyear-1)] / N6[1:(nyear-1)]
		h9_hw  <- C9_hw / N9
		
		
		N7[1] <- N6[1] * (1- h6_hw[1])
		N10[1] <- N9[1] * (1-h9_hw[1])  + Const_nimble$Stocking_2SW[1]
		
		for (t in 2:(nyear-1)){
			N7[t]  <- max(  ((N6[t] * (1- h6_hw[t])) * (1 - Const_nimble$theta6_delSp[t])) +
												((N6[t-1] * (1- h6_hw[t-1])) * Const_nimble$theta6_delSp[t-1]), 1)
		}
		for (t in 2:nyear){
			N10[t] <- max(((N9[t] * (1- h9_hw[t]))* (1 - Const_nimble$theta9_delSp[t])) +
											((N9[t-1] * (1- h9_hw[t-1]))* Const_nimble$theta9_delSp[t-1]) + Const_nimble$Stocking_2SW[t], 1)
		}
		
		
		
		### N1 ----------------------------------------------
		
		for (t in 1:(nyear-4)){
			N1[t] <- N7[t]*Const_nimble$eggs[1,t]*Const_nimble$prop_female[1,t] + N10[t]*Const_nimble$eggs[2,t]*Const_nimble$prop_female[2,t]
		}
		
		
		
		### N2 ----------------------------------------------
		for (t in 1:(nyear-4)){
			logN2 [t] <- log( Const_nimble$E_theta1[t]*N1[t] + 0.1)
			N2[t] <- exp(logN2[t])
		}
		
		
		
		### N3 ----------------------------------------------
		
		for (t in 1:(nyear-4)){
			for (k in 1:nSm){
				N3[t+1+k,k] <- p_smolt[t,k]*N2[t]
			}
		}
		
		
		### N1 1963-1970 ----------------------------------------------
		CV_theta1_pr <- Const_nimble$CV_theta1_pr
		sd_N1_pr <- sqrt(log(CV_theta1_pr^2 + 1))
		
		for(t in 1:(nSm+2)){
			log_N1_pr[t] <- rnorm(1,Const_nimble$mu_N1_pr,sd_N1_pr)
			N1_pr[t] <- exp(log_N1_pr[t])
		}
		
		
		
		### N2 1963-1970 ----------------------------------------------
		for(t in 1:(nSm+2)) {
			N2_pr[t] <- N1_pr[t]*Const_nimble$E_theta1[t]
		}  
		
		### N3 ----------------------------------------------
		
		# Reparition of smolt for years 63 to 70
		for (k in 1:nSm) {
			for (t in 1:(nSm+2)) {
				N3_pr[t+k+1,k] <- p_smolt_pr[t,k] * N2_pr[t]
			}
		}
		
		# The year 71 correspond to the 8th line in N3_pr (from 63 to 70)
		# The first year that have to be completed
		
		for(k in 1:nSm) {
			N3[1,k] <- N3_pr[8,k]
		}
		
		# step filling of years 72:77 from N3_pr of years 64 to 70
		for (k in 1:(nSm)) {
			for (kk in k:nSm) {
				N3[k+1,kk] <- N3_pr[k+8,kk]
			}
		}
		
		
		# N3_tot of all years 
		for (t in 1:(nyear-2)) { 
			N3_tot[t] <- sum(N3[t,1:nSm])
		}
		
		
		## From N9 to N4 (backward) ------------------------------------------------
		
		
		### N8_1 ----------------------------------------------
		
		for (t in 1:(nyear-1)){
			N8_1[t] <- N9[t+1] + C8[t]
		}
		
		
		### h8 ----------------------------------------------------------------------
		
		
		
		for (t in 1:(nyear-1)){
			h8[t] <- C8[t]/(N8_1[t]) 
		}
		
		
		
		### N8 ----------------------------------------------
		
		for (t in 1:(nyear-1)){  
			N8[t] <- N8_1[t] / theta8[t]
		}
		
		
		### N5_1 ----------------------------------------------
		
		for (t in 1:(nyear-1)){
			N5_1[t] <- N6[t] + C5[t]
		}
		
		
		for (t in 1:(nyear-1)){  
			h5[t] <- C5[t]/N5_1[t]
		}
		
		for (t in 1:(nyear-1)){
			N5[t] <- N5_1[t]/theta5[t]
		}
		
		### N4 ----------------------------------------------
		N4 <- N5+N8
		
		## Calculations theta3 & theta4 -----------------------
		
		### theta4 : ratio N5/(N5+N8) ----------------------------------------------
		for (t in 1:(nyear-1)){ 
				theta4[t] <- N5[t] / N4[t]
		}
		
		### theta3 - ratio N4 -> N3 ----------------------------------------------
		for (t in 1:(nyear-2)){ 
				theta3[t] <- min(N4[t+1]/N3_tot[t],0.8)
		}
		
		# Calculate logit_theta3 and logit_theta4
		logit_theta3 <- logit(theta3)
		logit_theta4 <- logit(theta4)
		sd_theta3 <- sd(logit_theta3) 
		sd_theta4 <- sd(logit_theta4)
		# Save variables ----------------------------------------------------------
		
		# Save variables that must be fixed to reasonnable values
		# when simulating from the Nimble model
		# Variables saved as a list
		log_N9_pr <- log(N9[1])
		thisfix <- list( N1_pr, N2_pr, log_N9_pr,  
										 N2, 
										 h6_hw, h9_hw,
										 h5, h8,
										 sd_theta3, sd_theta4,
										 logit_theta3, logit_theta4)
		
		
		
		
		thisnames <- 			c( "N1_pr", "N2_pr", "log_N9_pr",  
											"N2", 
											"h6_hw", "h9_hw",
											"h5", "h8",
											"sd_theta3", "sd_theta4",
											"logit_theta3", "logit_theta4")
		
		names(thisfix) <- thisnames
		Fixed_for_simul[[thischain]] <- thisfix
		
		sink()
		sink(type = "message")
		cli_process_done()
	}
	return(Fixed_for_simul)
}