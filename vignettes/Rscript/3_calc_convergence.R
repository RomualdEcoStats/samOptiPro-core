dfname <- get_df_varname(colnames(mcmc_nimble$chain1))
listvar <- unique(dfname$name)[unique(dfname$name) %in% var_to_look ]


logfile <- glue("log/convergence_{model_name}.log")
logfile_connection <- file(logfile, open = "wt")

cli_h2("Calculation for rhat and effective size, for model {model_name}")
cli_alert_info("log registered in {logfile}")
nvar <- length(listvar)
current_var <- listvar[1]
whichvar <- 1
totvar <- length(listvar)
cli_progress_step("Convergence calculations for variable {whichvar}/{totvar}: {current_var}")
for(current_var in listvar){
  whichvar <- which(listvar == current_var)
  cli_progress_update()
  filerhat <- glue("outputs_figures/Convergence/convergence_{current_var}_{model_name}.rds")
  sink(file = logfile_connection, append = TRUE)
  sink(file = logfile_connection,type = "message",append = TRUE)
  tmp <- dplyr::filter(dfname, name == current_var)
    fullrhat <- NULL
# if(whichvar == 8) sum("A","B")
    for(i in 1:nrow(tmp)){
      # message(i)
      
      niter.thinned <- nrow(mcmc_nimble[[1]])
      ch <- lapply(mcmc_nimble, function(tmpdf){
        return(coda::as.mcmc(tmpdf[(burnin_output+1):niter.thinned,tmp$fullname[i]]))
      })
      current_mcmc <- coda::as.mcmc.list(ch)

      # Gelman-Rubin Rhat
      current_rhat <- coda::gelman.diag(current_mcmc, confidence = 0.95, autoburnin = FALSE, multivariate = FALSE, transform = TRUE)
      # current_rhat <- rstan::rhat(current_matrix)
      current_rhat <- current_rhat$psrf
      current_rhat <- as.data.frame(current_rhat)
      current_rhat$fullname <- tmp$fullname[i]
      colnames(current_rhat) <- c("rhat","rhat_upperCI","fullname")
      tmp_rhat <- dplyr::left_join(current_rhat,tmp,by = "fullname")
      
      # Effective Size
      isnullvar <- lapply(ch, var)
      if(any(isnullvar == 0)){
        current_neff <- NA
      } else {
        current_neff <- coda::effectiveSize(current_mcmc)
      }
      current_neff <- data.frame("fullname"= tmp$fullname[i],
                                 "neff" = c(current_neff))
      tmp_rhat <- dplyr::left_join(current_neff,tmp_rhat,by = "fullname")
      if(is.null(fullrhat)){
        fullrhat <- tmp_rhat
      } else {
        fullrhat <- rbind(fullrhat, tmp_rhat)
      }
    }

    saveRDS(fullrhat, file = filerhat)
    sink()
    sink(type = "message")
  } 
current_var <- "All done"
cli_progress_update()
cli_progress_done()
