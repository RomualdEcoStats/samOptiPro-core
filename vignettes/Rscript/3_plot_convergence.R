cli_h2("Plotting Rhat and effective size")

current_var <- "logit_theta3"
# Loading Data ------------------------------------------------------------
totdf <- data.frame(fullname = character(),
                    neff = numeric(), 
                    rhat = numeric(),
                    rhat_upperCI = numeric(),
                    name = character(),
                    dim1 = numeric(), 
                    dim2 = numeric(), 
                    dim3 = numeric())
for(current_var in var_to_look){
  tmp <- readRDS(file = glue("outputs_figures/Convergence/convergence_{current_var}_{model_name}.rds"))
  totdf <- rbind(totdf, tmp)
}

# Reformatting data -------------------------------------------------------

cli_progress_step("Reformatting Rhat and Effective size data")

totdf$rhatfact <- cut(totdf$rhat, c(0,1, 1.05,1.1,1.5,5,10,Inf))
totdf$nefffact <- cut(totdf$neff, c(0,100,200,350,500,750,1000,5000,Inf))
totdf %>%
  group_by(name) %>% 
  summarize(neff_min =  min(neff, na.rm = TRUE),
            neff_q5 = quantile(neff, na.rm = TRUE, p = 0.05),
            median_rhat = median(rhat, na.rm = TRUE),
            prop_rhat = length(which(rhat > 1.05))/n()) %>% 
  mutate(neff.time = neff_q5/as.numeric(time_mcmc)*3600 ) -> totdf_summarized


totdf <- mutate(totdf, neff.time = neff/as.numeric(time_mcmc)*3600)
totdf$nefftime.fact <- cut(totdf$neff.time, c(0,5,10,20,30,40,50,Inf))





# Quantile 5% effective size ----------------------------------------------
min_neff <- quantile(totdf$neff, na.rm = TRUE, p = 0.05)
min_neff_round <- round(min_neff, digits = 2)

cli_progress_step("Plotting 5% quantiles in effective sample size")

g <- 
  ggplot(totdf_summarized) +
  geom_bar(aes(x = name, y = neff_q5), stat= "identity")+
  xlab("Variable")+
  ylab("Quantile 5% Effective sample size")+
  # scale_y_log10()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust= 1, vjust= 0.5))+
  geom_text(aes(x= name, y = neff_q5, label = paste0(round(neff_q5, digits = 0))), nudge_y = 100, size = 3)

cowplot::save_plot(glue("Figures/Convergence/quantile_0.05_neff_{model_name}.pdf"),
                   g, base_width = 20/cm(1), base_height = 12/cm(1))


# Distribution effective size ---------------------------------------------

cli_progress_step("Plotting distribution in effective sample size")

g <- ggplot(filter(totdf,!is.na(nefffact))) +
  geom_bar(aes(x = name, fill = nefffact), stat= "count", position = "stack")+
  xlab("Variable")+
  ylab("Parameter Number")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust= 1, vjust= 0.5))+
  scale_fill_brewer(type = "div", palette = "RdYlBu", drop = FALSE )
  # expand_limits(fill= mancut)+
  # geom_hline(yintercept = 0.95, linetype = 2, alpha = 0.5)

cowplot::save_plot(glue("Figures/Convergence/distrib_neff_{model_name}.pdf"),
                   g, base_width = 22/cm(1), base_height = 12/cm(1))


# Rhat Plot ---------------------------------------------------------------

cli_progress_step("Plotting median Gelman-Rubin Rhat")


# median rhat
g <- ggplot(totdf_summarized) +
  geom_bar(aes(x = name, y = median_rhat-1), stat= "identity")+
  xlab("Variable")+
  ylab("Median of Gelman-Rubin Rhat")+
  geom_hline(yintercept = 0.05, col = "red", linetype = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust= 1, vjust= 0.5)) +
  scale_y_continuous(labels = function(y) y + 1)

cowplot::save_plot(glue("Figures/Convergence/median_rhat_{model_name}.pdf"),
                   g, base_width = 20/cm(1), base_height = 12/cm(1))

cli_progress_step("Plotting Summarized Distribution Gelman-Rubin Rhat")

# proportion > 1.05
g <- ggplot(filter(totdf, !is.na(rhatfact)))+
  geom_bar(aes(x= name, fill = rhatfact), position = 'stack')+
  scale_fill_manual("Gelman-Rubin Rhat", values = c("#4575b4","#abd9e9", "#fee090","#fdae61","#f46d43","#d73027","#a50026"))+
  xlab("Variable")+
  ylab("Number of parameters in class")+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust= 1, vjust= 0.5))

cowplot::save_plot(glue("Figures/Convergence/distrib_rhat_{model_name}.pdf"),
                   g, base_width = 22/cm(1), base_height = 12/cm(1))


# full distrib rhat
cli_progress_step("Plotting Full Distribution Gelman-Rubin Rhat")
g <- ggplot(totdf) +
  geom_boxplot(aes(x = name, y = rhat-1))+
  xlab("Variable")+
  ylab("Gelman-Rubin Rhat")+
  geom_hline(yintercept = 0.05, col = "red", linetype = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust= 1, vjust= 0.5))+
  scale_y_continuous(labels = function(y) y + 1)

cowplot::save_plot(glue("Figures/Convergence/full_distrib_rhat_{model_name}.pdf"),
                   g, base_width = 20/cm(1), base_height = 12/cm(1))
cli_progress_done()
