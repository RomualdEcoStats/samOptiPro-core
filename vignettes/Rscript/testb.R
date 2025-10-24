# ---------------------------------------------------------------
# Script généré à partir du Rmd "Optimization model_2 Process, 'NUTS on alpha'"
# Auteur: Romuald H
# Date:  [date d'exécution]
# ---------------------------------------------------------------

# ---------------------------
# Chargement des packages et configuration globale
# ---------------------------
library(knitr)
library(coda)
library(ggrepel)
library(ggplot2)
library(magick)
library(MASS)  # Pour kde2d

knitr::opts_chunk$set(comment = NA, echo = FALSE, cache = TRUE, message = FALSE, warning = FALSE, error = FALSE, fig.width = 8, fig.height = 4, bg = "transparent", cache.lazy = FALSE)

# ---------------------------
# Step1: Model Configuration, Execution, and Diagnostics
# ---------------------------

# Step1-1: Scorff run model
source("Rscript/0_load_library.R")
source("Rscript/0_load_functions.R")

run_traceplot <- FALSE
run_WAIC <- TRUE
model_name <- "base_model"

modelfile <- "model/model_scorff_LCM_v1.R"
datafile <- "input_data/data_scorff_LCM_v2.rds"
Constfile <- "input_data/Const_scorff_LCM_v2.rds"
configfile <- "Rscript/1_configuration_base.R"
initfile <- "Rscript/0_generate_inits_base.R"
Data_nimble <- readRDS(datafile)
Const_nimble <- readRDS(Constfile)
ESS_threshold <- 1000  # For time to be reach to ESS =1000
n.iter <- 700e3        # total number of model iterations 
n.burnin <- n.iter * 0.3  # burnin defined in total number of model iterations
n.chains <- 3          # number of parallel chains
n.thin <- 250          # thin applied to keep a reasonable number of posteriors

burnin_output <- 250    # burnin defined relative to number of posterior used in outputs only
calculate_N <- function(n.iter, n.burnin, n.chains, n.thin) {
  N <- n.chains * floor((n.iter - n.burnin) / n.thin)
  return(N)
}

N <- calculate_N(n.iter, n.burnin, n.chains, n.thin)
print(N)  # Total number of posterior samples ie the number of iterations actually used for inference
np <- 0.1  # or 0.6/0.4/0.2/0.1 Proportion of nodes to be viewed
parallel_run <- TRUE

# Create log file
if (!dir.exists("log")) {
  dir.create("log")
}

# Names files to store results
source(configfile)  ## Check that all possible nodes are monitored before continuing.
source("Rscript/2_run_model.R")  # Run Model (will not run if name_file_mcmc exists)
source("Rscript/3_convergence.R")  # Convergence Analysis & plot
source("Rscript/4_parameters_summary_and_plot.R")  # Extract median & plot parameters
source("Rscript/5_fit_assessment.R")  # Fit assessment
if (run_WAIC) {
  source("Rscript/6_WAIC.R")
}
if (run_traceplot) {
  source("Rscript/7_traceplot.R")
}

# ---------------------------
# Step1-2: Check the convergence (Gelman−Rubin Rhat, neff) of model parameters
# ---------------------------
output_dir <- "output_images"
dir.create(output_dir, showWarnings = FALSE)

pdf_dir <- "~/Optimization_process/Toy_model/Scorff LCM_model2/Figures/Convergence"
#pdf_dir <- "C:/Users/hounyeme/Desktop/Scorff/Scorff LCM_model2/Scorff LCM_model2/Figures/Convergence"
pdf_files <- list.files(pdf_dir, pattern = "\\.pdf$", full.names = TRUE)

for (pdf in pdf_files) {
  img_list <- image_read(pdf)
  for (i in seq_along(img_list)) {
    img_path <- file.path(output_dir, paste0(basename(pdf), "_page", i, ".png"))
    image_write(img_list[i], path = img_path, format = "png")
  }
}

output_dir <- "output_images"
png_files <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
knitr::include_graphics(png_files)

## Once convergence has been verified, continue the analysis, otherwise resume running the model after making adjustments...

# ---------------------------
# Step2: Metrics calculation & visualisations
# ---------------------------

## Step2-1: Calculating ESS
results <- readRDS("C:/Users/hounyeme/Desktop/Scorff/Scorff LCM_model2/Scorff LCM_model2/outputs_hindcast/mcmc_results_base_model.rds")
niter_thinned <- nrow(results[[1]])

if (is.list(results) && length(dim(results[[1]])) == 2) {
  chains <- length(results)
  mcmc_list <- list()
  ess_full <- numeric()
  parameter_names <- colnames(results[[1]])
  
  for (parameter in parameter_names) {
    filtered_chains <- lapply(results, function(chain) {
      chain[(burnin_output + 1):niter_thinned, parameter, drop = FALSE]
    })
    parameter_mcmc_list <- lapply(filtered_chains, function(chain) {
      as.mcmc(as.matrix(chain))
    })
    mcmc_list[[parameter]] <- as.mcmc.list(parameter_mcmc_list)
    ess_full[parameter] <- effectiveSize(mcmc_list[[parameter]])
  }
  
  ess_values <- data.frame(
    Parameter = names(ess_full),
    ESS = as.numeric(ess_full),
    stringsAsFactors = FALSE
  )
  cat("Effective Sample Sizes (ESS) par nœud :\n")
  print(ess_values)
} else {
  stop("Erreur : 'results' n'est pas au format attendu pour le calcul de l'ESS.")
}

ess_values_filtered <- subset(ess_values, ESS > 0)
cat("ESS filtrés (ESS > 0) :\n")
print(ess_values_filtered)
removed_nodes <- setdiff(ess_values$Parameter, ess_values_filtered$Parameter)
cat("Nœuds supprimés (ESS = 0) :\n")
print(removed_nodes)

ess_values_filtered$Family <- sub("\\[.*", "", ess_values_filtered$Parameter)
ess_family <- ess_values_filtered %>%
  group_by(Family) %>%
  summarise(median_ESS = median(ESS)) %>%
  ungroup()
cat("Médiane ESS par famille :\n")
print(ess_family)
ess_quantile_5 <- quantile(ess_family$median_ESS, probs = 0.05)

p_bar <- ggplot(ess_family, aes(x = reorder(Family, median_ESS), y = median_ESS)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_text(aes(label = round(median_ESS, 1)), vjust = -0.5, size = 3) +
  geom_hline(yintercept = ess_quantile_5, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Quantile 5% Effective sample size",
       subtitle = paste("5th quantile =", round(ess_quantile_5, 1)),
       x = "Nodes",
       y = "Médian Effective Sample Size (ESS)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_bar)
ggsave("~/Optimization_process/Toy_model/Scorff LCM_model2/outputs_hindcast/ess_family_barplot.png",
       p_bar, width = 10, height = 6)

saveRDS(ess_family, "~/Optimization_process/Toy_model/Scorff LCM_model2/outputs_hindcast/sampler_results_by_family.rds")
cat("Résultats par famille sauvegardés dans '~/Optimization_process/Toy_model/Scorff LCM_model2/outputs_hindcast/sampler_results_by_family.rds'.\n")
saveRDS(ess_values_filtered, "~/Optimization_process/Toy_model/Scorff LCM_model2/outputs_hindcast/sampler_results.rds")
cat("Sampler results saved to '~/Optimization_process/Toy_model/Scorff LCM_model2/outputs_hindcast/sampler_results.rds'.\n")

## Step2-2: Remove nodes with NA and zeros values.
# ----------------------- Functions: Second Level Filtering --------------------
identify_nodes_with_na <- function(ess_values) {
  if (!is.data.frame(ess_values)) {
    stop("Error: `ess_values` must be a data.frame.")
  }
  nodes_with_na <- ess_values$Parameter[is.na(ess_values$ESS)]
  return(nodes_with_na)
}
cleaned_ess_values <- function(ess_values) {
  if (!is.data.frame(ess_values)) {
    stop("Error: `ess_values` must be a data.frame.")
  }
  nodes_with_na <- identify_nodes_with_na(ess_values)
  if (length(nodes_with_na) > 0) {
    cat("Nodes with NA ESS detected (to be removed):\n")
    print(nodes_with_na)
  } else {
    cat("No nodes with NA ESS detected.\n")
  }
  ess_values_cleaned <- ess_values[!is.na(ess_values$ESS), ]
  if (any(is.na(ess_values_cleaned$ESS))) {
    stop("Error: NA values persist after cleaning.")
  }
  return(ess_values_cleaned)
}
identify_nodes_with_zero <- function(ess_values) {
  if (!is.data.frame(ess_values)) {
    stop("Error: `ess_values` must be a data.frame.")
  }
  nodes_with_zero <- ess_values$Parameter[ess_values$ESS == 0]
  return(nodes_with_zero)
}
remove_zero_ess <- function(ess_values) {
  if (!is.data.frame(ess_values)) {
    stop("Error: `ess_values` must be a data.frame.")
  }
  nodes_with_zero <- identify_nodes_with_zero(ess_values)
  if (length(nodes_with_zero) > 0) {
    cat("Nodes with zero ESS detected (to be removed):\n")
    print(nodes_with_zero)
  } else {
    cat("No nodes with zero ESS detected.\n")
  }
  ess_values_cleaned_zero <- ess_values[ess_values$ESS != 0, ]
  if (any(ess_values_cleaned_zero$ESS == 0)) {
    stop("Error: Zero ESS values persist after cleaning.")
  }
  return(ess_values_cleaned_zero)
}
cat("Number of parameters before cleaning:", nrow(ess_values_filtered), "\n")
cat("Identifying nodes with NA ESS values...\n")
nodes_with_na <- identify_nodes_with_na(ess_values_filtered)
if (length(nodes_with_na) > 0) {
  cat("Nodes with NA ESS detected (to be removed):\n")
  print(nodes_with_na)
} else {
  cat("No nodes with NA ESS detected.\n")
}
cat("Cleaning ESS values by removing rows with NA ESS values...\n")
ess_no_na <- cleaned_ess_values(ess_values_filtered)
cat("Number of parameters after removing NA values:", nrow(ess_no_na), "\n")
cat("Identifying nodes with zero ESS values...\n")
nodes_with_zero <- identify_nodes_with_zero(ess_no_na)
if (length(nodes_with_zero) > 0) {
  cat("Nodes with zero ESS detected (to be removed):\n")
  print(nodes_with_zero)
} else {
  cat("No nodes with zero ESS detected.\n")
}
cat("Cleaning ESS values by removing rows with zero ESS values...\n")
ess_values_cleaned <- remove_zero_ess(ess_no_na)
cat("Cleaning complete. Remaining values:\n")
cat("Number of parameters after cleaning (NA and zeros removed):", nrow(ess_values_cleaned), "\n")
cat("Final cleaned ESS values (first few rows):\n")
print(head(ess_values_cleaned))
n_worst <- ceiling(length(ess_values_cleaned$Parameter) * np)  # Proportion of bottlenecks to be examined

## Step2-3: The three metrics after ESS
stopifnot("ess_values_cleaned must contain the 'Parameter' and 'ESS' columns" =
            all(c("Parameter", "ESS") %in% names(ess_values_cleaned)))
total_time <- as.numeric(time_mcmc)
Algorithmic_efficiency_vec <- ess_values_cleaned$ESS / N
Computational_efficiency_tot_vec <- 1 / ((ESS_threshold * total_time) / ess_values_cleaned$ESS)
Algorithmic_efficiency <- data.frame(
  Parameter = ess_values_cleaned$Parameter,
  Value = Algorithmic_efficiency_vec,
  stringsAsFactors = FALSE
)
Computational_efficiency_tot <- data.frame(
  Parameter = ess_values_cleaned$Parameter,
  Value = Computational_efficiency_tot_vec,
  stringsAsFactors = FALSE
)
bottleneck_nodes_by_Algorithmic_efficiency <- order(Algorithmic_efficiency$Value, decreasing = TRUE, na.last = NA)
bottleneck_nodes_by_Computational_efficiency_tot <- order(Computational_efficiency_tot$Value, decreasing = FALSE, na.last = NA)
all_param_names <- ess_values_cleaned$Parameter
valid_bottleneck_nodes_by_Algorithmic_efficiency <- bottleneck_nodes_by_Algorithmic_efficiency[
  bottleneck_nodes_by_Algorithmic_efficiency <= length(all_param_names)
]
valid_bottleneck_nodes_by_Computational_efficiency_tot <- bottleneck_nodes_by_Computational_efficiency_tot[
  bottleneck_nodes_by_Computational_efficiency_tot <= length(all_param_names)
]
Algorithmic_efficiency_bottleneck_node_names <- all_param_names[valid_bottleneck_nodes_by_Algorithmic_efficiency]
Computational_efficiency_bottleneck_node_names_tot <- all_param_names[valid_bottleneck_nodes_by_Computational_efficiency_tot]
worst_algorithmic_efficiency_nodes <- tail(Algorithmic_efficiency_bottleneck_node_names, n_worst)
worst_computational_efficiency_nodes_tot <- head(Computational_efficiency_bottleneck_node_names_tot, n_worst)
output_directory <- "~/Optimization_process/Toy_model/Scorff LCM_model2/outputs_hindcast"
algorithmic_efficiency_file <- file.path(output_directory, "Algorithmic_efficiency.rds")
computational_efficiency_tot_file <- file.path(output_directory, "Computational_efficiency_tot.rds")
saveRDS(Algorithmic_efficiency, file = algorithmic_efficiency_file)
saveRDS(Computational_efficiency_tot, file = computational_efficiency_tot_file)
cat("Efficiencies have been saved to:", output_directory, "\n")
combined_data <- data.frame(
  Algorithmic_Efficiency = Algorithmic_efficiency$Value,
  Computational_Efficiency_tot = Computational_efficiency_tot$Value,
  AssociatedNodes = ess_values_cleaned$Parameter,
  stringsAsFactors = FALSE
)
combined_data <- combined_data %>%
  mutate(Family = sub("\\[.*\\]", "", AssociatedNodes))
grouped_data <- combined_data %>%
  group_by(Family) %>%
  summarise(
    MedianAlgorithmicEfficiency = median(Algorithmic_Efficiency, na.rm = TRUE),
    MedianComputationalEfficiency_tot = median(Computational_Efficiency_tot, na.rm = TRUE),
    FamilySize = n()
  )
algorithmic_plot_raw <- ggplot(grouped_data, aes(x = reorder(Family, MedianAlgorithmicEfficiency),
                                                 y = MedianAlgorithmicEfficiency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Median Algorithmic Efficiency by Node Family",
       x = "Node Family", y = "Median Algorithmic Efficiency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
computational_plot_raw_tot <- ggplot(grouped_data, aes(x = reorder(Family, MedianComputationalEfficiency_tot),
                                                       y = MedianComputationalEfficiency_tot)) +
  geom_bar(stat = "identity", fill = "green") +
  labs(title = "Median Computational Efficiency (Total Time) by Node Family",
       x = "Node Family", y = "Median Computational Efficiency (Total Time)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(algorithmic_plot_raw)
print(computational_plot_raw_tot)
worst_nodes_raw_algorithmic <- grouped_data %>%
  arrange(MedianAlgorithmicEfficiency) %>%
  slice_head(n = n_worst)
worst_nodes_raw_computational_tot <- grouped_data %>%
  arrange(MedianComputationalEfficiency_tot) %>%
  slice_head(n = n_worst)
cat("Node families with the lowest median algorithmic efficiencies:\n")
print(worst_nodes_raw_algorithmic)
cat("Node families with the lowest median total computational efficiencies:\n")
print(worst_nodes_raw_computational_tot)
all_node_names <- setdiff(model.nimble$getNodeNames(), removed_nodes)
num_dependencies <- sapply(all_node_names, function(node) {
  deps <- model.nimble$getDependencies(nodes = node, self = FALSE)
  deps_clean <- deps[!grepl("^lifted_", deps)]
  length(deps_clean)
})
filtered_node_names <- names(num_dependencies)[num_dependencies > 0]
node_dependencies <- lapply(filtered_node_names, function(node) {
  deps <- model.nimble$getDependencies(nodes = node, self = FALSE)
  deps_clean <- deps[!grepl("^lifted_", deps)]
  deps_clean
})
names(node_dependencies) <- filtered_node_names
df_dependencies <- do.call(rbind, lapply(names(node_dependencies), function(node) {
  deps <- node_dependencies[[node]]
  if (length(deps) == 0) {
    data.frame(
      node = node,
      dependency = NA,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      node = node,
      dependency = deps,
      stringsAsFactors = FALSE
    )
  }
}))
head(df_dependencies)
# Optionnel: write.csv(df_dependencies, file = file.path(output_directory, "dependencies_per_node.csv"), row.names = FALSE)
df <- data.frame(
  nodename = filtered_node_names,
  n = as.numeric(num_dependencies[filtered_node_names]),
  stringsAsFactors = FALSE
)
df <- df %>%
  mutate(
    type = ifelse(grepl("^lifted_", nodename), "Lifted", "Regular"),
    family = sub("\\[.*\\]", "", nodename)
  )
df_summary <- df %>%
  group_by(type, family) %>%
  summarise(
    median_n = median(n),
    count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(median_n))
plot_regular <- ggplot(df_summary %>% filter(type == "Regular"), 
                       aes(x = reorder(family, median_n), y = median_n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = median_n), hjust = -0.2, size = 3.5) +
  coord_flip() +
  theme_minimal(base_size = 6) +
  labs(title = "Median Dependency by Family (Regular Nodes)",
       x = "Node Family", y = "Median Dependencies") +
  theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6, face = "bold"),
        axis.title.y = element_text(size = 6, face = "bold"),
        plot.margin = margin(4, 4, 4, 4))
print(plot_regular)
# Optionnel: ggsave("plot_regular.png", plot = plot_regular, width = 12, height = 8, dpi = 300)
node_of_interest <- "alpha"
dependencies_of_interest <- node_dependencies[[node_of_interest]]
if (length(dependencies_of_interest) > 0) {
  print(paste("Dependencies for node:", node_of_interest))
  print(dependencies_of_interest)
} else {
  print(paste("No dependencies for node:", node_of_interest))
}

# ---------------------------
# Step2-3-1: List the target nodes for each sampler
# ---------------------------
samplers <- conf.mcmc.model$getSamplers()
lapply(samplers, function(sampler) {
  list(
    Type = class(sampler),
    TargetNodes = sampler$target
  )
})
n.samplers <- length(samplers)
prop_worst <- ceiling(n.samplers * np)

# ---------------------------
# Step2-3-2: Extraction of nodes from samplers - Filtering of ESS values
# ---------------------------
sampler_runtime <- readRDS("C:/Users/hounyeme/Desktop/Scorff/Scorff LCM_model2/Scorff LCM_model2/outputs_hindcast/sampler_times_base_model_chain1.rds")
if (!exists("ess_values_cleaned")) {
  stop("Error: The cleaned ESS values (ess_values_cleaned) are not available in the workspace.")
}
if (!all(c("Parameter", "ESS") %in% names(ess_values_cleaned))) {
  stop("Error: 'ess_values_cleaned' must contain the columns 'Parameter' and 'ESS'.")
}
unique_sampler_nodes <- unique(unlist(lapply(samplers, function(s) {
  s$target
})))
filtered_ess_df <- ess_values_cleaned[ess_values_cleaned$Parameter %in% unique_sampler_nodes, ]
missing_nodes <- setdiff(unique_sampler_nodes, filtered_ess_df$Parameter)
ordered_sampler_nodes <- sort(unique_sampler_nodes[!(unique_sampler_nodes %in% missing_nodes)])
filtered_ess_df <- filtered_ess_df[order(filtered_ess_df$Parameter), ]
aligned_ess_df <- filtered_ess_df[match(ordered_sampler_nodes, filtered_ess_df$Parameter), ]
cat("\nFiltered and Aligned ESS values (aligned_ess_df):\n")
cat("\nSampler target nodes missing in ESS data:\n")
if (length(missing_nodes) > 0) {
  print(missing_nodes)
} else {
  cat("None\n")
}
combined_sampler_results_df <- data.frame(
  SamplerNode = ordered_sampler_nodes,
  SamplerRuntime = sampler_runtime,
  ESS = aligned_ess_df$ESS,
  stringsAsFactors = FALSE
)
max_length <- max(
  length(combined_sampler_results_df$SamplerNode),
  length(combined_sampler_results_df$SamplerRuntime),
  length(combined_sampler_results_df$ESS)
)
combined_sampler_results_df <- data.frame(
  SamplerNode = c(combined_sampler_results_df$SamplerNode, rep(NA, max_length - length(combined_sampler_results_df$SamplerNode))),
  SamplerRuntime = c(combined_sampler_results_df$SamplerRuntime, rep(NA, max_length - length(combined_sampler_results_df$SamplerRuntime))),
  ESS = c(combined_sampler_results_df$ESS, rep(NA, max_length - length(combined_sampler_results_df$ESS))),
  stringsAsFactors = FALSE
)
saveRDS(combined_sampler_results_df, "~/Optimization_process/Toy_model/Scorff LCM_model2/outputs_hindcast/sampler_results.rds")
cat("Sampler results saved to '~/Optimization_process/Toy_model/Scorff LCM_model2/outputs_hindcast/sampler_results.rds'.\n")

# ---------------------------
# Step3: Calculating node-specific metrics with samplers and visualisations
# ---------------------------
## Step3-1: Relationship between Rhat and ESS
if (!exists("mcmc_list") || !exists("ess_values_cleaned")) {
  stop("The 'mcmc_list' and/or 'ess_values_cleaned' objects are missing.")
}
cleaned_ESS <- ess_values_cleaned
mcmc_chains <- mcmc_list
node_names <- cleaned_ESS$Parameter
matching_nodes <- intersect(names(mcmc_chains), node_names)
if (length(matching_nodes) == 0) {
  stop("No matching node found between 'mcmc_chains' and 'cleaned_ESS'.")
}
filtered_mcmc <- mcmc_chains[matching_nodes]
mcmc_matrix <- do.call(cbind, lapply(filtered_mcmc, function(node) {
  do.call(cbind, lapply(node, function(chain) {
    chain_matrix <- as.matrix(chain)
    colnames(chain_matrix) <- attr(chain, "dimnames")[[2]]
    chain_matrix
  }))
}))
cat("Column names in mcmc_matrix:\n")
mcmc_df <- as.data.frame(mcmc_matrix)
cat("Overview of mcmc_df:\n")
rhat_list <- lapply(filtered_mcmc, function(node) {
  gelman.diag(as.mcmc.list(node), multivariate = FALSE)$psrf[, 1]
})
ess_list <- lapply(filtered_mcmc, function(node) {
  effectiveSize(as.mcmc.list(node))
})
family_names <- sub("\\[.*\\]", "", matching_nodes)
diagnostics_df <- data.frame(
  Node = matching_nodes,
  Family = family_names,
  ESS = unlist(ess_list),
  Rhat = unlist(rhat_list),
  stringsAsFactors = FALSE
)
aggregated_df <- diagnostics_df %>%
  group_by(Family) %>%
  summarise(
    MedianESS = median(ESS, na.rm = TRUE),
    MedianRhat = median(Rhat, na.rm = TRUE),
    FamilySize = n()
  )
cat("Aggregated data by family:\n")
rhat_ess_plot <- ggplot(aggregated_df, aes(x = MedianESS, y = MedianRhat)) +
  geom_point(size = 3, alpha = 0.7, color = "blue") +
  geom_text_repel(aes(label = Family), size = 3, max.overlaps = 10) +
  labs(
    title = "Rhat vs ESS by Node Family",
    x = "Median Effective Sample Size (ESS)",
    y = "Median Rhat"
  ) +
  theme_minimal()
print(rhat_ess_plot)

f.density.bivar <- function(x, y, nlevels, nb.points) {
  indices <- which(x <= quantile(x, prob = 0.995))
  x <- x[indices]
  y <- y[indices]
  xrange <- range(x); nbreaks.x <- 100
  yrange <- range(y); nbreaks.y <- 100
  xhist <- hist(x, breaks = seq(xrange[1], xrange[2], length.out = nbreaks.x), plot = FALSE)
  yhist <- hist(y, breaks = seq(yrange[1], yrange[2], length.out = nbreaks.y), plot = FALSE)
  nf <- layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), c(3, 1), c(1, 3), TRUE)
  layout.show(nf)
  par(mar = c(5, 5, 1, 1))
  plot(x[1:nb.points], y[1:nb.points], xlim = xrange, ylim = yrange,
       xlab = "", ylab = "", pch = ".", cex.lab = 1.5)
  dens2d <- kde2d(x = x, y = y, n = 100)
  contour(dens2d, nlevels = nlevels, drawlabels = FALSE, col = "red", lwd = 2, add = TRUE)
  mtext(text = expression(par1), side = 1, line = 3, cex = 1.3)
  mtext(text = expression(par2), side = 2, line = 3, cex = 1.3)
  par(mar = c(0, 5, 3, 1))
  barplot(xhist$density, axes = FALSE, space = 0, horiz = FALSE, col = "lightblue")
  mtext(text = expression(paste("Marginal pdf for ", par1)), side = 3, line = 1, cex = 1.3)
  par(mar = c(5, 0, 1, 3))
  barplot(yhist$density, axes = FALSE, space = 0, horiz = TRUE, col = "lightblue")
  mtext(text = expression(paste("Marginal pdf for ", par2)), side = 4, line = 1, las = 0, cex = 1.3)
}
panel.dens <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- density(x, adjust = 1)
  h$y <- h$y / max(h$y)
  lines(h, col = "black", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r_val <- cor(x, y)
  txt <- format(c(r_val, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex <- 0.8 / strwidth(txt)
  cex <- 2
  text(0.5, 0.5, txt, cex = cex)
}

if (!all(c("alpha", "k") %in% colnames(mcmc_df))) {
  stop("The 'alpha' and 'k' columns are not present in mcmc_df.")
}
alpha_vals <- mcmc_df$alpha
k_vals <- mcmc_df$k
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(density(alpha_vals), xlab = "alpha", ylab = "Density", main = "Posterior Density of alpha", col = "red", lwd = 2)
plot(density(k_vals), xlab = "k", ylab = "Density", main = "Posterior Density of k", col = "blue", lwd = 2)
par(mfrow = c(1, 1))
plot(k_vals, alpha_vals, xlab = "k", ylab = "alpha", main = "Joint Posterior Distribution of alpha and k", col = rgb(0, 0, 1, 0.5), pch = 16)
abline(lm(alpha_vals ~ k_vals), col = "red", lwd = 2)
pairs(cbind(alpha_vals, k_vals),
      labels = c("alpha", "k"), 
      lower.panel = panel.smooth,
      diag.panel = panel.dens, 
      upper.panel = panel.cor,
      cex.labels = 1.5, font.labels = 2)
kde_res <- kde2d(k_vals, alpha_vals, n = 200)
contour(kde_res, xlab = "k", ylab = "alpha", main = "Contour Plot of Joint Density")
points(k_vals, alpha_vals, pch = 16, col = rgb(0, 0, 1, 0.5))
plot(k_vals, alpha_vals, xlab = "k", ylab = "alpha", main = "Uncertainty in Regression Lines", col = "blue", pch = 16, cex = 0.6)
nb_lines <- 100
for (i in 1:nb_lines) {
  abline(a = alpha_vals[i], b = k_vals[i], col = rgb(0.7, 0.7, 0.7, 0.2), lwd = 0.5)
}
abline(median(alpha_vals), median(k_vals), col = "red", lwd = 2)
data_zoom <- data.frame(k = k_vals, alpha = alpha_vals)
data_zoom <- subset(data_zoom, k <= 130000 & alpha <= 0.05)
kde_zoom <- kde2d(data_zoom$k, data_zoom$alpha, n = 200)
ggplot(data_zoom, aes(x = k, y = alpha)) +
  geom_point(color = "blue", size = 0.5, alpha = 0.6) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.3) +
  scale_fill_viridis_c() +
  xlim(0.05, 130000) + ylim(0, 0.05) +
  labs(title = "Zoomed Scatter Plot with Density Contours", x = "k (zoomed)", y = "alpha (zoomed)") +
  theme_minimal()

# ---------------------------
# Step3-2: Sampler metrics calculation & visualisations
# ---------------------------
stopifnot("ess_values_cleaned must contain the 'Parameter' and 'ESS' columns" =
            all(c("Parameter", "ESS") %in% names(ess_values_cleaned)))
eff_min <- min(ess_values_cleaned$ESS, na.rm = TRUE)
Comp_Eff <- 1 / ((ESS_threshold * sampler_times) / ess_values_cleaned$ESS)
Computational_efficiency_df <- data.frame(
  Parameter = ess_values_cleaned$Parameter,
  value = Comp_Eff,
  stringsAsFactors = FALSE
)
bottleneck_indices <- order(Computational_efficiency_df$value, decreasing = FALSE, na.last = NA)
all_param_names <- ess_values_cleaned$Parameter
valid_bottleneck_nodes <- bottleneck_indices[bottleneck_indices <= length(all_param_names)]
CompEff_bottleneck_node_names <- all_param_names[valid_bottleneck_nodes]
worst_comp_eff_nodes <- head(CompEff_bottleneck_node_names, prop_worst)
output_directory <- "~/Optimization_process/Toy_model/Scorff LCM_model2/outputs_hindcast"
comp_efficiency_file <- file.path(output_directory, "Computational_efficiency.rds")
saveRDS(Computational_efficiency_df, file = comp_efficiency_file)
cat("Computational efficiency results have been saved to:", output_directory, "\n")
combined_data <- data.frame(
  Computational_Efficiency = Computational_efficiency_df$value,
  AssociatedNodes = Computational_efficiency_df$Parameter,
  stringsAsFactors = FALSE
)
combined_data <- combined_data %>%
  mutate(Family = sub("\\[.*\\]", "", AssociatedNodes))
grouped_data <- combined_data %>%
  group_by(Family) %>%
  summarise(
    MedianComputationalEfficiency = median(Computational_Efficiency, na.rm = TRUE),
    FamilySize = n()
  )
comp_eff_plot <- ggplot(grouped_data, aes(x = reorder(Family, MedianComputationalEfficiency),
                                          y = MedianComputationalEfficiency)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  labs(title = "Median Computational Efficiency by Node Family",
       x = "Node Family", y = "Median Computational Efficiency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(comp_eff_plot)
worst_eff_families <- grouped_data %>%
  arrange(MedianComputationalEfficiency) %>%
  slice_head(n = n_worst)
cat("Families with the lowest median computational efficiency:\n")
print(worst_eff_families)
sampler_info_list <- lapply(samplers, function(s) {
  data.frame(
    SamplerType = class(s),
    AssociatedNodes = paste(s$target, collapse = ", "),
    stringsAsFactors = FALSE
  )
})
sampler_info_df <- do.call(rbind, sampler_info_list)
colnames(sampler_info_df) <- c("SamplerType", "AssociatedNodes")
saveRDS(sampler_info_df, "sampler_info.rds")
sampler_times <- readRDS("C:/Users/hounyeme/Desktop/Scorff/Scorff LCM_model2/Scorff LCM_model2/outputs_hindcast/sampler_times_base_model_chain1.rds")
if (length(sampler_times) != nrow(sampler_info_df)) {
  stop("The number of sampler times does not match the number of rows in sampler_info_df.")
}
sampler_info_df$Time <- sampler_times
sampler_info_file <- file.path(output_directory, "combined_sampler_info.rds")
saveRDS(sampler_info_df, file = sampler_info_file)
cat("Sampler information saved to", sampler_info_file, "\n")
combined_sampler_data <- readRDS(sampler_info_file)
if (exists("sampler_nodes_filtered")) {
  combined_sampler_data <- combined_sampler_data %>%
    filter(AssociatedNodes %in% sampler_nodes_filtered)
} else {
  warning("sampler_nodes_filtered is not defined; using all sampler information.")
}
combined_sampler_data <- combined_sampler_data %>%
  mutate(Family = sub("\\[.*\\]", "", AssociatedNodes))
grouped_sampler_data <- combined_sampler_data %>%
  group_by(Family) %>%
  summarise(
    MedianTime = median(Time, na.rm = TRUE),
    FamilySize = n()
  )
time_plot <- ggplot(grouped_sampler_data, aes(x = reorder(Family, MedianTime),
                                              y = MedianTime)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Median Sampler Times by Node Family",
       x = "Node Family", y = "Median Time (seconds)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(time_plot)
slowest_families <- grouped_sampler_data %>%
  filter(!is.na(MedianTime)) %>%
  arrange(desc(MedianTime)) %>%
  slice_max(MedianTime, n = 5)
cat("Families with the highest raw median sampler times (slowest families):\n")
print(slowest_families)

# ---------------------------
# Step3: Calculating node-specific metrics with samplers and visualisations
# ---------------------------
if (!exists("mcmc_list") || !exists("ess_values_cleaned")) {
  stop("The 'mcmc_list' and/or 'ess_values_cleaned' objects are missing.")
}
cleaned_ESS <- ess_values_cleaned
mcmc_chains <- mcmc_list
node_names <- cleaned_ESS$Parameter
matching_nodes <- intersect(names(mcmc_chains), node_names)
if (length(matching_nodes) == 0) {
  stop("No matching node found between 'mcmc_chains' and 'cleaned_ESS'.")
}
filtered_mcmc <- mcmc_chains[matching_nodes]
mcmc_matrix <- do.call(cbind, lapply(filtered_mcmc, function(node) {
  do.call(cbind, lapply(node, function(chain) {
    chain_matrix <- as.matrix(chain)
    colnames(chain_matrix) <- attr(chain, "dimnames")[[2]]
    chain_matrix
  }))
}))
cat("Column names in mcmc_matrix:\n")
mcmc_df <- as.data.frame(mcmc_matrix)
cat("Overview of mcmc_df:\n")
rhat_list <- lapply(filtered_mcmc, function(node) {
  gelman.diag(as.mcmc.list(node), multivariate = FALSE)$psrf[, 1]
})
ess_list <- lapply(filtered_mcmc, function(node) {
  effectiveSize(as.mcmc.list(node))
})
family_names <- sub("\\[.*\\]", "", matching_nodes)
diagnostics_df <- data.frame(
  Node = matching_nodes,
  Family = family_names,
  ESS = unlist(ess_list),
  Rhat = unlist(rhat_list),
  stringsAsFactors = FALSE
)
aggregated_df <- diagnostics_df %>%
  group_by(Family) %>%
  summarise(
    MedianESS = median(ESS, na.rm = TRUE),
    MedianRhat = median(Rhat, na.rm = TRUE),
    FamilySize = n()
  )
cat("Aggregated data by family:\n")
rhat_ess_plot <- ggplot(aggregated_df, aes(x = MedianESS, y = MedianRhat)) +
  geom_point(size = 3, alpha = 0.7, color = "blue") +
  geom_text_repel(aes(label = Family), size = 3, max.overlaps = 10) +
  labs(
    title = "Rhat vs ESS by Node Family",
    x = "Median Effective Sample Size (ESS)",
    y = "Median Rhat"
  ) +
  theme_minimal()
print(rhat_ess_plot)

f.density.bivar <- function(x, y, nlevels, nb.points) {
  indices <- which(x <= quantile(x, prob = 0.995))
  x <- x[indices]
  y <- y[indices]
  xrange <- range(x); nbreaks.x <- 100
  yrange <- range(y); nbreaks.y <- 100
  xhist <- hist(x, breaks = seq(xrange[1], xrange[2], length.out = nbreaks.x), plot = FALSE)
  yhist <- hist(y, breaks = seq(yrange[1], yrange[2], length.out = nbreaks.y), plot = FALSE)
  nf <- layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), c(3, 1), c(1, 3), TRUE)
  layout.show(nf)
  par(mar = c(5, 5, 1, 1))
  plot(x[1:nb.points], y[1:nb.points], xlim = xrange, ylim = yrange,
       xlab = "", ylab = "", pch = ".", cex.lab = 1.5)
  dens2d <- kde2d(x = x, y = y, n = 100)
  contour(dens2d, nlevels = nlevels, drawlabels = FALSE, col = "red", lwd = 2, add = TRUE)
  mtext(text = expression(par1), side = 1, line = 3, cex = 1.3)
  mtext(text = expression(par2), side = 2, line = 3, cex = 1.3)
  par(mar = c(0, 5, 3, 1))
  barplot(xhist$density, axes = FALSE, space = 0, horiz = FALSE, col = "lightblue")
  mtext(text = expression(paste("Marginal pdf for ", par1)), side = 3, line = 1, cex = 1.3)
  par(mar = c(5, 0, 1, 3))
  barplot(yhist$density, axes = FALSE, space = 0, horiz = TRUE, col = "lightblue")
  mtext(text = expression(paste("Marginal pdf for ", par2)), side = 4, line = 1, las = 0, cex = 1.3)
}
panel.dens <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- density(x, adjust = 1)
  h$y <- h$y / max(h$y)
  lines(h, col = "black", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r_val <- cor(x, y)
  txt <- format(c(r_val, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex <- 0.8 / strwidth(txt)
  cex <- 2
  text(0.5, 0.5, txt, cex = cex)
}
if (!all(c("alpha", "k") %in% colnames(mcmc_df))) {
  stop("The 'alpha' and 'k' columns are not present in mcmc_df.")
}
alpha_vals <- mcmc_df$alpha
k_vals <- mcmc_df$k
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(density(alpha_vals), xlab = "alpha", ylab = "Density", main = "Posterior Density of alpha", col = "red", lwd = 2)
plot(density(k_vals), xlab = "k", ylab = "Density", main = "Posterior Density of k", col = "blue", lwd = 2)
par(mfrow = c(1, 1))
plot(k_vals, alpha_vals, xlab = "k", ylab = "alpha", main = "Joint Posterior Distribution of alpha and k", col = rgb(0, 0, 1, 0.5), pch = 16)
abline(lm(alpha_vals ~ k_vals), col = "red", lwd = 2)
pairs(cbind(alpha_vals, k_vals),
      labels = c("alpha", "k"), 
      lower.panel = panel.smooth,
      diag.panel = panel.dens, 
      upper.panel = panel.cor,
      cex.labels = 1.5, font.labels = 2)
kde_res <- kde2d(k_vals, alpha_vals, n = 200)
contour(kde_res, xlab = "k", ylab = "alpha", main = "Contour Plot of Joint Density")
points(k_vals, alpha_vals, pch = 16, col = rgb(0, 0, 1, 0.5))
plot(k_vals, alpha_vals, xlab = "k", ylab = "alpha", main = "Uncertainty in Regression Lines", col = "blue", pch = 16, cex = 0.6)
nb_lines <- 100
for (i in 1:nb_lines) {
  abline(a = alpha_vals[i], b = k_vals[i], col = rgb(0.7, 0.7, 0.7, 0.2), lwd = 0.5)
}
abline(median(alpha_vals), median(k_vals), col = "red", lwd = 2)
data_zoom <- data.frame(k = k_vals, alpha = alpha_vals)
data_zoom <- subset(data_zoom, k <= 130000 & alpha <= 0.05)
kde_zoom <- kde2d(data_zoom$k, data_zoom$alpha, n = 200)
ggplot(data_zoom, aes(x = k, y = alpha)) +
  geom_point(color = "blue", size = 0.5, alpha = 0.6) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.3) +
  scale_fill_viridis_c() +
  xlim(0.05, 130000) + ylim(0, 0.05) +
  labs(title = "Zoomed Scatter Plot with Density Contours", x = "k (zoomed)", y = "alpha (zoomed)") +
  theme_minimal()

# ---------------------------
# Step3-2: Sampler metrics calculation & visualisations
# ---------------------------
stopifnot("ess_values_cleaned must contain the 'Parameter' and 'ESS' columns" =
            all(c("Parameter", "ESS") %in% names(ess_values_cleaned)))
eff_min <- min(ess_values_cleaned$ESS, na.rm = TRUE)
Comp_Eff <- 1 / ((ESS_threshold * sampler_times) / ess_values_cleaned$ESS)
Computational_efficiency_df <- data.frame(
  Parameter = ess_values_cleaned$Parameter,
  value = Comp_Eff,
  stringsAsFactors = FALSE
)
bottleneck_indices <- order(Computational_efficiency_df$value, decreasing = FALSE, na.last = NA)
all_param_names <- ess_values_cleaned$Parameter
valid_bottleneck_indices <- bottleneck_indices[bottleneck_indices <= length(all_param_names)]
CompEff_bottleneck_node_names <- all_param_names[valid_bottleneck_indices]
worst_comp_eff_nodes <- head(CompEff_bottleneck_node_names, prop_worst)
output_directory <- "~/Optimization_process/Toy_model/Scorff LCM_model2/outputs_hindcast"
comp_efficiency_file <- file.path(output_directory, "Computational_efficiency.rds")
saveRDS(Computational_efficiency_df, file = comp_efficiency_file)
cat("Computational efficiency results have been saved to:", output_directory, "\n")
combined_data <- data.frame(
  Computational_Efficiency = Computational_efficiency_df$value,
  AssociatedNodes = Computational_efficiency_df$Parameter,
  stringsAsFactors = FALSE
)
combined_data <- combined_data %>%
  mutate(Family = sub("\\[.*\\]", "", AssociatedNodes))
grouped_data <- combined_data %>%
  group_by(Family) %>%
  summarise(
    MedianComputationalEfficiency = median(Computational_Efficiency, na.rm = TRUE),
    FamilySize = n()
  )
comp_eff_plot <- ggplot(grouped_data, aes(x = reorder(Family, MedianComputationalEfficiency),
                                          y = MedianComputationalEfficiency)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  labs(title = "Median Computational Efficiency by Node Family", x = "Node Family", y = "Median Computational Efficiency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(comp_eff_plot)
worst_eff_families <- grouped_data %>%
  arrange(MedianComputationalEfficiency) %>%
  slice_head(n = n_worst)
cat("Families with the lowest median computational efficiency:\n")
print(worst_eff_families)
sampler_info_list <- lapply(samplers, function(s) {
  data.frame(
    SamplerType = class(s),
    AssociatedNodes = paste(s$target, collapse = ", "),
    stringsAsFactors = FALSE
  )
})
sampler_info_df <- do.call(rbind, sampler_info_list)
colnames(sampler_info_df) <- c("SamplerType", "AssociatedNodes")
saveRDS(sampler_info_df, "sampler_info.rds")
sampler_times <- readRDS("C:/Users/hounyeme/Desktop/Scorff/Scorff LCM_model2/Scorff LCM_model2/outputs_hindcast/sampler_times_base_model_chain1.rds")
if (length(sampler_times) != nrow(sampler_info_df)) {
  stop("The number of sampler times does not match the number of rows in sampler_info_df.")
}
sampler_info_df$Time <- sampler_times
sampler_info_file <- file.path(output_directory, "combined_sampler_info.rds")
saveRDS(sampler_info_df, file = sampler_info_file)
cat("Sampler information saved to", sampler_info_file, "\n")
combined_sampler_data <- readRDS(sampler_info_file)
if (exists("sampler_nodes_filtered")) {
  combined_sampler_data <- combined_sampler_data %>%
    filter(AssociatedNodes %in% sampler_nodes_filtered)
} else {
  warning("sampler_nodes_filtered is not defined; using all sampler information.")
}
combined_sampler_data <- combined_sampler_data %>%
  mutate(Family = sub("\\[.*\\]", "", AssociatedNodes))
grouped_sampler_data <- combined_sampler_data %>%
  group_by(Family) %>%
  summarise(
    MedianTime = median(Time, na.rm = TRUE),
    FamilySize = n()
  )
time_plot <- ggplot(grouped_sampler_data, aes(x = reorder(Family, MedianTime),
                                              y = MedianTime)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Median Sampler Times by Node Family", x = "Node Family", y = "Median Time (seconds)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(time_plot)
slowest_families <- grouped_sampler_data %>%
  filter(!is.na(MedianTime)) %>%
  arrange(desc(MedianTime)) %>%
  slice_max(MedianTime, n = 5)
cat("Families with the highest raw median sampler times (slowest families):\n")
print(slowest_families)

# ---------------------------
# Step3: Calculating node-specific metrics with samplers and visualisations (fin)
# ---------------------------

# End of script.
