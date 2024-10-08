library(here)
library(dplyr)
library(tidyr)
library(forcats)
library(readr)
library(ggplot2)
library(rlang)

# Read Data
data_path <- here("Sim_Data", "CorrectedSE_10012024.rds")
df <- readRDS(data_path)

# Helper Function
process_data <- function(data) {
  results <- data %>%
    gather("var", "val", rb.joint.est:rmse_corrected.rel.est_corrected) %>%
    dplyr::select(-c(REPLICATIONS:WARNINGS)) %>%
    separate(col = var, into = c("stats", "method"), sep = "\\.") %>%
    spread(stats, val) %>%
    mutate(N = N_per_p*n_items,
           N_lab = as_factor(paste0("italic(N) == ", N_per_p*n_items)),
           path_lab = as_factor(paste0("\\beta_{xm} == ", path)),
           rel_lab = as_factor(paste0("rho == ", rel)))
  return(results)
}

# Plot function
plot_results <- function(data,
                         x_var, y_var, 
                         x_label, y_label, 
                         y_limits = NULL, 
                         color_var = "method", 
                         shape_var = "method", 
                         point_size = 3, 
                         add_lines = FALSE, 
                         remove_baseline = FALSE) {
  browser()
  # Optionally filter out "mmr" group
  if (remove_baseline) {
    data <- data %>% filter(!!sym(color_var) != "mmr")
  }
  
  plot <- data %>%
    ggplot(aes(x = factor(!!sym(x_var)), y = !!sym(y_var), color = !!sym(color_var), shape = !!sym(shape_var))) +
    geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0), size = point_size, alpha = 0.8) +
    facet_grid(path_lab ~ rel_lab, labeller = label_parsed) +
    scale_shape_manual(name = "Method",
                       values = if (remove_baseline) c("rapi" = 16, "upi" = 18, "tspa" = 15) else c("mmr" = 17, "rapi" = 16, "upi" = 18, "tspa" = 15),
                       labels = if (remove_baseline) c("rapi" = "RAPI", "upi" = "Matched-Pair UPI", "tspa" = "2S-PA-Int") else c("mmr" = "MMR", "rapi" = "RAPI", "upi" = "Matched-Pair UPI", "tspa" = "2S-PA-Int")) +
    scale_color_manual(name = "Method",
                       values = if (remove_baseline) c("rapi" = "#E69F00", "upi" = "#CC79A7", "tspa" = "#009E73") else c("mmr" = "pink", "rapi" = "#CC79A7", "upi" = "#d95f02", "tspa" = "#009E73"),
                       labels = if (remove_baseline) c("rapi" = "RAPI", "upi" = "Matched-Pair UPI", "tspa" = "2S-PA-Int") else c("mmr" = "MMR", "rapi" = "RAPI", "upi" = "Matched-Pair UPI", "tspa" = "2S-PA-Int")) +
    labs(x = x_label, y = y_label, color = "method", shape = "method") +
    theme_minimal() +
    theme(legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          legend.position = "right",
          strip.text = element_text(size = 12, face = "bold"),  # Adjust facet label text
          panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add black border to panels
          strip.background = element_rect(color = "black", size = 1),  # Add black border to facet labels
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray80"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),  # Remove vertical grid lines
          panel.grid.major.y = element_line(size = 0.5, linetype = 'solid', colour = "gray80"),  # Keep horizontal lines
          axis.line = element_line(color = "black"),  # Add black lines to x and y axes
          axis.ticks = element_line(color = "black"))
  
  if (add_lines) {
    plot <- plot + 
      geom_line(aes(group = !!sym(color_var), 
                    linetype = !!sym(shape_var)), 
                size = 0.7, alpha = 0.6) +  # Removed `position_dodge()`
      scale_linetype_manual(name = "Method",
                            values = if (remove_baseline) c("rapi" = "solid", "upi" = "dashed", "tspa" = "twodash") else c("mmr" = "dotted", "rapi" = "solid", "upi" = "dashed", "tspa" = "twodash"),
                            labels = if (remove_baseline) c("rapi" = "RAPI", "upi" = "Matched-Pair UPI", "tspa" = "2S-PA-Int") else c("mmr" = "MMR", "rapi" = "RAPI", "upi" = "Matched-Pair UPI", "tspa" = "2S-PA-Int"))
  }
  
  if (!is.null(y_limits)) {
    plot <- plot + ylim(y_limits)
  }
  
  return(plot)
}

# Split the data
type1data <- df %>% filter(path == 0)
powerdata_1 <- df %>% filter(path == 0.3)
powerdata_2 <- df %>% filter(path == 0.6)

type1_pd <- process_data(type1data)
power1_pd <- process_data(powerdata_1)
power2_pd <- process_data(powerdata_2)

# Write the processed data
write_csv(type1_pd, "Sim_Data/type1_10012024.csv")
write_csv(power1_pd, "Sim_Data/power1_10012024.csv")
write_csv(power2_pd, "Sim_Data/power2_10012024.csv")

# Plot the results
type1_plot <- read.csv("Sim_Data/type1_10012024.csv")
power1_plot <- read.csv("Sim_Data/power1_10012024.csv")
power2_plot <- read.csv("Sim_Data/power2_10012024.csv")

# Standard Bias
sd_bias_plot <- plot_results(power1_plot,
                             "N",
                             "sb",
                             "Sample Size (N)",
                             "Standardized Bias",
                             c(-0.5, 0.5),
                             add_lines = TRUE,
                             remove_baseline = TRUE)
sd_bias_plot <- sd_bias_plot +
  geom_hline(yintercept = c(-0.4, 0.4), linetype = "dashed", color = "red")

# Median-Mad Relative SE Bias
rse_bias_plot <- plot_results(power_plot,
                              "N",
                              "stdMed_rse_bias",
                              "Sample Size (N)",
                              "Relative SE Bias",
                              c(-0.3, 0.3),
                              remove_baseline = TRUE)
rse_bias_plot <- rse_bias_plot +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "red")

# Coverage rate
cov_bias_plot <- plot_results(power_plot,
                              "N",
                              "coverage_std",
                              "Sample Size (N)",
                              "Coverage Rate of 95% CI",
                              remove_baseline = TRUE)
cov_bias_plot <- cov_bias_plot +
  geom_hline(yintercept = c(0.91), linetype = "dashed", color = "red")

# RMSE
rmse_bias_plot <- plot_results(power_plot,
                               "N",
                               "rmse",
                               "Sample Size (N)",
                               "Root Mean Squre Error",
                               add_lines = TRUE)

# Type I error rate
type1_bias_plot <- plot_results(type1_plot,
                                "N",
                                "type1_lrt",
                                "Sample Size (N)",
                                "Type I Error Rate",
                                c(0, 0.1),
                                remove_baseline = TRUE)

# Statistical Power
power_bias_plot <- plot_results(power_plot,
                                "N",
                                "power_lrt",
                                "Sample Size (N)",
                                "Statistical Power",
                                remove_baseline = TRUE)
power_bias_plot <- power_bias_plot +
  geom_hline(yintercept = c(0.8), linetype = "dashed", color = "red")
