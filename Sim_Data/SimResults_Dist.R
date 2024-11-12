# Read Data
library(here)
dat_dir <- here("Sim_Data", "CorrectedSE_10212024_Raw")

# Wrong order: Define the name directly
# List all files matching the pattern
dat_list <- list.files(path = dat_dir, pattern = "results-row-.*\\.rds", full.names = TRUE)
file_numbers <- as.numeric(gsub(".*results-row-([0-9]+)\\.rds", "\\1", dat_list))
dat_list <- dat_list[order(file_numbers)]

results_list <- lapply(dat_list, function(file) {
  data <- readRDS(file)
})

data_list <- lapply(dat_list, function(file) {
  data <- readRDS(file)
  data[["results"]]  
})

condition_list <- do.call(rbind, lapply(dat_list, function(file) {
  data <- readRDS(file)
  data[["condition"]]  
}))

# Examine SE
extract_columns <- function(data, column_names, suffix) {
  columns <- setNames(
    lapply(column_names, function(name) {
      col_name <- paste0(name, suffix)
      if (col_name %in% names(data)) {
        return(data[[col_name]])
      } else {
        warning(paste("Column", col_name, "not found in the data. Returning NA."))
        return(rep(NA, nrow(data)))
      }
    }),
    paste0(column_names, suffix)
  )
  data.frame(columns)
}

remove_na_columns <- function(df_list) {
  lapply(df_list, function(df) {
    # Remove columns with all NA values
    df[, colSums(is.na(df)) < nrow(df)]
  })
}

column_names <- c("joint", "gsam", "lsam", "tspa", "rel")
est_list <- lapply(data_list, extract_columns, column_names = column_names, suffix = ".est")
se_list <- lapply(data_list, extract_columns, column_names = column_names, suffix = ".se")
se_corrected_list <- lapply(data_list, extract_columns, column_names = column_names, suffix = ".se_corrected")
se_corrected_list <- remove_na_columns(se_corrected_list)

# Examine Condition
# Condition 3
est <- est_list[[3]]
se <- se_list[[3]]
est_corrected <- est_list[[3]][c("tspa.est", "rel.est")]
se_corrected <- se_corrected_list[[3]]

# Density Plot
par(mfrow = c(2, 3))
for (i in 1:ncol(se)) {
  emp_sd <- sd(est[[i]], na.rm = TRUE)
  mean_se <- mean(se[[i]], na.rm = TRUE)  # Calculate the mean
  plot(density(se[[i]], na.rm = TRUE), 
       main = paste("Density Plot of SE -", colnames(se)[i]), 
       xlab = "SE", 
       col = "blue", 
       lwd = 2)
  abline(v = emp_sd, col = "red", lwd = 2, lty = 2)
  abline(v = mean_se, col = "green", lwd = 2, lty = 2)
}
par(mfrow = c(1, 1))

# Mean and Median
apply(est, 2, sd, na.rm = TRUE) # Empirical SD
apply(se, 2, mean, na.rm = TRUE) # Mean SE
apply(se_corrected, 2, mean, na.rm = TRUE) # Mean SE
apply(est, 2, mean, na.rm = TRUE) # Mean Estimate

# Boxplot
par(mfrow = c(2, 3))
for (i in 1:ncol(se)) {
  boxplot(se[[i]], 
          main = paste("Boxplot of SE -", colnames(se)[i]), 
          ylab = "SE", 
          col = "lightgreen",
          ylim = c(0, 0.25))
}
par(mfrow = c(1, 1))

# Outliers
Q1 <- quantile(se, 0.25, na.rm = TRUE)
Q3 <- quantile(se, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1
lower_bound <- (Q1 - 1.5 * IQR)
upper_bound <- (Q3 + 1.5 * IQR)
outliers <- se[se < lower_bound | se > upper_bound]
percentage <- length(outliers) / sum(!is.na(se)) * 100




