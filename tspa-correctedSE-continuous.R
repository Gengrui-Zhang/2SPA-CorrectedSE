library(lavaan)
library(SimDesign)
library(dplyr)
library(here)

# Source Function
r_scripts <- list.files(here("R"), pattern = "\\.R$", full.names = TRUE)
lapply(r_scripts, source)

# ========================================= Simulation Conditions ========================================= #

# Design Factor
DESIGNFACTOR <- createDesign(
  N = c(100, 500),
  rel = c(0.7, 0.9),
  reg = c(0, 0.3)
  # num_items = c(3, 10) # Need to change. Now it's hardcoding
)

# Fixed parameters
FIXED_PARAMETER <- list(
  joint_mod = "fy =~ y1 + y2 + y3
                fx =~ x1 + x2 + x3
                fy ~ fx",
  
  cfa_mod = "fy =~ y1 + y2 + y3
              fx =~ x1 + x2 + x3",
  
  rel_mod = "fy =~ NA * fs_fy + ly * fs_fy
              fx =~ NA * fs_fx + lx * fs_fx
              fs_fy ~~ ev1 * fs_fy
              fs_fx ~~ ev2 * fs_fx
              fx ~~ 1 * fx
              fy ~ beta * fx
              fy ~~ dvy * fy
              # constraints
              ev2 == ep2 * lx^2
              1 == beta^2 + dvy
              ev1 == ep1 * ly^2"
)

# ========================================= Data Generation ========================================= #

sim_sem_dat <- function(num_obs, Lambda, Psi, Theta) {
  eta_delta <- MASS::mvrnorm(
    num_obs,
    mu = rep(0, sum(dim(Lambda))),  
    Sigma = Matrix::bdiag(list(Psi, Theta))  # Block diagonal covariance matrix
  )
  
  ind <- eta_delta[, 1:2] %*% t(Lambda) + eta_delta[, -(1:2)] # Observed Indicators
  
  dat <- data.frame(ind)
  names(dat) <- c(paste0("x", 1:3), paste0("y", 1:3))
  return(dat)
}

# Data Generation Function
generate_dat <- function (condition, fixed_objects = NULL) {
  N <- condition$N
  rel = condition$rel
  reg = condition$reg
  # num_items = condition$num_items
  
  # Model Parameters
  Alpha <- c(0, 0)
  Lambda_1 <- c(.8, .5, .6)
  Lambda_2 <- c(.6, .7, .5)
  Lambda <- cbind(c(Lambda_1, rep(0, 3)),
                  c(rep(0, 3), Lambda_2))
  Psi <- matrix(c(1, reg, 
                  reg, 1), 
                nrow = 2) # reg = 0.5
  # Compute error variance
  Error_1 <- sum(Lambda_1)^2*(1 - rel)/rel*c(0.44, 0.33, 0.23)
  Error_2 <- sum(Lambda_2)^2*(1 - rel)/rel*c(0.44, 0.33, 0.23)
  Theta <- diag(c(Error_1, Error_2)) 
  
  # Simulate Data
  sim_sem_dat(num_obs = N, Lambda = Lambda, Psi = Psi, Theta = Theta)
}


# ========================================= Data Analysis ========================================= #
# ====== Helper Function ========= #

# Helper Function 1
# Adjusting a parameter in a model updating the model with the new parameter values and computing the objective function
ff <- function(par, par_id, obj) {
  est1 <- partable(obj)$est
  est1[par_id] <- par
  glist1 <- obj@Model@GLIST
  pt1 <- lavInspect(obj, "partable")
  for (l in seq_along(glist1)) {
    to_fill <- which(pt1[[l]] != 0)
    glist1[[l]][to_fill] <- est1[pt1[[l]][to_fill]] # Updating parameters
  }
  lavaan:::lav_model_objective(
    obj@Model,
    GLIST = glist1,
    lavsamplestats = obj@SampleStats,
    lavdata = obj@Data
  )
}

# Helper Function 2
# Convert unstandardized parameter estimates into standardized estimates
unstd2std <- function(par, par_id, obj) {
  est1 <- partable(obj)$est
  est1[par_id] <- par
  glist1 <- obj@Model@GLIST
  pt1 <- lavInspect(obj, "partable")
  for (l in seq_along(glist1)) {
    to_fill <- which(pt1[[l]] != 0)
    glist1[[l]][to_fill] <- est1[pt1[[l]][to_fill]]
  }
  standardizedSolution(obj,
                       type = "std.lv",
                       est = est1,
                       GLIST = glist1
  )$est[par_id]
}

# Helper Function 3
# Computing standardized parameters for tspa
std_est <- function(ev, obj, fs) {
  fsT1 <- attr(obj, which = "fsT")
  diag(fsT1) <- ev
  mobj <- tspa(" fy ~ fx",
               data = fs,
               fsT = fsT1, fsL = attr(fs, "fsL")
  )
  standardizedSolution(mobj, type = "std.lv")$est.std[8]
}

# Helper Function 4
compute_std_est <- function(ep, mod, fs) {
  mod_tmp <- gsub("ep2", ep[2], mod)
  mod_tmp <- gsub("ep1", ep[1], mod_tmp)
  standardizedSolution(sem(mod_tmp, data = fs))$est[6:7]
}

# ====== Analyses Function ========= #

# Joint Model
analyze_joint <- function(condition, dat, fixed_objects) {
  
  # Extract Model
  joint_mod <- FIXED_PARAMETER$joint_mod
  
  # Fit Model
  mjoint <- sem(joint_mod, data = dat)
  joint_est <- standardizedSolution(mjoint)[7, ]
  
  # Check convergence
  converge <- ifelse(lavInspect(mjoint, "converged"), 1, 0)
  
  # Extract Parameter
  out <- c(joint_est[["est.std"]], joint_est[["se"]], converge)
  names(out) <- c("est", "se", "converge")
  
  # Check if any NA exists in the output
  if (anyNA(out)) {
    stop("The model did not obtain SE")
  }
  return(out)
}

# SAM Global Model
analyze_gsam <- function(condition, dat, fixed_objects) {
  
  # Extract Model
  joint_mod <- FIXED_PARAMETER$joint_mod
  
  # Fit Model
  mgsam <- sam(joint_mod, data = dat, sam.method = "global")
  gsam_est <- standardizedSolution(mgsam)[7, ]
  
  # Check convergence
  converge <- ifelse(lavInspect(mgsam, "converged"), 1, 0)
  
  # Extract Parameter
  out <- c(gsam_est[["est.std"]], gsam_est[["se"]], converge)
  names(out) <- c("est", "se", "converge")
  
  # Check if any NA exists in the output
  if (anyNA(out)) {
    stop("The model did not obtain SE")
  }
  return(out)
}

# SAM Local Model
analyze_lsam <- function(condition, dat, fixed_objects) {
  
  # Extract Model
  joint_mod <- FIXED_PARAMETER$joint_mod
  
  # Fit Model
  mlsam <- sam(joint_mod, data = dat)
  lsam_est <- standardizedSolution(mlsam)[7, ]
  
  # Check convergence
  converge <- ifelse(lavInspect(mlsam, "converged"), 1, 0)
  
  # Extract Parameter
  out <- c(lsam_est[["est.std"]], lsam_est[["se"]], converge)
  names(out) <- c("est", "se", "converge")
  
  # Check if any NA exists in the output
  if (anyNA(out)) {
    stop("The model did not obtain SE")
  }
  return(out)
}

# 2S-PA Model with Corrected SE
analyze_tspa <- function(condition, dat, fixed_objects) {

  # Extract Model
  cfa_mod <- FIXED_PARAMETER$cfa_mod
  
  # Get Factor Score
  cfa_fit <- cfa(cfa_mod,
                 std.lv = TRUE, 
                 data = dat,
                 bounds = "pos.ov.var"
  )
  fs <- get_fs_lavaan(cfa_fit, method = "Bartlett", vfsLT = TRUE)
  
  # Fit Model
  m2spa1 <- tspa("fy ~ fx",
                 data = fs,
                 fsT = attr(fs, "fsT"), fsL = attr(fs, "fsL")
  )
  tspa1_est <- standardizedSolution(m2spa1)[8, ]

  # Corrected Standard Errors
  v1 <- attr(fs, "vfsLT")[c(5, 7), c(5, 7)]
  # Gradient
  grad_std <- numDeriv::grad(std_est,
                             x = diag(attr(fs, "fsT")), 
                             obj = m2spa1,
                             fs = fs)
  vc_corrected2 <- lavInspect(m2spa1, what = "vcov.std")[1, 1] +
    t(grad_std) %*% v1 %*% grad_std
  
  # Check convergence
  converge <- ifelse(lavInspect(m2spa1, "converged"), 1, 0)

  # Extract Parameter
  out <- c(tspa1_est[["est.std"]], tspa1_est[["se"]],
           tspa1_est[["est.std"]], sqrt(vc_corrected2[1, 1]),
           converge)
  names(out) <- c("est", "se",
                  "est_corrected", "se_corrected",
                  "converge")
  
  # Check if any NA exists in the output
  if (anyNA(out)) {
    stop("The model did not obtain SE")
  }
  return(out)
}

# Reliability Model with Corrected Standard Error
analyze_rel <- function(condition, dat, fixed_objects) {

  # Extract Model
  cfa_mod <- FIXED_PARAMETER$cfa_mod
  rel_mod <- FIXED_PARAMETER$rel_mod
  
  # Get Factor Score
  cfa_fit <- cfa(cfa_mod,
                 std.lv = TRUE, 
                 data = dat,
                 bounds = "pos.ov.var"
  )
  fs <- get_fs_lavaan(cfa_fit, method = "Bartlett", vfsLT = TRUE)
  
  # Fit Model: Reliability
  ep <- diag(attr(fs, "fsT"))
  mod_tmp <- gsub("ep2", ep[2], rel_mod)
  mod_tmp <- gsub("ep1", ep[1], mod_tmp)
  m2spa2 <- sem(mod_tmp, data = fs)
  tspa2_est <- standardizedSolution(m2spa2)[6, ]
  
  # Corrected SE
  v1 <- attr(fs, "vfsLT")[c(5, 7), c(5, 7)]
  
  jac_std <- numDeriv::jacobian(compute_std_est,
                                x = diag(attr(fs, "fsT")), 
                                mod = rel_mod,
                                fs = fs)
  vc_corrected3 <- lavInspect(m2spa2, what = "vcov.std")[5:6, 5:6] +
    jac_std %*% v1 %*% t(jac_std)
  
  # Check convergence
  converge <- ifelse(lavInspect(m2spa2, "converged"), 1, 0)
  
  # Extract Parameter
  out <- c(tspa2_est[["est.std"]], tspa2_est[["se"]],
           tspa2_est[["est.std"]], sqrt(vc_corrected3[1, 1]),
           converge)
  names(out) <- c("est", "se",
                  "est_corrected", "se_corrected",
                  "converge")
  
  # Check if any NA exists in the output
  if (anyNA(out)) {
    stop("The model did not obtain SE")
  }
  return(out)
}

# ========================================= Results Summary ========================================= #

# Helper function: robust bias
robust_bias <- function(est, se, par, trim = 0, type = NULL) {
  output <- numeric(ncol(est))
  for (i in seq_len(ncol(est))) {
    if (type == "raw") {
      output[i] <- mean((est[,i] - par), na.rm = TRUE)
    } else if (type == "standardized") {
      output[i] <- (mean(est[,i], na.rm = TRUE) - par)/sd(est[,i], na.rm = TRUE)
    } else if (type == "trim") {
      output[i] <- mean(est[,i], trim = trim, na.rm = TRUE) - par
    } else if (type == "median") {
      output[i] <- (median(est[,i], na.rm = TRUE) - par) / mad(est[,i], na.rm = TRUE)
    } else {
      output[i] <- (mean(est[,i], trim = trim, na.rm = TRUE) - par) / sd(est[,i], na.rm = TRUE)
    }
  }
  names(output) <- colnames(est)
  return(output)
}

# Helper function: relative SE bias
rse_bias <- function(est, se, trim = 0, type = "raw") {
  if (type == "raw") {
    se_mean <- apply(se, 2, mean, na.rm = T)
    se_sd <- apply(est, 2L, sd, na.rm = T)
    rse_bias <- se_mean / se_sd - 1
  } else if (type == "median") {
    se_median <- apply(se, 2, median, na.rm = TRUE)
    mad_sd <- apply(est, 2, function(x) mad(x, na.rm = TRUE))
    rse_bias <- se_median / mad_sd - 1
  } else if (type == "trim") {
    se_mean <- apply(se, 2, mean, trim = trim, na.rm = TRUE)
    se_sd <- apply(est, 2L, sd, na.rm = T)
    rse_bias <- se_mean / se_sd - 1
  }
  return(rse_bias)
}

# Helper function: detecting outliers for SE
outlier_se <- function(se) {
  results <- c()
  for(column in colnames(se)) {
    # Calculate Q1, Q3, and IQR
    Q1 <- quantile(se[,column], 0.25, na.rm = TRUE)
    Q3 <- quantile(se[,column], 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    # Determine outliers
    lower_bound <- (Q1 - 1.5 * IQR)
    upper_bound <- (Q3 + 1.5 * IQR)
    outliers <- se[,column][se[,column] < lower_bound | se[,column] > upper_bound]
    # Calculate the percentage of outliers
    percentage <- length(outliers) / sum(!is.na(se[,column])) * 100
    results[column] <- percentage
  }
  return(results)
}

# Helper function for convergence rate
convergence_rate <- function(converge) {
  apply(converge, 2, function(x) 1-(sum(is.na(x)) / length(x)))
}

# Helper function for calculating coverage rate, Type I error rate, and power
ci_stats <- function(est, se, par, stats_type) {
  
  # Calculate the confidence intervals (std)
  lo_95 <- est - qnorm(.975) * se
  up_95 <- est + qnorm(.975) * se
  ci_est <- vector("list", length = ncol(est))
  names(ci_est) <- colnames(est)
  
  # Construct confidence intervals for each method
  for (i in seq_len(ncol(est))) {
    ci_est[[i]] <- cbind(lo_95[,i], up_95[,i])
  }

  # Determine which statistic to calculate
  if (stats_type == "Coverage") {
    return(sapply(ci_est, function(ci) mean(ci[,1] <= par & ci[,2] >= par)))
  } else if (stats_type == "TypeI") {
    return(sapply(ci_est, function(ci) mean(ci[,1] > 0 | ci[,2] < 0)))
  } else if (stats_type == "Power") {
    return(sapply(ci_est, function(ci) (1 - mean(ci[,1] < 0 & ci[,2] > 0))))
  } else {
    stop("Invalid stats_type specified. Please choose from 'Coverage', 'TypeI', or 'Power'.")
  }
}

# Evaluation Function
evaluate_res <- function (condition, results, fixed_objects = NULL) {
  
  # Population parameter
  pop_par <- condition$reg
  
  # Parameter estimates
  est <- results[, grep(".est$", colnames(results))]
  se <- results[, grep(".se$", colnames(results))]
  converge <- results[, grep(".converge$", colnames(results))]
  # Check if corrected estimates exist
  est_corrected <- results[, grep(".est_corrected$", colnames(results)), drop = FALSE]
  se_corrected <- results[, grep(".se_corrected$", colnames(results)), drop = FALSE]

  c(rb = robust_bias(est,
                     se,
                     pop_par,
                     type = "raw"),
    rb_corrected = robust_bias(est_corrected,
                               se_corrected,
                               pop_par,
                               type = "raw"),
    sb = robust_bias(est,
                     se,
                     pop_par,
                     type = "standardized"),
    sb_corrected = robust_bias(est_corrected,
                               se_corrected,
                               pop_par,
                               type = "standardized"),
    # trim_bias = robust_bias(est,
    #                         se,
    #                         pop_par,
    #                         trim = 0.2,
    #                         type = "trim"), # 20% trimmed mean
    # stdMed_bias = robust_bias(est,
    #                           se,
    #                           pop_par,
    #                           type = "median"),
    mean_se = colMeans(se, na.rm = TRUE),
    mean_se_corrected = colMeans(se_corrected, na.rm = TRUE),
    raw_rsb = rse_bias(est,
                       se,
                       type = "raw"),
    raw_rsb_corrected = rse_bias(est_corrected,
                                 se_corrected,
                                 type = "raw"),
    # stdMed_rse_bias = rse_bias(est,
    #                            se,
    #                            type = "median"),
    # trim_rse_bias = rse_bias(est,
    #                          se,
    #                          trim = 0.2,
    #                          type = "trim"),
    # outlier_se = outlier_se(se),
    coverage = ci_stats(est,
                        se,
                        pop_par,
                        "Coverage"),
    coverage_corrected = ci_stats(est_corrected,
                                  se_corrected,
                                  pop_par,
                                  "Coverage"),
    type1 = ci_stats(est,
                     se,
                     pop_par,
                     "TypeI"),
    type1_corrected = ci_stats(est_corrected,
                               se_corrected,
                               pop_par,
                               "TypeI"),
    power = ci_stats(est,
                     se,
                     pop_par,
                     "Power"),
    power_corrected = ci_stats(est_corrected,
                               se_corrected,
                               pop_par,
                               "Power"),
    rmse = RMSE(na.omit(est),
                parameter = pop_par),
    rmse_corrected = RMSE(na.omit(est_corrected),
                          parameter = pop_par),
    convergence = convergence_rate(converge)
  )
}

# ========================================= Run Experiment ========================================= #

res <- runSimulation(design = DESIGNFACTOR,
                     replications = 200,
                     generate = generate_dat,
                     analyse = list(joint = analyze_joint,
                                    gsam = analyze_gsam,
                                    lsam = analyze_lsam,
                                    tspa = analyze_tspa,
                                    rel = analyze_rel),
                     summarise = evaluate_res,
                     fixed_objects = FIXED_PARAMETER,
                     seed = rep(66330, nrow(DESIGNFACTOR)),
                     packages = "lavaan",
                     filename = "CorrectedSE_09242024",
                     parallel = TRUE,
                     ncores = 6,
                     save = TRUE,
                     save_results = TRUE)


# Notes:
# 1. Warning to error
# 2. Number of items
# 3. Try regression coefficients: larger regression coefficients could give more window of bias
# 4. Add more comments (especially for analysis functions)

