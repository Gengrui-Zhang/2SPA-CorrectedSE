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
  rel = c(0.7, 0.9),
  N_per_p = c(6000, 25000, 100000), # Sample Size per Indicator
  path = list(c(0, 0.7),
              c(0.3, 0.4),
              c(0.6, 0.1)), # May include negative paths
  n_items = c(5, 10, 20), # Number of Indicator
  cor = c(0, 0.3, 0.6) # Correlation between the two factors
) 

# Fixed parameters
FIXED_PARAMETER <- list(
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

sim_sem_dat <- function(num_obs, Lambda, Psi, Theta, n_items) {
  eta_delta <- MASS::mvrnorm(
    num_obs,
    mu = rep(0, sum(dim(Lambda))),  
    Sigma = Matrix::bdiag(list(Psi, Theta))  # Block diagonal covariance matrix
  )
  ind <- eta_delta[, 1:ncol(Psi)] %*% t(Lambda) + eta_delta[, -(1:ncol(Psi))] # Observed Indicators
  dat <- data.frame(ind)
  names(dat) <- c(paste0("x", 1:n_items), paste0("m", 1:n_items), paste0("y", 1:n_items))
  return(dat)
}

# Data Generation Function
generate_dat <- function (condition, fixed_objects = NULL) {
  n_items <- condition$n_items
  N <- condition$N_per_p * n_items # Sample Size
  rel = condition$rel
  cor = condition$cor
  path = unlist(condition$path)
  
  # Model Parameters
  # Alpha <- rep(0, 2)
  Lambda_1 <- rep(c(0.4, 0.5, 0.6, 0.7, 0.8), each = n_items / 5)
  Lambda_2 <- rep(c(0.45, 0.55, 0.65, 0.75, 0.85), each = n_items / 5)
  Lambda_3 <- rep(c(0.35, 0.45, 0.55, 0.65, 0.75), each = n_items / 5)
  Lambda <- cbind(c(Lambda_1, rep(0, n_items), rep(0, n_items)), 
                  c(rep(0, n_items), Lambda_2, rep(0, n_items)),
                  c(rep(0, n_items), rep(0, n_items), Lambda_3))
  Psi <- matrix(c(1, cor, path[1],
                  cor, 1, path[2],
                  path[2], path[2], 1),
                nrow = 3)
  
  # Compute error variance
  err_prop <- 0.9^(0:(n_items - 1))
  err_prop <- err_prop / sum(err_prop)
  Error_1 <- sum(Lambda_1)^2*(1 - rel)/rel*err_prop
  Error_2 <- sum(Lambda_2)^2*(1 - rel)/rel*err_prop
  Error_3 <- sum(Lambda_3)^2*(1 - rel)/rel*err_prop
  Theta <- diag(c(Error_1, Error_2, Error_3)) 
  
  # Simulate Data
  sim_sem_dat(num_obs = N, Lambda = Lambda, 
              Psi = Psi, Theta = Theta,
              n_items = n_items)
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
std_est <- function(ev, obj, fs, target = c("fx", "fm")) {
  target <- match.arg(target)
  fsT1 <- attr(obj, which = "fsT")
  diag(fsT1) <- ev
  mobj <- tspa("fy ~ fx
                fy ~ fm",
               data = fs,
               fsT = fsT1, fsL = attr(fs, "fsL"))
  std_est <- standardizedSolution(mobj, type = "std.lv") %>%
    filter(lhs == "fy", op == "~", rhs == target)
  unname(std_est$est.std)  
}


# Helper Function 4
compute_std_est <- function(ep, mod, fs) {
  mod_tmp <- gsub("ep2", ep[2], mod)
  mod_tmp <- gsub("ep1", ep[1], mod_tmp)
  standardizedSolution(sem(mod_tmp, data = fs))$est[6:7]
}

# Helper Function 5
syntax_update <- function(n_items, path_include = TRUE) {
  fy_loadings <- paste0("y", 1:n_items, collapse = " + ")
  fx_loadings <- paste0("x", 1:n_items, collapse = " + ")
  fm_loadings <- paste0("m", 1:n_items, collapse = " + ")
  lavaan_syntax <- paste0(
    "fy =~ ", fy_loadings, "\n",
    "fx =~ ", fx_loadings, "\n",
    "fm =~ ", fm_loadings
  )
  if (path_include) {
    lavaan_syntax <- paste0(lavaan_syntax, 
                            "\nfy ~ fx",
                            "\nfy ~ fm")
  }
  return(lavaan_syntax)
}

# ====== Analyses Function ========= #

# Joint Model
analyze_joint <- function(condition, dat, fixed_objects) {
  
  # Initialize a local warning counter for this replication
  local_warning_counter <- 0
  
  # Fit the model and capture warnings
  result <- withCallingHandlers(
    {
      # Extract Model
      joint_mod <- syntax_update(condition$n_items, path_include = TRUE)
      
      # Fit Model
      mjoint <- sem(joint_mod, data = dat)
      mjoint_est <- standardizedSolution(mjoint)
      joint_est_yx <-  mjoint_est %>%
        filter(lhs == "fy", op == "~", rhs == "fx")
      joint_est_ym <-  mjoint_est %>%
        filter(lhs == "fy", op == "~", rhs == "fm")
      
      # Convergence Check
      converge <- ifelse(lavInspect(mjoint, "converged"), 1, 0)
      
      # Extract Parameter
      out <- c(joint_est_yx[["est.std"]], joint_est_yx[["se"]], 
               joint_est_ym[["est.std"]], joint_est_ym[["se"]], 
               local_warning_counter, converge)
      names(out) <- c("est_yx", "se_yx", "est_ym", "se_ym", 
                      "warnings_count", "convergence")
      
      # Check if any NA exists in the output
      if (anyNA(out)) {
        stop("The model did not obtain SE")
      }
      return(out)
    },
    warning = function(w) {
      # Increment the local warning counter
      local_warning_counter <<- local_warning_counter + 1
    }
  )
  
  return(result)
}

# SAM Global Model
analyze_gsam <- function(condition, dat, fixed_objects) {
  
  # Initialize a local warning counter for this replication
  local_warning_counter <- 0
  
  # Fit the model and capture warnings
  result <- withCallingHandlers(
    {
      # Extract Model
      joint_mod <- syntax_update(condition$n_items, path_include = TRUE)
      
      # Fit Model
      mgsam <- sam(joint_mod, data = dat, sam.method = "global")
      mgsam_est <- standardizedSolution(mgsam)
      gsam_est_yx <-  mgsam_est %>%
        filter(lhs == "fy", op == "~", rhs == "fx")
      gsam_est_ym <-  mgsam_est %>%
        filter(lhs == "fy", op == "~", rhs == "fm")
      
      # Convergence Check
      converge <- ifelse(lavInspect(mgsam, "converged"), 1, 0)
      
      # Extract Parameter
      out <- c(gsam_est_yx[["est.std"]], gsam_est_yx[["se"]], 
               gsam_est_ym[["est.std"]], gsam_est_ym[["se"]], 
               local_warning_counter, converge)
      names(out) <- c("est_yx", "se_yx", "est_ym", "se_ym", 
                      "warnings_count", "convergence")
      
      # Check if any NA exists in the output
      if (anyNA(out)) {
        stop("The model did not obtain SE")
      }
      return(out)
    },
    warning = function(w) {
      # Increment the local warning counter
      local_warning_counter <<- local_warning_counter + 1
    }
  )
  
  return(result)
}

# SAM Local Model
analyze_lsam <- function(condition, dat, fixed_objects) {
  
  # Initialize a local warning counter for this replication
  local_warning_counter <- 0
  
  # Fit the model and capture warnings
  result <- withCallingHandlers(
    {
      # Extract Model
      joint_mod <- syntax_update(condition$n_items, path_include = TRUE)
      
      # Fit Model
      mlsam <- sam(joint_mod, data = dat)
      mlsam_est <- standardizedSolution(mlsam)
      lsam_est_yx <-  mlsam_est %>%
        filter(lhs == "fy", op == "~", rhs == "fx")
      lsam_est_ym <-  mlsam_est %>%
        filter(lhs == "fy", op == "~", rhs == "fm")
      
      # Convergence Check
      converge <- ifelse(lavInspect(mlsam, "converged"), 1, 0)
      
      # Extract Parameter
      out <- c(lsam_est_yx[["est.std"]], lsam_est_yx[["se"]], 
               lsam_est_ym[["est.std"]], lsam_est_ym[["se"]], 
               local_warning_counter, converge)
      names(out) <- c("est_yx", "se_yx", "est_ym", "se_ym", 
                      "warnings_count", "convergence")
      
      # Check if any NA exists in the output
      if (anyNA(out)) {
        stop("The model did not obtain SE")
      }
      return(out)
    },
    warning = function(w) {
      # Increment the local warning counter
      local_warning_counter <<- local_warning_counter + 1
    }
  )
  
  return(result)
}

# 2S-PA Model with Corrected SE
analyze_tspa <- function(condition, dat, fixed_objects) {
  
  # Initialize a local warning counter for this replication
  local_warning_counter <- 0
  
  # Fit the model and capture warnings
  result <- withCallingHandlers(
    {
      # Extract Model
      cfa_mod <- syntax_update(condition$n_items, path_include = FALSE)
      
      # Get Factor Score
      cfa_fit <- cfa(cfa_mod,
                     std.lv = TRUE, 
                     data = dat,
                     bounds = "pos.ov.var") 
        # Sets boundaries for parameter estimation to ensure positive variances for observed variables
      
      fs <- get_fs_lavaan(cfa_fit, method = "Bartlett", vfsLT = TRUE)
      
      # Fit Model
      m2spa1 <- tspa("fy ~ fx
                      fy ~ fm",
                     data = fs,
                     fsT = attr(fs, "fsT"), fsL = attr(fs, "fsL"))
      
      m2spa1_est <- standardizedSolution(m2spa1)
      tspa1_est_fx <-  m2spa1_est %>%
        filter(lhs == "fy", op == "~", rhs == "fx")
      tspa1_est_fm <-  m2spa1_est %>%
        filter(lhs == "fy", op == "~", rhs == "fm")
      
      # Corrected Standard Errors
      v1 <- attr(fs, "vfsLT")[c(10, 13, 15), c(10, 13, 15)] # fy, fx, fm
      
      # Gradient
      grad_std_x <- numDeriv::grad(std_est,
                                   x = diag(attr(fs, "fsT")), 
                                   obj = m2spa1,
                                   fs = fs,
                                   target = "fx")
      vc_corrected2_x <- lavInspect(m2spa1, what = "vcov.std")[1, 1] +
        t(grad_std_x) %*% v1 %*% grad_std_x
      
      grad_std_m <- numDeriv::grad(std_est,
                                   x = diag(attr(fs, "fsT")), 
                                   obj = m2spa1,
                                   fs = fs,
                                   target = "fm")
      vc_corrected2_m <- lavInspect(m2spa1, what = "vcov.std")[1, 1] +
        t(grad_std_m) %*% v1 %*% grad_std_m
      
      # Convergence Check
      converge <- ifelse(lavInspect(m2spa1, "converged"), 1, 0)
      
      # Extract Parameter
      out <- c(tspa1_est[["est.std"]], tspa1_est[["se"]],
               tspa1_est[["est.std"]], sqrt(vc_corrected2[1, 1]), 
               local_warning_counter,
               converge)
      names(out) <- c("est", "se", "est_corrected", "se_corrected", 
                      "warnings_count", "convergence")
      
      # Check if any NA exists in the output
      if (anyNA(out)) {
        stop("The model did not obtain SE")
      }
      return(out)
    },
    warning = function(w) {
      # Increment the local warning counter
      local_warning_counter <<- local_warning_counter + 1
    }
  )
  
  return(result)
}

# Reliability Model with Corrected Standard Error
analyze_rel <- function(condition, dat, fixed_objects) {
  
  # Initialize a local warning counter for this replication
  local_warning_counter <- 0
  
  # Fit the model and capture warnings
  result <- withCallingHandlers(
    {
      # Extract Model
      cfa_mod <- syntax_update(condition$n_items, path_include = FALSE)
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
      m2spa2_est <- standardizedSolution(m2spa2)
      tspa2_est <- m2spa2_est %>%
        filter(lhs == "fy", op == "~", rhs == "fx")
      
      # Corrected SE
      v1 <- attr(fs, "vfsLT")[c(5, 7), c(5, 7)]
      
      jac_std <- numDeriv::jacobian(compute_std_est,
                                    x = diag(attr(fs, "fsT")), 
                                    mod = rel_mod,
                                    fs = fs)
      vc_corrected3 <- lavInspect(m2spa2, what = "vcov.std")[5:6, 5:6] +
        jac_std %*% v1 %*% t(jac_std)
      
      # Convergence Check
      converge <- ifelse(lavInspect(m2spa2, "converged"), 1, 0)
      
      # Extract Parameter
      out <- c(tspa2_est[["est.std"]], tspa2_est[["se"]],
               tspa2_est[["est.std"]], sqrt(vc_corrected3[1, 1]), 
               local_warning_counter, converge)
      names(out) <- c("est", "se", "est_corrected", "se_corrected", 
                      "warnings_count", "convergence")
      
      # Check if any NA exists in the output
      if (anyNA(out)) {
        stop("The model did not obtain SE")
      }
      return(out)
    },
    warning = function(w) {
      # Increment the local warning counter
      local_warning_counter <<- local_warning_counter + 1
    }
  )
  
  return(result)
}

# ========================================= Results Summary ========================================= #
# Helper function for convergence rate
convergence_rate <- function(converge) {
  apply(converge, 2, function(x) 1-(sum(is.na(x)) / length(x)))
}

# Helper function for warning sum
warning_sum <- function(count) {
  apply(count, 2, sum, na.rm = TRUE)
}

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
  pop_par <- condition$path
  
  # Parameter estimates
  est <- results[, grep(".est$", colnames(results))]
  se <- results[, grep(".se$", colnames(results))]
  # Check if corrected estimates exist
  est_corrected <- results[, grep(".est_corrected$", colnames(results)), drop = FALSE]
  se_corrected <- results[, grep(".se_corrected$", colnames(results)), drop = FALSE]
  # Convergence and Warning
  convergences <- results[, grep(".convergence$", colnames(results))]
  warnings <- results[, grep(".warnings_count$", colnames(results))]
  
  c(rbias = robust_bias(est,
                        se,
                        pop_par,
                        type = "raw"),
    rbias_corrected = robust_bias(est_corrected,
                                  se_corrected,
                                  pop_par,
                                  type = "raw"),
    sbias = robust_bias(est,
                        se,
                        pop_par,
                        type = "standardized"),
    sbias_corrected = robust_bias(est_corrected,
                                  se_corrected,
                                  pop_par,
                                  type = "standardized"),
    mean_se = colMeans(se, na.rm = TRUE),
    mean_se_corrected = colMeans(se_corrected, na.rm = TRUE),
    raw_rsb = rse_bias(est,
                       se,
                       type = "raw"),
    raw_rsb_corrected = rse_bias(est_corrected,
                                 se_corrected,
                                 type = "raw"),
    stdMed_rsb = rse_bias(est,
                          se,
                          type = "median"),
    stdMed_rsb_corrected = rse_bias(est_corrected,
                                    se_corrected,
                                    type = "median"),
    out_se = outlier_se(se),
    out_se_corrected = outlier_se(se_corrected),
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
    convergence = convergence_rate(convergences),
    warning_total = warning_sum(warnings)
  )
}

# ========================================= Run Experiment ========================================= #

res <- runSimulation(design = DESIGNFACTOR,
                     replications = 2000,
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
                     filename = "CorrectedSE_11112024",
                     parallel = TRUE,
                     ncores = 30,
                     save = TRUE,
                     save_results = TRUE)

# runSimulation(design = DESIGNFACTOR[1:2, ],
#               replications = 5,
#               generate = generate_dat,
#               analyse = list(joint = analyze_joint,
#                              gsam = analyze_gsam,
#                              lsam = analyze_lsam,
#                              tspa = analyze_tspa,
#                              rel = analyze_rel),
#               summarise = evaluate_res,
#               fixed_objects = FIXED_PARAMETER)
