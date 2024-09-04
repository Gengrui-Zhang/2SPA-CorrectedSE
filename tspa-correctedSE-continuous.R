seed_value <- 66234
set.seed(seed_value)

library(lavaan)
library(SimDesign)
library(dplyr)
library(here)

# Source Function
r_scripts <- list.files(here("R"), pattern = "\\.R$", full.names = TRUE)
lapply(r_scripts, source)

# Helper Function 1
# Adjusting a parameter in a model
# Updating the model with the new parameter values
# Computing the objective function
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
std_est <- function(ev, obj) {
  fsT1 <- attr(obj, which = "fsT")
  diag(fsT1) <- ev
  mobj <- tspa(" fy ~ fx",
               data = fs,
               fsT = fsT1, fsL = attr(fs, "fsL")
  )
  standardizedSolution(mobj, type = "std.lv")$est.std[8]
}

# Helper Function 4
compute_std_est <- function(ep, mod) {
  mod_tmp <- gsub("ep2", ep[2], mod)
  mod_tmp <- gsub("ep1", ep[1], mod_tmp)
  standardizedSolution(sem(mod_tmp, data = fs))$est[6:7]
}

# Helper Function 5
# Data Generation
generate_data <- function(num_obs, Lambda, Psi, Theta) {
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

# Example
# Lambda <- matrix(c(.8, .5, .6, rep(0, 6), .6, .7, .5), ncol = 2)
# Psi <- diag(c(1, 1))  
# Theta <- diag(c(1.33, 1.22, 1.11, 1.30, 1.20, 1.10)) 
# dat <- generate_data(num_obs = 100, Lambda = Lambda, Psi = Psi, Theta = Theta)
# head(dat)

# Design Factor
DESIGNFACTOR <- createDesign(
  N = c(100, 500),
  rel = c(0.7, 0.9)
  # num_items = c(3, 10) # Need to change. Now it's hardcoding
)

# Data Generation
GenData <- function (condition) {
  N <- condition$N
  rel = condition$rel
  # num_items = condition$num_items
  
  # Model Parameters
  Alpha <- c(0, 0)
  Lambda_1 <- c(.8, .5, .6)
  Lambda_2 <- c(.6, .7, .5)
  Lambda <- cbind(c(Lambda_1, rep(0, 3)),
                  c(rep(0, 3), Lambda_2))
  Psi <- diag(c(1, 1))
  # Compute error variance
  Error_1 <- sum(Lambda_1)^2*(1 - rel)/rel*c(0.44, 0.33, 0.23)
  Error_2 <- sum(Lambda_2)^2*(1 - rel)/rel*c(0.44, 0.33, 0.23)
  Theta <- diag(c(Error_1, Error_2)) 
  
  # Simulate Data
  generate_data(num_obs = N, Lambda = Lambda, Psi = Psi, Theta = Theta)
}

# Example
# GenData(condition = DESIGNFACTOR[1,])

# Parameter Extraction
extract_res <- function(condition, dat) {
  # Joint Model
  mod <- "
  fy =~ y1 + y2 + y3
  fx =~ x1 + x2 + x3
  fy ~ fx
    "
  mjoint <- sem(mod, data = dat)
  joint_est <- standardizedSolution(mjoint)[7, ]
  
  # SAM
  ## Global
  mgsam <- sam(mod, data = dat, sam.method = "global")
  gsam_est <- standardizedSolution(mgsam)[7, ]
  ## Local
  mlsam <- sam(mod, data = dat)
  lsam_est <- standardizedSolution(mlsam)[7, ]
  
  # 2S-PA
  ## Error
  cfa_mod <- "
  fy =~ y1 + y2 + y3
  fx =~ x1 + x2 + x3
    "
  cfa_fit <- cfa(cfa_mod,
                 std.lv = TRUE, data = dat,
                 bounds = "pos.ov.var"
  )
  fs <- get_fs_lavaan(cfa_fit, method = "Bartlett", vfsLT = TRUE)
  m2spa1 <- tspa(" fy ~ fx",
                 data = fs,
                 fsT = attr(fs, "fsT"), fsL = attr(fs, "fsL")
  )
  tspa1_est <- standardizedSolution(m2spa1)[8, ]
  ### Corrected Standard Errors
  v1 <- attr(fs, "vfsLT")[c(5, 7), c(5, 7)]
  # Gradient
  grad_std <- numDeriv::grad(std_est,
                             x = diag(attr(fs, "fsT")), obj = m2spa1
  )
  vc_corrected2 <- lavInspect(m2spa1, what = "vcov.std")[1, 1] +
    t(grad_std) %*% v1 %*% grad_std
  
  ## Reliability
  mod2 <- "
  fy =~ NA * fs_fy + ly * fs_fy
  fx =~ NA * fs_fx + lx * fs_fx
  fs_fy ~~ ev1 * fs_fy
  fs_fx ~~ ev2 * fs_fx
  fx ~~ 1 * fx
  fy ~ beta * fx
  fy ~~ dvy * fy
  # constraints
  ev2 == ep2 * lx^2
  1 == beta^2 + dvy
  ev1 == ep1 * ly^2
    "
  ep <- diag(attr(fs, "fsT"))
  mod_tmp <- gsub("ep2", ep[2], mod2)
  mod_tmp <- gsub("ep1", ep[1], mod_tmp)
  m2spa2 <- sem(mod_tmp, data = fs)
  tspa2_est <- standardizedSolution(m2spa2)[6, ]
  
  ## Corrected standard errors
  # Jacobian
  jac_std <- numDeriv::jacobian(compute_std_est,
                                x = diag(attr(fs, "fsT")), mod = mod2
  )
  vc_corrected3 <- lavInspect(m2spa2, what = "vcov.std")[5:6, 5:6] +
    jac_std %*% v1 %*% t(jac_std)
}

# Example
# extract_res(condition = DESIGNFACTOR[1,], dat = GenData(condition = DESIGNFACTOR[1,]))
