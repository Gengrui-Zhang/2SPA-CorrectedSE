library(lavaan) # for the PoliticalDemocracy data
library(blavaan)

dat <- generate_dat(condition)
cfa_mod <- syntax_update(condition$n_items, path_include = FALSE)
fs <- get_fs(dat, 
             model = cfa_mod,
             method = "Bartlett", 
             vfsLT = TRUE)

model <- '
  # Measurement model
  x =~ 1*fs_fx
  y =~ 1*fs_fy
  
  # Error variance variance of fx_fs
  fs_fx ~~ 0.01117451 * fs_fx
  fs_fy ~~ 0.0649131 * fs_fy
  
  # Structural model
  y ~ x
'

# Blavaan
summary(bsem(model,
    data = fs))

# lavaan
summary(sem(model,
            data = fs))

