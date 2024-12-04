//
// This Stan program defines a simple regression model
// with known measurement error in 'x'. 
//
data {
    int<lower=0> N;  // number of observations
    vector[N] fs_fy;  // latent outcome
    vector[N] fs_fx;  // latent predictor
    real<lower=0> fs_fx_se;  // reliability of fs_fx [Should be standard error]
    real<lower=0> fs_fy_se;  // reliability of fs_fy
}
parameters {
    real beta0;  // intercept for regression
    real beta1;  // slope for regression
    real<lower=0> sigma;  // residual standard deviation for fy
    real<lower=0> sd_fx;    // SD of fx
    vector[N] zfx;  // true latent predictor (standardized)
    vector[N] fy;  // true latent outcome
}
transformed parameters {
    vector[N] fx = sd_fx .* zfx;  // true latent predictor
}
model {
    zfx ~ std_normal();  // Prior for latent predictor
    fy ~ normal(beta0 + beta1 * fx, sigma);
    fs_fx ~ normal(fx, fs_fx_se);
    fs_fy ~ normal(fy, fs_fy_se);
}
generated quantities {
    real vfy = beta1^2 * sd_fx^2 + sigma^2;
    real std_beta1 = beta1 * sd_fx / sqrt(vfy);
}
