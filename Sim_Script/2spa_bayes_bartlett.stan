//
// This Stan program defines a simple regression model
// with known measurement error in 'x'. 
//
data {
    int<lower=0> N;  // number of observations
    vector[N] fs_fy;  // latent outcome
    vector[N] fs_fx;  // latent predictor
    vector<lower=0>[N] fs_fx_se;  // reliability of fs_fx [Should be standard error]
    vector<lower=0>[N] fs_fy_se;  // reliability of fs_fy
}
parameters {
    real beta0;  // intercept for regression
    real beta1;  // slope for regression
    real<lower=0> sigma;  // residual standard deviation for fy
    vector[N] fx;  // true latent predictor
    vector[N] fy;  // true latent outcome
}
model {
    fx ~ std_normal();  // Prior for latent predictor
    fs_fx ~ normal(fx, fs_fx_se);
    fs_fy ~ normal(fy, fs_fy_se);
    fy ~ normal(beta0 + beta1 * fx, sigma);
}
