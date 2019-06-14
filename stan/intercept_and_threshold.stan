data {
  // number of observations and group sizes
  int<lower=0> n_obs;
  int<lower=0> n_fish;
  int<lower=0> n_year;
  // indexing vectors
  int<lower=1, upper=n_fish> id_fish[n_obs];
  int<lower=1, upper=n_year> id_year[n_obs];
  // indicator variables
  int<lower=0, upper=1> is_f[n_obs];
  int<lower=0, upper=1> is_m[n_obs];
  // successive otolith sizes
  vector<lower=0>[n_obs] z0;
  vector<lower=0>[n_obs] z1;
}

parameters {
  // fixed intercept term (females, pre-maturation)
  real alpha;
  // intercept contrasts
  real d_alpha_male;
  real d_alpha_mature;
  // slope contrats (relative to slope of 1)
  real<lower= 0.0, upper=1.0> beta_1;
  real<lower=-1.0, upper=0.0> beta_2;
  // ****
  real<lower=-10.0, upper=10.0> gamma;
  // sd + random intercept deviation due to individual fish
  real<lower=0.00, upper=1.00> sigma_fish;
  vector[n_fish] u_fish;
  // sd + random intercept deviation due to observation year
  real<lower=0.00, upper=1.00> sigma_year;
  vector[n_year] u_year;
  // fixed threshold for growth switch
  real<lower=0.25, upper=1.75> eta; 
  // residual standard deviation
  real<lower=0.00, upper=1.00> sigma;
}

transformed parameters {
    vector[n_obs] tau;
    vector[n_obs] z1hat;
    vector[n_obs] i1_f;
    vector[n_obs] i1_m;
    vector[n_obs] i2_f;
    vector[n_obs] i2_m;

    for (i in 1:n_obs) {
      // thesholding function: -->1 when past the tau threshold, -->0 otherwise
      tau[i] = 1 / (1 + exp(-20 * (z0[i] - eta - gamma * u_fish[id_fish[i]])));
      // fish/year-specific intercepts
      i1_f[i] = alpha + u_fish[id_fish[i]] + u_year[id_year[i]];
      i1_m[i] = i1_f[i] + d_alpha_male; 
      i2_f[i] = i1_f[i] + d_alpha_mature;
      i2_m[i] = i1_f[i] + d_alpha_male + d_alpha_mature;
      // predicted otolith size next year
      z1hat[i] = (i1_f[i]*is_f[i] + i1_m[i]*is_m[i] + z0[i] + beta_1*z0[i]) * (1 - tau[i]) + 
                 (i2_f[i]*is_f[i] + i2_m[i]*is_m[i] + z0[i] + beta_2*z0[i]) * tau[i];
    }
}

model { 
  // intercept terms
  alpha ~ normal(0.25, 0.1);
  d_alpha_male ~ normal(0.0, 0.1);
  d_alpha_mature ~ normal(0.0, 0.1);
  // slope terms
  beta_1 ~ normal(0, 0.1);
  beta_2 ~ normal(0, 0.1);
  // ****
  gamma ~ normal(0, 3.0);
  // random intercept deviation due to individual
  sigma_fish ~ cauchy(0.0, 0.25);
  u_fish ~ normal(0.0, sigma_fish);
  // random intercept deviation due to year
  sigma_year ~ cauchy(0.0, 0.25);
  u_year ~ normal(0.0, sigma_year);
  // threshold for growth switch (weakly informative, data-derived prior)
  eta ~ normal(1.0, 0.2);
  // data likelihood
  sigma ~ cauchy(0, 0.25);
  z1 ~ normal(z1hat, sigma);
}

generated quantities {
    vector[n_obs] sim_z1;
    vector[n_obs] log_lik;
    
    for (i in 1:n_obs) {
      sim_z1[i] = normal_rng(z1hat[i], sigma);
    }

    for (i in 1:n_obs) {
      log_lik[i] = normal_lpdf(z1[i] | z1hat[i], sigma); 
    }
}
