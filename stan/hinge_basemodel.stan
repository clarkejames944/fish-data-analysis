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
  // fixed intercept term (females)
  real<lower= 0.0, upper=2.0> alpha;
  // fixed threshold term (females)
  real<lower= 0.4, upper=2.0> eta;
  // intercept contrasts
  real d_alpha_male;
  // threshold contrasts
  real d_eta_male;
  // slope (relative to slope of 1)
  real<lower= 0.0, upper=1.0> beta_1;
  real<lower=-1.0, upper=0.0> beta_2;
  // sd + random intercept deviation due to individual fish
  real<lower=0.00, upper=1.00> sigma_fish;
  vector[n_fish] u_fish;
  // sd + random intercept deviation due to observation year
  real<lower=0.00, upper=1.00> sigma_year;
  vector[n_year] u_year;
  // residual standard deviation
  real<lower=0.00, upper=1.00> sigma;
}

transformed parameters {
    vector[n_obs] z1hat;
    vector[n_obs] IN;
    vector[n_obs] TR;

    for (i in 1:n_obs) {
      // intercepts
      IN[i] = alpha + d_alpha_male*is_m[i] + u_fish[id_fish[i]] + u_year[id_year[i]] ;
      // thresholds 
      TR[i] = eta + d_eta_male*is_m[i];
      // predicted otolith size next year
      z1hat[i] = (IN[i] + (1 + beta_1) * (z0[i] - TR[i])) / (1 + exp(+(z0[i] - TR[i]) / 0.02)) + 
                 (IN[i] + (1 + beta_2) * (z0[i] - TR[i])) / (1 + exp(-(z0[i] - TR[i]) / 0.02));
    }
}

model { 
  // intercept terms
  alpha ~ normal(1.0, 0.25);
  d_alpha_male ~ normal(0.0, 0.25);
  // threshold terms
  eta ~ normal(1.0, 0.5);
  d_eta_male ~ normal(0.0, 0.25);
  // slope terms
  beta_1 ~ normal(0, 0.25);
  beta_2 ~ normal(0, 0.25);
  // random intercept deviation due to individual
  sigma_fish ~ cauchy(0.0, 0.25);
  u_fish ~ normal(0.0, sigma_fish);
  // random intercept deviation due to year
  sigma_year ~ cauchy(0.0, 0.25);
  u_year ~ normal(0.0, sigma_year);
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
