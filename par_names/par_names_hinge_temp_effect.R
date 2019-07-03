# vector of names of required parameters
par_names <- 
  c(# intercept parameters
    "alpha", "d_alpha_male", "d_alpha_temp",
    # threshold size
    "eta", "d_eta_male", "d_eta_temp",
    # pre-threshold slope
    "beta_1", "d_beta_1_male", "d_beta_1_temp",
    # post-threshold slope
    "beta_2", "d_beta_2_male", "d_beta_2_temp",
    # variance components
    "sigma_fish", "sigma_year", "sigma",
    # 'random' year effects
    "u_year"
  )