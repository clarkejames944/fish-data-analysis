# vector of names of required parameters
par_names <- 
  c(# intercept parameters
    "alpha", "d_alpha_male", "d_alpha_temp", "d_alpha_zone",
    # threshold size
    "eta", "d_eta_male", "d_eta_temp", "d_eta_zone",
    # pre-threshold slope
    "beta_1", "d_beta_1_male", "d_beta_1_temp", "d_beta_1_zone",
    # post-threshold slope
    "beta_2", "d_beta_2_male", "d_beta_2_temp", "d_beta_2_zone",
    # variance components
    "sigma_fish", "sigma_year", "sigma",
    # 'random' year effects
    "u_year"
  )