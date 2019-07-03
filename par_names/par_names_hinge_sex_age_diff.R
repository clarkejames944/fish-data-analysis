# vector of names of required parameters
par_names <- 
  c(# intercept parameters
    "alpha", "d_alpha_male", "d_alpha_age",
    # threshold size
    "eta", "d_eta_male", "d_eta_age",
    # pre-threshold slope
    "beta_1", "d_beta_1_male", "d_beta_1_age",
    # post-threshold slope
    "beta_2", "d_beta_2_male", "d_beta_2_age",
    # variance components
    "sigma_fish", "sigma_year", "sigma",
    # 'random' year effects
    "u_year"
  )