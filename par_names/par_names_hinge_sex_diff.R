# vector of names of required parameters
par_names <- 
  c(# intercept parameters
    "alpha", "d_alpha_male",
    # threshold size
    "eta", "d_eta_male",
    # pre-threshold slope
    "beta_1", "d_beta_1_male",
    # post-threshold slope
    "beta_2", "d_beta_2_male",
    # variance components
    "sigma_fish", "sigma_year", "sigma",
    # 'random' year effects
    "u_year"
  )