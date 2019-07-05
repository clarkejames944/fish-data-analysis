# vector of names of required parameters
par_names <- 
  c(# intercept parameters
    "alpha", "d_alpha_male", "d_alpha_EBS", "d_alpha_ETAS", "d_alpha_WTAS", "d_alpha_temp",
    # threshold size
    "eta", "d_eta_male", "d_eta_EBS", "d_eta_ETAS", "d_eta_WTAS", "d_eta_temp",
    # pre-threshold slope
    "beta_1", "d_beta_1_male", "d_beta_1_EBS", "d_beta_1_ETAS", "d_beta_1_WTAS", "d_beta_1_temp",
    # post-threshold slope
    "beta_2", "d_beta_2_male", "d_beta_2_EBS", "d_beta_2_ETAS", "d_beta_2_WTAS", "d_beta_2_temp",
    # variance components
    "sigma_fish", "sigma_year", "sigma",
    # 'random' year effects
    "u_year"
  )