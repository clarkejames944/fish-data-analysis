## Define the vital rate functions ----

# growth function: *pre*-threshold
# --> given you are size z and age a now, returns the  pdf of size next
g_z1z_1 <- function(z1, z0, u, a, is_m, m_par) {
  with(m_par, {
    # intercept
    IN = alpha + d_alpha_male * is_m + d_alpha_EBS * is_EBS + d_alpha_ETAS * is_ETAS + d_alpha_WTAS * is_WTAS + d_alpha_temp * temp  + u + u_year
    # threshold
    TR = eta + d_eta_male * is_m + d_eta_EBS * is_EBS + d_eta_ETAS * is_ETAS + d_eta_WTAS * is_WTAS + d_eta_temp * temp 
    # pre-threshold slope
    B1 = beta_1 + d_beta_1_male * is_m + d_beta_1_EBS * is_EBS + d_beta_1_ETAS * is_ETAS + d_beta_1_WTAS * is_WTAS + d_beta_1_temp * temp
    # expected otolith size next year
    z1_hat = IN + (1 + B1) * (z0 - TR)
    # return the density function
    dnorm(z1, mean = z1_hat, sd = sigma)
  })
}

# growth function, *post*-threshold
# --> given you are size z and age a now, returns the  pdf of size next
g_z1z_2 <- function(z1, z0, u, a, is_m, m_par) {
  with(m_par, {
    # intercept
    IN = alpha + d_alpha_male * is_m + d_alpha_EBS * is_EBS + d_alpha_ETAS * is_ETAS + d_alpha_WTAS * is_WTAS + d_alpha_temp * temp + u + u_year
    # threshold
    TR = eta + d_eta_male * is_m + d_eta_EBS * is_EBS + d_eta_ETAS * is_ETAS + d_eta_WTAS * is_WTAS + d_eta_temp * temp 
    # post-threshold slope
    B2 = beta_2 + d_beta_2_male * is_m + d_beta_2_EBS * is_EBS + d_beta_2_ETAS * is_ETAS + d_beta_2_WTAS * is_WTAS + d_beta_2_temp * temp
    # expected otolith size next year
    z1_hat = IN + (1 + B2) * (z0 - TR)
    # return the density function
    dnorm(z1, mean = z1_hat, sd = sigma)
  })
}

# switch function --> given that you are size z and age a now, returns the 
# probability you will switch growth trajectories
p_sw <- function(z1, z0, u, a, is_m, m_par) {
  with(m_par, {
    # threshold
    TR = eta + d_eta_male * is_m + d_eta_EBS * is_EBS + d_eta_ETAS * is_ETAS + d_eta_WTAS * is_WTAS + d_eta_temp * temp 
    # return the switch probability
    1 / (1 + exp(-(z1 - TR) / 0.02))
  })
}



