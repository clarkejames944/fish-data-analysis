## Define the vital rate functions ----

# growth function: *pre*-threshold
# --> given you are size z and age a now, returns the  pdf of size next
g_z1z_1 <- function(z1, z0, u, a, is_m, m_par) {
  with(m_par, {
    # intercept
    IN = alpha + d_alpha_male * is_m + d_alpha_temp * (temp - mean(temp)) + d_alpha_zone * zone  + u + u_year
    # threshold
    TR = eta + d_eta_male * is_m + d_eta_temp * (temp - mean(temp)) + d_eta_zone * zone
    # pre-threshold slope
    B1 = beta_1 + d_beta_1_male * is_m + d_beta_1_temp * (temp - mean(temp)) + d_beta_1_zone * zone
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
    IN = alpha + d_alpha_male * is_m + d_alpha_temp * (temp - mean(temp)) + d_alpha_zone * zone + u + u_year
    # threshold
    TR = eta + d_eta_male * is_m + d_eta_temp * (temp - mean(temp)) + d_eta_zone * zone
    # post-threshold slope
    B2 = beta_2 + d_beta_2_male * is_m + d_beta_2_temp * (temp - mean(temp)) + d_beta_2_zone * zone
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
    TR = eta + d_eta_male * is_m + d_eta_temp * (temp - mean(temp))+ d_eta_zone * zone
    # return the switch probability
    1 / (1 + exp(-(z1 - TR) / 0.02))
  })
}