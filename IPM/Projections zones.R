source("analysis/setup.R")

###############################################################################
## Read in the model: ----

# read in the model we'll use
stan_model <- readRDS(file = "models/hinge_zone_temp_sex_diff.rds")

# vector of names of required parameters
# source according to which model is being analysed
source("par_names/par_names_hinge_zone_temp_sex_diff.R")

# extract posterior draws of required parameters
par_posterior <- rstan::extract(stan_model, par_names) # GLOBAL

# calculate posterior mean of individual fish effects 
#this is extracting all the iterations and then the mean for each individual is calculated
#try and get a set of 1000 iterations by looping across

#In this the mean of the 4000 iterations for each individual is taken
#leaves a vector of mean values
u_fish <- apply(rstan::extract(stan_model, "u_fish")[[1]], 2, mean) 


#try <- rstan::extract(stan_model, "u_fish")[[1]]

#try <- as.data.frame(try)
#View(try)
#View(u_fish)
#glimpse(u_fish)
#summary(try)
# cleanup (remove the massive model object)
rm(stan_model); gc()


###############################################################################
## get all the state data into one place ----

state_var_names <- c("id_fish", "z0", "a", "is_f", "is_m", "is_EBS","is_ETAS","is_WTAS","is_NSW")
#this prep_stan_data is a function created in the setup.R document
#this is grabbing the appropriate data from fishdat_cut
#then it is taking a subset based on the variable names we have 
#just specified
all_states <- prep_stan_data(fishdat_cut)[state_var_names]
all_states <- as.data.frame(all_states)
#u_fish corresponds with the id_fish column in this dataset
all_states$u_fish <- u_fish[all_states$id_fish]


###############################################################################
## construct initial distribution ----

# correlation between size and individual effect across ages
u_z_cor <- 
  all_states %>% 
  group_by(a, is_f, is_m) %>% 
  summarise(u_z_cor = cor(u_fish, z0))
# inspect - should be +ve and go up with age...
u_z_cor


##### Change the zone initial distribution here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

all_states <- filter(all_states, is_EBS=="1")





# grad the sizes at initial age + the individual effects
z_init <- filter(all_states, a == min(a))$z0
u_init <- filter(all_states, a == min(a))$u_fish
# log size is roughly normally distributed--i.e. has approx. log-normal dist
car::qqPlot(log(z_init))
# the rf's should be normally distributed--right tail is too fat (hard to fix)
# we will just leave this the way it is--shouldn't make too much of a difference
car::qqPlot(u_init)
# raw relationship
plot(z_init, u_init, pch = 20, cex = 0.5)
# evidence for non-linear effect?
#polynomial regression
init_mod_p4 <- lm(u_init ~ poly(z_init, 4))
summary(init_mod_p4) # no!
# p-value only significant for the simple linear relationship
# so just use the simple regression
init_mod <- lm(u_init ~ z_init)
# residual distribution is a bit better (still not great)
car::qqPlot(resid(init_mod))

# GLOBAL -- used below
# list holding parameters describing initial distribution
init_dens_par <- list(
  # parameters of the initial size distribution
  meanlog_z = mean(log(z_init)),
  sdlog_z   = sd(log(z_init)),
  # parameters of the conditional distribution of u
  intercept    = coef(init_mod)[1],
  slope        = coef(init_mod)[2],
  sd_u         = sigma(init_mod)
)


###############################################################################
## Define functions to implement the model ----


# function to select parameter vector from posterior 
mk_m_par <- function(par_posterior, samp) {
  # extractor functions
  get_mean <- function(x) 
    if (length(dim(x)) == 1) mean(x) else apply(x, 2, mean)
  get_samp <- function(x) 
    if (length(dim(x)) == 1) x[samp] else x[samp,]
  # safety first
  if(samp > length(par_posterior[[1]]))
    stop("Insufficient number of samples")
  # return posterior mean (`samp` = 0) or one draw indexed by `samp`
  if (samp == 0) {
    return( lapply(par_posterior, get_mean) )
  } else {
    return( lapply(par_posterior, get_samp) )
  }
}

# function to draw the year effect
draw_u_year <- function(m_par, method = "constant") {
  if (method == "constant") {
    m_par$u_year <- 0.0
  } else stop ("method not recognised")
  return(m_par)
}
# for now, just sets the year effect to 0, i.e. constant environment case
# if I wanted to discover the effect of environmental stochatisity I would 
# change this

###############################################################################
# obtain the appropriate vital rate functions according to 
# the model that is being analysed
source("vitals/vital_rate_hinge_zone_temp_sex_diff.R")

###############################################################################
## Define the growth-switch kernel functions ----

# kernel when below the threshold
G_z1z_1_1 <- function (z1, z0, u, a, is_m, m_par) {
  (1 - p_sw(z1, z0, u, a, is_m, m_par)) * g_z1z_1(z1, z0, u, a, is_m, m_par)
}

# kernel after *first* crossing the threshold
G_z1z_2_1 <- function (z1, z0, u, a, is_m, m_par) {
  p_sw(z1, z0, u, a, is_m, m_par) * g_z1z_1(z1, z0, u, a, is_m, m_par)
}

# kernel after crossing the threshold
G_z1z_2_2 <- function (z1, z0, u, a, is_m, m_par) {
  g_z1z_2(z1, z0, u, a, is_m, m_par)
}

###############################################################################
## Define functions to implement the model ----

# calculate the mesh points, mesh width and store with their upper / lower 
# bounds of size and random effect + upper / lower bounds of age
mk_i_par <- function(N_z, L_z, U_z, N_u, L_u, U_u, L_a, U_a, sex) {
  # otolith size
  h_z <- (U_z - L_z) / N_z 
  m_z <-  L_z + (seq_len(N_z) - 0.5) * h_z
  # individual effect
  h_u <- (U_u - L_u) / N_u
  m_u <-  L_u + (seq_len(N_u) - 0.5) * h_u
  # iteration kernel nodes - !!! order of vars matters here !!! 
  nodes <- expand.grid(z1 = m_z, z0 = m_z, u = m_u)
  # return it all
  list(
    # size
    m_z = m_z, h_z = h_z, N_z = N_z, L_z = L_z, U_z = U_z,
    # individual
    m_u = m_u, h_u = h_u, N_u = N_u, L_u = L_u, U_u = U_u,
    # age
    L_a = L_a, U_a = U_a,
    # iteration kernel nodes
    z1 = nodes$z1, z0 = nodes$z0, u = nodes$u, 
    # sex indicator (for males)
    is_m = ifelse(sex == "M", 1, 0)
  )
}

# construct the three iteration kernels needed to project cohorts
mk_K <- function(i_par, m_par) {
  with(i_par, {
    # store the kernels in a list
    K_list <- list()
    # below the threshold
    K_list$K_z1z_1_1 <- G_z1z_1_1(z1, z0, u, a, is_m, m_par) * h_z
    # after *first* crossing the threshold
    K_list$K_z1z_2_1 <- G_z1z_2_1(z1, z0, u, a, is_m, m_par) * h_z
    # after crossing the threshold
    K_list$K_z1z_2_2 <- G_z1z_2_2(z1, z0, u, a, is_m, m_par) * h_z
    #
    return(K_list)
  })
}

# function to construct initial distribution
mk_init_dist <- function(i_par, m_par, init_dens_par) {
  # nodes on which we constrcut the initial distribution
  nodes <- expand.grid(z1 = i_par$m_z, z0 = NA, 
                       u = i_par$m_u, a = i_par$L_a, is_m = i_par$is_m)
  with(nodes, {
    # size * individual effect distribution
    init_dens <- 
      with(init_dens_par, {
        u_hat <- intercept + slope * z1
        dlnorm(z1, meanlog_z, sdlog_z) * dnorm(u, u_hat, sd_u)
      })
    # probability that the state lies pre-/post- threshold
    p <- p_sw(z1, z0, u, a, is_m, m_par)
    # return the distribution split into pre-/post- components
    return( list( (1-p) * init_dens, p * init_dens ) )
  })
}


# convenience function to construct ij indices to index into the block
# diagonal (sparse) iteration matrices
construct_K_index <- function(N_z, N_u) {
  expand.grid(
    i = seq_len(N_z), # row
    j = seq_len(N_z), # column
    b = seq_len(N_u)  # block
  ) %>% 
    mutate(
      i = i + (b-1) * N_z, 
      j = j + (b-1) * N_z
    )
} 

# function to iterate a cohort 
iterate <- function(i_par, m_par, init_dens_par) {
  # 
  index_set <- construct_K_index(i_par$N_z, i_par$N_u)
  # storage vectors
  nt_1 <- nt_2 <- list()
  # initial state
  nt_0 <- mk_init_dist(i_par, m_par, init_dens_par)
  nt_1[[1]] <- nt_0[[1]]
  nt_2[[1]] <- nt_0[[2]]
  # counters
  a <- i_par$L_a 
  i <- 1
  # keep going until we reach the maximum allowed age
  while(a < i_par$U_a) {
    # update the model parameters for current year
    m_par_use <- draw_u_year(m_par, method = "constant")
    # build the iteratiopn matrices
    K_set <- mk_K(i_par, m_par_use)
    # pre-threshold growth
    K_1_1 <- sparseMatrix(index_set$i, index_set$j, x = K_set$K_z1z_1_1)
    nt_1[[i+1]] <- (K_1_1 %*% nt_1[[i]])[,1]
    # post-threshold growth
    K_2_1 <- sparseMatrix(index_set$i, index_set$j, x = K_set$K_z1z_2_1)
    K_2_2 <- sparseMatrix(index_set$i, index_set$j, x = K_set$K_z1z_2_2)
    nt_2[[i+1]] <- (K_2_1 %*% nt_1[[i]])[,1] + (K_2_2 %*% nt_2[[i]])[,1]
    # increment counters
    a <- a + 1; i <- i + 1
  }
  list(nt_1 = nt_1, nt_2 = nt_2)
}


###############################################################################
## Using model ----

# grab the posterior mean model parameters
m_par <- mk_m_par(par_posterior, 0)
m_par$temp <- 1
m_par$zone <- "WTAS"
m_par$is_EBS = ifelse(m_par$zone == "EBS", 1, 0)
m_par$is_ETAS = ifelse(m_par$zone == "ETAS", 1, 0)
m_par$is_WTAS = ifelse(m_par$zone == "WTAS", 1, 0)
m_par$is_NSW = ifelse(m_par$zone == "NSW", 1, 0)


# set up the  parameters to control the numerics
i_par <- mk_i_par(
  N_z = 80, L_z =  0.10, U_z =  3.30, 
  N_u = 20,  L_u = -0.15, U_u = +0.15, 
  L_a = min(all_states$a), U_a = max(all_states$a),
  sex = "F" # <- code works on one sex at a time so must be scalar
)
# run the cohort dynamics
system.time(
  cohort_dynamics <- iterate(i_par, m_par, init_dens_par)
) # not super fast w/ above params 

# quick plot to check results...
plt_2d <- function(nt_1, nt_2) {
  plt <- function(n_t, i_par, zmax, col) {
    dim(n_t) <- c(i_par$N_z, i_par$N_u)
    image(i_par$m_z, i_par$m_u, n_t, col = col)
  }
  mz <- max(c(sapply(nt_1, max), sapply(nt_2, max)))
  par(mfcol = c(i_par$U_a - i_par$L_a + 1, 2), 
      mar = c(2, 2, 1, 1), oma = c(1, 1, 3, 1))
  col <- viridis(1000)
  sapply(nt_1, plt, i_par, mz, col)
  sapply(nt_2, plt, i_par, mz, col)
  mtext("  pre-threshold \n", side = 3, line = 0, outer = TRUE, adj = 0) 
  mtext(" post-threshold \n", side = 3, line = 0, outer = TRUE, adj = 1) 
}
plt_2d(cohort_dynamics$nt_1, cohort_dynamics$nt_2)   

# ... simpler to work with the data in 'tidy' format

# function to convert to 'tidy' format
tidy_output <- function(i_par, x) {
  # state vars
  out <- expand.grid(
    z = i_par$m_z,                     # size
    u = i_par$m_u,                     # individual effect
    a = seq.int(i_par$L_a, i_par$U_a), # age
    s = 1:2                            # stage
  )
  # density at state
  out$d <- unlist(c(x$nt_1, x$nt_2))
  return(out)
}

# tidy up the data
res <- tidy_output(i_par, cohort_dynamics)

# sanity check proportion in each age class should be = ~1
res %>% 
  group_by(a) %>% 
  summarise(p_stage = sum(d) * i_par$h_z * i_par$h_u)

# proportion by age / stage
res %>% 
  group_by(a, s) %>% 
  summarise(p_stage = sum(d) * i_par$h_z * i_par$h_u) 

# mean otolith size in each stage
res %>% 
  group_by(a) %>% 
  summarise(mean_z =  sum(d * z) * i_par$h_z * i_par$h_u) 

# or do all the first two multivariate moments
dd <- i_par$h_z * i_par$h_u 
res %>% 
  group_by(a) %>%
  summarise(
    # normalisation constant
    n      = sum(d) * dd,
    # moments
    mean_z = sum(d * z) * dd / n,
    mean_u = sum(d * u) * dd / n,
    var_z  = sum(d * (z - mean_z)^2) * dd / n,
    var_u  = sum(d * (u - mean_u)^2) * dd / n,
    cov_zu = sum(d * (u - mean_u) * (z - mean_z)) * dd / n,
    # derived quantities
    sd_z   = sqrt(var_z),
    sd_u   = sqrt(var_u),
    cor_zu = cov_zu / (sd_z * sd_u)
  ) 

# size distribution by age / stage
plt_data <- res %>% 
  group_by(a, s, z) %>% 
  summarise(d = sum(d) * i_par$h_u) 
ggplot(plt_data, aes(x = z, y = d, group = a)) +
  geom_line() + facet_wrap(~ s, nrow = 2) + 
  theme_minimal() 

# size distribution by age 
plt_data <- res %>% 
  group_by(a, z) %>% 
  summarise(d = sum(d) * i_par$h_u) 
ggplot(plt_data, aes(x = z, y = d, group = a)) +
  geom_line() + theme_minimal() 

