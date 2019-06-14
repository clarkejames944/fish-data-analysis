source("analysis/setup.R")

#######################################################################
## Exploratory models (again) ----

fit_grow_mod <- function(gam_form, mod_data) {
  # sex-specific datasets
  mod_data_f <- filter(mod_data, sex == "F")
  mod_data_m <- filter(mod_data, sex == "M")
  # sex-specific models
  mod_f <- gam(gam_form, data = mod_data_f)
  mod_m <- gam(gam_form, data = mod_data_m)
  # return the models
  list(mod_f = mod_f, mod_m = mod_m)
}

plot_models <- function(models, xlim = c(0, 6), ylim = c(0.0, 0.8)) {
  par(mfrow = c(1,2))
  # plot predictions from female model
  z1_pred <- predict(models$mod_f, exclude = c("s(FishID)", "s(Year_f)" ))
  plot(
    z1_pred - z0 ~ z0, data = models$mod_f$model,
    col = "red", pch = 20, cex = 0.5, xlim = xlim, ylim = ylim,
    main = "Females", xlab = "Size (z)", ylab = "Increment (z' - z)"
  )
  # plot predictions from male model
  z1_pred <- predict(models$mod_m, exclude = c("s(FishID)", "s(Year_f)" ))
  plot(
    z1_pred - z0 ~ z0, data = models$mod_m$model,
    col = "blue", pch = 20, cex = 0.5, xlim = xlim, ylim = ylim,
    main = "Males", xlab = "Size (z)", ylab = "Increment (z' - z)"
  )
}


a <- 2.0
mod_data <- fish %>% 
  mutate(z0 = z0^a, z1 = z1^a) %>% 
  filter(Age >= 2, Age <= 99)

gam_form_1 <- z1 ~ s(z0, k = 25) + s(Age) + s(z0, k = 25, by = Age)
g_mods_1 <- fit_grow_mod(gam_form_1, mod_data)
# save(gam_mods_1, file = "g_mods_1.rda")

gam_form_2.1 <- z1 ~ s(z0, k = 25) + s(Age) + s(z0, k = 25, by = Age) + s(FishID, bs = "re") + s(Year_f, bs = "re")
g_mods_2.1 <- fit_grow_mod(gam_form_2.1, mod_data)
# save(gam_mods_2.1, file = "g_mods_2.1.rda")

gam_form_2 <- z1 ~ s(z0, k = 25, by = Age_f) + s(FishID, bs = "re") + s(Year_f, bs = "re")
g_mods_2 <- fit_grow_mod(gam_form_2, mod_data)
# save(gam_mods_2, file = "g_mods_2.rda")

gam_form_3 <- z1 ~ s(z0, k = 25) + s(FishID, bs = "re") + s(Year_f, bs = "re")
g_mods_3 <- fit_grow_mod(gam_form_3, mod_data)
# save(gam_mods_3, file = "g_mods_3.rda")

gam_form_4 <- z1 ~ s(z0, k = 25)
g_mods_4 <- fit_grow_mod(gam_form_4, mod_data)

load("models/gam_mods_1.rda")
load("models/gam_mods_2.rda")
load("models/gam_mods_2.1.rda")
load("models/gam_mods_3.rda")

plot_models(g_mods_1, ylim = c(0.1, 0.45))
plot_models(g_mods_2, ylim = c(0.1, 0.45))
plot_models(g_mods_3, ylim = c(0.1, 0.45))
plot_models(g_mods_4, ylim = c(0.1, 0.45))


