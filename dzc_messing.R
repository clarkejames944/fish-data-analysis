library(tidyverse)

## 0. preamble ----

# box-cox transformation with translation
bct <- function(y, l1, l2 = 0) {
  if (l1 == 0)
    log(y + l2)
  else
    ((y - l2) ^ l1 - 1) / l1
}
# pick the subset we are working with
site_sex_data <- EBSf
# common transformations
trans_data <- site_sex_data %>% 
  mutate(year_f = factor(Year))

## 2. experiment with b-c transformations ----

lam <- 1.6

trans_data <- site_sex_data %>%
  mutate(
    oto_size = bct(oto_size, l1 = lam, l2 = 0),
    prev = bct(prev, l1 = lam, l2 = 0)
  )

# fit a flexible growth model *by maturation status 
# N.B. this is just an arbitrary age threshold (>5 == 'adult')
grow_mod <-
  gam(
    oto_size ~ prev + s(prev, k = 50, by = maturity),
    #+ s(year_f, bs = "re"),
    family = gaussian(link = "identity"),
    data = trans_data
  )
# check the choice of k (needs to be higher than the default)
gam.check(grow_mod)
# quick look at the terms and R^2
summary(grow_mod)
# extract fitted values and (r)esponse (i.e. raw) residuals
trans_data <- trans_data %>%
  mutate(
    preds   =   predict(grow_mod, type = "response"),
    resid_r = residuals(grow_mod, type = "response")
  )
# now examine the empirical mean-variance relationship
ggplot(trans_data, aes(x = preds, y = log(resid_r ^ 2))) +
  geom_point(aes(colour = maturity), size = 1, alpha = 0.1) +
  geom_smooth(method = 'gam',
              formula = y ~ s(x, k = 100), se = FALSE) +
  geom_smooth(method = 'lm',
              se = FALSE, linetype = 2)
# calculate scaled residuals based on m-v relationship
var_mod <- gam(log(resid_r ^ 2) ~ s(preds), data = trans_data)
trans_data <- mutate(trans_data,
                     sc_resid = resid_r / sqrt(exp(predict(var_mod))))
# look at the distributional assumps. w/ scaled resids
with(trans_data, car::qqp(sc_resid[maturity == "juvenile"])) # less 'normal'
with(trans_data, car::qqp(sc_resid[maturity == "adult"])) # more 'normal'
# plot the growth data and the fitted values
ggplot(trans_data, aes(x = prev)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(y = oto_size, colour = maturity),
             size = 1, alpha = 0.1) +
  geom_point(aes(y = preds), size = 0.25, colour = "steelblue") +
  ylab("s'") + xlab("s")
