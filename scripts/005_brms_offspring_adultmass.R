###############################################################################
##  Script: Effect of parental age on Offspring adult mass
##  Note:   Analysis restricted to offspring that survived to adulthood
###############################################################################


#===========================================================================
# 1. SETUP------------------------------------------------------------------
# ===========================================================================

## Libraries ------------------------------------------------------------------
library(dplyr) #v.1.1.4
library(tidyr) #v.1.3.1
library(ggplot2) #v.4.0.0
library(patchwork) #v.1.3.0
library(ggpubr) #v.0.6.1
library(ggdist) #v.3.3.3
library(emmeans) #v.1.10.1
library(brms) #v.2.21.0
library(bayesplot) #v.1.11.1
library(bayestestR) #v.0.17.0
library(tidybayes) #v.3.0.7
library(parameters) #v.0.28.2
library(sjPlot)  #v.2.8.16
library(car) #v.3.1.2
library(loo)  #v.2.8.0
library(lme4) #v.1.1.35.2
library(gt) #v.10.1
library(gtsummary) #v.1.7.2

## Stan backend ---------------------------------------------------------------
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

## Colour scheme for bayesplot ------------------------------------------------
color_scheme_set("teal")


# ===========================================================================
# 2. DATA ORGANISATION--------------------------------------------------------
# ===========================================================================

## Load and filter data -------------------------------------------------------
F1data <- readRDS("./raw data/F1_filtered_data06022025.RDS")
length(unique(F1data$F1_ID))  # 1316 offspring

# Retain only adults with a recorded mass
data1 <- F1data %>%
  filter(Include_in_adult == "Y") %>%
  filter(adult_surv == 1) %>%
  filter(!is.na(F1_adultmass))

## Scale continuous predictors ------------------------------------------------
# All predictors scaled to mean = 0, SD = 1; a 1-unit change reflects 1 SD
data1$avg_age <- as.numeric(scale(data1$avg.age, center = TRUE, scale = TRUE))

# Within-individual (delta) age: scaled but NOT mean-centred (as variable is already centre)
# the within-subject structure; NAs replaced with 0
data1 <- data1 %>%
  group_by(PairID) %>%
  mutate(
    delta_age_scaled = as.numeric(scale(within_subject_age,
                                        center = FALSE,
                                        scale = TRUE)),
    delta_age_scaled = replace(delta_age_scaled, is.na(delta_age_scaled), 0)
  ) %>%
  ungroup()

data1$Mother_bodymass_scaled <- as.numeric(
  scale(data1$Mother_bodymass, center = TRUE, scale = TRUE)
)

# Rename to match model code throughout
data1 <- data1 %>% rename(delta.age = delta_age_scaled)

## Sample sizes ---------------------------------------------------------------
length(unique(data1$Mother_ID))  # number of mothers/parent pairs (77)
length(unique(data1$F1_ID))      # number of F1 offspring (937)


# ===========================================================================
# 3. SUMMARY STATISTICS------------------------------------------------------
# ===========================================================================

## Distribution of F1 adult mass (all individuals) ---------------------------
ggplot(F1data, aes(x = F1_adultmass)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", colour = "black") +
  labs(x = "F1 adult body mass (g)", y = "Count") +
  theme_classic()

## Overall mean and SD (used to build intercept priors later)-----------------------
sum_overall <- data1 %>%
  summarise(
    F1_mass     = mean(F1_adultmass),
    sd_total    = sd(F1_adultmass),
    n_offspring = n()
  )

range(data1$F1_adultmass)  # 0.277–1.539 g

## Mean F1 mass by parental temperature treatment ----------------------------
sum_by_temp <- data1 %>%
  group_by(Temp) %>%
  summarise(
    F1_mass     = mean(F1_adultmass),
    sd_total    = sd(F1_adultmass),
    n_offspring = n(),
    n_parents   = n_distinct(PairID)
  )

## Mean F1 mass by parental age timepoint ------------------------------------
sum_by_timepoint <- data1 %>%
  group_by(Timepoint) %>%
  summarise(
    F1_mass     = mean(F1_adultmass),
    sd_total    = sd(F1_adultmass),
    n_offspring = n(),
    n_parents   = n_distinct(PairID)
  )

## Mean F1 mass by delta age -------------------------------------------------
sum_by_deltaage <- data1 %>%
  group_by(delta.age) %>%
  summarise(
    F1_mass     = mean(F1_adultmass),
    sd_total    = sd(F1_adultmass),
    n_offspring = n(),
    n_parents   = n_distinct(PairID)
  )
## Collinearity checks --------------------------------------------------------
X <- model.matrix(
  ~ avg.age + delta.age + Mother_bodymass_scaled + Temp + F1_sex,
  data = data.frame(
    avg.age                       = data1$avg.age,
    delta.age                     = data1$delta.age,
    Mother_bodymass_scaled        = data1$Mother_bodymass_scaled,
    Temp                          = data1$Temp,
    F1_sex                        = data1$F1_sex
  )
)

cor(X[, -1])  # correlation matrix (excluding intercept column)

vif(lm(F1_adultmass ~ avg_age + delta.age +
        + Mother_bodymass_scaled + Temp + F1_sex,
       data = data1))


# ===========================================================================
# 4. MODEL STRUCTURE SELECTION-----------------------------------------------
# ===========================================================================
# Compare Gaussian family models varying in random-effects structure and
# whether a quadratic delta-age term is needed (or if this is too complex)
#In these models, using default priors. Prior specification occurs later

## 4.1 Baseline: random intercept only, linear delta age --------------------
mod1.1_basegaussian <- brm(
  F1_adultmass ~ avg_age + delta.age + Mother_bodymass_scaled + Temp + F1_sex +
    (1 | PairID),
  family  = gaussian,
  data    = data1,
  iter    = 5000,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE),
  cores   = 4
)
saveRDS(mod1.1_basegaussian,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_basegaussian.rda")
mod1.1_basegaussian <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_basegaussian.rda")

#model diagnostics
loo_basegaussian <- loo(mod1.1_basegaussian, save_psis = TRUE, moment_match = TRUE)
plot(loo_basegaussian)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_basegaussian, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(mod1.1_basegaussian, type ="loo_pit_overlay", ndraws=100) 


## 4.2 Quadratic delta-age term (orthogonal quadratic)-----------------------------
mod1.1_quadraticage <- brm(
  F1_adultmass ~ avg_age + 
    poly(delta.age,2) +
    Mother_bodymass_scaled +
    Temp + 
    F1_sex + 
    (1 | PairID),
  family  = gaussian,
  data    = data1,
  iter    = 5000,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE),
  cores   = 4
)
saveRDS(mod1.1_quadraticage,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_quadraticage.rda")
mod1.1_quadraticage <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_quadraticage.rda")

#model diagnostics
loo_quadraticage <- loo(mod1.1_quadraticage, save_psis = TRUE, moment_match = TRUE)
plot(loo_quadraticage)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_quadraticage, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(mod1.1_quadraticage, type ="loo_pit_overlay", ndraws=100) 


## 4.3 Random slopes (intercept + delta age) --------------------------------
#random slopes fitted with linear delta age term
mod1.1_randomslopes <- brm(
  F1_adultmass ~ avg_age + 
    delta.age + 
    Temp +
    Mother_bodymass_scaled + 
    F1_sex +
    (1 + delta.age | PairID),
  family  = gaussian,
  data    = data1,
  iter    = 5000,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE),
  cores   = 4
)
saveRDS(mod1.1_randomslopes,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_randomslopes.rda")
mod1.1_randomslopes <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_randomslopes.rda")

#model diagnostics
loo_randomslopes <- loo(mod1.1_randomslopes, save_psis = TRUE, moment_match = TRUE)
plot(loo_randomslopes)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_randomslopes, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(mod1.1_randomslopes, type ="loo_pit_overlay", ndraws=100) 


## 4.4 Random slopes + quadratic delta-age ----------------------------------
#random slopes fitted with quadratic population level delta age effect
#random slopes themselves are not allowed to change as a quadratic function. Assumes a shared shape of trajectory across individuals
mod1.1_quadraticslopes <- brm(
  F1_adultmass ~ avg_age + 
    poly(delta.age,2)+
    Temp + Mother_bodymass_scaled +
    F1_sex +
    (1 + delta.age | PairID),
  family  = gaussian,
  data    = data1,
  iter    = 5000,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE),
  cores   = 4
)
saveRDS(mod1.1_quadraticslopes,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_quadraticslopes.rda")
mod1.1_quadraticslopes <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_quadraticslopes.rda")

#model diagnostics
loo_quadraticslopes <- loo(mod1.1_quadraticslopes, save_psis = TRUE, moment_match = TRUE)
plot(loo_quadraticslopes)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_quadraticslopes, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(mod1.1_quadraticslopes, type ="loo_pit_overlay", ndraws=100) 


## 4.5 LOO comparison: random-effects and quadratic structure ---------------
loo_family <- loo(
  mod1.1_basegaussian, mod1.1_randomslopes,
  mod1.1_quadraticage, mod1.1_quadraticslopes,
  moment_match = TRUE
)
saveRDS(loo_family,
        "scripts/model_outputs/Offspring Trait Models/adult mass/loo_family.rda")

## 4.6 Distributional models: evidence for sigma submodel ------------------
# Does parental age (and/or temperature) predict residual variance?
mod1.1_sigma <- brm(
  bf(F1_adultmass ~ avg_age + 
       delta.age + 
       Temp + 
       F1_sex + 
       Mother_bodymass_scaled +
       (1 + delta.age | PairID),
     sigma ~ delta.age),
  family  = gaussian,
  data    = data1,
  iter    = 5000,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE),
  cores   = 4
)
saveRDS(mod1.1_sigma,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_sigma.rda")
mod1.1_sigma <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_sigma.rda")

#model diagnostics
loo_sigma <- loo(mod1.1_sigma, save_psis = TRUE, moment_match = TRUE)
plot(loo_sigma)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_sigma, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(mod1.1_sigma, type ="loo_pit_overlay", ndraws=100) 

#4.7. parental age and temperature on sigma---------------------
mod1.1_sigma_temp <- brm(
  bf(F1_adultmass ~ avg_age + 
       delta.age + 
       Temp +
       F1_sex + 
       Mother_bodymass_scaled +
       (1 + delta.age | PairID),
     sigma ~ delta.age + Temp),
  family  = gaussian,
  data    = data1,
  iter    = 5000,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE),
  cores   = 4
)
saveRDS(mod1.1_sigma_temp,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_sigma_temp.rda")
mod1.1_sigma_temp <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_sigma_temp.rda")

#model diagnostics
loo_sigma_temp <- loo(mod1.1_sigma_temp, save_psis = TRUE, moment_match = TRUE)
plot(loo_sigma_temp)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_sigma_temp, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(mod1.1_sigma_temp, type ="loo_pit_overlay", ndraws=100) 


#Temperature effects are very weakly identified; broad and flat posteriors

## LOO comparison across distributional model variants ----------------------
loo_dpars <- loo(mod1.1_randomslopes, mod1.1_sigma, mod1.1_sigma_temp,
                 moment_match = TRUE)
saveRDS(loo_dpars,
        "scripts/model_outputs/Offspring Trait Models/adult mass/loo_dpars.rda")
#Selecting model with just delta age on sigma moving forward

# ===========================================================================
# 5. PRIOR SPECIFICATION-----------------------------------------------------
# ===========================================================================
# Priors assume all continuous predictors are mean-centred and SD-scaled.
# Intercept prior informed by observed mean (0.764 g) and SD (0.188 g).

range(data1$F1_adultmass)   # 0.277–1.539 g
hist(rnorm(1000, 0.76, 0.188))  # visualise intercept prior
log(0.19)                   # sigma intercept on log scale ≈ -1.66

## Diffuse priors -------------------------------------------------------------
#Flat priors that are very similar to the brms defaults
diffuse_priors <- c(
  prior(normal(0.764, 0.188),    class = "Intercept"),
  prior(normal(log(0.19), 0.5),  class = "Intercept", dpar = "sigma"),
  prior(normal(0, 1),            class = "b",          dpar = "sigma"),
  prior(normal(0, 1),            class = "b"),
  prior(student_t(3, 0, 2.5),   class = "sd",         lb = 0),
  prior(lkj(2),                  class = "cor")
)

## Weak priors [SELECTED] ------------------------------------------------------
#Slightly regularising, assume most effects  are <1SD of the response,
#rarely culminate in a >1.96SD change in the response (following the Gaussian distribution)
weak_priors <- c(
  prior(normal(0.764, 0.188),    class = "Intercept"),
  prior(normal(log(0.19), 0.5),  class = "Intercept", dpar = "sigma"),
  prior(normal(0, 0.5),          class = "b",          dpar = "sigma"),
  prior(normal(0, 0.188),        class = "b"),
  prior(exponential(10),         class = "sd",         lb = 0), #tight to stop adult mass from inflating beyond biological limits (i.e., 3g)
  prior(lkj(2),                  class = "cor")
)

## Moderate priors ------------------------------------------------------------
#Similar to weak priors, but assume even smaller effect sizes
moderate_priors <- c(
  prior(normal(0.764, 0.188),    class = "Intercept"),
  prior(normal(log(0.19), 0.5),  class = "Intercept", dpar = "sigma"),
  prior(normal(0, 0.5),          class = "b",          dpar = "sigma"),
  prior(normal(0, 0.09),         class = "b"),
  prior(exponential(10),         class = "sd",         lb = 0),
  prior(lkj(2),                  class = "cor")
)


# ===========================================================================
# 6. PRIOR PREDICTIVE CHECKS-------------------------------------------------
# ===========================================================================

## 6.1. Prior-only model: diffuse -------------------------------------------------
diffuse_prior_model <- brm(
  bf(F1_adultmass ~ avg_age + 
       delta.age +
       Temp + 
       F1_sex +
       Mother_bodymass_scaled +
       (1 + delta.age | PairID),
     sigma ~ delta.age),
  family       = gaussian,
  data         = data1,
  prior        = diffuse_priors,
  sample_prior = "only",
  iter         = 5000,
  cores        = 4
)
saveRDS(diffuse_prior_model,
        "scripts/model_outputs/Offspring Trait Models/adult mass/diffuse_prior_model.rda")
diffuse_prior_model <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/diffuse_prior_model.rda")

#prior predictive cumulative density
diffuseprior_cumvdis <- pp_check(diffuse_prior_model, ndraws = 100,
                                 type = "ecdf_overlay") +
  scale_x_continuous(limits = c(0, 2)) +
  theme_classic()

## 6.2. Prior-only model: weak ----------------------------------------------------
weak_prior_model <- brm(
  bf(F1_adultmass ~ avg_age + 
       delta.age + 
       Temp + 
       F1_sex + 
       Mother_bodymass_scaled +
       (1 + delta.age | PairID),
     sigma ~ delta.age),
  family       = gaussian,
  data         = data1,
  prior        = weak_priors,
  sample_prior = "only",
  iter         = 5000,
  cores        = 4
)
saveRDS(weak_prior_model,
        "scripts/model_outputs/Offspring Trait Models/adult mass/weak_prior_model.rda")
weak_prior_model <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/weak_prior_model.rda")

#weak prior cumulative density
weakprior_cumvdis <- pp_check(weak_prior_model, ndraws = 100,
                              type = "ecdf_overlay", discrete = FALSE) +
  xlim(c(0, 2)) +
  theme_classic()

## 6.3. Prior-only model: moderate ------------------------------------------------
moderate_prior_model <- brm(
  bf(F1_adultmass ~ avg_age + delta.age + Temp + F1_sex + Mother_bodymass_scaled +
       (1 + delta.age | PairID),
     sigma ~ delta.age),
  family       = gaussian,
  data         = data1,
  prior        = moderate_priors,
  sample_prior = "only",
  iter         = 5000,
  cores        = 4
)
saveRDS(moderate_prior_model,
        "scripts/model_outputs/Offspring Trait Models/adult mass/moderate_prior_model.rda")
moderate_prior_model <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/moderate_prior_model.rda")

#moderate prior cumulative density
moderateprior_cumvdis <- pp_check(moderate_prior_model, ndraws = 100,
                                  type = "ecdf_overlay", discrete = FALSE) +
  xlim(c(0, 2)) +
  theme_classic()


# ===========================================================================
# 7. POSTERIOR MODELS---------------------------------------------------------
# ===========================================================================

## 7.1 Diffuse priors ---------------------------------------------------------
mod1.1_diffusepriors <- brm(
  bf(F1_adultmass ~ avg_age + 
       delta.age + 
       Temp + 
       F1_sex + 
       Mother_bodymass_scaled +
       (1 + delta.age | PairID),
     sigma ~ delta.age),
  family  = gaussian,
  data    = data1,
  prior   = diffuse_priors,
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  save_pars = save_pars(all = TRUE),
  cores   = 4
)
saveRDS(mod1.1_diffusepriors,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_diffusepriors.rda")
mod1.1_diffusepriors <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_diffusepriors.rda")

#Changing the colour scheme to pink here
color_scheme_set("pink")

#posterior cumulative density
diffuseprior_ecdf <- pp_check(mod1.1_diffusepriors, ndraws = 100,
                              type = "ecdf_overlay") +
  xlim(0, 2) + theme_classic()

## 7.2 Weak priors  [SELECTED MODEL] -----------------------------------------
mod1.1_weakpriors <- brm(
  bf(F1_adultmass ~ avg_age + 
       delta.age + 
       Temp + 
       F1_sex + 
       Mother_bodymass_scaled +
       (1 + delta.age | PairID),
     sigma ~ delta.age),
  family  = gaussian,
  data    = data1,
  prior   = weak_priors,
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  save_pars = save_pars(all = TRUE),
  cores   = 4
)
saveRDS(mod1.1_weakpriors,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_weakpriors.rda")
mod1.1_weakpriors <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_weakpriors.rda")

#posterior cumulative density
weakprior_ecdf <- pp_check(mod1.1_weakpriors, ndraws = 100,
                              type = "ecdf_overlay") +
  xlim(0, 2) + theme_classic()


## 7.3 Moderate priors --------------------------------------------------------
mod1.1_moderatepriors <- brm(
  bf(F1_adultmass ~ avg_age +
       delta.age + 
       Temp + 
       Mother_bodymass_scaled + 
       F1_sex +
       (1 + delta.age | PairID),
     sigma ~ delta.age),
  family  = gaussian,
  data    = data1,
  prior   = moderate_priors,
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  save_pars = save_pars(all = TRUE),
  cores   = 4
)
saveRDS(mod1.1_moderatepriors,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_moderatepriors.rda")
mod1.1_moderatepriors <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_moderatepriors.rda")

#moderate prior cumulative density
moderateprior_ecdf <- pp_check(mod1.1_moderatepriors, ndraws = 100,
                               type = "ecdf_overlay") +
  xlim(0, 2) + theme_classic()

## 7.4 Default (flat) priors --------------------------------------------------
mod1.1_defaultpriors <- brm(
  bf(F1_adultmass ~ 
       avg_age + 
       delta.age + 
       Temp + 
       F1_sex + 
       Mother_bodymass_scaled +
       (1 + delta.age | PairID),
     sigma ~ delta.age),
  family  = gaussian,
  data    = data1,
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  cores   = 4
)
saveRDS(mod1.1_defaultpriors,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_defaultpriors.rda")
mod1.1_defaultpriors <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod1.1_defaultpriors.rda")

## 7.5 Prior sensitivity: LOO comparison across prior specifications ---------
loo_prior_sensitivity <- loo(mod1.1_moderatepriors, mod1.1_weakpriors,
                             mod1.1_diffusepriors, moment_match = TRUE)
saveRDS(loo_prior_sensitivity,
        "scripts/model_outputs/Offspring Trait Models/adult mass/loo_prior_sensitivity.rda")


# ===========================================================================
# 8. SELECTED MODEL DIAGNOSTICS  (mod1.1_weakpriors)
# ===========================================================================

## LOO-CV --------------------------------------------------------------------
loo_weakpriors <- loo(mod1.1_weakpriors, save_psis = TRUE)
plot(loo_weakpriors)
psis_weakpriors <- loo_weakpriors$psis_object
psis_weights    <- weights(psis_weakpriors)

## Posterior predictive checks -----------------------------------------------
#Posterior simiulated and empirical (observed) means
weakprior_mean     <- ppc_stat(yrep = posterior_predict(mod1.1_weakpriors),
                               y    = data1$F1_adultmass,
                               stat = "mean") + theme_classic()
#posterior cumulative density
weakprior_ecdf     <- pp_check(mod1.1_weakpriors, ndraws = 100,
                               type = "ecdf_overlay") +
  xlim(0, 2) + theme_classic()

#posterior probability density
weakprior_pdf      <- pp_check(mod1.1_weakpriors, ndraws = 100) +
  xlim(0, 2) + theme_classic()

#posterior LOO Q-Q plot
weakprior_loo_qq   <- pp_check(mod1.1_weakpriors, type = "loo_pit_qq",
                               ndraws = 100) + theme_classic()
#LOO-PIT values fall neatly along the diagonal line (i.e. expected uniform distribution)

#LOO uniformity plot
weakprior_loo_unif <- pp_check(mod1.1_weakpriors, type = "loo_pit_overlay",
                               ndraws = 100) + theme_classic()
#predictive interval plot
weakprior_intervals <- ppc_loo_intervals(
  y    = data1$F1_adultmass,
  yrep = posterior_predict(mod1.1_weakpriors),
  psis_weakpriors
) + theme_classic()

## Combined fit plot ---------------------------------------------------------
fit_plots <- ggarrange(
  weakprior_mean, weakprior_ecdf,     weakprior_pdf,
  weakprior_loo_qq, weakprior_loo_unif, weakprior_intervals,
  nrow = 2, ncol = 3,
  labels = c("A", "B", "C", "D", "E", "F")
)

ggsave("./bayesian_plots/model fit plots/adult mass/selectedmodelfit.png",
       plot   = fit_plots,
       device = "png",
       width  = 500, height = 400, units = "mm")

## Bayes R² ------------------------------------------------------------------
bayes_R2(mod1.1_weakpriors, re.form = NA)   # fixed effects only (i.e. marginal)
bayes_R2(mod1.1_weakpriors, re.form = NULL)  # including random effects (i.e., conditional)

## MAP estimates and pd ------------------------------------------------------
MAP_adultmass <- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  test        = "pd",
  ci_method   = "HDI",
  centrality  = "MAP",
  component   = "all",
  effects     = "full",
  ci          = 0.95
)
saveRDS(MAP_adultmass,
        "scripts/model_outputs/Offspring Trait Models/adult mass/MAP_adultmass.rda")

## ROPE -----------------------------------------------------------------------
rope_adultmass <- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1 #ROPE estimates need to use the full posterior distribution, hence why they're calculated seperately
)
saveRDS(rope_adultmass,
        "scripts/model_outputs/Offspring Trait Models/adult mass/rope_adultmass.rda")


## Selective disappearance test -----------------------------------------------
hypothesis(mod1.1_weakpriors, "avg_age - delta.age > 0")
hypothesis(mod1.1_weakpriors, "avg_age - delta.age < 0")

# ===========================================================================
# 9. PRIOR vs. POSTERIOR SENSITIVITY PLOTS
# ===========================================================================

## Combined prior–posterior ECDF panel ---------------------------------------
priorpost_ecdf_panel <- ggarrange(
  diffuseprior_cumvdis, weakprior_cumvdis,  moderateprior_cumvdis,
  diffuseprior_ecdf,    weakprior_ecdf,     moderateprior_ecdf,
  nrow = 2, ncol = 3,
  labels = c("A", "B", "C", "D", "E", "F")
)

ggsave("./bayesian_plots/model fit plots/adult mass/priorvpost_ecdf.png",
       plot   = priorpost_ecdf_panel,
       device = "png",
       dpi    = 300,
       width  = 490, height = 340, units = "mm")


# ===========================================================================
# 10. HYPOTHESIS TESTING: INTERACTION MODELS
# ===========================================================================

## 10.1 Three-way interaction: delta age × temperature × offspring sex -------
#Asks whether any temperature dependent parental age effects are conditional on offspring sex
mod2.1 <- brm(
  bf(F1_adultmass ~ avg_age + delta.age + Temp + Mother_bodymass_scaled +
       F1_sex + F1_sex*delta.age*Temp + (1 + delta.age | PairID),
     sigma ~ delta.age),
  family  = gaussian,
  data    = data1,
  prior   = weak_priors,
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  cores   = 4
)
saveRDS(mod2.1,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod2.1.rda")
mod2.1 <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod2.1.rda")

## 10.2 Two-way interaction: delta age × temperature (mu only) ---------------
#Asnwers a core question of our paper: are parental age effects conditional on environmental temperatures?
mod3.1 <- brm(
  bf(F1_adultmass ~ avg_age + delta.age + Temp + F1_sex +
       Mother_bodymass_scaled + delta.age:Temp + (1 + delta.age | PairID),
     sigma ~ delta.age),
  family  = gaussian,
  data    = data1,
  prior   = weak_priors,
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  cores   = 4
)
saveRDS(mod3.1,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod3.1.rda")
mod3.1 <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod3.1.rda")

## MAP estimates for interaction model (mod3.1) ------------------------------
#For reporting in the main manuscript
MAP_adultmass_mod3.1 <- bayestestR::describe_posterior(
  mod3.1,
  test       = "pd",
  ci_method  = "HDI",
  centrality = "MAP",
  effects    = "full",
  ci         = 0.95
)
saveRDS(MAP_adultmass_mod3.1,
        "scripts/model_outputs/Offspring Trait Models/adult mass/MAP_adultmass_mod3.1.rda")

#ROPE estimates for the interactive effect (using the full posterior distribution)
rope_adultmass_mod3.1 <- bayestestR::describe_posterior(
  mod3.1,
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1
)
saveRDS(rope_adultmass_mod3.1,
        "scripts/model_outputs/Offspring Trait Models/adult mass/rope_adultmass_mod3.1.rda")

## 10.3 Two-way interaction: delta age × temperature (mu and sigma) ----------
#Same as above, but also test the interaction on trait variance
#I.e. does any age-dependent change in residual variability depend on environmental temperatures
mod3.2 <- brm(
  bf(F1_adultmass ~ avg_age + delta.age + Temp + F1_sex +
       Mother_bodymass_scaled + delta.age:Temp + (1 + delta.age | PairID),
     sigma ~ delta.age + Temp + Temp:delta.age),
  family  = gaussian,
  data    = data1,
  prior   = weak_priors,
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  cores   = 4
)
saveRDS(mod3.2,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod3.2.rda")
mod3.2 <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod3.2.rda")

## 10.4 Two-way interaction: delta age × offspring sex (mu only) -------------
#Asks whether parental age effects are dependent on offspring sex
mod4.1 <- brm(
  bf(F1_adultmass ~ avg_age + delta.age + Temp + F1_sex +
       Mother_bodymass_scaled + delta.age:F1_sex + (1 + delta.age | PairID),
     sigma ~ delta.age),
  family  = gaussian,
  data    = data1,
  prior   = weak_priors,
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  cores   = 4
)
saveRDS(mod4.1,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod4.1.rda")
mod4.1 <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod4.1.rda")

## 10.5 Two-way interaction: delta age × offspring sex (mu and sigma) --------
#Same as above, but the effects are also placed on sigma
mod4.2 <- brm(
  bf(F1_adultmass ~ avg_age + delta.age + Temp + F1_sex +
       Mother_bodymass_scaled + delta.age:F1_sex + (1 + delta.age | PairID),
     sigma ~ delta.age + F1_sex + delta.age:F1_sex),
  family  = gaussian,
  data    = data1,
  prior   = weak_priors,
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  cores   = 4
)
saveRDS(mod4.2,
        "scripts/model_outputs/Offspring Trait Models/adult mass/mod4.2.rda")
mod4.2 <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult mass/mod4.2.rda")


# ===========================================================================
# 11. MODEL SELECTION---------------------------------------------------------
# ===========================================================================
loo_hypothesis_models <- loo(
  mod1.1_weakpriors, mod2.1, mod3.1, mod4.1, mod3.2, mod4.2,
  moment_match = TRUE
)
saveRDS(loo_hypothesis_models,
        "scripts/model_outputs/Offspring Trait Models/adult mass/loo_hypothesis_models.rda")

# ===========================================================================
# 12. POSTERIOR DISTRIBUTION PLOTS  (mod1.1_weakpriors)-----------------------
# ===========================================================================
# Plots styled using stat_halfeye (ggdist)

post1 <- as_draws_df(mod1.1_weakpriors) #the posetriors of the single-effect predictors
post2<-as_draws_df(mod3.1) #for the posteriors of the core interactive effect (delta age x temperature)
post3<-as_draws_df(mod4.1) #posteriors for the offspring sex interaction

## Extract posterior draws ---------------------------------------------------
posterior_df_mod1.1 <-
  data.frame(
    "μ: Parents' Δage"              = post1$b_delta.age,
    "μ: Parents' average age"       = post1$b_avg_age,
    "μ: Parent Temperature (28.0°C)"       = post1$b_Temp28,
    "μ: Parent Temperature (30.5°C)"       = post1$b_Temp30.5,
    "μ: Offspring sex (Male)"       = post1$b_F1_sexM,
    "μ: Mother's adult mass"        = post1$b_Mother_bodymass_scaled,
    "μ: Δage × Temperature (28.0°C)"            = post2$`b_delta.age:Temp28`,
    "μ: Δage × Temperature (30.5°C)"            = post2$`b_delta.age:Temp30.5`,
    "μ: Δage × Offspring sex (Male)" =  post3$`b_delta.age:F1_sexM`,
    "σ: Parents' Δage"              = post1$b_sigma_delta.age,
    check.names = FALSE
    ) %>%

  pivot_longer(everything(),
               names_to  = "parameter",
               values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = c(
    "μ: Parents' Δage",
    "μ: Parents' average age",
    "μ: Parent Temperature (28.0°C)",
    "μ: Parent Temperature (30.5°C)",
    "μ: Offspring sex (Male)",
    "μ: Mother's adult mass",
    "μ: Δage × Temperature (28.0°C)",
    "μ: Δage × Temperature (30.5°C)",
    "μ: Δage × Offspring sex (Male)",
    "σ: Parents' Δage"
  )))

y_levels <- rev(c(
  "μ: Parents' Δage",
  "μ: Parents' average age",
  "μ: Parent Temperature (28.0°C)",
  "μ: Parent Temperature (30.5°C)",
  "μ: Offspring sex (Male)",
  "μ: Mother's adult mass",
  "μ: Δage × Temperature (28.0°C)",
  "μ: Δage × Temperature (30.5°C)",
  "μ: Δage × Offspring sex (Male)",
  "σ: Parents' Δage" # this is the one we want unshaded in the plot below (since we didn't use ROPE for sigma values)
))



## Halfeye posterior plot ----------------------------------------------------
posterior_plot_mod1.1 <- ggplot(
  posterior_df_mod1.1,
  aes(x = value, y = parameter, fill = parameter)
) +
  annotate("rect",
           xmin = -0.02, xmax = 0.02,
           ymin = 1 - 0.5,  ymax = length(y_levels) - 1 + 0.5,
           alpha = 0.10, fill = "grey20")+
  stat_halfeye(
    scale         = 0.8,
    adjust        = 1,
    justification = -0.1,
    .width        = 0.95,
    slab_colour   = "grey85",
    slab_linewidth = 4,
    slab_alpha    = 0.8,
    linewidth     = 10,
    point_size    = 15,
    point_colour  = "white",
    shape         = 21,
    stroke        = 1.5
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey21") +
  scale_fill_manual(values = c(
    "#4A6479", "#A9BAC2", "#5f9289", "#e653a1",
    "#f16122", "#d4a5a5", "#e5a9f5", "#d2e9f5", "#8C2F4B", "#9b7fc7"
  )) +
  xlim(c(-0.20, 0.10)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 30)
  ) +
  labs(x = "Fixed-effect posterior estimates\n(on μ = mean adult mass (g); on σ = residual SD, log scale)",
       y = NULL)

ggsave("./bayesian_plots/posterior plots/adult mass/001_hypothesis1_halfeye.png",
       plot   = posterior_plot_mod1.1,
       device = "png",
       width  = 580, height = 400, units = "mm")


## Marginal means via emmeans ------------------------------------------------
pairwise_age  <- emmeans(mod1.1_weakpriors, "delta.age", type = "response",
                         at = list(delta.age = c(-1.96, 1.96)))
pairs(pairwise_age)

## Post-hoc standardised mean difference (SMD) --------------------------------
# Compares extreme ends of within-individual parental age; standardised by
# observed SD of F1 adult mass. Descriptive only — not used for inference.
sd_adultmass <- sd(data1$F1_adultmass)

mass_SMD <- summary(pairwise_age) %>%
  summarise(
    SMD       = (emmean[2]    - emmean[1])    / sd_adultmass,
    SMD_lower = (lower.HPD[2] - upper.HPD[1]) / sd_adultmass,
    SMD_upper = (upper.HPD[2] - lower.HPD[1]) / sd_adultmass
  )

saveRDS(mass_SMD,
        "scripts/model_outputs/Offspring Trait Models/mass_SMD.rda")
#This is exported to create our posterior effect size plot

# ===========================================================================
# 13. MODEL PREDICTIONS AND RESULTS PLOTS--------------------------------------
# ===========================================================================

## 13.1 Prediction grids: within- and between-individual age effects ---------
# Within-individual: vary delta.age, hold avg_age constant
df_predict_within <- expand.grid(
  avg_age                = mean(data1$avg_age),
  delta.age              = unique(data1$delta.age),
  Mother_bodymass_scaled = mean(data1$Mother_bodymass_scaled),
  Temp                   = unique(data1$Temp),
  F1_sex                 = unique(data1$F1_sex),
  PairID                 = unique(data1$PairID)[1]
)
df_predict_within$F1_sex <- factor(df_predict_within$F1_sex, levels = c("F", "M"))

pred_within <- as.data.frame(
  fitted(mod1.1_weakpriors, df_predict_within, re_formula = NA)
)
df_predict_within <- df_predict_within %>%
  mutate(prediction = pred_within$Estimate,
         lower      = pred_within$Q2.5,
         upper      = pred_within$Q97.5)

within_age_newdat <- df_predict_within %>%
  group_by(delta.age) %>%
  summarise(
    mean  = mean(prediction),
    lower = mean(lower),
    upper = mean(upper)
  ) %>%
  ungroup() %>%
  mutate(within_subject_age = delta.age * sd(data1$within_subject_age))

# 13.2 Between-individual: vary avg_age, hold delta.age constant
df_predict_between <- expand.grid(
  avg_age                = unique(data1$avg_age),
  delta.age              = mean(data1$delta.age, na.rm = TRUE),
  Mother_bodymass_scaled = mean(data1$Mother_bodymass_scaled),
  Temp                   = unique(data1$Temp),
  F1_sex                 = unique(data1$F1_sex),
  PairID                 = unique(data1$PairID)[1]
)
df_predict_between$F1_sex <- factor(df_predict_between$F1_sex, levels = c("F", "M"))

pred_between <- as.data.frame(
  fitted(mod1.1_weakpriors, df_predict_between, re_formula = NA)
)
df_predict_between <- df_predict_between %>%
  mutate(prediction = pred_between$Estimate,
         lower      = pred_between$Q2.5,
         upper      = pred_between$Q97.5)

between_age_newdat <- df_predict_between %>%
  group_by(avg_age) %>%
  summarise(
    prediction = mean(prediction),
    lower      = mean(lower),
    upper      = mean(upper)
  ) %>%
  mutate(
    avg_age_raw     = avg_age * sd(data1$avg.age) + mean(data1$avg.age),
    avg_age_centred = avg_age_raw - mean(avg_age_raw)
  )

## 13.3 Raw within-individual means (for overlaying on plot) ----------------
breaks_wk <- seq(-4.5, 4.5, by = 1)
labels_wk  <- -4:4

data1 <- data1 %>%
  mutate(delta_age_bin = cut(within_subject_age, breaks = breaks_wk,
                             labels = labels_wk))

raw_deltaage_means <- data1 %>%
  select(-within_subject_age) %>%
  filter(!is.na(delta_age_bin)) %>%
  rename(within_subject_age = delta_age_bin) %>%
  mutate(within_subject_age = as.numeric(as.character(within_subject_age))) %>%
  group_by(within_subject_age) %>%
  summarise(
    n         = sum(!is.na(F1_adultmass)),
    mean_size = mean(F1_adultmass, na.rm = TRUE),
    se_size   = ifelse(n > 1, sd(F1_adultmass, na.rm = TRUE) / sqrt(n), NA_real_),
    .groups   = "drop"
  )

## 13.4 Main parental age effect plot ----------------------------------------
adultmass_plot <- ggplot(data1,
                         aes(x = within_subject_age, y = F1_adultmass)) +
  geom_point(position = position_jitter(width = 0.2, height = 0.001),
             shape = 21, size = 6, stroke = 1.8, alpha = 0.7,
             colour = "white", fill = "grey") +
  geom_line(data = within_age_newdat,
            aes(x = within_subject_age, y = mean,
                linetype = "Within-individual",
                colour   = "Within-individual"),
            linewidth = 5) +
  geom_ribbon(data = within_age_newdat,
              aes(y = NULL, ymin = lower, ymax = upper,
                  fill   = "Within-individual",
                  colour = "Within-individual"),
              alpha = 0.2, linewidth = 2, show.legend = FALSE) +
  geom_line(data = between_age_newdat,
            aes(x = avg_age_centred, y = prediction,
                linetype = "Between-individual",
                colour   = "Between-individual"),
            linewidth = 4) +
  geom_ribbon(data = between_age_newdat,
              aes(x = avg_age_centred, y = NULL,
                  ymin   = lower, ymax = upper,
                  fill   = "Between-individual",
                  colour = "Between-individual"),
              alpha = 0.1, linewidth = 1, show.legend = FALSE) +
  geom_linerange(data = raw_deltaage_means,
                 aes(y    = mean_size,
                     ymin = mean_size - se_size,
                     ymax = mean_size + se_size,
                     colour = "Within-individual"),
                 linewidth = 3, show.legend = FALSE) +
  geom_point(data = raw_deltaage_means,
             aes(x = within_subject_age, y = mean_size,
                 fill = "Within-individual"),
             shape = 21, stroke = 1.8, size = 18, alpha = 1,
             colour = "white", show.legend = FALSE) +
  scale_x_continuous(breaks = -4:4) +
  scale_y_continuous(breaks = c(0.0, 0.3, 0.6, 0.9, 1.2, 1.5)) +
  scale_fill_manual(name   = "Parental age effect",
                    values = c("Within-individual"  = "#4A6479",
                               "Between-individual" = "#8C2F4B")) +
  scale_colour_manual(name  = "Parental age effect",
                      values = c("Within-individual"  = "#4A6479",
                                 "Between-individual" = "#8C2F4B")) +
  scale_linetype_manual(name   = "Parental age effect",
                        values = c("Within-individual"  = "solid",
                                   "Between-individual" = "dashed")) +
  guides(
    linetype = guide_legend(order = 1, title.position = "top",
                            title.hjust = 0.5,
                            override.aes = list(size = 15)),
    fill     = guide_legend(order = 1, title.position = "top",
                            title.hjust = 0.5),
    colour   = guide_legend(order = 1, title.position = "top",
                            title.hjust = 0.5)
  ) +
  theme_classic() +
  theme(
    legend.position    = "top",
    legend.box         = "horizontal",
    legend.title       = element_text(size = 50),
    legend.text        = element_text(size = 50),
    axis.title         = element_text(size = 50),
    axis.text          = element_text(size = 50),
    strip.text         = element_text(size = 50),
    panel.background   = element_rect(fill = "transparent", colour = NA),
    plot.background    = element_rect(fill = "transparent", colour = NA),
    legend.background  = element_rect(fill = "transparent", colour = NA)
  ) +
  labs(x = "Parents' adult age at reproduction (weeks; mean-centred)",
       y = "Offspring adult mass (g)")

ggsave("./bayesian_plots/offspring trait plots/F1 mass/001_adultmass_plot_tidied.png",
       plot = adultmass_plot, bg = "transparent",
       device = "png", width = 450, height = 700, units = "mm")

## 13.5 Probability density function (PDF) for F1 adult mass -----------------
#function to estimate the porbability density of a Gaussian distribution
normal_pdf <- function(mu, sigma, x) dnorm(x, mean = mu, sd = sigma)

#Data frame to estimate the probability densities from
newdat_pdf <- expand.grid(
  avg_age                = mean(data1$avg_age),
  delta.age              = c(-1.95, 0, 1.85),
  Mother_bodymass_scaled = mean(data1$Mother_bodymass_scaled),
  Temp                   = unique(data1$Temp),
  F1_sex                 = unique(data1$F1_sex),
  PairID                 = unique(data1$PairID)[1]
)
#sequence of bodymass to be estimated (rather than estimating every possible value in the whole range)
bodymass_seq <- seq(min(data1$F1_adultmass, na.rm = TRUE),
                    max(data1$F1_adultmass, na.rm = TRUE),
                    length.out = 250)

#Estimating the probability density
posterior_pdf_draws <- mod1.1_weakpriors %>%
  linpred_draws(
    newdat_pdf,
    value            = "mu",
    allow_new_levels = TRUE,
    transform        = TRUE,
    re_formula       = NA,
    dpar             = "sigma",
    ndraws           = 3000,
    seed             = 123
  ) %>%
  ungroup() %>%
  crossing(x = bodymass_seq) %>%
  mutate(f = normal_pdf(mu, sigma, x)) %>%
  group_by(x, .draw, delta.age) %>% #carrying across posterior uncertainty, and estimating how the densities differ by parental age
  summarise(f = median(f), .groups = "drop")

#summarising the full posteriors 
pdf_summary <- posterior_pdf_draws %>%
  group_by(x, delta.age) %>%
  summarise(
    median = median(f),
    lower  = quantile(f, 0.025),
    upper  = quantile(f, 0.975),
    .groups = "drop"
  ) %>%
  mutate(timepoint = case_when(
    delta.age == -1.95 ~ "Early-Aged",
    delta.age ==  0    ~ "Middle-Aged",
    delta.age ==  1.85 ~ "Late-Aged"
  )) %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("Early-Aged", "Middle-Aged", "Late-Aged")))

#estimating the posterior median (for plotting purposes)
median_mass_pdf <- mod1.1_weakpriors %>%
  linpred_draws(
    newdat_pdf,
    value            = "mu",
    allow_new_levels = TRUE,
    transform        = TRUE,
    re_formula       = NA,
    dpar             = "sigma",
    ndraws           = 5000,
    seed             = 123
  ) %>%
  ungroup() %>%
  group_by(.draw, delta.age) %>%
  summarise(median_mass = median(mu), .groups = "drop") %>%
  group_by(delta.age) %>%
  median_hdi(median_mass) %>%
  mutate(timepoint = case_when(
    delta.age == -1.95 ~ "Early-Aged",
    delta.age ==  0    ~ "Middle-Aged",
    delta.age ==  1.85 ~ "Late-Aged"
  ))

#putting the pdf plot together
mass_pdf_plot <- ggplot(pdf_summary,
                        aes(x = x, y = median, colour = timepoint)) +
  geom_line(aes(colour = timepoint), linewidth = 4, alpha = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper,
                  fill = timepoint, colour = timepoint, alpha = timepoint),
              linetype = "dashed", linewidth = 1) +
  geom_segment(data = median_mass_pdf,
               aes(x = median_mass, xend = median_mass,
                   y = 0, yend = 3.19, colour = timepoint),
               linetype = "dashed", linewidth = 3, alpha = 0.9) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) +
  scale_y_continuous(breaks = c(0, 0.80, 1.60, 2.40, 3.20),
                     limits = c(0, 3.3)) +
  scale_colour_manual(name   = "Parents' adult age at reproduction",
                      values = c("#f96161", "#66b2b2", "#066594")) +
  scale_fill_manual(name     = "Parents' adult age at reproduction",
                    values   = c("#f96161", "#66b2b2", "#066594")) +
  scale_alpha_manual(name    = "Parents' adult age at reproduction",
                     values  = c(0.6, 0.1, 0.6)) +
  guides(
    alpha  = guide_legend(order = 1, title.position = "top",
                          title.hjust = 0.5,
                          override.aes = list(size = 10)),
    fill   = guide_legend(order = 1, title.position = "top",
                          title.hjust = 0.5),
    colour = guide_legend(order = 1, title.position = "top",
                          title.hjust = 0.5)
  ) +
  theme_classic() +
  theme(
    legend.position   = "top",
    legend.box        = "horizontal",
    legend.title      = element_text(size = 50),
    legend.text       = element_text(size = 50),
    axis.title        = element_text(size = 50),
    axis.text         = element_text(size = 50),
    strip.text        = element_text(size = 50, face = "bold"),
    panel.background  = element_rect(fill = "transparent", colour = NA),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)
  ) +
  labs(x = "Offspring adult mass (g)", y = "Probability density, f(x)")

ggsave("./bayesian_plots/offspring trait plots/F1 mass/003_adultmass_pdf.png",
       plot = mass_pdf_plot, bg = "transparent",
       device = "png", width = 520, height = 400, units = "mm")

## 14.6 Combined inference panel (main effect + PDF) -------------------------
inference_panel <- ggarrange(
  adultmass_plot, mass_pdf_plot,
  ncol = 2, nrow = 1,
  labels    = c("A", "B"),
  font.label = list(size = 50, face = "bold"),
  label.x   = c(0.05, 0.05),
  widths    = c(0.7, 1)
)

ggsave("./bayesian_plots/offspring trait plots/F1 mass/inference_plots.png",
       plot   = inference_panel,
       device = "png",
       width  = 980, height = 570, units = "mm")

#Creating a table to export------------------------------------------------------------------------------
# ── helper ────────────────────────────────────────────────────────────────────
#Model draws from selected models
base_model_draws<-as_draws_df(mod1.1_weakpriors) #just single effects, and the final selected model
sigma_dropped<-as_draws_df(mod1.1_sigma_temp)
interaction_draws<-as_draws_df(mod2.1) 
interaction_draws2<-as_draws_df(mod3.1) 
interaction_draws2.1<-as_draws_df(mod3.2) 
interaction_draws3<-as_draws_df(mod4.1)
interaction_draws3.1<-as_draws_df(mod4.2)


#estimating the difference between average age and delta age
differences<-data.frame(
  selectivedis = base_model_draws$b_avg_age - base_model_draws$b_delta.age
)

SDy<-sd(data1$F1_adultmass)
#function for generating estimates
summarise_param <- function(draws, param, section,
                            rope = TRUE, #using default ROPE range automatically calculated from bayestest
                            rope_range = c(-0.1*SDy, +0.1*SDy),
                            pd = TRUE) {
  draws <- as.numeric(draws) 
  h <- bayestestR::hdi(draws, ci = 0.95)
  tibble(
    Section    = section,
    Parameter  = param,
    MAP        = round(as.numeric(map_estimate(draws)), 3),
    `95% HDI`  = paste0("[", round(h$CI_low, 3), ", ", round(h$CI_high, 3), "]"),
    `% in ROPE` = if (rope) paste0(round(as.numeric(rope(draws, range = rope_range, ci = 1)$ROPE_Percentage) * 100, 1), "%") else "—",
    pd         = if (pd)   paste0(round(as.numeric(p_direction(draws)) * 100, 1), "%") else "—"
  )
}



# ── 1. Fixed effects (location) — ROPE + pd ───────────────────────────────────
fe<- bind_rows(
  summarise_param(base_model_draws$b_Intercept,        "Intercept",            "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws$b_avg_age,   "Parents' average age (z-scaled)", "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws$b_delta.age, "Parents' Δage (z-scaled)",        "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws$b_Temp28, "Temperature: 28.0°C",        "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws$b_Temp30.5, "Temperature: 30.5°C",        "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws$b_F1_sexM, "Offspring sex (Male)", "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws$b_Mother_bodymass_scaled, "Mother's adult mass (g); z-scaled", "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(differences$selectivedis, "Average age - Δage", "Location submodel (Fixed effects)", rope = FALSE)
)


# ── 1. dropped interaction (location) — ROPE + pd ───────────────────────────────────
fe_interaction2 <- bind_rows(summarise_param(interaction_draws$`b_delta.age:Temp28:F1_sexM`, "Parents' Δage x  Temperature: 28.0°C x Offspring sex (Male)","Location submodel (Dropped interactions)", rope = TRUE),
                             summarise_param(interaction_draws$`b_delta.age:Temp30.5:F1_sexM`, "Parents' Δage x  Temperature: 30.5°C x Offspring sex (Male)","Location submodel (Dropped interactions)", rope = TRUE),
                             summarise_param(interaction_draws2$`b_delta.age:Temp28`, "Parents' Δage x  Temperature: 28.0°C", "Location submodel (Dropped interactions)", rope = TRUE),
                             summarise_param(interaction_draws2$`b_delta.age:Temp30.5`, "Parents' Δage x  Temperature: 30.5°C", "Location submodel (Dropped interactions)", rope = TRUE),
                             summarise_param(interaction_draws3$`b_delta.age:F1_sexM`, "Parents' Δage x Offspring sex (Male)","Location submodel (Dropped interactions)", rope = TRUE)
)

# ── 3.Random effects ──────────────────────────────────────────────────
# SDs 
re_sd_explore<- bind_rows(
  summarise_param(base_model_draws$sd_PairID__Intercept, "σ intercept",    "Random effects (Location)", pd = FALSE),
  summarise_param(base_model_draws$sd_PairID__delta.age, "σ slope Δage",    "Random effects (Location)", pd = FALSE),
  summarise_param(base_model_draws$cor_PairID__Intercept__delta.age, "r intercept ~ slope Δage",    "Random effects (Location)", pd = TRUE)
)

#---4. Sigma parameter──────────────────────────────────────────────────
sigma<- bind_rows(
  summarise_param(base_model_draws$b_sigma_Intercept, "Sigma intercept (σ)", "Distributional parameters (scale component)",  rope = FALSE),
  summarise_param(base_model_draws$b_sigma_delta.age, "Parents' Δage (z-scaled)", "Distributional parameters (scale component)", rope= FALSE)
  )

#---5. Sigma dropped terms----------------------------------------------
sigma_dropped2<-bind_rows(
  summarise_param(sigma_dropped$b_sigma_Temp28, "Temperature: 28.0°C", "Distributional parameters (Dropped terms)", rope = FALSE),
  summarise_param(sigma_dropped$b_sigma_Temp30.5, "Temperature: 30.5°C", "Distributional parameters (Dropped terms)", rope = FALSE),
  summarise_param(interaction_draws2.1$`b_sigma_delta.age:Temp28`, "Parents' Δage x Temperature: 28.0°C", "Distributional parameters (Dropped terms)",  rope = FALSE),
  summarise_param(interaction_draws2.1$`b_sigma_delta.age:Temp30.5`, "Parents' Δage x Temperature: 30.5°C", "Distributional parameters (Dropped terms)",  rope = FALSE),
  summarise_param(interaction_draws3.1$`b_sigma_delta.age:F1_sexM`, "Parents' Δage x Offspring sex (Male)", "Distributional parameters (Dropped terms)", rope = FALSE))
  

#---5. Conditional and marginal Bayes R² ------------------------------------------------------------------
marginal<-bayes_R2(mod1.1_weakpriors, re.form = NA, summary = FALSE)   # fixed effects only (i.e., marginal)
conditional<-bayes_R2(mod1.1_weakpriors, re.form = NULL, summary = FALSE)  # including random effects (i.e., conditional)

bayes<-bind_rows(
  summarise_param(marginal, "Marginal R²", "Bayes R²", pd= FALSE, rope= FALSE),
  summarise_param(conditional, "Conditional R²", "Bayes R²", pd= FALSE, rope =FALSE)
)


# ── 5. Combine and render ─────────────────────────────────────────────────────
bind_rows(fe, fe_interaction2, re_sd_explore, sigma, sigma_dropped2, bayes) %>%
  gt(groupname_col = "Section") %>%
  cols_label(Parameter = "Parameter", MAP = "MAP",
             `95% HDI` = "95% HDI", `% in ROPE` = "% in ROPE", pd = "pd") %>%
  tab_header(
    title    = "Adult mass model summary: Gaussian model"
  ) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) %>%
  tab_footnote("ROPE = [−0.019, 0.019] on unit scale (grams); applied to location fixed effects only.",
               locations = cells_column_labels(`% in ROPE`)) %>%
  tab_footnote("pd = probability of direction.",
               locations = cells_column_labels(pd)) %>%
  opt_stylize(style = 1)%>%
  gtsave("./tables/adultmass_updated.docx")


###############################################################################
##  END OF SCRIPT-------------------------------------------------------------------
###############################################################################
