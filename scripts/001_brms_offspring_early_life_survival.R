###############################################################################
##  Script: Effect of parental age on F1 early-life survival
##  Note:   Early-life survival = survival to week 4 of life 
#(early-life survival coded as binary 0/1)
##          All 1317 F1 offspring included in this analysis
###############################################################################

# ===========================================================================
# 1. SETUP-------------------------------------------------------------------
# ===========================================================================

## Libraries ------------------------------------------------------------------
library(dplyr) #v.1.1.4
library(tidyr)  #v.1.3.1
library(ggplot2) #v.4.0.0
library(patchwork) #v.1.3.0
library(ggpubr) #v.0.6.1
library(ggdist) #v.3.3.3
library(ggridges) #v.0.5.6
library(emmeans) #v.1.10.1
library(brms) #v.2.21.0
library(bayesplot) #v.1.11.1
library(bayestestR) #v.0.17.0
library(tidybayes) #v.3.0.7
library(parameters) #v.0.28.2
library(sjPlot) #v.2.8.16
library(loo) #v.2.8.0
library(car) #v.3.1.2
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

## Load data ------------------------------------------------------------------
F1data <- readRDS('./raw data/F1_filtered_data06022025.RDS')

## Sample sizes ---------------------------------------------------------------
length(unique(F1data$PairID))  # 78 parent pairs
length(F1data$F1_ID)           # 1317 offspring

## Create early-life survival variable ----------------------------------------
# early_life_surv = 1 if offspring survived past week 4, 0 otherwise
F1data <- F1data %>%
  mutate(early_life_surv = ifelse(total_lifespan <= 4, 0, 1))

# 25.05% of animals died in their first four weeks of life
mean(F1data$early_life_surv == 0, na.rm = TRUE)

## Scale continuous predictors ------------------------------------------------
# All predictors scaled to mean = 0, SD = 1; a 1-unit change = 1 SD
F1data$avg_age <- as.numeric(scale(F1data$avg.age, center = TRUE, scale = TRUE))

# Within-individual (delta) age: scaled but NOT mean-centred; NAs replaced with 0
F1data <- F1data %>%
  group_by(PairID) %>%
  mutate(
    delta.age = as.numeric(scale(within_subject_age,
                                 center = FALSE,
                                 scale = TRUE)),
    delta.age = replace(delta.age, is.na(delta.age), 0)
  ) %>%
  ungroup()

## Sample sizes ---------------------------------------------------------------
length(unique(F1data$Mother_ID))  # number of mothers (N = 78)
length(unique(F1data$F1_ID))      # number of F1 offspring (n = 1317)


# ===========================================================================
# 3. SUMMARY STATISTICS--------------------------------------------------------
# ===========================================================================
ggplot(F1data, aes(x = early_life_surv)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +  # Create histogram
  facet_wrap(~Timepoint) +  
  labs(
    x = "Early-life survival",
    y = "Count"
  ) +
  theme_classic()


## Overall mean and SD --------------------------------------------------------
sum_overall <- F1data %>%
  summarise(
    mean_surv   = mean(early_life_surv),
    sd_total    = sd(early_life_surv),
    n_offspring = n(),
    n_parents   = n_distinct(PairID)
  )
# mean_surv = 0.749, SD = 0.434, n = 1317, n_parents = 78

## Mean early survival by parental age timepoint ------------------------------
sum_by_timepoint <- F1data %>%
  group_by(Timepoint) %>%
  summarise(
    mean_surv   = mean(early_life_surv),
    sd_total    = sd(early_life_surv),
    n_offspring = n(),
    n_parents   = n_distinct(PairID)
  )

## Mean early survival by parental temperature treatment ----------------------
sum_by_temp <- F1data %>%
  group_by(Temp) %>%
  summarise(
    mean_surv   = mean(early_life_surv),
    sd_total    = sd(early_life_surv),
    n_offspring = n(),
    n_parents   = n_distinct(PairID)
  )

#Collinearity checks
X <- model.matrix(
  ~ avg_age + delta.age + cum_succesful_matings +
    Temp + F1_sex,
  data = data.frame(
    avg_age = F1data $avg_age,
    delta.age = F1data $delta.age,
    Temp = F1data$Temp,
    cum_succesful_matings = F1data$cum_successful_matings,
    F1_sex = F1data$F1_sex
  )
)
# Correlation matrix
cor(X[, -1])  # correlation matrix (excluding intercept column)


# Variance inflation factors (VIFs)
vif(lm(early_life_surv ~ avg_age + delta.age + cum_successful_matings+
         Temp, data = F1data))

#removing cum_succesful matings
vif(lm(early_life_surv ~ avg_age + delta.age +
         Temp, data = F1data))


# ===========================================================================
# 4. MODEL STRUCTURE SELECTION------------------------------------------------
# ===========================================================================
# Compare Bernoulli models varying in random-effects structure and
# whether a quadratic delta-age term is supported by the data.

## 4.1 Baseline: random intercept only, linear delta age ---------------------
mod1.1_basebernoulli <- brm(
  early_life_surv ~ avg_age + delta.age + Temp + (1 | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  save_pars = save_pars(all = TRUE),
  cores     = 4
)
saveRDS(mod1.1_basebernoulli,
        "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_basebernoulli.rda")
mod1.1_basebernoulli <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_basebernoulli.rda")

#model diagnostics
loo_basebernoulli <- loo(mod1.1_basebernoulli, save_psis = TRUE, moment_match = TRUE)
plot(loo_basebernoulli)

## 4.2 Quadratic delta-age term ---------------------------------------------
mod1.1_quadratic <- brm(
  early_life_surv ~ avg_age + 
    poly(delta.age,2) + 
    Temp + 
    (1 | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  save_pars = save_pars(all = TRUE),
  cores     = 4
)
saveRDS(mod1.1_quadratic,
        "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_quadratic.rda")
mod1.1_quadratic <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_quadratic.rda")

loo_quadratic <- loo(mod1.1_quadratic, save_psis = TRUE, moment_match = TRUE)

## 4.3 Random slopes (intercept + delta age), quadratic term -----------------
#females are not allowed to vary in the shape of their slopes here (assumes a shared quadratic trajectory across all individuals)
mod1.1_randomslopes <- brm(
  early_life_surv ~ avg_age + 
    poly(delta.age, 2) + 
    Temp +
    (1 + delta.age | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  save_pars = save_pars(all = TRUE),
  cores     = 4
)
saveRDS(mod1.1_randomslopes,
        "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_randomslopes.rda")
mod1.1_randomslopes <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_randomslopes.rda")

loo_randomslopes <- loo(mod1.1_randomslopes, save_psis = TRUE, moment_match = TRUE)

## 4.4 Random slopes, linear delta age [REFERENCE MODEL STRUCTURE] ----------
mod1.1_randomslope_linear <- brm(
  early_life_surv ~ avg_age + 
    delta.age +
    Temp + 
    (1 + delta.age | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  save_pars = save_pars(all = TRUE),
  cores     = 4
)
saveRDS(mod1.1_randomslope_linear,
        "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_randomslope_linear.rda")
mod1.1_randomslope_linear <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_randomslope_linear.rda")

loo_randomslope_linear <- loo(mod1.1_randomslope_linear, save_psis = TRUE, moment_match = TRUE)

## 4.5 LOO comparison: random-effects structure and quadratic term -----------
loo_structure <- loo_compare(
  loo_basebernoulli, loo_randomslopes,
  loo_quadratic,     loo_randomslope_linear
)
saveRDS(loo_structure,
        "scripts/model_outputs/Offspring Trait Models/early life surv/modelfamily2.rda")

# Selected structure: random slopes, linear delta age (no convincing evidence for a quadratci age effect)


# ===========================================================================
# 5. PRIOR SPECIFICATION------------------------------------------------------
# ===========================================================================
# Priors assume all continuous predictors are mean-centred and SD-scaled.
# Intercept prior informed by observed mean early-life survival (0.749):
# qlogis(0.749) = 1.093286

## Diffuse priors -------------------------------------------------------------
diffuse_priors <- c(
  prior(normal(1.093286, 1.5), class = "Intercept"),
  prior(normal(0, 1.5),        class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd", lb = 0),
  prior(lkj(2),                class = "cor")
)

## Weak priors [SELECTED] -----------------------------------------------------
weak_priors <- c(
  prior(normal(1.093286, 1), class = "Intercept"),
  prior(normal(0, 1),        class = "b"),
  prior(exponential(1),      class = "sd", lb = 0),
  prior(lkj(2),              class = "cor")
)

## Moderate priors ------------------------------------------------------------
moderate_priors <- c(
  prior(normal(1.093286, 0.5), class = "Intercept"),
  prior(normal(0, 0.5),        class = "b"),
  prior(exponential(3),        class = "sd", lb = 0),
  prior(lkj(2),                class = "cor")
)


# ===========================================================================
# 6. PRIOR PREDICTIVE CHECKS--------------------------------------------------
# ===========================================================================

##6.1. Prior-only model: diffuse -------------------------------------------------
diffuse_prior_model <- brm(
  early_life_surv ~ avg_age +
    delta.age + 
    Temp + 
    (1 + delta.age | PairID),
  family       = bernoulli(link = "logit"),
  data         = F1data,
  prior        = diffuse_priors,
  sample_prior = "only",
  iter         = 5000,
  control      = list(adapt_delta = 0.98),
  save_pars    = save_pars(all = TRUE),
  cores        = 4
)
saveRDS(diffuse_prior_model,
        "scripts/model_outputs/Offspring Trait Models/early life surv/diffuse_prior_model.rda")
diffuse_prior_model <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/early life surv/diffuse_prior_model.rda")

color_scheme_set("teal")

#Prior draws (I.e., the prior cumulative density function)
diffuse_prior_cumvdis <- pp_check(diffuse_prior_model, ndraws = 100) +
  xlim(-0.5, 1.5) + theme_classic()

## Prior-only model: weak ----------------------------------------------------
weak_prior_model <- brm(
  early_life_surv ~ avg_age +
    delta.age + 
    Temp + 
    (1 + delta.age | PairID),
  family       = bernoulli(link = "logit"),
  data         = F1data,
  prior        = weak_priors,
  sample_prior = "only",
  iter         = 5000,
  control      = list(adapt_delta = 0.98),
  save_pars    = save_pars(all = TRUE),
  cores        = 4
)
saveRDS(weak_prior_model,
        "scripts/model_outputs/Offspring Trait Models/early life surv/weak_prior_model.rda")
weak_prior_model <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/early life surv/weak_prior_model.rda")

#Prior draws (I.e., the prior cumulative density function)
weak_prior_cumvdis <- pp_check(weak_prior_model, ndraws = 100,
                               type = "dens_overlay") +
  xlim(-0.5, 1.5) + theme_classic()

## Prior-only model: moderate ------------------------------------------------
moderate_prior_model <- brm(
  early_life_surv ~ avg_age + 
    delta.age + 
    Temp + 
    (1 + delta.age | PairID),
  family       = bernoulli(link = "logit"),
  data         = F1data,
  prior        = moderate_priors,
  sample_prior = "only",
  iter         = 5000,
  control      = list(adapt_delta = 0.98),
  save_pars    = save_pars(all = TRUE),
  cores        = 4
)
saveRDS(moderate_prior_model,
        "scripts/model_outputs/Offspring Trait Models/early life surv/moderate_prior_model.rda")
moderate_prior_model <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/early life surv/moderate_prior_model.rda")

#Prior draws (I.e., the prior cumulative density function)
moderate_prior_cumvdis <- pp_check(moderate_prior_model, ndraws = 100) +
  xlim(-0.5, 1.5) + theme_classic()


# ===========================================================================
# 7. POSTERIOR MODELS---------------------------------------------------------
# ===========================================================================

## 7.1 Diffuse priors ---------------------------------------------------------
mod1.1_diffusepriors <- brm(
  early_life_surv ~ avg_age + 
    delta.age + 
    Temp + 
    (1 + delta.age | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  prior     = diffuse_priors,
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  save_pars = save_pars(all = TRUE),
  cores     = 4
)
saveRDS(mod1.1_diffusepriors,
        "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_diffusepriors.rda")
mod1.1_diffusepriors <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_diffusepriors.rda")

color_scheme_set("pink")

#Posterior draws (I.e., the posterior cumulative density function)
diffuse_prior_pdf <- pp_check(mod1.1_diffusepriors, ndraws = 100) +
  xlim(-0.5, 1.5) + theme_classic()

## 7.2 Weak priors  [SELECTED MODEL] -----------------------------------------
mod1.1_weakpriors <- brm(
  early_life_surv ~ avg_age + 
    delta.age + 
    Temp + 
    (1 + delta.age | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  prior     = weak_priors,
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  save_pars = save_pars(all = TRUE),
  cores     = 4
)
saveRDS(mod1.1_weakpriors,
        "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_weakpriors.rda")
mod1.1_weakpriors <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_weakpriors.rda")

#Posterior draws (I.e., the posterior cumulative density function)
weak_prior_pdf <- pp_check(mod1.1_weakpriors, ndraws = 100,
                           type = "dens_overlay") +
  xlim(-0.5, 1.5) + theme_classic()

## 7.3 Moderate priors --------------------------------------------------------
mod1.1_moderatepriors <- brm(
  early_life_surv ~ avg_age + 
    delta.age + 
    Temp +
    (1 + delta.age | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  prior     = moderate_priors,
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  save_pars = save_pars(all = TRUE),
  cores     = 4
)
saveRDS(mod1.1_moderatepriors,
        "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_moderatepriors.rda")
mod1.1_moderatepriors <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/early life surv/mod1.1_moderatepriors.rda")

#Posterior draws (I.e., the posterior cumulative density function)
moderate_prior_pdf <- pp_check(mod1.1_moderatepriors, ndraws = 100,
                               type = "dens_overlay") +
  xlim(-0.5, 1.5) + theme_classic()

## 7.4 Prior sensitivity: LOO comparison across prior specifications ----------
loo_prior_sensitivity <- loo_compare(
  loo(mod1.1_weakpriors),
  loo(mod1.1_diffusepriors),
  loo(mod1.1_moderatepriors)
)
saveRDS(loo_prior_sensitivity,
        "scripts/model_outputs/Offspring Trait Models/early life surv/modelpriors2.rda")
#Selecting the model with weakly-informative priors moving forwards

# ===========================================================================
# 8. SELECTED MODEL DIAGNOSTICS  (mod1.1_weakpriors)-------------------------
# ===========================================================================

## LOO-CV --------------------------------------------------------------------
color_scheme_set("blue")
loo_weak    <- loo(mod1.1_weakpriors, save_psis = TRUE)
plot(loo_weak)
psis_weak   <- loo_weak$psis_object
psis_weights <- weights(psis_weak)

## Posterior predictive checks -----------------------------------------------

#posterior simulated and empirical mean
weakprior_mean     <- ppc_stat(yrep = posterior_predict(mod1.1_weakpriors),
                               y    = F1data$early_life_surv,
                               stat = "mean") + theme_classic()

#posterior simulated and empirical bar chart
weakprior_bars     <- pp_check(mod1.1_weakpriors, ndraws = 100,
                               type = "bars") + theme_classic()

#posterior simulated and empirical rootogram
weakprior_rootogram <- pp_check(mod1.1_weakpriors, ndraws = 100,
                                type = "rootogram", style = "hanging") +
  theme_classic()

## Bayes R² ------------------------------------------------------------------
bayes_R2(mod1.1_weakpriors, re.form = NA)   # marginal (fixed effects only)
bayes_R2(mod1.1_weakpriors, re.form = NULL)  # conditional (incl. random effects)
#both fixed and random effects do a poor job at explaining variance in early life survival

## MAP estimates and pd -------------------------------------------------------
MAP_earlylifesurv <- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  test       = "pd",
  ci_method  = "HDI",
  centrality = "MAP",
  effects    = "full",
  ci         = 0.95
)
saveRDS(MAP_earlylifesurv,
        "scripts/model_outputs/Offspring Trait Models/early life surv/MAP_earlylifesurv.rda")

## ROPE -----------------------------------------------------------------------
rope_earlylifesurv <- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  rope_range = c(-0.18, 0.18), 
  test      = "rope",
  ci_method = "HDI",
  ci        = 1 #using the full posterior distribution to estimate the % in ROPE
)
saveRDS(rope_earlylifesurv,
        "scripts/model_outputs/Offspring Trait Models/early life surv/earlylifesurvrope1.rda")

## Selective disappearance test -----------------------------------------------
hypothesis(mod1.1_weakpriors, "avg_age - delta.age > 0")
hypothesis(mod1.1_weakpriors, "avg_age - delta.age < 0")


# ===========================================================================
# 9. PRIOR vs. POSTERIOR SENSITIVITY PLOTS-----------------------------------
# ===========================================================================
## Combined prior–posterior panel ---------------------------------------------
priorpost_panel <- ggarrange(
  diffuse_prior_cumvdis, weak_prior_cumvdis, moderate_prior_cumvdis,
  diffuse_prior_pdf,     weak_prior_pdf,     moderate_prior_pdf,
  nrow = 2, ncol = 3,
  labels = c("A", "B", "C", "D", "E", "F")
)

ggsave("./bayesian_plots/model fit plots/early life surv/priorvpost_pdf.png",
       plot   = priorpost_panel,
       device = "png",
       dpi    = 300,
       width  = 290, height = 140, units = "mm")


# ===========================================================================
# 10. HYPOTHESIS TESTING: INTERACTION MODEL-----------------------------------
# ===========================================================================
#Interaction answers a core question: Does the parents environmental temperature mediate the observed parental age effects?

## 10.1 Two-way interaction: delta age × temperature -------------------------
mod2.1 <- brm(
  early_life_surv ~ avg_age + 
    delta.age +
    Temp +
    delta.age:Temp +
    (1 + delta.age | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  prior     = weak_priors,
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  cores     = 4
)
saveRDS(mod2.1,
        "scripts/model_outputs/Offspring Trait Models/early life surv/mod2.1.rda")
mod2.1 <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/early life surv/mod2.1.rda")

## MAP and ROPE for interaction model -----------------------------------------
MAP_earlylifesurv_mod2.1 <- bayestestR::describe_posterior(
  mod2.1,
  test       = "pd",
  ci_method  = "HDI",
  centrality = "MAP",
  effects    = "full",
  ci         = 0.95
)
saveRDS(MAP_earlylifesurv_mod2.1,
        "scripts/model_outputs/Offspring Trait Models/early life surv/earlylifesurv2rope2.rda")

rope_earlylifesurv_mod2.1 <- bayestestR::describe_posterior(
  mod2.1,
  test      = "rope",
  ci_method = "HDI",
  ci        = 1
)
saveRDS(rope_earlylifesurv_mod2.1,
        "scripts/model_outputs/Offspring Trait Models/early life surv/earlylifesurv2rope.rda")


# ===========================================================================
# 11. MODEL SELECTION---------------------------------------------------------
# ===========================================================================

loo_selected   <- loo(mod1.1_weakpriors)
loo_interaction <- loo(mod2.1)

early_life_fit <- loo_compare(loo_selected, loo_interaction)
saveRDS(early_life_fit,
        "scripts/model_outputs/Offspring Trait Models/early life surv/early_life_fit.rda")

# Selected model: mod1.1_weakpriors (no interaction)
#still reporting the posteriors for the interaction (they answer a core question of the manuscript)

# ===========================================================================
# 12. POSTERIOR DISTRIBUTION PLOTS---------------------------------------------
# ===========================================================================
# Draws from the SELECTED model (mod1.1_weakpriors) for main effects.
# Interaction term draws sourced from the hypothesis model (mod2.1) and
# included in the posterior plot to show the magnitude of the removed terms.

post1 <- as_draws_df(mod1.1_weakpriors)
post2 <- as_draws_df(mod2.1)

## Build combined posterior data frame ----------------------------------------
posterior_df_mod1.1 <- data.frame(
  "μ: Parents' Δage"                   = post1$b_delta.age,
  "μ: Parents' average age"            = post1$b_avg_age,
  "μ: Parent Temperature (28.0°C)"     = post1$b_Temp28,
  "μ: Parent Temperature (30.5°C)"     = post1$b_Temp30.5,
  "μ: Δage × Temperature (28.0°C)"    = post2$`b_delta.age:Temp28`,
  "μ: Δage × Temperature (30.5°C)"    = post2$`b_delta.age:Temp30.5`,
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
    "μ: Δage × Temperature (28.0°C)",
    "μ: Δage × Temperature (30.5°C)"
  )))

y_levels <- rev(c(
  "μ: Parents' Δage",
  "μ: Parents' average age",
  "μ: Parent Temperature (28.0°C)",
  "μ: Parent Temperature (30.5°C)",
  "μ: Δage × Temperature (28.0°C)",
  "μ: Δage × Temperature (30.5°C)"
))

## Halfeye posterior plot ------------------------------------------------------
posterior_plot_earlylifesurv <- ggplot(
  posterior_df_mod1.1,
  aes(x = value, y = parameter, fill = parameter)
) +
  annotate("rect",
           xmin = -0.18, xmax = 0.18,
           ymin = -Inf, ymax = Inf,
           alpha = 0.08, fill = "grey20") +
  stat_halfeye(
    scale          = 0.8,
    adjust         = 1,
    justification  = -0.1,
    .width         = 0.95,
    slab_colour    = "grey85",
    slab_linewidth = 4,
    slab_alpha     = 0.8,
    linewidth      = 10,
    point_size     = 15,
    point_colour   = "white",
    shape          = 21,
    stroke         = 1.5
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey21") +
  scale_fill_manual(values = c(
    "#4A6479", "#A9BAC2", "#5f9289", "#e653a1",
    "#e5a9f5", "#d2e9f5"
  )) +
  xlim(c(-1.0, 1.2)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title      = element_text(size = 30),
    axis.text       = element_text(size = 30)
  ) +
  labs(x = "Fixed-effect posterior estimates\n(μ = mean log-odds of early-life survival)",
       y = NULL)

ggsave("./bayesian_plots/posterior plots/early life surv/001_hypothesis1_halfeye.png",
       plot   = posterior_plot_earlylifesurv,
       device = "png",
       width  = 580, height = 400, units = "mm")


# ===========================================================================
# 13. MODEL PREDICTIONS AND RESULTS PLOTS------------------------------------
# ===========================================================================

## 13.1 Prediction grids: within- and between-individual age effects ---------

# Within-individual: vary delta.age, hold avg_age constant
df_predict_within <- expand.grid(
  avg_age   = mean(F1data$avg_age),
  delta.age = unique(F1data$delta.age),
  Temp      = unique(F1data$Temp),
  PairID    = unique(F1data$PairID)[1]
)

pred_within <- as.data.frame(fitted(mod1.1_weakpriors, df_predict_within,
                                    re_formula = NA))
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
  mutate(within_subject_age = delta.age * sd(F1data$within_subject_age))

# Between-individual: vary avg_age, hold delta.age constant
df_predict_between <- expand.grid(
  avg_age   = unique(F1data$avg_age),
  delta.age = mean(F1data$delta.age),
  Temp      = unique(F1data$Temp),
  PairID    = unique(F1data$PairID)[1]
)

pred_between <- as.data.frame(fitted(mod1.1_weakpriors, df_predict_between,
                                     re_formula = NA))
df_predict_between <- df_predict_between %>%
  mutate(prediction = pred_between$Estimate,
         lower      = pred_between$Q2.5,
         upper      = pred_between$Q97.5)

avg_age_slopes <- df_predict_between %>%
  group_by(avg_age) %>%
  summarise(prediction = mean(prediction),
            lower      = mean(lower),
            upper      = mean(upper)) %>%
  mutate(
    avg_age_raw     = avg_age * sd(F1data$avg.age) + mean(F1data$avg.age),
    avg_age_centred = avg_age_raw - mean(avg_age_raw)
  )

## 13.2 Raw within-individual means ------------------------------------------
breaks_wk <- seq(-4.2, 4.2, by = 1)
labels_wk  <- -4:3

F1data <- F1data %>%
  mutate(delta_age_bin = cut(within_subject_age, breaks = breaks_wk,
                             labels = labels_wk))

raw_deltaage <- F1data %>%
  select(-within_subject_age) %>%
  filter(!is.na(delta_age_bin)) %>%
  rename(within_subject_age = delta_age_bin) %>%
  mutate(within_subject_age = as.numeric(as.character(within_subject_age))) %>%
  group_by(within_subject_age) %>%
  summarise(
    n        = sum(!is.na(early_life_surv)),
    mean_surv = mean(early_life_surv, na.rm = TRUE),
    se_surv  = ifelse(n > 1, sd(early_life_surv, na.rm = TRUE) / sqrt(n),
                      NA_real_),
    .groups  = "drop"
  )

## 13.3 Main parental age effect plot -----------------------------------------
early_surv_plot <- ggplot(F1data,
                          aes(x = within_subject_age, y = early_life_surv)) +
  geom_point(position = position_jitter(width = 0.2, height = 0.05),
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
  geom_line(data = avg_age_slopes,
            aes(x = avg_age_centred, y = prediction,
                linetype = "Between-individual",
                colour   = "Between-individual"),
            linewidth = 4) +
  geom_ribbon(data = avg_age_slopes,
              aes(x = avg_age_centred, y = NULL,
                  ymin   = lower, ymax = upper,
                  fill   = "Between-individual",
                  colour = "Between-individual"),
              alpha = 0.1, linewidth = 1, show.legend = FALSE) +
  geom_linerange(data = raw_deltaage,
                 aes(y = mean_surv,
                     ymin = mean_surv - se_surv,
                     ymax = mean_surv + se_surv,
                     colour = "Within-individual"),
                 linewidth = 3, show.legend = FALSE) +
  geom_point(data = raw_deltaage,
             aes(x = within_subject_age, y = mean_surv,
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
    legend.position   = "top",
    legend.box        = "horizontal",
    legend.title      = element_text(size = 50),
    legend.text       = element_text(size = 50),
    axis.title        = element_text(size = 50),
    axis.text         = element_text(size = 50),
    strip.text        = element_text(size = 50),
    panel.background  = element_rect(fill = "transparent", colour = NA),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)
  ) +
  labs(x = "Parents' adult age at reproduction (weeks; mean-centred)",
       y = "Early-life survival probability")

ggsave("./bayesian_plots/offspring trait plots/F1 early life survival/001_early_life_surv_plot.png",
       plot = early_surv_plot, bg = "transparent",
       device = "png", width = 600, height = 500, units = "mm")

## 13.4 Posterior ridge density (PDF) plot ------------------------------------
newdat_pdf <- expand.grid(
  avg_age                   = mean(F1data$avg_age),
  delta.age                 = c(-1.95, 0, 1.85),
  Temp                      = unique(F1data$Temp),
  PairID                    = unique(F1data$PairID)[1]
)

proportional_posterior <- mod1.1_weakpriors %>%
  linpred_draws(
    newdat_pdf,
    value            = "mu",
    allow_new_levels = TRUE,
    transform        = TRUE,
    re_formula       = NA,
    ndraws           = 3000,
    seed             = 123
  )

summary_probability <- proportional_posterior %>%
  ungroup() %>%
  group_by(delta.age, .draw) %>%
  summarise(mu = mean(mu), .groups = "drop") %>%
  mutate(timepoint = case_when(
    delta.age == -1.95 ~ "Early-Aged",
    delta.age ==  0    ~ "Middle-Aged",
    delta.age ==  1.85 ~ "Late-Aged"
  )) %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("Early-Aged", "Middle-Aged", "Late-Aged")))

surv_density <- ggplot(summary_probability,
                       aes(x = mu, y = timepoint,
                           colour = timepoint, fill = timepoint,
                           alpha  = timepoint)) +
  geom_density_ridges(linewidth = 2, scale = 0.8,
                      rel_min_height = 0) +
  scale_fill_manual(values  = c("Early-Aged"  = "#f96161",
                                "Middle-Aged" = "#66b2b2",
                                "Late-Aged"   = "#066594")) +
  scale_colour_manual(values = c("Early-Aged"  = "#f96161",
                                 "Middle-Aged" = "#66b2b2",
                                 "Late-Aged"   = "#066594")) +
  scale_alpha_manual(values  = c("Early-Aged"  = 0.6,
                                 "Middle-Aged" = 0.1,
                                 "Late-Aged"   = 0.6)) +
  scale_x_continuous(limits = c(0.65, 0.85)) +
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
    legend.title      = element_text(size = 50),
    legend.text       = element_text(size = 50),
    axis.title        = element_text(size = 50),
    axis.text         = element_text(size = 50),
    axis.line.y       = element_blank(),
    panel.background  = element_rect(fill = "transparent", colour = NA),
    plot.background   = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA)
  ) +
  labs(x     = "Early-life survival probability",
       y     = "Probability density, f(x)",
       fill  = "Parents' age at reproduction",
       color = "Parents' age at reproduction",
       alpha = "Parents' age at reproduction")

ggsave("./bayesian_plots/offspring trait plots/F1 early life survival/003_pdf.png",
       plot = surv_density, bg = "transparent",
       device = "png", width = 520, height = 400, units = "mm")

## 13.5 Combined inference panel (main effect + PDF) --------------------------
inference_panel <- ggarrange(
  surv_density, early_surv_plot,
  ncol = 2, nrow = 1,
  labels     = c("A", "B"),
  font.label = list(size = 50, face = "bold"),
  label.x    = c(0.05, 0.05),
  widths     = c(0.9, 1.2)
)

ggsave("./bayesian_plots/offspring trait plots/F1 early life survival/inference_plots.png",
       plot   = inference_panel,
       device = "png",
       width  = 980, height = 570, units = "mm")


# ===========================================================================
# 14. STANDARDISED MEAN DIFFERENCE (SMD)
# ===========================================================================
#for post-hoc inference, creates a standardised effect size

pairwise_estimates <- emmeans(mod1.1_weakpriors, "delta.age", type = "response",
                              at = list(delta.age = c(-1.96, 1.96)))
pairs(pairwise_estimates)

sd_surv <- sd(F1data$early_life_surv)

early_life_surv_SMD <- summary(pairwise_estimates) %>%
  summarise(
    SMD       = (response[2]    - response[1])    / sd_surv,
    SMD_lower = (lower.HPD[2]   - upper.HPD[1])   / sd_surv,
    SMD_upper = (upper.HPD[2]   - lower.HPD[1])   / sd_surv
  )

saveRDS(early_life_surv_SMD,
        "scripts/model_outputs/Offspring Trait Models/early_life_surv_SMD")


#Creating a table to export------------------------------------------------------------------------------

# ── helper ────────────────────────────────────────────────────────────────────
#Model draws from selected models
base_model_draws<-as_draws_df(mod1.1_weakpriors) #just single effects, and the final selected model
interaction_draws<-as_draws_df(mod2.1) 

#estimating the difference between average age and delta age
differences<-data.frame(
  selectivedis = base_model_draws$b_avg_age - base_model_draws$b_delta.age
)

#function for generating estimates
summarise_param <- function(draws, param, section,
                            rope = TRUE, 
                            rope_range = c(-0.18, +0.18),
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
  summarise_param(differences$selectivedis, "Average age - Δage", "Location submodel (Fixed effects)", rope = FALSE)
)


# ── 1. dropped interaction (location) — ROPE + pd ───────────────────────────────────
fe_interaction <- bind_rows(summarise_param(interaction_draws$`b_delta.age:Temp28`, "Parents' Δage x  Temperature: 28.0°C", "Location submodel (Dropped interactions)", rope = TRUE),
                             summarise_param(interaction_draws$`b_delta.age:Temp30.5`, "Parents' Δage x  Temperature: 30.5°C", "Location submodel (Dropped interactions)", rope = TRUE),
)

# ── 3.Random effects ──────────────────────────────────────────────────
# SDs 
re_sd_explore<- bind_rows(
  summarise_param(base_model_draws$sd_PairID__Intercept, "σ intercept",    "Random effects (Location)", pd = FALSE),
  summarise_param(base_model_draws$sd_PairID__delta.age, "σ slope Δage",    "Random effects (Location)", pd = FALSE),
  summarise_param(base_model_draws$cor_PairID__Intercept__delta.age, "r intercept ~ slope Δage",    "Random effects (Location)", pd = TRUE)
)


#---5. Conditional and marginal Bayes R² ------------------------------------------------------------------
marginal<-bayes_R2(mod1.1_weakpriors, re.form = NA, summary = FALSE)   # fixed effects only (i.e., marginal)
conditional<-bayes_R2(mod1.1_weakpriors, re.form = NULL, summary = FALSE)  # including random effects (i.e., conditional)

bayes<-bind_rows(
  summarise_param(marginal, "Marginal R²", "Bayes R²", pd= FALSE, rope= FALSE),
  summarise_param(conditional, "Conditional R²", "Bayes R²", pd= FALSE, rope =FALSE)
)


# ── 5. Combine and render ─────────────────────────────────────────────────────
bind_rows(fe, fe_interaction, re_sd_explore, bayes) %>%
  gt(groupname_col = "Section") %>%
  cols_label(Parameter = "Parameter", MAP = "MAP",
             `95% HDI` = "95% HDI", `% in ROPE` = "% in ROPE", pd = "pd") %>%
  tab_header(
    title    = "Early-life survival model summary: Bernoulli model"
  ) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) %>%
  tab_footnote("ROPE = [−0.18, 0.18] on log-odds scale; applied to location fixed effects only.",
               locations = cells_column_labels(`% in ROPE`)) %>%
  tab_footnote("pd = probability of direction.",
               locations = cells_column_labels(pd)) %>%
  opt_stylize(style = 1)%>%
  gtsave("./tables/earlysurvival_updated.docx")


###############################################################################
##  END OF SCRIPT----------------------------------------------------------------
###############################################################################
