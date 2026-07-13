###############################################################################
##  Script: Effect of parental age on F1 juvenile survival (survival from week four to adulthood)
##  Note:   Juvenile survival = surviving to adulthood, conditional on
##          having survived the first four weeks of life (i.e., early-life survival)
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

## Stan backend ---------------------------------------------------------------
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

## Colour scheme for bayesplot ------------------------------------------------
color_scheme_set("teal")


# ===========================================================================
# 2. DATA ORGANISATION
# ===========================================================================

## Load data ------------------------------------------------------------------
F1data <- readRDS('./raw data/F1_filtered_data06022025.RDS')

## Create early-life survival indicator and filter to early-life survivors ----
# Offspring that died in weeks 1–4 are excluded from this analysis;
F1data <- F1data %>%
  mutate(early_life_surv = ifelse(total_lifespan <= 4, 0, 1)) %>%
  filter(!early_life_surv == 0) #remove animals that fell out from the study in the early-life analysis

#Juvenile survival coded as the variable "Adult survival" column already coded into raw spreadsheet
#This column represents whether an individual survived to adulthood (1 = yes, 0 = no)

## Scaling continuous predictors ------------------------------------------------
#avg_age: between-individual (mean and SD centred)
F1data$avg_age <- as.numeric(scale(F1data$avg.age, center = TRUE, scale = TRUE))

# delta.age: within-individual (scaled by SD but NOT mean-centred to preserve
# the within-subject contrast); NAs replaced with 0
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
length(unique(F1data$PairID))  # 77 parent pairs
length(F1data$F1_ID)           # 987 offspring


# ===========================================================================
# 3. SUMMARY STATISTICS--------------------------------------------------------
# ===========================================================================
ggplot(F1data, aes(x = adult_surv)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +  # Create histogram
  facet_wrap(~Timepoint) +  
  labs(
    x = "Juvenile survival",
    y = "Count"
  ) +
  theme_classic() #Very little chance of not making it to adulthood if you survived your first four weeks
#Most individuals that survived their first month make it to adulthood

## Overall juvenile survival rate ---------------------------------------------
sum_dat <- F1data %>%
  summarise(
    mean_surv    = mean(adult_surv, na.rm = TRUE),
    sd_surv      = sd(adult_surv, na.rm = TRUE),
    n_offspring  = n(),
    n_pairs      = n_distinct(PairID)
  )
# mean_surv ≈ 0.959; SD ≈ 0.197; n = 987; pairs = 77

## Juvenile survival by parental age timepoint --------------------------------
sum_dat2 <- F1data %>%
  group_by(Timepoint) %>%
  summarise(
    mean_surv   = mean(adult_surv, na.rm = TRUE),
    sd_surv     = sd(adult_surv, na.rm = TRUE),
    n_offspring = n(),
    n_pairs     = n_distinct(PairID)
  )

# Logit-scale intercept used to set the prior
qlogis(0.959)  # ≈ 3.152319


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
vif(lm(adult_surv ~ avg_age + delta.age + cum_successful_matings+
         Temp, data = F1data))

#removing cum_succesful matings (cumulative number of succesful mating attempts)
vif(lm(adult_surv ~ avg_age + delta.age +
         Temp, data = F1data))


# ===========================================================================
# 4. MODEL FAMILY SELECTION
# ===========================================================================
# Compare: random intercept only vs. random slopes for delta.age;
# linear vs. quadratic delta.age term — all with default priors (for now)

## 4.1. Base model: random intercept -------------------------------------------
mod1.1_basebernoulli <- brm(
  adult_surv ~ avg_age + 
    delta.age + 
    Temp + 
    (1 | PairID),
  family  = bernoulli(link = "logit"),
  data    = F1data,
  save_pars = save_pars(all = TRUE),
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  cores   = 4
)
saveRDS(mod1.1_basebernoulli,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_basebernoulli.rda")
mod1.1_basebernoulli <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_basebernoulli.rda")

## 4.2. Quadratic delta.age term -----------------------------------------------
mod1.1_quadratic <- brm(
  adult_surv ~ avg_age + 
    poly(delta.age, 2) +
    Temp +
    (1 | PairID),
  family  = bernoulli(link = "logit"),
  data    = F1data,
  save_pars = save_pars(all = TRUE),
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  cores   = 4
)
saveRDS(mod1.1_quadratic,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_quadratic.rda")
mod1.1_quadratic <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_quadratic.rda")


## 4.3. Random slopes + quadratic ----------------------------------------------
#females aren't allowed to vary in the shape of their random slopes (assumes a shared shape across females)
mod1.1_randomslopes <- brm(
  adult_surv ~ avg_age + 
    poly(delta.age,2) +
    Temp +
    (1 + delta.age | PairID),
  family  = bernoulli(link = "logit"),
  data    = F1data,
  save_pars = save_pars(all = TRUE),
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  cores   = 4
)
saveRDS(mod1.1_randomslopes,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_randomslopes.rda")
mod1.1_randomslopes <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_randomslopes.rda")

## 4.4. Random slopes + linear delta.age (selected structure) ------------------
mod1.1_randomslope_linear <- brm(
  adult_surv ~ avg_age + 
    delta.age + 
    Temp + 
    (1 + delta.age | PairID),
  family  = bernoulli(link = "logit"),
  data    = F1data,
  save_pars = save_pars(all = TRUE),
  iter    = 5000,
  control = list(adapt_delta = 0.98),
  cores   = 4
)
saveRDS(mod1.1_randomslope_linear,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_randomslope_linear.rda")
mod1.1_randomslope_linear <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_randomslope_linear.rda")

## LOO comparison of model structures -----------------------------------------
modelfamily2 <- loo_compare(
  loo(mod1.1_randomslope_linear, save_psis = TRUE),
  loo(mod1.1_randomslopes,       save_psis = TRUE),
  loo(mod1.1_quadratic,          save_psis = TRUE),
  loo(mod1.1_basebernoulli,      save_psis = TRUE)
)
print(modelfamily2, simplify = FALSE)
saveRDS(modelfamily2,
        "scripts/model_outputs/Offspring Trait Models/adult surv/modelfamily2.rda")
# Selected structure: linear delta.age + random slopes → mod1.1_randomslope_linear


# ===========================================================================
# 5. PRIOR SPECIFICATION-------------------------------------------------------
# ===========================================================================
# Intercept centred at qlogis(0.959) ≈ 3.152319 (empirical juvenile survival rate)
# All predictors are z-scored; priors on log-odds scale.

## Diffuse priors -------------------------------------------------------------
diffusepriors <- c(
  prior(normal(3.152319, 1.5), class = "Intercept"),
  prior(normal(0, 1.5),        class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd", lb = 0),
  prior(lkj(2),                class = "cor")
)

## Weakly informative priors (selected) ---------------------------------------
weakpriors <- c(
  prior(normal(3.152319, 1), class = "Intercept"),
  prior(normal(0, 1),        class = "b"),
  prior(exponential(1),      class = "sd", lb = 0),
  prior(lkj(2),              class = "cor")
)

## Moderate priors ------------------------------------------------------------
moderatepriors <- c(
  prior(normal(3.152319, 0.5), class = "Intercept"),
  prior(normal(0, 0.5),        class = "b"),
  prior(exponential(3),        class = "sd", lb = 0),
  prior(lkj(2),                class = "cor")
)


# ===========================================================================
# 6. PRIOR PREDICTIVE CHECKS--------------------------------------------------
# ===========================================================================

## 6.1. Diffuse prior-only model -----------------------------------------------
diffuse_prior_model <- brm(
  adult_surv ~ avg_age + delta.age + Temp + (1 + delta.age | PairID),
  family        = bernoulli(link = "logit"),
  data          = F1data,
  prior         = diffusepriors,
  sample_prior  = "only",
  save_pars     = save_pars(all = TRUE),
  iter          = 5000,
  control       = list(adapt_delta = 0.98),
  cores         = 4
)
saveRDS(diffuse_prior_model,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/diffuse_prior_model.rda")
diffuse_prior_model <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult surv/diffuse_prior_model.rda")

#Prior draws (I.e., the prior cumulative density function)
diffuse_prior_ppc <- pp_check(diffuse_prior_model, ndraws = 100) +
  xlim(-0.5, 1.5) + theme_classic()

## 6.2. Weak prior-only model --------------------------------------------------
weak_prior_model <- brm(
  adult_surv ~ avg_age + delta.age + Temp + (1 + delta.age | PairID),
  family        = bernoulli(link = "logit"),
  data          = F1data,
  prior         = weakpriors,
  sample_prior  = "only",
  save_pars     = save_pars(all = TRUE),
  iter          = 5000,
  control       = list(adapt_delta = 0.98),
  cores         = 4
)
saveRDS(weak_prior_model,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/weak_prior_model.rda")
weak_prior_model <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult surv/weak_prior_model.rda")

#Prior draws (I.e., the prior cumulative density function)
weak_prior_ppc <- pp_check(weak_prior_model, ndraws = 100, type = "dens_overlay") +
  xlim(-0.5, 1.5) + theme_classic()

## 6.3. Moderate prior-only model ----------------------------------------------
moderate_prior_model <- brm(
  adult_surv ~ avg_age + delta.age + Temp + (1 + delta.age | PairID),
  family        = bernoulli(link = "logit"),
  data          = F1data,
  prior         = moderatepriors,
  sample_prior  = "only",
  save_pars     = save_pars(all = TRUE),
  iter          = 5000,
  control       = list(adapt_delta = 0.98),
  cores         = 4
)
saveRDS(moderate_prior_model,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/moderate_prior_model.rda")
moderate_prior_model <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult surv/moderate_prior_model.rda")

#Prior draws (I.e., the prior cumulative density function)
moderate_prior_ppc <- pp_check(moderate_prior_model, ndraws = 100) +
  xlim(-0.5, 1.5) + theme_classic()

# ===========================================================================
# 7.POSTERIOR DRAWS-----------------------------------------------------------
# ===========================================================================

## 7.1. Diffuse priors ---------------------------------------------------------
mod1.1_diffusepriors <- brm(
  adult_surv ~ avg_age + delta.age + Temp + (1 + delta.age | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  prior     = diffusepriors,
  save_pars = save_pars(all = TRUE),
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  cores     = 4
)
saveRDS(mod1.1_diffusepriors,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_diffusepriors.rda")
mod1.1_diffusepriors <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_diffusepriors.rda")

color_scheme_set("pink")

#Posterior draws (I.e., the posterior cumulative density function)
diffuse_prior_pdf <- pp_check(mod1.1_diffusepriors, ndraws = 100) +
  xlim(-0.5, 1.5) + theme_classic()



## 7.2. Weakly informative priors (SELECTED MODEL) -----------------------------
mod1.1_weakpriors <- brm(
  adult_surv ~ avg_age + delta.age + Temp + (1 + delta.age | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  prior     = weakpriors,
  save_pars = save_pars(all = TRUE),
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  cores     = 4
)
saveRDS(mod1.1_weakpriors,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_weakpriors.rda")
mod1.1_weakpriors <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_weakpriors.rda")

#Posterior draws (I.e., the posterior cumulative density function)
weak_prior_pdf <- pp_check(mod1.1_weakpriors, ndraws = 100) +
  xlim(-0.5, 1.5) + theme_classic()


## 7.3. Moderate priors --------------------------------------------------------
mod1.1_moderatepriors <- brm(
  adult_surv ~ avg_age + delta.age + Temp + (1 + delta.age | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  prior     = moderatepriors,
  save_pars = save_pars(all = TRUE),
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  cores     = 4
)
saveRDS(mod1.1_moderatepriors,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_moderatepriors.rda")
mod1.1_moderatepriors <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult surv/mod1.1_moderatepriors.rda")

moderate_prior_pdf <- pp_check(mod1.1_moderatepriors, ndraws = 100) +
  xlim(-0.5, 1.5) + theme_classic()


# ===========================================================================
# 8. PRIOR SENSITIVITY ANALYSIS-----------------------------------------------
# ===========================================================================

## LOO comparison across prior specifications ---------------------------------
loomodweak     <- loo(mod1.1_weakpriors,    save_psis = TRUE)
loodiffuse     <- loo(mod1.1_diffusepriors, save_psis = TRUE)
loomodmoderate <- loo(mod1.1_moderatepriors, save_psis = TRUE)
loomodstruct   <- loo(mod1.1_randomslope_linear, save_psis = TRUE)

modelpriors_compare <- loo_compare(loomodweak, loodiffuse,
                                   loomodmoderate, loomodstruct)
print(modelpriors_compare, simplify = FALSE)
saveRDS(modelpriors_compare,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/modelpriors_compare.rda")

# ===========================================================================
# PRIOR vs. POSTERIOR SENSITIVITY PLOTS-----------------------------------
# ===========================================================================
## Combined prior–posterior panel ---------------------------------------------
priorpost_panel <- ggarrange(
  diffuse_prior_ppc , weak_prior_ppc, moderate_prior_ppc,
  diffuse_prior_pdf,     weak_prior_pdf,     moderate_prior_pdf,
  nrow = 2, ncol = 3,
  labels = c("A", "B", "C", "D", "E", "F")
)

ggsave("./bayesian_plots/model fit plots/adult surv/priorvpost_pdf.png",
       plot   = priorpost_panel,
       device = "png",
       dpi    = 300,
       width  = 290, height = 140, units = "mm")


# ===========================================================================
# 9. SELECTED MODEL DIAGNOSTICS  (mod1.1_weakpriors)-------------------------
# ===========================================================================

## LOO-CV --------------------------------------------------------------------
color_scheme_set("blue")
loo_weak    <- loo(mod1.1_weakpriors, save_psis = TRUE)
plot(loo_weak)
psis_weak   <- loo_weak$psis_object
psis_weights <- weights(psis_weak)

#Posterior probability density
weak_post_ppc  <- pp_check(mod1.1_weakpriors, ndraws = 100, type = "dens_overlay") +
  xlim(-0.5, 1.5) + theme_classic()

#posterior simulated and empirical mean
weakprior_mean<- ppc_stat(yrep = posterior_predict(mod1.1_weakpriors),
                               y    = F1data$adult_surv,
                               stat = "mean") + theme_classic()

#posterior simulated and empirical bar chart
weak_post_bars <- pp_check(mod1.1_weakpriors, ndraws = 100, type = "bars") +
  theme_classic()

#posterior simulated and empirical rootogram
weakprior_rootogram <- pp_check(mod1.1_weakpriors, ndraws = 100,
                                type = "rootogram", style = "hanging") +
  theme_classic()

## Bayes R² -------------------------------------------------------------------
bayes_R2(mod1.1_weakpriors, re.form = NA)   # population-level (i.e., marginal)
bayes_R2(mod1.1_weakpriors, re.form = NULL) # including random effects (i.e., conditional)
#Very little of the total variation in juvenile survival explained by the fixed and random effects

# ===========================================================================
# 9. HYPOTHESIS TESTING — TEMPERATURE × PARENTAL AGE INTERACTION
# ===========================================================================
#The core interaction: Does temperature mediate the average parental age effect?

## mod2.1: delta.age × Temperature interaction --------------------------------
mod2.1 <- brm(
  adult_surv ~ avg_age + delta.age + Temp + delta.age:Temp + (1 + delta.age | PairID),
  family    = bernoulli(link = "logit"),
  data      = F1data,
  prior     = weakpriors,
  save_pars = save_pars(all = TRUE),
  iter      = 5000,
  control   = list(adapt_delta = 0.98),
  cores     = 4
)
saveRDS(mod2.1,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/mod2.1.rda")
mod2.1 <- readRDS(
  "scripts/model_outputs/Offspring Trait Models/adult surv/mod2.1.rda")

summary(mod2.1)
parameters::parameters(mod2.1)

## MAP + ROPE for interaction model -------------------------------------------
bayestestR::describe_posterior(
  mod2.1,
  test       = "pd",
  centrality = "MAP",
  ci         = 0.95,
  ci_method  = "HDI",
  effects    = "full"
)
bayestestR::describe_posterior(
  mod2.1,
  rope_range = c(-0.18, 0.18),
  test      = "rope",
  ci        = 1,
  ci_method = "HDI"
)


# ===========================================================================
# 10. MODEL SELECTION — LOO-CV
# ===========================================================================

## Compare selected structure (mod1.1_weakpriors) vs. interaction (mod2.1) ----
loomodweak        <- loo(mod1.1_weakpriors, save_psis = TRUE)
loomodinteraction <- loo(mod2.1,            save_psis = TRUE)

hypothesisloo <- loo_compare(loomodweak, loomodinteraction)
print(hypothesisloo, simplify = FALSE)
saveRDS(hypothesisloo,
        file = "scripts/model_outputs/Offspring Trait Models/adult surv/adult_fit.rda")
# mod1.1_weakpriors selected: simpler model with equivalent predictive fit


# ===========================================================================
# 11. MAP ESTIMATES AND ROPE
# ===========================================================================

## MAP estimates with 95 % HDI (selected model) --------------------------------
MAPS <- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  test       = "pd",
  centrality = "MAP",
  ci         = 0.95,
  ci_method  = "HDI",
  effects    = "full"
)
print(MAPS)

## ROPE (100 % HDI — full posterior in ROPE) ----------------------------------
# ROPE range derives from the default ±0.18 SD rule on the logit scale
Inrope <- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  test      = "rope",
  rope_range = c(-0.18, 0.18),
  ci        = 1,
  ci_method = "HDI"
)
print(Inrope)


# ===========================================================================
# 12. POSTERIOR DISTRIBUTIONS--------------------------------------------------
# ===========================================================================
# The stat_halfeye plot below combines:
#   • Main-effect posteriors from mod1.1_weakpriors (selected model)
#   • Interaction posteriors from mod2.1 (hypothesis model, not selected)

## Extract posterior draws -----------------------------------------------------
post1 <- as_draws_df(mod1.1_weakpriors)  # selected model
post2 <- as_draws_df(mod2.1)             # interaction model

## Assemble long-format data frame --------------------------------------------
posterior_df <- data.frame(
  "μ: Parents' average age"             = post1$b_avg_age,
  "μ: Parents' Δage"                    = post1$b_delta.age,
  "μ: Temperature (28.0°C)"             = post1$b_Temp28,
  "μ: Temperature (30.5°C)"             = post1$b_Temp30.5,
  "μ: Δage × Temperature (28.0°C)"      = post2$`b_delta.age:Temp28`,
  "μ: Δage × Temperature (30.5°C)"      = post2$`b_delta.age:Temp30.5`,
  check.names = FALSE
)

posterior_long <- posterior_df %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "parameter",
    values_to = "draw"
  ) %>%
  mutate(
    parameter = factor(parameter, levels = rev(c(
      "μ: Δage × Temperature (28.0°C)",
      "μ: Δage × Temperature (30.5°C)",
      "μ: Temperature (28.0°C)",
      "μ: Temperature (30.5°C)",
      "μ: Parents' average age",
      "μ: Parents' Δage"
    ))))

y_levels <- rev(c(
  "μ: Parents' Δage",
  "μ: Parents' average age",
  "μ: Parent Temperature (28.0°C)",
  "μ: Parent Temperature (30.5°C)",
  "μ: Δage × Temperature (28.0°C)",
  "μ: Δage × Temperature (30.5°C)"
))

## stat_halfeye posterior plot -------------------------------------------------
posterior_plot_juvsurv <- ggplot(
  posterior_long,
  aes(x = draw, y = parameter, fill = parameter)
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
  xlim(c(-2.0, 1.5)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title      = element_text(size = 30),
    axis.text       = element_text(size = 30)
  ) +
  labs(x = "Fixed-effect posterior estimates\n(μ = mean log-odds of juvenile survival)",
       y = NULL)

ggsave("./bayesian_plots/posterior plots/adult surv/001_hypothesis1_halfeye.png",
       plot   = posterior_plot_juvsurv,
       device = "png",
       width  = 580, height = 380, units = "mm")




# ===========================================================================
# 13. MODEL PREDICTIONS-------------------------------------------------------
# ===========================================================================

## 13.1. Within-individual predictions (varying delta.age, fixed avg_age) ------
df_predict_within <- expand.grid(
  avg_age = mean(F1data$avg_age),
  delta.age = unique(F1data$delta.age),
  Temp    = unique(F1data$Temp),
  PairID  = unique(F1data$PairID)[1]
)

pred_within <- fitted(mod1.1_weakpriors, df_predict_within, re_formula = NA)
df_predict_within$prediction <- pred_within[, "Estimate"]
df_predict_within$lower      <- pred_within[, "Q2.5"]
df_predict_within$upper      <- pred_within[, "Q97.5"]

newdat_within <- df_predict_within %>%
  group_by(delta.age) %>%
  summarise(
    mean  = mean(prediction),
    lower = mean(lower),
    upper = mean(upper),
    .groups = "drop"
  ) %>%
  mutate(within_subject_age = delta.age * sd(F1data$within_subject_age))

## 13.2. Between-individual predictions (varying avg_age, fixed delta.age) ------
df_predict_between <- expand.grid(
  avg_age   = unique(F1data$avg_age),
  delta.age = mean(F1data$delta.age),
  Temp      = unique(F1data$Temp),
  PairID    = unique(F1data$PairID)[1]
)

pred_between <- fitted(mod1.1_weakpriors, df_predict_between, re_formula = NA)
df_predict_between$prediction <- pred_between[, "Estimate"]
df_predict_between$lower      <- pred_between[, "Q2.5"]
df_predict_between$upper      <- pred_between[, "Q97.5"]

avg_age_slopes <- df_predict_between %>%
  group_by(avg_age) %>%
  summarise(
    prediction = mean(prediction),
    lower      = mean(lower),
    upper      = mean(upper),
    .groups    = "drop"
  ) %>%
  mutate(
    avg.age2               = avg_age * sd(F1data$avg.age) + mean(F1data$avg.age),
    avg_timepoint_centered = avg.age2 - mean(avg.age2)
  )

## 13.3. Raw data summaries for overlay ----------------------------------------
breaks <- seq(-4.2, 4.2, by = 1)
labels <- -4:3.8

F1data <- F1data %>%
  mutate(delta_age_bin = cut(within_subject_age, breaks = breaks, labels = labels))

raw_deltaage <- F1data %>%
  filter(!is.na(delta_age_bin)) %>%
  mutate(within_subject_age = as.numeric(as.character(delta_age_bin))) %>%
  group_by(within_subject_age) %>%
  summarise(
    n         = sum(!is.na(adult_surv)),
    mean_surv = mean(adult_surv, na.rm = TRUE),
    se_surv   = ifelse(n > 1, sd(adult_surv, na.rm = TRUE) / sqrt(n), NA_real_),
    .groups   = "drop"
  )

## 13.4. Prediction plot --------------------------------------------------------
adult_surv_plot <- ggplot(F1data,
                          aes(x = within_subject_age, y = adult_surv)) +
  geom_point(
    position = position_jitter(width = 0.2, height = 0.05),
    shape = 21, size = 6, stroke = 1.8,
    alpha = 0.7, colour = "white", fill = "grey"
  ) +
  geom_line(
    data = newdat_within,
    aes(x = within_subject_age, y = mean,
        colour = "Within-individual", linetype = "Within-individual"),
    linewidth = 5
  ) +
  geom_ribbon(
    data = newdat_within,
    aes(x = within_subject_age, y = mean, ymin = lower, ymax = upper,
        fill = "Within-individual"),
    alpha = 0.2, show.legend = FALSE
  ) +
  geom_line(
    data = avg_age_slopes,
    aes(x = avg_timepoint_centered, y = prediction,
        colour = "Between-individual", linetype = "Between-individual"),
    linewidth = 4
  ) +
  geom_ribbon(
    data = avg_age_slopes,
    aes(x = avg_timepoint_centered,  y = prediction, ymin = lower, ymax = upper,
        fill = "Between-individual"),
    alpha = 0.1, show.legend = FALSE
  ) +
  geom_linerange(
    data = raw_deltaage,
    aes(x = within_subject_age, y = mean_surv,
        ymin = mean_surv - se_surv, ymax = mean_surv + se_surv,
        colour = "Within-individual"),
    linewidth = 3, show.legend = FALSE
  ) +
  geom_point(
    data = raw_deltaage,
    aes(x = within_subject_age, y = mean_surv, fill = "Within-individual"),
    shape = 21, stroke = 1.8, size = 18, colour = "white", show.legend = FALSE
  ) +
  scale_fill_manual(
    name   = "Parental age effect",
    values = c("Within-individual" = "#4A6479", "Between-individual" = "#8C2F4B")
  ) +
  scale_colour_manual(
    name   = "Parental age effect",
    values = c("Within-individual" = "#4A6479", "Between-individual" = "#8C2F4B")
  ) +
  scale_linetype_manual(
    name   = "Parental age effect",
    values = c("Within-individual" = "solid", "Between-individual" = "dashed")
  ) +
  scale_x_continuous(breaks = -4:4) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(
    x = "Parents' adult age at reproduction (weeks; mean-centred)",
    y = "Juvenile survival probability"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    legend.box      = "horizontal"
  )

ggsave(
  filename = "./bayesian_plots/offspring trait plots/F1 adult survival/001_adult_surv_plot.png",
  plot     = adult_surv_plot,
  bg       = "transparent",
  device   = "png", dpi = 300,
  width = 600, height = 500, units = "mm"
)


# ===========================================================================
# 14. PROBABILITY DENSITY — POSTERIOR PREDICTIONS
# ===========================================================================
# Ridge density plots showing posterior-predicted juvenile survival probability
# for offspring of early-, middle-, and late-aged parents.

## Generate posterior predictions ----------------------------------------------
newdat_pdf <- expand.grid(
  avg_age   = mean(F1data$avg_age),
  delta.age = c(-1.95, 0, 1.85),
  Temp      = unique(F1data$Temp),
  PairID    = unique(F1data$PairID)[1]
)

proportional_posterior <- mod1.1_weakpriors %>%
  linpred_draws(
    newdat_pdf,
    value              = "mu",
    allow_new_levels   = TRUE,
    transform          = TRUE,
    re_formula         = NA,
    ndraws             = 3000,
    seed               = 123
  )

## Summarise across temperature -----------------------------------------------
posterior_draws_pdf <- proportional_posterior %>%
  ungroup() %>%
  group_by(delta.age, .draw) %>%
  summarise(mu = mean(mu), .groups = "drop")

summary_probability <- posterior_draws_pdf %>%
  mutate(
    timepoint = case_when(
      delta.age == -1.95 ~ "Early-Aged",
      delta.age ==  0    ~ "Middle-Aged",
      delta.age ==  1.85 ~ "Late-Aged"
    ),
    timepoint = factor(timepoint, levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))
  )

## Ridge density plot ----------------------------------------------------------
juvsurv_density <- ggplot(
  summary_probability,
  aes(x = mu, y = timepoint,
      colour = timepoint, fill = timepoint, alpha = timepoint)
) +
  geom_density_ridges(linewidth = 2, scale = 0.8, rel_min_height = 0) +
  scale_fill_manual(
    values = c("Early-Aged" = "#f96161", "Middle-Aged" = "#66b2b2", "Late-Aged" = "#066594")
  ) +
  scale_colour_manual(
    values = c("Early-Aged" = "#f96161", "Middle-Aged" = "#66b2b2", "Late-Aged" = "#066594")
  ) +
  scale_alpha_manual(
    values = c("Early-Aged" = 0.6, "Middle-Aged" = 0.1, "Late-Aged" = 0.6)
  ) +
  labs(
    x      = "Juvenile survival probability",
    y      = "Probability density, f(x)",
    fill   = "Parents' age at reproduction",
    colour = "Parents' age at reproduction",
    alpha  = "Parents' age at reproduction"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    axis.line.y     = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )

ggsave(
  filename = "./bayesian_plots/offspring trait plots/F1 adult survival/003_juvsurv_pdf.png",
  plot     = juvsurv_density,
  bg       = "transparent",
  device   = "png", dpi = 300,
  width = 520, height = 400, units = "mm"
)

## Combined inference panel ---------------------------------------------------
inference_panel <- ggarrange(
  juvsurv_density, adult_surv_plot,
  ncol       = 2, nrow = 1,
  labels     = c("A", "B"),
  font.label = list(size = 50, face = "bold"),
  widths     = c(0.9, 1.2)
)

ggsave(
  filename = "./bayesian_plots/offspring trait plots/F1 adult survival/inference_plots.png",
  plot     = inference_panel,
  device   = "png", dpi = 300,
  width = 980, height = 570, units = "mm"
)


# ===========================================================================
# 15. EFFECT SIZES
# ===========================================================================

## Standardised mean difference (unit scale) ----------------------------------
# Compare offspring from early-aged vs. late-aged parents (±1.96 SD of delta.age)
pairwise_estimates <- emmeans(
  mod1.1_weakpriors, "delta.age",
  type = "response",
  at   = list(delta.age = c(-1.96, 1.96))
)

sd_juvsurv <- sd(F1data$adult_surv, na.rm = TRUE)

juvsurv_SMD <- summary(pairwise_estimates) %>%
  summarise(
    SMD       = (response[2] - response[1]) / sd_juvsurv,
    SMD_lower = (lower.HPD[2] - upper.HPD[1]) / sd_juvsurv,
    SMD_upper = (upper.HPD[2] - lower.HPD[1]) / sd_juvsurv
  )

print(juvsurv_SMD)

## Save SMD for combined forest plot ------------------------------------------
saveRDS(juvsurv_SMD,
        file = "scripts/model_outputs/Offspring Trait Models/juvsurv_SMD")

## Between- vs. within-individual contrast ------------------------------------
# Test whether average vs. within-individual parental age effects differ
hypothesis(mod1.1_weakpriors, "avg_age - delta.age > 0")
hypothesis(mod1.1_weakpriors, "avg_age - delta.age < 0")

## Temperature pairwise contrasts ---------------------------------------------
pairwise_temp <- emmeans(mod1.1_weakpriors, "Temp", type = "response")
pairs(pairwise_temp)


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
    title    = "Juvenile survival model summary: Bernoulli model"
  ) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) %>%
  tab_footnote("ROPE = [−0.18, 0.18] on log-odds scale; applied to location fixed effects only.",
               locations = cells_column_labels(`% in ROPE`)) %>%
  tab_footnote("pd = probability of direction.",
               locations = cells_column_labels(pd)) %>%
  opt_stylize(style = 1)%>%
  gtsave("./tables/juvenilesurvival_updated.docx")




###############################################################################
### END OF SCRIPT ###
###############################################################################
