###################################################################################
##  Script: Effect of parental age on Offspring hatching success
##  Note:   Analysis restricted to a small subset of offspring were mated at adulthood
##################################################################################

#===========================================================================
# 1. SETUP------------------------------------------------------------------
# ===========================================================================

## Libraries ------------------------------------------------------------------
library(dplyr) #v.1.1.4
library(tidyr) #v.1.3.1
library(ggplot2) #v.4.0.0
library(tidyverse) #v.2.0.0
library(brms) #v.2.21.0
library(tidybayes) #v.3.0.7
library(lme4) #v.1.1.35.2
library(bayesplot)  #v.1.11.1
library(bayestestR) #v.0.17.0
library(car) #v.3.1.2
library(loo)  #v.2.8.0
library(emmeans) #v.1.10.1
library(patchwork) #v.1.3.0
library(ggpubr) #v.0.6.1
library(parameters) #v.28.2
library(posterior) #v.1.6.0
library(ggdist) #v.3.3.3
library(sjPlot) #v.2.8.16
library(ggh4x) #v.0.3.1
library(VGAM) #v.1.1.13
library(gt) #v.10.1
library(gtsummary) #v.1.7.2

## Stan backend ---------------------------------------------------------------
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

# ===========================================================================
# 2. DATA ORGANISATION--------------------------------------------------------
# ===========================================================================
eggdata<- readRDS('./raw data/F1_Fecundity_16092025.RDS')

#Creating a proportion that hatched variable
eggdata<-eggdata %>% 
  filter(!Totaleggcount == 0) %>% 
  mutate(prophatched=Totalhatchcount/Totaleggcount) %>% 
  filter(!(F1_ID == "LES-2-1489" & mating_date == as.Date("2024-06-18"))) #removing this individual as it had unreliable hatching count data

#sample sizes for this analyses
length(unique(eggdata$F1_ID))#82 females
length(eggdata$Totalhatchcount)#82 batches of eggs (1 per F1 female)
length(unique(eggdata$PairID)) #from 45 parent pairs (F0 animals)

# ===========================================================================
# 3. SUMMARY STATISTICS----------------------------------------------------------
# ===========================================================================

#Visualising and summarising the data ------------------
hatch_count_graph <- ggplot(eggdata, aes(x = F0_timepoint_binned, y = prophatched, group = F0_timepoint_binned, fill = Temp)) +
  geom_jitter(width = 0.1, col = 'grey70') +
  geom_boxplot(alpha = 0.3) + 
  xlab("Parents' age at reproduction (weeks)") +
  ylab('Number of eggs laid') +
  facet_wrap(~Temp,
             labeller = as_labeller(c("25.5"='25.5°C', "28"='28.0°C', "30.5"='30.5°C'))) +
  scale_fill_manual(values = c('#2f4b7c', 'orange2', "#8B475D"))+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 35),
        axis.text = element_text(size = 35),
        strip.text.x = element_text(size = 25, face = "bold"))
#should be noted that hatching success is generally low overall...

#Distribution of F1 hatching success
hist_hatch_count <- ggplot(eggdata, aes(x =  Totalhatchcount)) +
  geom_histogram() +
  theme_classic() +
  ylab('Count') +
  xlab('Hatchling Count') #Data consists mostly of zeroes --> zero-inflated beta-binomial?

#What is the proportion of zeroes
mean(eggdata$Totalhatchcount == 0, na.rm = TRUE) #Around 34% of the data consist of zero hatchling counts

## Overall mean and SD (used to build intercept priors later)-----------------------
sum_dat1<-eggdata %>% 
  summarise(mean_prophatched = mean (prophatched),
            sd_total= sd(prophatched),
            n_eggs=n())
#Mean proportion that hatched = 0.163, SD = 0.192, n = 82 valid egg batches

#Raw mean effect of parental age on hatching success
sum_dat2<-eggdata %>% 
  group_by(F0_timepoint_binned) %>% 
  summarise(mean_prophatched = mean (prophatched),
            sd_total= sd(prophatched),
            n = n(),
            n_parents=n_distinct(PairID))

#F0_timepoint_binned mean_prophatched sd_total     n n_females
#early                        0.116    0.160    47        32
#late                         0.226    0.218    35        26


# ===========================================================================
# 4. MODEL STRUCTURE SELECTION-----------------------------------------------
# ===========================================================================
# Model is trial data. Compared binomial to beta-binomial families parameterised in brms. 
#Also tested a series of models that varied in their random-effects structure 
#In these models, using default priors. Prior specification occurs later

#4.1. binomial model------------------------------------------
default_binomial<-brm(Totalhatchcount|trials(Totaleggcount)~
                        F0_timepoint_binned+
                        Temp+
                        (1|PairID),
                      family=binomial(link="logit"),
                      data=eggdata,
                      iter=5000,
                      control=list(adapt_delta=0.95, max_treedepth = 15),
                      save_pars = save_pars(all = TRUE),
                      core=4)
summary(default_binomial)
saveRDS(default_binomial, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/default_binomial.rda")
default_binomial<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 hatch/default_binomial.rda")

#MODEL DIAGNOSTICS
loo_binom<-loo(default_binomial,  save_psis = TRUE)
plot(loo_binom) #56 problematic Pareto K observations --> these are incredibly high, not a good model

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(default_binomial, type="loo_pit", ndraws=100) #Model is very overdispersed
pp_check(default_binomial, type ="loo_pit_overlay", ndraws=100)


#4.2. beta binomial---------------------------------------------
default_betabinomial<-brm(Totalhatchcount|trials(Totaleggcount)~
                            F0_timepoint_binned+
                            Temp+
                            max_attempt_round+
                            (1|PairID),
                          family=beta_binomial(link="logit"),
                          data=eggdata,
                          iter=5000,
                          control=list(adapt_delta=0.95),
                          save_pars = save_pars(all = TRUE),
                          core=4)
summary(default_betabinomial)
saveRDS(default_betabinomial, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/default_betabinomial.rda")
default_betabinomial<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 hatch/default_betabinomial.rda")

#MODEL DIAGNOSTICS
loo_betabinom<-loo(default_betabinomial,  save_psis = TRUE)
plot(loo_betabinom) #1 problematic K observations, these are just above 0.7

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(default_betabinomial, type="loo_pit_qq", ndraws=100) #Model is underdispersed in the mid-range --> model expects more extreme values
pp_check(default_betabinomial, type ="loo_pit_overlay", ndraws=100) #Not a fantastic fit

#4.3. Zero inflated beta binomial-----------------------------------------------

#-------------------------Zero-inflated beta-binomial---------------------------
default_zeroinfbetabinomial<-brm(Totalhatchcount|trials(Totaleggcount)~
                                   F0_timepoint_binned+
                                   Temp+
                                   (1|PairID),
                                 family=zero_inflated_beta_binomial(link="logit"),
                                 data=eggdata,
                                 iter=5000,
                                 control=list(adapt_delta=0.98),
                                 save_pars = save_pars(all = TRUE),
                                 core=4)
summary(default_zeroinfbetabinomial)
saveRDS(default_zeroinfbetabinomial, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/default_zeroinfbetabinomial.rda")
default_zeroinfbetabinomial<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 hatch/default_zeroinfbetabinomial.rda")

#MODEL DIAGNOSTICS
loo_zeroinfbetabinom<-loo(default_zeroinfbetabinomial,  save_psis = TRUE)
plot(loo_zeroinfbetabinom) #2 problematic K observations, these are just above 0.7

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(default_zeroinfbetabinomial, type="loo_pit_qq", ndraws=100) #Model is underdispersed in the mid-range --> model expects more extreme values
pp_check(default_zeroinfbetabinomial, type ="loo_pit_overlay", ndraws=100) #Not a fantastic fit

#-------------------COMPARING MODEL FIT-----------------------------------------
#Binomial model is very obviously overdispersed, with very unreliable elpd estimates
#Comparison of model families using K-fold CV

kfold_binom<-kfold(default_binomial, k =10, cores=4)
kfold_betabinom<-kfold(default_betabinomial, k =10, cores=4)
kfold_zerobetabinom<-kfold(default_zeroinfbetabinomial, k =10, cores=4)

#modelsetup
cv_list <- list(
  binomial  = kfold_binom,
  betabinomial = kfold_betabinom,
  zeroinflatedbetabinomial = kfold_zerobetabinom
)
saveRDS(cv_list, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/loo_modelsetup.rda")
cv_list<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 hatch/loo_modelsetup.rda")
loo_modelfit<-loo_compare(cv_list)

#4.4. MODELLING DISTRIBUTIONAL PARAMETERS---------------------------------------

#4.5. Effects of parental age+temp on both zi and phi-------------
mod1.1_dpar_full<-brm(bf(Totalhatchcount|trials(Totaleggcount)~
                           F0_timepoint_binned+
                           Temp+
                           (1|PairID), 
                         zi~F0_timepoint_binned+Temp, #effect of age and temperature on zero-inflation
                         phi~F0_timepoint_binned+Temp), #effect of age and temperature on precision
                      family=zero_inflated_beta_binomial,
                      data=eggdata,
                      #With default priors added
                      iter=5000,
                      control=list(adapt_delta=0.999),
                      save_pars = save_pars(all = TRUE),
                      core=4)

summary(mod1.1_dpar_full)
saveRDS(mod1.1_dpar_full, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/mod1.1_dpar_full.rda")
mod1.1_dpar_full<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 hatch/mod1.1_dpar_full.rda")

#MODEL DIAGNOSTICS
loo_dpar_full<-loo(mod1.1_dpar_full,  save_psis = TRUE)
plot(loo_dpar_full) 

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(mod1.1_dpar_full, type="loo_pit_qq", ndraws=100) #Model with fixed effects on dpars performs better
pp_check(mod1.1_dpar_full, type ="loo_pit_overlay", ndraws=100)



#4.6. Effects of JUST parental age on both zi and phi-------------
mod1.1_dpar_age<-brm(bf(Totalhatchcount|trials(Totaleggcount)~
                          F0_timepoint_binned+
                          Temp+
                          (1|PairID), 
                        zi~F0_timepoint_binned,
                        phi~F0_timepoint_binned),
                     family=zero_inflated_beta_binomial,
                     data=eggdata,
                     #With default priors added
                     iter=5000,
                     control=list(adapt_delta=0.99),
                     save_pars = save_pars(all = TRUE),
                     core=4)
summary(mod1.1_dpar_age)
saveRDS(mod1.1_dpar_age, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/mod1.1_dpar_age.rda")
mod1.1_dpar_age<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 hatch/mod1.1_dpar_age.rda")

#MODEL DIAGNOSTICS
loo_dpar_age-loo(mod1.1_dpar_age,  save_psis = TRUE)
plot(mod1.1_dpar_age) 

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(mod1.1_dpar_age, type="loo_pit_qq", ndraws=100) 
pp_check(mod1.1_dpar_age, type ="loo_pit_overlay", ndraws=100)
#From LOO PIT Q-Q PLOT, this model appears to follow the expected uniform distribution the best


#------------Comparison of distributional models--------------------------------
Loo_distributional<-loo(default_zeroinfbetabinomial, mod1.1_dpar_full, mod1.1_dpar_age, moment_match=TRUE)
saveRDS(Loo_distributional, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/Loo_distributional.rda")
Loo_distributional<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 hatch/Loo_distributional.rda")

# ===========================================================================
# 5. PRIOR SPECIFICATION-----------------------------------------------------
# ===========================================================================

#---------------------Setting up the priors-------------------------------------
#mean proportion of eggs that hatched
mean(eggdata$prophatched) #0.1628852

#mu intercept
qlogis(0.1628852) #Intercept is -1.49 in logit space

#zi intercept
mean(eggdata$Totalhatchcount == 0, na.rm = TRUE) #all zeroes is 0.34
qlogis(0.34)#on the logit scale, this is -0.663

#phi intercept --> up until these models, have use data-informed priors
#uncertain on the value it'll take based on the data -->but want to allow it to be very broad
#a normal(log(10),1) prior can be used here
log(10) #centered on 2.303 --> weakly informative normal prior on log scale
log(5)

#Can't back-calculate the SD onto the logit scale as it is a non-linear transformation
#If you permit B parameters to have effects >+/-1.5 on the plogis scale, 
#it starts to place probability mass on the lower and upper bounds
#Using a general normal 0,1 prior here as the response is bounded by 0 and 1


#------------------------SETTING THE PRIORS-------------------------------------
#need diffuse priors, weakly regularising priors, moderate priors

#mu and zi on logit scale, phi on log scale

## Diffuse priors -------------------------------------------------------------
#Flat priors that are very similar to the brms defaults
diffusepriors <- c(
  prior(normal(-1.636074, 1.5), class = "Intercept"), 
  prior(normal(log(10), 1.5), dpar = "phi", class = "Intercept"),
  prior(normal(-0.663, 1.5), dpar = "zi", class = "Intercept"),
  prior(normal(0,1.5), class = "b"),                
  prior(normal(0,1.5), dpar="phi", class = "b"),
  prior(normal(0,1.5), dpar="zi", class = "b"),
  prior(student_t(3, 0, 2.5), lb =0, class = "sd")             
)

#Weak priors [selected]------------------------------------------------------- 
weakpriors <- c( 
  prior(normal(-1.636074, 1), class = "Intercept"), 
  prior(normal(log(10), 1), dpar = "phi", class = "Intercept"),
  prior(normal(-0.663, 1), dpar = "zi", class = "Intercept"),
  prior(normal(0,1), class = "b"), 
  prior(normal(0,1), dpar="phi", class = "b"),
  prior(normal(0,1), dpar="zi", class = "b"),
  prior(exponential(1), lb =0, class = "sd")           
)

#Moderate priors --------------------------------------------------------------
moderatepriors <- c(
  prior(normal(-1.636074, 0.5), class = "Intercept"), 
  prior(normal(log(10), 0.5), dpar = "phi", class = "Intercept"),
  prior(normal(-0.663, 0.5),dpar = "zi", class = "Intercept"),
  prior(normal(0,0.5), class = "b"), 
  prior(normal(0,0.5), dpar="phi", class = "b"),
  prior(normal(0,0.5), dpar="zi", class = "b"),
  prior(exponential(2), lb =0, class = "sd")            
)


# ===========================================================================
# 6. PRIOR PREDICTIVE CHECKS-------------------------------------------------
# ===========================================================================

#6.1. Diffuse priors-------------------------------------
diffuse_prior_model <- brm(bf(Totalhatchcount|trials(Totaleggcount)~
                                F0_timepoint_binned+
                                Temp+
                                (1|PairID),
                              phi~F0_timepoint_binned,
                              zi~F0_timepoint_binned),
                           family=zero_inflated_beta_binomial(link="logit"),
                           data=eggdata,
                           #With specified priors added
                           prior = diffusepriors,
                           sample_prior = "only",
                           control=list(adapt_delta=0.98),
                           save_pars = save_pars(all = TRUE),
                           iter=5000,
                           cores = 4)
summary(diffuse_prior_model)
saveRDS(diffuse_prior_model, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/diffuse_prior_model.rda")
diffuse_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/diffuse_prior_model.rda")


#setting the colour scheme for the plots
color_scheme_set("teal")

#Prior draws
diffuse_prior_cumvdis<-pp_check(diffuse_prior_model, 
                                ndraws=100, 
                                type = "ecdf_overlay")+xlim(0,800)+theme_classic()
#Diffuse prior places too much mass on zero hatching success

#6.2. Weak priors-------------------------------------------------------------------
weak_prior_model <- brm(bf(Totalhatchcount|trials(Totaleggcount)~
                             F0_timepoint_binned+
                             Temp+
                             (1|PairID),
                           phi~F0_timepoint_binned,
                           zi~F0_timepoint_binned),
                        family=zero_inflated_beta_binomial(link="logit"),
                        data=eggdata,
                        #With specified priors added
                        prior = weakpriors,
                        sample_prior = "only",
                        control=list(adapt_delta=0.98),
                        save_pars = save_pars(all = TRUE),
                        iter=5000,
                        cores = 4)
summary(weak_prior_model)
saveRDS(weak_prior_model, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/weak_prior_model.rda")
weak_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/weak_prior_model.rda")

#Prior draws
weak_prior_cumvdis<-pp_check(weak_prior_model, 
                             ndraws=100, 
                             type = "ecdf_overlay")+xlim(0,800)+theme_classic()
#Pulls estimates slightly away from zero, better than the diffuse prior

#6.3. Moderate priors------------------------------------------------------------------------
moderate_prior_model <- brm(bf(Totalhatchcount|trials(Totaleggcount)~
                                 F0_timepoint_binned+
                                 Temp+
                                 (1|PairID),
                               phi~F0_timepoint_binned,
                               zi~F0_timepoint_binned),
                            family=zero_inflated_beta_binomial(link="logit"),
                            data=eggdata,
                            #With specified priors added
                            prior = moderatepriors,
                            sample_prior = "only",
                            control=list(adapt_delta=0.98),
                            save_pars = save_pars(all = TRUE),
                            iter=5000,
                            cores = 4)
summary(moderate_prior_model)
saveRDS(moderate_prior_model  , file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/moderate_prior_model.rda")
moderate_prior_model <-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/moderate_prior_model.rda")

#Prior draws
moderate_prior_cumvdis<-pp_check(moderate_prior_model, 
                                 ndraws=100, 
                                 type = "ecdf_overlay")+xlim(0,800)+theme_classic()
#Follows the data too tightly...

# ===========================================================================
# 7. POSTERIOR MODELS---------------------------------------------------------
# ===========================================================================

#7.1. Diffuse priors-------------------------------------
mod1.1_diffusepriors <- brm(bf(Totalhatchcount|trials(Totaleggcount)~
                                 F0_timepoint_binned+
                                 Temp+
                                 (1|PairID),
                               phi~F0_timepoint_binned,
                               zi~F0_timepoint_binned),
                            family=zero_inflated_beta_binomial(link="logit"),
                            data=eggdata,
                            #With specified priors added
                            prior = diffusepriors,
                            control=list(adapt_delta=0.98),
                            save_pars = save_pars(all = TRUE),
                            iter=5000,
                            cores = 4)
summary(mod1.1_diffusepriors)
saveRDS(mod1.1_diffusepriors, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/mod1.1_diffusepriors.rda")
mod1.1_diffusepriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/mod1.1_diffusepriors.rda")

#setting the colour scheme for the plots
color_scheme_set("pink")

#Posterior draws
mod1.1_diffusepriors_cumvdis<-pp_check(mod1.1_diffusepriors, 
                                       ndraws=100, 
                                       type = "ecdf_overlay")+xlim(0,800)+theme_classic()


#7.2. Weak priors [SELECTED]-------------------------------------
mod1.1_weakprior <- brm(bf(Totalhatchcount|trials(Totaleggcount)~
                             F0_timepoint_binned+
                             Temp+
                             (1|PairID),
                           phi~F0_timepoint_binned,
                           zi~F0_timepoint_binned),
                        family=zero_inflated_beta_binomial(link="logit"),
                        data=eggdata,
                        #With specified priors added
                        prior = weakpriors,
                        control=list(adapt_delta=0.99),
                        save_pars = save_pars(all = TRUE),
                        iter=5000,
                        cores = 4)
summary(mod1.1_weakprior) #Bulk of parental age effect sits above zero
plot(mod1.1_weakprior)
saveRDS(mod1.1_weakprior, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/mod1.1_weakprior.rda")
mod1.1_weakprior<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/mod1.1_weakprior.rda")

#Posterior draws
mod1.1_weakpriors_cumvdis<-pp_check(mod1.1_weakprior,
                                    ndraws=100, 
                                    type = "ecdf_overlay")+xlim(0,800)+theme_classic()

#7.3. Moderate priors--------------------------------------------
mod1.1_moderateprior <- brm(bf(Totalhatchcount|trials(Totaleggcount)~
                                 F0_timepoint_binned+
                                 Temp+
                                 (1|PairID),
                               phi~F0_timepoint_binned,
                               zi~F0_timepoint_binned),
                            family=zero_inflated_beta_binomial(link="logit"),
                            data=eggdata,
                            #With specified priors added
                            prior = moderatepriors,
                            control=list(adapt_delta=0.95),
                            save_pars = save_pars(all = TRUE),
                            iter=5000,
                            cores = 4)
summary(mod1.1_moderateprior)
saveRDS(mod1.1_moderateprior  , file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/mod1.1_moderateprior.rda")
mod1.1_moderateprior <-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/mod1.1_moderateprior.rda")

#Posterior draws
mod1.1_moderateprior_cumvdis<-pp_check(mod1.1_moderateprior, 
                                       ndraws=100, 
                                       type = "ecdf_overlay")+xlim(0,800)+theme_classic()


##---------------------Prior sensitivity analysis---------------------------------------
Loo_prior_performance<-LOO(mod1.1_diffusepriors,
                           mod1.1_weakprior, 
                           mod1.1_moderateprior, 
                           moment_match = TRUE)
saveRDS(Loo_prior_performance, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/priorfitsummary.rda")


# ===========================================================================
# 8. SELECTED MODEL DIAGNOSTICS  (mod1.1_weakprior)--------------------------
# ===========================================================================

## LOO-CV --------------------------------------------------------------------
loo_weakpriors <- loo(mod1.1_weakprior, save_psis = TRUE)
plot(loo_weakpriors)
psis_weakpriors <- loo_weakpriors$psis_object
psis_weights    <- weights(psis_weakpriors)


## Posterior predictive checks -----------------------------------------------
#Posterior simiulated and empirical (observed) means
weakprior_mean     <- ppc_stat(yrep = posterior_predict(mod1.1_weakprior),
                               y    = eggdata$Totalhatchcount,
                               stat = "mean") + theme_classic() #Mean is slightly higher than estimated

#posterior cumulative density
weakprior_ecdf     <- pp_check(mod1.1_weakprior, ndraws = 100,
                               type = "ecdf_overlay") +
  xlim(0, 600) + theme_classic()

#posterior probability density
weakprior_pdf      <- pp_check(mod1.1_weakprior, ndraws = 100) +
  xlim(0, 600) + theme_classic()

#posterior LOO Q-Q plot
weakprior_loo_qq   <- pp_check(mod1.1_weakprior, type = "loo_pit_qq",
                               ndraws = 100) + theme_classic()
#LOO-PIT values fall neatly along the diagonal line (i.e. expected uniform distribution)

#LOO uniformity plot
weakprior_loo_unif <- pp_check(mod1.1_weakprior, type = "loo_pit_overlay",
                               ndraws = 100) + theme_classic()
#predictive interval plot
weakprior_intervals <- ppc_loo_intervals(
  y    = eggdata$Totalhatchcount,
  yrep = posterior_predict(mod1.1_weakprior),
  psis_weakpriors
) + theme_classic()

## Combined fit plot ---------------------------------------------------------
fit_plots <- ggarrange(
  weakprior_mean, weakprior_ecdf, weakprior_pdf,
  weakprior_loo_qq, weakprior_loo_unif, weakprior_intervals,
  nrow = 2, ncol = 3,
  labels = c("A", "B", "C", "D", "E", "F")
)

ggsave("./bayesian_plots/model fit plots/F1 hatch/selectedmodelfit.png",
       plot   = fit_plots,
       device = "png",
       width  = 500, height = 400, units = "mm")

## Bayes R² ------------------------------------------------------------------
bayes_R2(mod1.1_weakprior, re.form = NA)   # fixed effects only (i.e. marginal)
bayes_R2(mod1.1_weakprior, re.form = NULL)  # including random effects (i.e., conditional)

## MAP estimates and pd ------------------------------------------------------
MAP_hatch<- bayestestR::describe_posterior(
  mod1.1_weakprior,
  test        = "pd",
  ci_method   = "HDI",
  centrality  = "MAP",
  component   = "all",
  effects     = "full",
  ci          = 0.95
)
saveRDS(MAP_hatch,
        "scripts/model_outputs/Offspring Trait Models/F1 hatch/MAP_hatch.rda")

## ROPE -----------------------------------------------------------------------
rope_hatch <- bayestestR::describe_posterior(
  mod1.1_weakprior,
  rope_range=c(-0.18, 0.18), #the reccomended default for models with a logitsic function
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1 #ROPE estimates need to use the full posterior distribution, hence why they're calculated seperately
)
saveRDS(rope_hatch,
        "scripts/model_outputs/Offspring Trait Models/F1 hatch/rope_hatch.rda")


# ===========================================================================
# 9. PRIOR vs. POSTERIOR SENSITIVITY PLOTS
# ===========================================================================
priorpostplots<-ggarrange(diffuse_prior_cumvdis, weak_prior_cumvdis, 
                          moderate_prior_cumvdis,
                          mod1.1_diffusepriors_cumvdis, 
                          mod1.1_weakpriors_cumvdis, 
                          mod1.1_moderateprior_cumvdis, nrow = 2, ncol = 3,
                          labels = c("A", "B", "c", "D", "E", "F"))

ggsave(filename = "./bayesian_plots/model fit plots/F1 hatch/priorvpost_ecdf.png",
       plot = priorpostplots, 
       device = "png", 
       width = 420, 
       height = 250, 
       units = "mm")

## ===========================================================================
# 10. HYPOTHESIS TESTING: INTERACTION MODELS
#Models test the core interactions of the paper: are parental age effects environmentally dependent?
# ===========================================================================

#10.1. Is there a temperature*parental age interaction?-------
mod2.1<-brm(bf(Totalhatchcount|trials(Totaleggcount)~
                 F0_timepoint_binned+
                 Temp+
                 F0_timepoint_binned:Temp+
                 (1|PairID),
               phi~F0_timepoint_binned,
               zi~F0_timepoint_binned),
            data=eggdata,
            family = zero_inflated_beta_binomial(link="logit"),
            #With specified priors added
            prior = weakpriors,
            iter=5000,
            control=list(adapt_delta=0.98),
            save_pars = save_pars(all = TRUE),
            core=4)
summary(mod2.1) #No evidence for an effect of the interaction
saveRDS(mod2.1, file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/mod2.1.rda")
mod2.1<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 hatch/mod2.1.rda")

## MAP estimates for interaction model (mod3.1) ------------------------------
#For reporting in the main manuscript
MAP_hatch_mod2.1 <- bayestestR::describe_posterior(
  mod2.1,
  test       = "pd",
  ci_method  = "HDI",
  centrality = "MAP",
  effects    = "full",
  ci         = 0.95
)
saveRDS(MAP_hatch_mod2.1,
        "scripts/model_outputs/Offspring Trait Models/F1 hatch/MAP_hatch_mod2.1.rda")

#ROPE estimates for the interactive effect (using the full posterior distribution)
rope_hatch_mod2.1 <- bayestestR::describe_posterior(
  mod2.1,
  rope_range= c(-0.18, 0.18),
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1
)
saveRDS(rope_hatch_mod2.1,
        "scripts/model_outputs/Offspring Trait Models/F1 hatch/rope_hatch_mod2.1.rda")

#No evidence for an interaction between temperature and parental age

#------------------------Model comparison between model with single-effects and interaction term----------------------------------
hypothesis_loo<-loo(mod1.1_weakprior, 
                    mod2.1, 
                    moment_match=TRUE)
saveRDS(hypothesis_loo , file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/hypothesis_loo.rda")
hypothesis_loo<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 fecundity/hypothesis_loo.rda")
#No difference in the predictive abilities between model 1 and model 2

#---------------------ESTIMATED MARGINAL MEANS---------------------------------
pairwise_estimates<- emmeans(mod1.1_weakprior, "F0_timepoint_binned", type="response")
pairs(pairwise_estimates)


# ===========================================================================
# 11. POSTERIOR DISTRIBUTION PLOTS  (mod1.1_weakprior)-----------------------
# ===========================================================================
# Plots styled using stat_halfeye (ggdist) for consistency with other traits

post1 <- as_draws_df(mod1.1_weakprior) #for single effect predictors
post2<-as_draws_df(mod2.1) #for interaction term

## Build combined posterior data frame ----------------------------------------
posterior_df_hatch <- data.frame(
  "µ: Parents' age at reproduction (Late-Aged)" = post1$b_F0_timepoint_binnedlate,
  "µ: Parents' temperature treatment (28.0°C)" = post1$b_Temp28,
  "µ: Parents' temperature treatment (30.5°C)" = post1$b_Temp30.5,
  "ϕ: Parents' age at reproduction (Late-Aged)" = post1$b_phi_F0_timepoint_binnedlate,
  "z: Parents' age at reproduction (Late-Aged)" = post1$b_zi_F0_timepoint_binnedlate,
  "μ: Parents' age × Temperature (28.0°C)"    = post2$`b_F0_timepoint_binnedlate:Temp28`,
  "μ: Parents' age × Temperature (30.5°C)"    = post2$`b_F0_timepoint_binnedlate:Temp30.5`,
  check.names = FALSE
) %>%
  pivot_longer(everything(),
               names_to  = "parameter",
               values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = c(
    "µ: Parents' age at reproduction (Late-Aged)",
    "µ: Parents' temperature treatment (28.0°C)",
    "µ: Parents' temperature treatment (30.5°C)",
    "μ: Parents' age × Temperature (28.0°C)",
    "μ: Parents' age × Temperature (30.5°C)" ,
    "ϕ: Parents' age at reproduction (Late-Aged)",
    "z: Parents' age at reproduction (Late-Aged)"
  )))

y_levels <- rev(c(
  "µ: Parents' age at reproduction (Late-Aged)",
  "µ: Parents' temperature treatment (28.0°C)",
  "µ: Parents' temperature treatment (30.5°C)",
  "μ: Parents' age × Temperature (28.0°C)",
  "μ: Parents' age × Temperature (30.5°C)" ,
  "ϕ: Parents' age at reproduction (Late-Aged)",
  "z: Parents' age at reproduction (Late-Aged)"
))

## Halfeye posterior plot ------------------------------------------------------
posterior_plot_hatch<- ggplot(
  posterior_df_hatch,
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
    "#e5a9f5", "#4c9279", "#d7c9a3"
  )) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title      = element_text(size = 30),
    axis.text       = element_text(size = 30)
  ) +
  labs(x = "Fixed-effect posterior estimates\n(on logit scale: μ = mean log-odds of hatching, zi = structural zero inflation, \non log scale: ϕ = precision)",
       y = NULL)

ggsave("./bayesian_plots/posterior plots/F1 hatch/001_hypothesis_halfeye.png",
       plot   = posterior_plot_hatch,
       device = "png",
       width  = 580, height = 400, units = "mm")

## Post-hoc standardised mean difference (SMD) --------------------------------

#Post-hoc standardised mean difference (on unit scale)
pairwise_estimates<- emmeans(mod1.1_weakprior, "F0_timepoint_binned", type = "response")
#This is now converted to the raw scale where the extreme ends of parental age are compared
sd_life<-sd(eggdata$prophatched)

#Then standardising by SD
hatching_SMD <- summary(pairwise_estimates) %>% 
  summarise(
    SMD = (prob[2] - prob[1])/sd_life,
    SMD_lower = (lower.HPD[2] - upper.HPD[1])/sd_life,
    SMD_upper = (upper.HPD[2] - lower.HPD[1])/sd_life
  ) 

#Saving for use in combined forest plot of all traits
saveRDS(hatching_SMD, file = "scripts/model_outputs/Offspring Trait Models/hatching_SMD")


# ===========================================================================
# 12. MODEL PREDICTIONS AND RESULTS PLOTS-------------------------------------
# ===========================================================================

#12.1: Prediction grid
df_predict_within<-expand.grid(
  Totaleggcount = unique(as.numeric(eggdata$Totaleggcount)),
  F0_timepoint_binned=unique(eggdata$F0_timepoint_binned),
  Temp=unique(eggdata$Temp),
  PairID=unique(eggdata$PairID)[1])


#calculating model predictions and extracting standard errors
pred_within <- fitted(mod1.1_weakprior, df_predict_within, re_formula=NA)
pred_within<-as.data.frame(pred_within)
df_predict_within$prediction <- pred_within$Estimate
df_predict_within$lower <- pred_within$Q2.5
df_predict_within$upper<-pred_within$Q97.5

# Calculate mean predictions and SE for plotting
newdat <- df_predict_within %>%
  group_by(F0_timepoint_binned) %>%
  summarise(
    mean  = sum(prediction) / sum(Totaleggcount),
    lower = sum(lower) / sum(Totaleggcount),
    upper = sum(upper) / sum(Totaleggcount)
  )  %>% 
  ungroup()



#12.2. Probability of hatching success plot-----
hatchingsuccess_plot <- ggplot(data = eggdata,
                               aes(x = F0_timepoint_binned, 
                                   y = prophatched,
                                   fill =F0_timepoint_binned,
                                   colour=F0_timepoint_binned)) +
  geom_point(position = position_jitter(width = 0.2, height=0.01),
             shape=21,
             size=8,
             stroke=1.8,
             alpha = 0.7,
             colour="white",
             fill= "grey")+
  geom_linerange(data =newdat,
                 aes(y=mean,
                     ymin = lower,
                     ymax = upper),
                 linewidth = 3,
                 show.legend = FALSE)+
  geom_point(data = newdat,
             aes(x = F0_timepoint_binned,
                 y = mean),
             shape=21, 
             stroke=1.8,
             size= 25,
             alpha=1,
             colour="white",
             show.legend = FALSE)+
  theme_classic() + 
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  #axis.title.x = element_blank(),  
  #axis.text.x = element_blank(),   
  #axis.ticks.x = element_blank())+
  theme(#legend.position = "none",
    legend.title=element_text(size=50),
    legend.text=element_text(size=50),
    strip.text = element_text(size=50),
    panel.background = element_rect(fill = "transparent", color = NA),  # Make panel background transparent
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA), # Make plot background transparent
  )  +
  scale_color_manual(name = "Parents' age at reproduction", labels = c("Early-Aged", "Late-Aged"), values = c("#f96161", "#066594"))+
  scale_fill_manual(name = "Parents' age at reproduction",labels = c("Early-Aged", "Late-Aged"), values = c("#f96161", "#066594"))+
  scale_x_discrete(labels=c("Early-Aged", "Late-Aged"))+
  labs(x = "Parents' adult age at reproduction", y = "Proportion of offspring eggs that hatched")


# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 hatching success/001_F1_prophatched2.png",
       plot = hatchingsuccess_plot,
       bg="transparent",
       device = "png", 
       width = 355, 
       height = 650, 
       units = "mm")


#12.3. Number of eggs that hatched plot-----------------------------------------------
newdat_number <- df_predict_within %>%
  group_by(F0_timepoint_binned) %>%
  summarise(
    mean  = mean(prediction),
    lower = mean(lower),
    upper = mean(upper)
  )  %>% 
  ungroup()


# plotting the predicted values onto a ggplot######
hatchingsuccess_number <- ggplot(data = eggdata,
                                 aes(x = F0_timepoint_binned, 
                                     y = Totalhatchcount,
                                     fill = F0_timepoint_binned,
                                     colour=F0_timepoint_binned)) +
  geom_point(position = position_jitter(width = 0.2, height=0.1),
             shape=21,
             size=8,
             stroke=1.8,
             alpha = 0.7,
             colour="white",
             fill= "grey")+
  geom_linerange(data =newdat_number,
                 aes(y=mean,
                     ymin = lower,
                     ymax = upper),
                 linewidth = 3,
                 show.legend = FALSE)+
  geom_point(data = newdat_number,
             aes(x = F0_timepoint_binned,
                 y = mean),
             shape=21, 
             stroke=1.8,
             size= 25,
             alpha=1,
             colour="white",
             show.legend = FALSE)+
  theme_classic() + 
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(#legend.position = "none",
    legend.title=element_text(size=50),
    legend.text=element_text(size=50),
    strip.text = element_text(size=50),
    panel.background = element_rect(fill = "transparent", color = NA),  # Make panel background transparent
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA), # Make plot background transparent
  )+
  scale_color_manual(name = "Parents' age at reproduction", labels = c("Early-Aged", "Late-Aged"), values = c("#f96161", "#066594"))+
  scale_fill_manual(name = "Parents' age at reproduction",labels = c("Early-Aged", "Late-Aged"), values = c("#f96161", "#066594"))+
  scale_x_discrete(labels=c("Early-Aged", "Late-Aged"))+
  labs(x = "Parents' age at reproduction", y = "Number of offspring eggs that hatched")

#Number of eggs that hatched
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 hatching success/001_F1_numberhatched2.png",
       plot = hatchingsuccess_number,
       bg="transparent",
       device = "png", 
       width = 355, 
       height = 650, 
       units = "mm")



#12.4. Probability mass function for number of hatchlings-------------

#------Holding the number of eggs laid constant (hatching success estimated at median fecundity)

#Proxy data
newdat <- expand.grid(
  Totaleggcount = unique(eggdata$Totaleggcount),
  F0_timepoint_binned = unique(eggdata$F0_timepoint_binned),
  Temp= unique(eggdata$Temp),
  PairID=unique(eggdata$PairID[1]))

#posterior draws
proportional_posterior <- mod1.1_weakprior %>% 
  linpred_draws(
    newdat,
    transform =TRUE,
    value = "mu",
    allow_new_levels = TRUE, 
    re_formula = NA,  # population-level predictions
    dpar = c("mu","phi", "zi"),
    ndraws = 3000, 
    seed = 123
  )%>%
  mutate(alpha = mu * phi,
         beta  = (1 - mu) * phi)

#Thinning the posterior draw number
posterior_df <- proportional_posterior %>% slice_sample(n = 1000)


# Population median clutch size
n_med <- round(median(eggdata$Totaleggcount, na.rm = TRUE))

# Compute posterior PMF for number of eggs that hatched
summary_pmf <- map_dfr(0:n_med, function(x_val) {
  posterior_df %>%
    mutate(
      f = dzoibetabinom.ab(
        x      = as.integer(x_val),
        size   = as.integer(n_med),
        shape1 = alpha,
        shape2 = beta,
        pstr0  = zi
      )
    ) %>%
    group_by(F0_timepoint_binned) %>%
    summarise(
      x      = x_val,      # number of eggs hatched
      prop = x_val/n_med,
      median = median(f, na.rm = TRUE),
      lower  = quantile(f, 0.025, na.rm = TRUE),
      upper  = quantile(f, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
})


#Getting the median development times
median_hatching<- mod1.1_weakprior %>% 
  epred_draws(
    newdat,
    value = "mu",
    allow_new_levels=TRUE,
    transform = TRUE,
    re_formula=NA,#population level predictions
    dpar = c("phi", "zi"),
    ndraws = 5000, 
    seed = 123
  ) %>% 
  ungroup() %>%
  group_by(.draw, F0_timepoint_binned) %>% 
  summarise(median_hatching= median(mu)) %>% 
  group_by(F0_timepoint_binned) %>% 
  median_hdi(median_hatching)
#This is the posterior median across all eggs...


#--------------------------Probability density for hatching success------------------
summary_probability_binned <- summary_pmf %>%
  mutate(x_bin = floor(x/10)*10) %>%  # bin width = 10 hatchlings
  group_by(x_bin, F0_timepoint_binned) %>%
  summarise(median = sum(median),
            lower = sum(lower),
            upper = sum(upper), .groups = "drop")


#Probability density function---------------------------------------------------
f1hatching_pdf<-ggplot(data=summary_probability_binned,
                       aes(x = x_bin, 
                           y = median, 
                           colour= F0_timepoint_binned)) +
  geom_col(position = position_dodge(width = 0.95),    # discrete bars
           fill="grey75",
           color = "grey75", width = 2.2, alpha = 1) +
  # Transparent ribbons for 95% credible intervals
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 10), 
                width = 3, linewidth = 6) +
  geom_segment(data = median_hatching,
               aes(x = median_hatching, xend = median_hatching, y = 0, yend = 0.85, color = F0_timepoint_binned),
               linetype = "dashed", linewidth = 4, alpha = 0.7)+
  theme_classic()+
  scale_y_continuous(limits=c(0, 0.89))+
  coord_cartesian(xlim = c(0, 178))+
  theme(axis.title=element_text(size=45),
        axis.text=element_text(size=45))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=45),
        strip.text = element_text(size=45, face="bold"),
        panel.background = element_rect(fill = "transparent", color = NA),  # Make panel background transparent
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), # Make plot background transparent
  ) +
  scale_color_manual(name = "Parents' age at reproduction", labels = c("Early-Aged", "Late-Aged"), values = c("#f96161", "#066594"))+
  scale_alpha_manual(name = "Parents' age at reproduction",labels = c("Early-Aged", "Late-Aged"), values = c(0.4,0.4))+
  scale_fill_manual(name = "Parents' age at reproduction",labels = c("Early-Aged", "Late-Aged"), values = c("#f96161", "#066594"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  theme(
    legend.position = "top", 
    legend.box = "horizontal"
  ) +
  labs(x = "Number of offspring eggs that hatched", y = "Probability mass, f(x)")

# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 hatching success/003_fecundity_pmf.png",
       plot = f1hatching_pdf, 
       bg="transparent",
       device = "png", 
       width = 400, 
       height = 400, 
       units = "mm")

#------------------------Saving model plots into one panel----------------------

#Making sure the plots have the appropraite scale design etc.
new_offspring_num<- hatchingsuccess_number +
  theme(axis.title=element_text(size=45),
        axis.text=element_text(size=45))+
  theme(legend.position = "none")

new_offspring_prop<-hatchingsuccess_plot +
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=45),
        axis.text=element_text(size=45))+
  theme(legend.position = "none")

#Saving the prior and poosterior plots for default versus moderate priors
inference<-ggarrange(
  ggarrange(new_offspring_prop, new_offspring_num, nrow = 2, labels = c("A", "B"),
            font.label = list(size = 50, face = "bold"),
            label.x = 0.16,
            heights = c(0.94, 1),align = "v"),
  f1hatching_pdf,
  ncol = 2,
  nrow=1,
  labels = c("", "C"),
  font.label = list(size = 50, face = "bold"),
  label.x=0.1,
  widths = c(0.85, 1.6))


ggsave(filename = "./bayesian_plots/offspring trait plots/F1 hatching success/inference_plots.png",
       plot = inference, 
       device = "png", 
       width = 850, 
       height = 630, 
       units = "mm")


#Creating a table to export------------------------------------------------------------------------------
# ── helper ────────────────────────────────────────────────────────────────────
#Model draws from selected models
base_model_draws<-as_draws_df(mod1.1_weakprior)
interaction_draws<-as_draws_df(mod2.1) 

#function for generating estimates
summarise_param <- function(draws, param, section,
                            rope = TRUE, #using default ROPE range automatically calculated from bayestest
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
  summarise_param(base_model_draws$b_F0_timepoint_binnedlate,   "Parents' age category: Late-aged", "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws$b_Temp28, "Temperature: 28.0°C",        "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws$b_Temp30.5, "Temperature: 30.5°C",        "Location submodel (Fixed effects)", rope = TRUE)
)


# ── 1. dropped interaction (location) — ROPE + pd ───────────────────────────────────
fe_interaction2 <- bind_rows(summarise_param(interaction_draws2$`b_delta.age:Temp28`, "Parents' Δage x  Temperature: 28.0°C", "Location submodel (Dropped interactions)", rope = TRUE),
                             summarise_param(interaction_draws2$`b_delta.age:Temp30.5`, "Parents' Δage x  Temperature: 30.5°C", "Location submodel (Dropped interactions)", rope = TRUE)
)

# ── 3.Random effects ──────────────────────────────────────────────────
# SDs 
re_sd_explore<- bind_rows(
  summarise_param(base_model_draws$sd_PairID__Intercept, "σ intercept",    "Random effects (Location)", pd = FALSE)
)

#---4. Shape parameter──────────────────────────────────────────────────
precision<- bind_rows(
  summarise_param(base_model_draws$b_phi_Intercept, "Precision intercept (ϕ)", "Distributional parameters (scale component)",  rope = FALSE),
  summarise_param(base_model_draws$b_phi_F0_timepoint_binnedlate, "Parents' age category: Late-aged", "Distributional parameters (scale component)", rope= FALSE)
)

#---5. zero-inflation parameter──────────────────────────────────────────────────
zeroinf<- bind_rows(
  summarise_param(base_model_draws$b_zi_Intercept, "Zero-inflation intercept (z)", "Distributional parameters (zero-inflation component)",  rope = FALSE),
  summarise_param(base_model_draws$b_zi_F0_timepoint_binnedlate, "Parents' age category: Late-aged", "Distributional parameters (zero-inflation component)", rope= FALSE)
)


#---6. Conditional and marginal Bayes R² ------------------------------------------------------------------
marginal<-bayes_R2(mod1.1_weakprior, re.form = NA, summary = FALSE)   # fixed effects only (i.e., marginal)
conditional<-bayes_R2(mod1.1_weakprior, re.form = NULL, summary = FALSE)  # including random effects (i.e., conditional)

bayes<-bind_rows(
  summarise_param(marginal, "Marginal R²", "Bayes R²", pd= FALSE, rope= FALSE),
  summarise_param(conditional, "Conditional R²", "Bayes R²", pd= FALSE, rope =FALSE)
)


# ── 5. Combine and render ─────────────────────────────────────────────────────
bind_rows(fe, fe_interaction2, re_sd_explore, precision, zeroinf, bayes) %>%
  gt(groupname_col = "Section") %>%
  cols_label(Parameter = "Parameter", MAP = "MAP",
             `95% HDI` = "95% HDI", `% in ROPE` = "% in ROPE", pd = "pd") %>%
  tab_header(
    title    = "Offspring hatching success model summary: zero-inflated beta-binomial model"
  ) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) %>%
  tab_footnote("ROPE = [−0.18, 0.18] on log-odds scale; applied to location fixed effects only.",
               locations = cells_column_labels(`% in ROPE`)) %>%
  tab_footnote("pd = probability of direction.",
               locations = cells_column_labels(pd)) %>%
  opt_stylize(style = 1)%>%
  gtsave("./tables/hatchingsuccess_updated.docx")


###############################################################################
##END OF SCRIPT----------------------------------------------------------------
###############################################################################

