###################################################################################
##  Script: Effect of parental age on Offspring fecundity
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
library(gt) #v.10.1
library(gtsummary) #v.1.7.2

## Stan backend ---------------------------------------------------------------
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores())

# ===========================================================================
# 2. DATA ORGANISATION--------------------------------------------------------
# ===========================================================================
data <- readRDS('./raw data/F1_Fecundity_16092025.RDS')

#sample sizes for this analyses
length(unique(data$F1_ID))#83 females
length(data$Totaleggcount)#83 batches of eggs (1 per F1 female)
length(unique(data$PairID)) #from 45 parent pairs (F0 animals)

# ===========================================================================
# 3. SUMMARY STATISTICS----------------------------------------------------------
# ===========================================================================

#Visualising and summarising the data ------------------
egg_count_graph <- ggplot(data, aes(x = F0_timepoint_binned, y = Totaleggcount, group = F0_timepoint_binned, fill = Temp)) +
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


## Distribution of F1 fecundity ---------------------------------------------------
hist_egg_count <- ggplot(data, aes(x =  Totaleggcount)) +
  geom_histogram() +
  theme_classic() +
  ylab('Count') +
  xlab('Egg Count') 

## Overall mean and SD (used to build intercept priors later)-----------------------
sum_dat1<-eggdata %>% 
  summarise(mean_fecundity = mean (Totaleggcount),
            sd_total= sd(Totaleggcount),
            n_eggs=n())
#Mean number of eggs = 288, SD = 224, n = 83 egg batches

#2. Mean effect of parental age on fecundity
sum_dat2<-eggdata %>% 
  group_by(F0_timepoint_binned) %>% 
  summarise(mean_fecundity = mean (Totaleggcount),
            sd_total= sd(Totaleggcount),
            n = n(),
            n_females=n_distinct(PairID))

# F0_timepoint_binned mean_fecundity sd_total  n_F1    n_F0

#early                         237.     197.    47       32
#late                          354.     242.    36       26


# ===========================================================================
# 4. MODEL STRUCTURE SELECTION-----------------------------------------------
# ===========================================================================
# Model is count data. Compared Poisson to negative binomial families. 
#Also tested a series of models that varied in their random-effects structure 
#In these models, using default priors. Prior specification occurs later.

#4.1. Poisson model----------------------------------------------------------
mod1.1_poisson<-brm(Totaleggcount~
                      F0_timepoint_binned+
                      Temp+
                      (1|PairID), 
                    family=poisson,
                    data=data,
                    iter=5000,
                    control=list(adapt_delta=0.98),
                    save_pars = save_pars(all = TRUE),
                    core=4)
summary(mod1.1_poisson)
saveRDS(mod1.1_poisson, file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_poisson.rda")
mod1.1_poisson<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_poisson.rda")

#MODEL DIAGNOSTICS
loo_poisson<-loo(mod1.1_poisson,  save_psis = TRUE, moment_match=TRUE)
plot(mod1.1_poisson) 

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(mod1.1_poisson, type="loo_pit_qq", ndraws=100) #Model is very overdispersed
pp_check(mod1.1_poisson, type ="loo_pit_overlay", ndraws=100) #Not a fantastic fit


#4.2. Negative binomial model----------------------------------------------------------
mod1.1_negbinom<-brm(Totaleggcount~F0_timepoint_binned+
                       Temp+
                       (1|PairID), 
                     family=brms::negbinomial,
                     data=data,
                     iter=5000,
                     control=list(adapt_delta=0.98),
                     save_pars = save_pars(all = TRUE),
                     core=4)
summary(mod1.1_negbinom) 
saveRDS(mod1.1_negbinom, 
        file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_negbinom.rda")
mod1.1_negbinom<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_negbinom.rda")

#MODEL DIAGNOSTICS
loo_negbinom<-loo(mod1.1_negbinom,  save_psis = TRUE, moment_match=TRUE)
plot(loo_negbinom) #No problematic loo diagnostics

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(mod1.1_negbinom, type="loo_pit_qq", ndraws=100) #performs exceptionally well
pp_check(mod1.1_negbinom, type ="loo_pit_overlay", ndraws=100) #a fantastic fit

#Negative binomial looks ideal

#-------------------COMPARING MODEL FIT-----------------------------------------
#Model setup Kfold values
#Poisson model has probelmatic pareto k estimates...using kfold instead here to get ELPD values
kfold_poisson<-kfold(mod1.1_poisson, k =10, cores=1)
kfold_negbinom<-kfold(mod1.1_negbinom, k =10, cores=1)

#modelsetup
cv_list <- list(
  poisson  = kfold_poisson,
  negbinom = kfold_negbinom
)
saveRDS(cv_list, file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/loo_modelsetup.rda")
cv_list<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 fecundity/loo_modelsetup.rda")

#Loo comparison of negative binomial and Poisson models
loo_modelfit<-loo_compare(cv_list)


#Evidence for fixed effects on the distributional parameters--------
#4.3. Both temperature and age on the shape parameter
mod1.1_allonshape<-brm(bf(Totaleggcount~
                            F0_timepoint_binned+
                            Temp+
                            (1|PairID),
                          shape~F0_timepoint_binned+
                            Temp),
                       family=brms::negbinomial,
                       data=data,
                       #With default priors added
                       iter=5000,
                       control=list(adapt_delta=0.99),
                       save_pars = save_pars(all = TRUE),
                       core=4)
summary(mod1.1_allonshape) 
saveRDS(mod1.1_allonshape, file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_allonshape.rda")
mod1.1_allonshape<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_allonshape.rda")

#Model diganostics
loo_shapetemp<-loo(mod1.1_allonshape)
plot(loo_shapetemp)

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(mod1.1_allonshape, type="loo_pit_qq", ndraws=100)
pp_check(mod1.1_allonshape, type ="loo_pit_overlay", ndraws=100) 


#4.4. Testing the effect of just parental age on the shape parameter
mod1.1_shape<-brm(bf(Totaleggcount~
                       F0_timepoint_binned+
                       Temp+
                       (1|PairID),
                     shape~F0_timepoint_binned),
                  family=brms::negbinomial,
                  data=data,
                  iter=5000,
                  control=list(adapt_delta=0.99),
                  save_pars = save_pars(all = TRUE),
                  core=4)
summary(mod1.1_shape) 
saveRDS(mod1.1_shape, 
        file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_shape.rda")
mod1.1_shape<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_shape.rda")

#MODEL DIAGNOSTICS
loo_positive<-loo(mod1.1_shape,  save_psis = TRUE, moment_match=TRUE)
plot(mod1.1_shape) #No problematic loo diagnostics

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(mod1.1_shape, type="loo_pit_qq", ndraws=100) #performs very well
pp_check(mod1.1_shape, type ="loo_pit_overlay", ndraws=100) 

#LOO comparisons of models with distributional parameters-------------------
loo_dpars<-loo(mod1.1_negbinom, 
               mod1.1_allonshape, 
               mod1.1_shape, 
               moment_match=TRUE)
saveRDS(loo_dpars, file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/loo_dpars.rda")
loo_dpars<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 fecundity/loo_dpars.rda")
#No difference in model fit between the three models.
#Selecting the model with just parental age on the shape parameter

# ===========================================================================
# 5. PRIOR SPECIFICATION-----------------------------------------------------
# ===========================================================================

#---------------------Setting up the priors-------------------------------------
var(data$Totaleggcount) #variance of 50304.03
mean(data$Totaleggcount) #mean of 287.739
sd(data$Totaleggcount) #sd of 224.286 eggs

#To get an idea of the ideal data informed shape parameter
#shape = 
287.739^2/(50304.03 -287.739)
#=1.653
#Model uses the log link for the shape parameter
log(1.655)
#Setting prior for shape prior intercept at 0.503 (with an SD of 1)

#mean egg count on log scale
log(287.739) #5.662054

#SD on log scale
sqrt(log(1 + (224.286^2 /287.739^2))) #0.689 on log scale


#------------------------SETTING THE PRIORS-------------------------------------
#need diffuse priors, weakly regularising priors, moderate priors

#specifying the priors for effects
hist(rnorm(1000,0,1)) #Allows large effects for a chnage in the predictor
#Allows up to an 80% increase or 45% decrease in egg counts per SD increase in the predictor

## Diffuse priors -------------------------------------------------------------
#Flat priors that are very similar to the brms defaults
diffusepriors <- c(
  prior(normal(5.66, 0.689), class = "Intercept"),               
  prior(normal(0.503, 2), class = "Intercept", dpar="shape"),
  prior(normal(0,1), class = "b"),    
  prior(normal(0,2), class = "b", dpar = "shape"),
  prior(student_t(3, 0, 2.5), class = "sd", lb = 0)
)


## Weak priors [SELECTED] ------------------------------------------------------
weakpriors <- c(
  prior(normal(5.66, 0.689), class = "Intercept"),  
  prior(normal(0,0.689), class = "b"),     
  prior(normal(0.503, 1), class = "Intercept", dpar="shape"),
  prior(normal(0,1), class = "b", dpar = "shape"),
  prior(exponential(1), class = "sd", lb = 0)
)

## Moderate priors ------------------------------------------------------------
moderatepriors <- c(
  prior(normal(5.66, 0.345), class = "Intercept"), 
  prior(normal(0.503, 1), class = "Intercept", dpar="shape"),
  prior(normal(0,0.345), class = "b"),  
  prior(normal(0,1), class = "b", dpar = "shape"),
  prior(exponential(2), class = "sd", lb = 0)
)

# ===========================================================================
# 6. PRIOR PREDICTIVE CHECKS-------------------------------------------------
# ===========================================================================

#6.1. Diffuse priors-------------------------------------
diffuse_prior_model <- brm(bf(Totaleggcount~
                                F0_timepoint_binned+
                                Temp+
                                (1|PairID),
                              shape~F0_timepoint_binned),
                           data=data,
                           family=brms::negbinomial(),
                           #With specified priors added
                           prior = diffusepriors,
                           sample_prior = "only",
                           control=list(adapt_delta=0.95),
                           iter=5000,
                           cores = 4)
summary(diffuse_prior_model)
saveRDS(diffuse_prior_model, file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/diffuse_prior_model.rda")
diffuse_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/diffuse_prior_model.rda")


#Prior draws
color_scheme_set("teal")
diffuse_prior_cumvdis<-pp_check(diffuse_prior_model, ndraws=100, type = "ecdf_overlay")+
  theme_classic()+
  xlim(0,2000) 

#The diffuse prior overestimates the number of zeroes within the data
#Expects more zeroes than what are actually present 

#6.2. Weak priors----------------------------------------------------
weak_prior_model <- brm(bf(Totaleggcount~
                             F0_timepoint_binned+
                             Temp+
                             (1|PairID),
                           shape~F0_timepoint_binned),
                        data=data,
                        family=brms::negbinomial(),
                        #With specified priors added
                        prior = weakpriors,
                        sample_prior = "only",
                        control=list(adapt_delta=0.95),
                        iter=5000,
                        cores = 4)
summary(weak_prior_model)
saveRDS(weak_prior_model, file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/weak_prior_model.rda")
weak_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/weak_prior_model.rda")

#Prior draws
weak_prior_cumvdis<-pp_check(weak_prior_model, 
                             ndraws=100, 
                             type = "ecdf_overlay")+theme_classic()+xlim(0,2000) 
#These priors are slightlye better looking than the diffuse ones

#6.3. Moderate priors----------------------------------------------------
moderate_prior_model <- brm(bf(Totaleggcount~
                                 F0_timepoint_binned+
                                 Temp+
                                 (1|PairID),
                               shape~F0_timepoint_binned),
                            data=data,
                            family=brms::negbinomial(),
                            #With specified priors added
                            prior = moderatepriors,
                            sample_prior = "only",
                            iter=5000,
                            cores = 4)
summary(moderate_prior_model)
saveRDS(moderate_prior_model  , file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/moderate_prior_model.rda")
moderate_prior_model <-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/moderate_prior_model.rda")


#Prior draws
moderate_prior_cumvdis<-pp_check(moderate_prior_model, 
                                 ndraws=100, 
                                 type = "ecdf_overlay")+theme_classic()+xlim(0,2000) 

#Good, but may follow the observed data too tightly...


# ===========================================================================
# 7. POSTERIOR MODELS---------------------------------------------------------
# ===========================================================================

#7.1. Diffuse priors-----------------------------------------------------------
mod1.1_diffusepriors <- brm(bf(Totaleggcount~
                                 F0_timepoint_binned+
                                 Temp+
                                 (1|PairID),
                               shape~F0_timepoint_binned),
                            data=data,
                            family=brms::negbinomial(),
                            #With specified priors added
                            prior = diffusepriors,
                            control=list(adapt_delta=0.95),
                            save_pars = save_pars(all = TRUE),
                            iter=5000,
                            cores = 4)
summary(mod1.1_diffusepriors)
saveRDS(mod1.1_diffusepriors, file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_diffusepriors.rda")
mod1.1_diffusepriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_diffusepriors.rda")

color_scheme_set("pink")

#Posterior draws
mod1.1_diffusepriors_cumvdis<-pp_check(mod1.1_diffusepriors, 
                                       ndraws=100, 
                                       type = "ecdf_overlay")+theme_classic()+xlim(0,2000) 


#7.2. Weak priors-[SELECTED MODEL]------------------------------------------------------------
mod1.1_weakprior <- brm(bf(Totaleggcount~
                             F0_timepoint_binned+
                             Temp+
                             (1|PairID),
                           shape~F0_timepoint_binned),
                        data=data,
                        save_warmup =TRUE,
                        family=brms::negbinomial(),
                        #With specified priors added
                        prior = weakpriors,
                        control=list(adapt_delta=0.98),
                        save_pars = save_pars(all = TRUE),
                        iter=5000,
                        cores = 4)
summary(mod1.1_weakprior)
saveRDS(mod1.1_weakprior, file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_weakprior.rda")
mod1.1_weakprior<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_weakprior.rda")

#Posterior draws
mod1.1_diffusepriors_cumvdis<-pp_check(mod1.1_weakprior, 
                                       ndraws=100, 
                                       type = "ecdf_overlay")+theme_classic()+xlim(0,2000) 

#7.3.Moderate priors--------------------------------------------

#moderate prior model
mod1.1_moderateprior <- brm(bf(Totaleggcount~
                                 F0_timepoint_binned+
                                 Temp+
                                 (1|PairID),
                               shape~F0_timepoint_binned),
                            data=data,
                            family=brms::negbinomial(),
                            #With specified priors added
                            prior = moderatepriors,
                            control=list(adapt_delta=0.95),
                            save_pars = save_pars(all = TRUE),
                            iter=5000,
                            cores = 4)
summary(mod1.1_moderateprior)
saveRDS(mod1.1_moderateprior  , file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_moderateprior.rda")
mod1.1_moderateprior <-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_moderateprior.rda")

#Posterior draws
mod1.1_moderateprior_cumvdis<-pp_check(mod1.1_moderateprior,
                                       ndraws=100, 
                                       type = "ecdf_overlay")+theme_classic()+xlim(0,2000) 


#7.4. Default model (no specified priors)----------------------------------------------------
mod1.1_defaultprior <- brm(bf(Totaleggcount~
                                F0_timepoint_binned+
                                Temp+
                                (1|PairID),
                              shape~F0_timepoint_binned),
                           data=data,
                           family=brms::negbinomial(),
                           #No specified priors added
                           control=list(adapt_delta=0.95),
                           save_pars = save_pars(all = TRUE),
                           iter=5000,
                           cores = 4)
summary(mod1.1_defaultprior)
saveRDS(mod1.1_defaultprior, file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_defaultprior.rda")
mod1.1_defaultprior<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod1.1_defaultprior.rda")

#---------------------Prior sensitivity analysis---------------------------------------
Loo_prior_performance<-LOO(mod1.1_defaultprior, 
                           mod1.1_diffusepriors, 
                           mod1.1_weakprior, 
                           mod1.1_moderateprior, 
                           moment_match=TRUE)
saveRDS(Loo_prior_performance, file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/priorfitsummary.rda")
Loo_prior_perfromance<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 fecundity/priorfitsummary.rda")


# ===========================================================================
# 8. SELECTED MODEL DIAGNOSTICS  (mod1.1_weakprior)-------------------------
# ===========================================================================

## LOO-CV --------------------------------------------------------------------
loo_weakpriors <- loo(mod1.1_weakprior, save_psis = TRUE)
plot(loo_weakpriors)
psis_weakpriors <- loo_weakpriors$psis_object
psis_weights    <- weights(psis_weakpriors)

## Posterior predictive checks -----------------------------------------------
#Posterior simiulated and empirical (observed) means
weakprior_mean     <- ppc_stat(yrep = posterior_predict(mod1.1_weakprior),
                               y    = data$Totaleggcount,
                               stat = "mean") + theme_classic()

#posterior cumulative density
weakprior_ecdf     <- pp_check(mod1.1_weakprior, ndraws = 100,
                               type = "ecdf_overlay") +
  xlim(0, 1000) + theme_classic()

#posterior probability density
weakprior_pdf      <- pp_check(mod1.1_weakprior, ndraws = 100) +
  xlim(0, 1000) + theme_classic()

#posterior LOO Q-Q plot
weakprior_loo_qq   <- pp_check(mod1.1_weakprior, type = "loo_pit_qq",
                               ndraws = 100) + theme_classic()
#LOO-PIT values fall neatly along the diagonal line (i.e. expected uniform distribution)

#LOO uniformity plot
weakprior_loo_unif <- pp_check(mod1.1_weakprior, type = "loo_pit_overlay",
                               ndraws = 100) + theme_classic()
#predictive interval plot
weakprior_intervals <- ppc_loo_intervals(
  y    = data$Totaleggcount,
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

ggsave("./bayesian_plots/model fit plots/F1 fecundity/selectedmodelfit.png",
       plot   = fit_plots,
       device = "png",
       width  = 500, height = 400, units = "mm")

## Bayes R² ------------------------------------------------------------------
bayes_R2(mod1.1_weakprior, re.form = NA)   # fixed effects only (i.e. marginal)
bayes_R2(mod1.1_weakprior, re.form = NULL)  # including random effects (i.e., conditional)

## MAP estimates and pd ------------------------------------------------------
MAP_fecundity <- bayestestR::describe_posterior(
  mod1.1_weakprior,
  test        = "pd",
  ci_method   = "HDI",
  centrality  = "MAP",
  component   = "all",
  effects     = "full",
  ci          = 0.95
)
saveRDS(MAP_fecundity,
        "scripts/model_outputs/Offspring Trait Models/F1 fecundity/MAP_fecundity.rda")

#manually estimating the ROPE range (as 10% of the response SD)
mean_y<-mean(data$Totaleggcount) #mean = 287.73

#Standard deviation on log scale
sd_y <- data %>%
  summarise(sd_log10 = sd(log10(Totaleggcount), na.rm = TRUE)*2.302585) %>%
  pull(sd_log10)

#SD = 1.342141

rope_response <- c(-0.1 * sd_y, 0.1 * sd_y)
#Rope = -0.13, 0.13 (very similar to default -0.1, +0.1)

## ROPE -----------------------------------------------------------------------
rope_fecundity <- bayestestR::describe_posterior(
  mod1.1_weakprior,
  rope_range=c(-0.134, 0.134),
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1 #ROPE estimates need to use the full posterior distribution, hence why they're calculated seperately
)
saveRDS(rope_adultmass,
        "scripts/model_outputs/Offspring Trait Models/adult mass/rope_adultmass.rda")

# ===========================================================================
# 9. PRIOR vs. POSTERIOR SENSITIVITY PLOTS------------------------------------
# ===========================================================================

## Combined prior–posterior ECDF panel ---------------------------------------
priorpostplots<-ggarrange(diffuse_prior_cumvdis, weak_prior_cumvdis, 
                          moderate_prior_cumvdis, mod1.1_diffusepriors_cumvdis, 
                          weakprior_ecdf, mod1.1_moderateprior_cumvdis, nrow = 2, ncol = 3,
                          labels = c("A", "B", "c", "D", "E", "F"))

ggsave(filename = "./bayesian_plots/model fit plots/F1 fecundity/priorvpost_ecdf.png",
       plot = priorpostplots, 
       device = "png", 
       width = 390, 
       height = 220, 
       units = "mm")

# ===========================================================================
# 10. HYPOTHESIS TESTING: INTERACTION MODELS
#Models test the core interactions of the paper: are parental age effects environmentally dependent?
# ===========================================================================

#10.1. Fitting an age: temperature interaction
mod2.1<-brm(bf(Totaleggcount~
                 F0_timepoint_binned+
                 Temp+
                 Temp:F0_timepoint_binned+
                 (1|PairID),
               shape~F0_timepoint_binned),
            data=data,
            family=brms::negbinomial(),
            #With specified priors added
            prior = weakpriors,
            control=list(adapt_delta=0.98),
            save_pars = save_pars(all = TRUE),
            iter=5000,
            cores = 4)
summary(mod2.1) 
saveRDS(mod2.1 , file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod2.1.rda")
mod2.1<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 fecundity/mod2.1.rda")


## MAP estimates for interaction model (mod3.1) ------------------------------
#For reporting in the main manuscript
MAP_fecundity_mod2.1 <- bayestestR::describe_posterior(
  mod2.1,
  test       = "pd",
  ci_method  = "HDI",
  centrality = "MAP",
  effects    = "full",
  ci         = 0.95
)
saveRDS(MAP_fecundity_mod2.1,
        "scripts/model_outputs/Offspring Trait Models/F1 fecundity/MAP_fecundity_mod2.1.rda")

#ROPE estimates for the interactive effect (using the full posterior distribution)
rope_fecundity_mod2.1 <- bayestestR::describe_posterior(
  mod2.1,
  rope_range= c(-0.134, 0.134),
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1
)
saveRDS(rope_fecundity_mod2.1,
        "scripts/model_outputs/Offspring Trait Models/F1 fecundity/rope_fecundity_mod2.1.rda")

#No evidence for an interaction between temperature and parental age

#------------------------Model selection for hypothesis series----------------------------------
hypothesis_loo<-loo(mod1.1_weakprior, mod2.1, moment_match=TRUE)
saveRDS(hypothesis_loo , file = "scripts/model_outputs/Offspring Trait Models/F1 fecundity/hypothesis_loo.rda")
hypothesis_loo<-readRDS("scripts/model_outputs/Offspring Trait Models/F1 fecundity/hypothesis_loo.rda")
#No difference in the predictive abilities between model 1 and model 2

#----------------------ESTIMATED MARGINAL MEANS FROM BEST-FITTING MODEL---------
#Marginal effect of parental age on offspring fecundity
pairwise_estimates<- emmeans(mod1.1_weakprior, "F0_timepoint_binned", type="response")
pairs(pairwise_estimates) 

# ===========================================================================
# 11. POSTERIOR DISTRIBUTION PLOTS  (mod1.1_weakpriors)-----------------------
# ===========================================================================
# Plots styled using stat_halfeye (ggdist) for consistency with other traits

post1 <- as_draws_df(mod1.1_weakprior) #for single effect predictors
post2<-as_draws_df(mod2.1) #for interaction term

## Build combined posterior data frame ----------------------------------------
posterior_df <- data.frame(
  "µ: Parents' age at reproduction (Late-Aged)" = post1$b_F0_timepoint_binnedlate,
  "µ: Parents' temperature treatment (28.0°C)" = post1$b_Temp28,
  "µ: Parents' temperature treatment (30.5°C)" = post1$b_Temp30.5,
  "ϕ: Parents' age at reproduction (Late-Aged)" = post1$b_shape_F0_timepoint_binnedlate,
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
    "ϕ: Parents' age at reproduction (Late-Aged)"
  )))

y_levels <- rev(c(
  "µ: Parents' age at reproduction (Late-Aged)",
  "µ: Parents' temperature treatment (28.0°C)",
  "µ: Parents' temperature treatment (30.5°C)",
  "μ: Parents' age × Temperature (28.0°C)",
  "μ: Parents' age × Temperature (30.5°C)" ,
  "ϕ: Parents' age at reproduction (Late-Aged)"
))

## Halfeye posterior plot ------------------------------------------------------
posterior_plot_fecundity<- ggplot(
  posterior_df,
  aes(x = value, y = parameter, fill = parameter)
) +
  annotate("rect",
           xmin = -0.134, xmax = 0.134,
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
    "#e5a9f5", "#4c9279"
  )) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title      = element_text(size = 30),
    axis.text       = element_text(size = 30)
  ) +
  labs(x = "Fixed-effect posterior estimates\n(on μ = mean offspring fecundity, log-scale; on ϕ = shape parameter, log scale)",
       y = NULL)

ggsave("./bayesian_plots/posterior plots/F1 fecundity/001_hypothesis_halfeye.png",
       plot   = posterior_plot_fecundity,
       device = "png",
       width  = 580, height = 400, units = "mm")

## Post-hoc standardised mean difference (SMD) --------------------------------
pairwise_estimates<- emmeans(mod1.1_weakprior, "F0_timepoint_binned", type = "response")
#This is now converted to the raw scale where the extreme ends of parental age are compared
sd_life<-sd(data$Totaleggcount)

#standardised emmeans --> estimating difference across the entire range of parental age

#Then standardising by SD
emm_life_scaled <- summary(pairwise_estimates)  %>%
  mutate(emmean_sd = prob/ sd_life,
         lower_sd = lower.HPD / sd_life,
         upper_sd = upper.HPD / sd_life)

# Calculate SMD and its CI
SMD <- with(emm_life_scaled, {
  emmean_sd[2] - emmean_sd[1]              # SMD
})
SMD_lower <- with(emm_life_scaled, {
  lower_sd[2] - upper_sd[1]                # lower bound of SMD
})
SMD_upper <- with(emm_life_scaled, {
  upper_sd[2] - lower_sd[1]                # upper bound of SMD
})

# Combine results
fecundity_SMD<-data.frame(
  SMD = SMD,
  SMD_lower = SMD_lower,
  SMD_upper = SMD_upper
)

#Saving for use in combined forest plot of all traits
saveRDS(fecundity_SMD, file = "scripts/model_outputs/Offspring Trait Models/fecundity_SMD")


# ===========================================================================
# 13. MODEL PREDICTIONS AND RESULTS PLOTS-------------------------------------
# ===========================================================================

#13.2: creating a prediction data frame
df_predict_within<-expand.grid(
  F0_timepoint_binned=unique(data$F0_timepoint_binned),
  Temp=unique(data$Temp),
  PairID=unique(data$PairID)[1])


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
    mean = mean(prediction),
    lower= mean(lower),
    upper=mean(upper)) %>% 
  ungroup()


# plotting the predicted values onto a ggplot######
Fecundity_plot <- ggplot(data = data,
                         aes(x = F0_timepoint_binned, 
                             y = Totaleggcount,
                             fill = F0_timepoint_binned,
                             colour=F0_timepoint_binned)) +
  geom_point(position = position_jitter(width = 0.2, height=0.1),
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
  theme(#legend.position = "none",
    legend.title=element_text(size=50),
    legend.text=element_text(size=50),
    strip.text = element_text(size=50),
    panel.background = element_rect(fill = "transparent", color = NA),  # Make panel background transparent
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA), # Make plot background transparent
  )  +
  scale_color_manual(name = "Parents' age at reproduction (F0)", labels = c("Early-Aged", "Late-Aged"), values = c("#f96161", "#066594"))+
  scale_fill_manual(name = "Parents' age at reproduction (F0)",labels = c("Early-Aged", "Late-Aged"), values = c("#f96161", "#066594"))+
  scale_x_discrete(labels=c("Early-Aged", "Late-Aged"))+
  labs(x = "Parents' age at reproduction", y = "Number of eggs laid by offspring")


# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 fecundity/001_F1_meanfecundity_updated2.png",
       plot = Fecundity_plot,
       bg="transparent",
       device = "png", 
       width = 355, 
       height = 650, 
       units = "mm")


#Plotting the Posterior mass function for F1 fecundity---------------------------------------------

negbinom_pdf <- function(mu, shape, x) {
  dnbinom(x, mu=mu, size=shape)
} #Calculates the probability density function based on the models mu and sigma parameters



#Proxy data
newdat <- expand.grid(
  F0_timepoint_binned = unique(data$F0_timepoint_binned),
  Temp= unique(data$Temp),
  PairID=unique(data$PairID[1]))


# Create a sequence over which to evaluate the PDF
fecundity_seq <- seq(
  floor(min(eggdata$Totaleggcount, na.rm = TRUE)),
  ceiling(max(eggdata$Totaleggcount, na.rm = TRUE)),
  by = 1
)

#Creating the PMF for F1 fecundity
proportional_posterior <- mod1.1_weakprior %>% 
  linpred_draws(
    newdat,
    value = "mu",
    allow_new_levels = TRUE, 
    transform = TRUE,
    re_formula = NA,  # population-level predictions
    dpar = "shape",
    ndraws = 3000, 
    seed = 123
  ) 

#Calculating the PMF for each draw
posterior_draws<-proportional_posterior %>% 
  crossing(x = fecundity_seq) %>%        
  mutate(f = negbinom_pdf(mu, shape, x)) %>% 
  group_by(x, .draw, F0_timepoint_binned) %>% 
  summarise(f = median(f), .groups = "drop")

#summarizing across draws
summary_probability <-posterior_draws %>%
  group_by(x, F0_timepoint_binned) %>%
  summarise(
    median = median(f),  
    lower  = quantile(f, 0.025),
    upper  = quantile(f, 0.975),
    .groups = "drop"
  )


#Getting the median fecundity
median_fecundity<- mod1.1_weakprior %>% 
  linpred_draws(
    newdat,
    value = "mu",
    allow_new_levels=TRUE,
    transform = TRUE,
    re_formula=NA,#population level predictions
    dpar = "shape",
    ndraws = 5000, 
    seed = 123
  ) %>% 
  ungroup() %>%
  group_by(.draw, F0_timepoint_binned) %>% 
  summarise(median_fecundity= median(mu)) %>% 
  group_by(F0_timepoint_binned) %>% 
  median_hdi(median_fecundity)

#--------------------------PLOTTING probability mass------------------

#Transforming from the PDF to the PMF (condensing counts into bins of 50 eggs)
summary_probability_binned <- summary_probability %>%
  mutate(x_bin = floor(x/50)*50) %>% 
  group_by(x_bin, F0_timepoint_binned) %>%
  summarise(median = sum(median),
            lower = sum(lower),
            upper = sum(upper), .groups = "drop")


#Probability density function---------------------------------------------------
f1fecundity_pdf<-ggplot(data=summary_probability_binned,
                        aes(x = x_bin, 
                            y = median, 
                            colour= F0_timepoint_binned)) +
  geom_col(position = position_dodge(width = 0.95),  
           fill = "grey75",
           color = "grey75",
           width = 7.8, alpha = 1) +
  # Transparent ribbons for 95% credible intervals
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 50), 
                width = 15, linewidth = 4) +
  geom_segment(data = median_fecundity,
               aes(x = median_fecundity, xend = median_fecundity, y = 0, yend = 0.3, color = F0_timepoint_binned),
               linetype = "dashed", linewidth = 4, alpha = 0.7)+
  theme_classic()+
  scale_y_continuous(limits=c(0, 0.30))+
  theme(axis.title=element_text(size=45),
        axis.text=element_text(size=45))+
  theme(legend.position = "top",
        legend.title=element_text(size=45),
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
  labs(x = "Number of eggs laid by offspring", y = "Probability mass, f(x)")

# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 fecundity/003_fecundity_pmf.png",
       plot = f1fecundity_pdf, 
       bg="transparent",
       device = "png", 
       width = 400, 
       height = 400, 
       units = "mm")



#----------------Combined plot--------------------------------------------------

#Saving the prior and poosterior plots for default versus moderate priors
inference2<-ggarrange(
  Fecundity_plot,
  f1fecundity_pdf,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B"),
  font.label = list(size = 50, face = "bold"),
  widths = c(1, 1.5) 
)

ggsave(filename = "./bayesian_plots/offspring trait plots/F1 fecundity/inference_plots.png",
       plot = inference2, 
       device = "png", 
       width = 740, 
       height = 540, 
       units = "mm")


#Creating a table to export------------------------------------------------------------------------------
# ── helper ────────────────────────────────────────────────────────────────────
#Model draws from selected models
base_model_draws<-as_draws_df(mod1.1_weakprior)
interaction_draws<-as_draws_df(mod2.1) 

#function for generating estimates
summarise_param <- function(draws, param, section,
                            rope = TRUE, #using default ROPE range automatically calculated from bayestest
                            rope_range = c(-0.134, +0.134),
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
shape<- bind_rows(
  summarise_param(base_model_draws$b_shape_Intercept, "Shape intercept (ϕ)", "Distributional parameters (scale component)",  rope = FALSE),
  summarise_param(base_model_draws$b_shape_F0_timepoint_binnedlate, "Parents' age category: Late-aged", "Distributional parameters (scale component)", rope= FALSE)
)

#---6. Conditional and marginal Bayes R² ------------------------------------------------------------------
marginal<-bayes_R2(mod1.1_weakprior, re.form = NA, summary = FALSE)   # fixed effects only (i.e., marginal)
conditional<-bayes_R2(mod1.1_weakprior, re.form = NULL, summary = FALSE)  # including random effects (i.e., conditional)

bayes<-bind_rows(
  summarise_param(marginal, "Marginal R²", "Bayes R²", pd= FALSE, rope= FALSE),
  summarise_param(conditional, "Conditional R²", "Bayes R²", pd= FALSE, rope =FALSE)
)


# ── 5. Combine and render ─────────────────────────────────────────────────────
bind_rows(fe, fe_interaction2, re_sd_explore, shape, bayes) %>%
  gt(groupname_col = "Section") %>%
  cols_label(Parameter = "Parameter", MAP = "MAP",
             `95% HDI` = "95% HDI", `% in ROPE` = "% in ROPE", pd = "pd") %>%
  tab_header(
    title    = "Offspring fecundity model summary: Negative binomial model"
  ) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) %>%
  tab_footnote("ROPE = [−0.134, 0.134] on log scale; applied to location fixed effects only.",
               locations = cells_column_labels(`% in ROPE`)) %>%
  tab_footnote("pd = probability of direction.",
               locations = cells_column_labels(pd)) %>%
  opt_stylize(style = 1)%>%
  gtsave("./tables/fecundity_updated.docx")


###############################################################################
##END OF SCRIPT----------------------------------------------------------------
###############################################################################
