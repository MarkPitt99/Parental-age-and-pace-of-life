###############################################################################
##  Script: Effect of parental age on F1 total lifespan
##  Note:   Analysis restricted to offspring that survived  past their first four weeks
#So mean survival estimates are conditional on having already survived your first month
#Script provides an insight into average longevity only - see BaSTA scripts for how mortality
#parameters are affected by parental age
###############################################################################


# ===========================================================================
# 1. SETUP--------------------------------------------------------------------
# ===========================================================================

## Libraries ------------------------------------------------------------------
library(dplyr) #v.1.1.4
library(tidyr) #v.1.3.1
library(ggplot2) #v.4.0.0
library(ggpubr) #v.0.6.1
library(ggdist) #v.3.3.3
library(patchwork) #v.1.3.0
library(brms) #v.2.21.0
library(bayesplot) #v.1.11.1
library(bayestestR) #v.0.17.0
library(posterior) #v.1.6.0
library(tidybayes) #v.3.0.7
library(emmeans) #v.1.10.1
library(parameters) #v.0.28.2
library(loo) #v.2.8.0
library(sjPlot) #v.2.8.16
library(lme4) #v.1.1.35.2
library(ggdist) #v.3.3.3
library(car) #v.3.1.2
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
length(unique(F1data$F1_ID))  # 1316 offspring (in total)

#Only retaining animals that survived their first four weeks
data1<-F1data%>%
  filter(!is.na(total_lifespan)) %>% 
  filter(total_lifespan>4) 

#Creating an "event" column to signal when the offspring died
#required for time to event models (e.g., Weibull, lognormal, gamma, exponential)
data1$event<-ifelse(is.na(data1$adult_lifespan), 0,1)

#Creating a lifespan variable
data1$total_lifespan <- as.numeric(difftime(data1$F1_death_date, data1$F1_hatch, units = "weeks"))

## Scaling continuous predictors ------------------------------------------------

# All predictors scaled to mean = 0, SD = 1; a 1-unit change reflects 1 SD
data1$avg_age <- as.numeric(scale(data1$avg.age, center = TRUE, scale = TRUE))


# Within-individual (delta) age: scaled but NOT mean-centred to preserve
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
length(unique(data1$Mother_ID))  # number of mothers (n = 77)
length(unique(data1$F1_ID))      # number of F1 offspring (987 entered the model, so survived their first month)

# ===========================================================================
# 3. SUMMARY STATISTICS--------------------------------------------------------
# ===========================================================================

## Distribution of F1 total lifespan ------------------------------------------
ggplot(data1, aes(x = total_lifespan)) +
  geom_histogram(binwidth = 1, fill = "skyblue", colour = "black") +
  labs(x = "F1 total lifespan", y = "Count") +
  theme_classic()

## Overall mean and SD (used to build intercept priors) ----------------------
sum_overall <- data1 %>%
  summarise(
    F1_mass     = mean(total_lifespan),
    sd_total    = sd(total_lifespan),
    n_offspring = n()
  )

range(data1$total_lifespan)  #4.14 weeks – 30.57 weeks

#mean lifespan of 18.3 weeks, SD of 3.97 weeks, n =987, no.parents = 77

## Collinearity checks --------------------------------------------------------
X <- model.matrix(
  ~ avg.age + delta.age + Temp + F1_sex,
  data = data.frame(
    avg.age                       = data1$avg.age,
    delta.age                     = data1$delta.age,
    Temp                          = data1$Temp,
    F1_sex                        = data1$F1_sex,
    cumulative_matings            = data1$cum_successful_matings #correlated with delta age, removed from inference
  )
)

cor(X[, -1])  # correlation matrix (excluding intercept column)

vif(lm(total_lifespan ~ avg.age + delta.age + Temp +cum_successful_matings,
       data = data1)) #number of succesful matings tightly correlated with delta age

# ===========================================================================
# 4. MODEL STRUCTURE SELECTION------------------------------------------------
# ===========================================================================

#Deciding on a family to use for the models----------------------
#Data is survival data (discrete time-to-event)

## 4.1: Weibull model-----------------------------------------------------------
default_weibull<-brm(total_lifespan|cens(1-event)~
                       avg_age+
                       delta.age+
                       Temp+
                       (1|PairID),
                     family=weibull(),
                     data=data1,
                     iter=5000,
                     save_pars = save_pars(all = TRUE),
                     cores = 4)
saveRDS(default_weibull, 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/defaultweibull.rda")
default_weibull<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/defaultweibull.rda")

#diagnostics
loo_Weibull<- loo(default_weibull, save_psis = TRUE, moment_match = TRUE)
plot(loo_Weibull)

#Loo probability integral transform (PIT) plots
pp_check(default_weibull, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(default_weibull, type ="loo_pit_overlay", ndraws=100) 

#Inspecting the fit with a KM curve
default_weibull_KM<- pp_check(default_weibull,
                              status_y=data1$event,
                              type="km_overlay", ndraws=100) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  ylab("Cumulative Survival Probability, S(x)") + coord_cartesian(xlim=c(0, 20))+
  xlab("Offspring adult age (weeks)") + coord_cartesian(xlim=c(0, 20))

## 4.2: Lognormal model---------------------------------------------------------
default_lognormal<-brm(total_lifespan|cens(1-event)~
                         avg_age+
                         delta.age+
                         Temp+
                         (1|PairID),
                       family=brms::lognormal(),
                       data=data1,
                       iter=5000,
                       save_pars = save_pars(all = TRUE),
                       cores = 4)
saveRDS(default_lognormal, 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/defaultlognormal.rda")
default_lognormal<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/defaultlognormal.rda")

#diagnostics
loo_lognormal<- loo(default_lognormal, save_psis = TRUE, moment_match = TRUE)
plot(loo_lognormal)

#Loo probability integral transform (PIT) plots
pp_check(default_lognormal, type="loo_pit_qq", ndraws=100)
pp_check(default_lognormal, type ="loo_pit_overlay", ndraws=100) 


#Inspecting the fit with a KM curve
default_lognormal_KM<- pp_check(default_lognormal,
                              status_y=data1$event,
                              type="km_overlay", ndraws=100) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  ylab("Cumulative Survival Probability, S(x)") + coord_cartesian(xlim=c(0, 20))+
  xlab("Offspring adult age (weeks)") + coord_cartesian(xlim=c(0, 20))

#4.3: Default exponential-------------------------------------------------------
#model assumes a constant mortality risk through time
default_exponential<-brm(total_lifespan|cens(1-event)~
                           avg_age+
                           delta.age+
                           Temp+
                           (1|PairID),
                         family=exponential(),
                         data=data1,
                         iter=5000,
                         save_pars = save_pars(all = TRUE),
                         cores = 4)
saveRDS(default_exponential, 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/defaultexponential.rda")
default_exponential<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/defaultexponential.rda")

#Diagnostics
loo_exponential<- loo(default_exponential, save_psis = TRUE, moment_match = TRUE)
plot(loo_exponential)

#Loo probability integral transform (PIT) plots
pp_check(default_exponential, type="loo_pit_qq", ndraws=100) 
pp_check(default_exponential, type ="loo_pit_overlay", ndraws=100) 


#Inspecting the fit with a KM curve
default_exponential_KM<- pp_check(default_exponential,
                                status_y=data1$event,
                                type="km_overlay", ndraws=100) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  ylab("Cumulative Survival Probability, S(x)") + coord_cartesian(xlim=c(0, 20))+
  xlab("Offspring adult age (weeks)") + coord_cartesian(xlim=c(0, 20))


#--------------------Comparison of the time to event models------------------------------------------
LOO_allsurvivalmodels<-loo(default_weibull, default_lognormal, 
                           default_exponential,
                           moment_match = TRUE, reloo = TRUE)
saveRDS(LOO_allsurvivalmodels, 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/LOO_allsurvivalmodels.rda")

loomod1<-loo(default_weibull, reloo = TRUE)
loomod2<-loo(default_lognormal, reloo=TRUE)
loomod3<-loo(default_exponential, reloo=TRUE)

loo_allsurvivalmodels2<-loo_compare(loomod1,
                                    loomod2,
                                    loomod3)
saveRDS(loo_allsurvivalmodels2, file = "scripts/model_outputs/Offspring Trait Models/total lifespan/LOO_allsurvivalmodels2.rda")
#Choosing the Weibull model as the best fitting

#Deciding on whether the effect of delta age should be a linear or a quadratic term-----
## 4.4 Quadratic delta-age term ------------------------------------------------
mod1.1_quadraticage<-brm(total_lifespan|cens(1-event)~
                           avg_age+
                           poly(delta.age,2)+ #No substantial evidence for the quadratic term
                           Temp+
                           (1|PairID),
                         family=weibull,
                         data=data1,
                         iter=5000,
                         control=list(adapt_delta=0.98),
                         save_pars = save_pars(all = TRUE),
                         core=4)
summary(mod1.1_quadraticage) 
saveRDS(mod1.1_quadraticage, 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_quadraticage.rda")
mod1.1_quadraticage<-readRDS("scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_quadraticage.rda")

#MODEL DIAGNOSTICS
loo_quadratic<-loo(mod1.1_quadraticage, save_psis = TRUE, moment_match = TRUE)
plot(loo_quadratic)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_quadraticage, type="loo_pit_qq", ndraws=100) 
pp_check(mod1.1_quadraticage, type ="loo_pit_overlay", ndraws=100) 


## 4.5 Random slopes (intercept + delta age) -----------------------------------
#where population-level slope of delta age is a linear slope
mod1.1_randomslopes<-brm(total_lifespan|cens(1-event)~
                           avg_age+ 
                           delta.age + 
                           Temp + 
                           (1+delta.age|PairID),
                         family=weibull,
                         data=data1,
                         iter=5000,
                         save_pars = save_pars(all = TRUE),
                         control=list(adapt_delta=0.95),
                         core=4)
summary(mod1.1_randomslopes)
saveRDS(mod1.1_randomslopes , file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_randomslopes.rda")
mod1.1_randomslopes<-readRDS("scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_randomslopes.rda")

#diagnostics
loo_randomslopes<-loo(mod1.1_randomslopes, save_psis=TRUE, moment_match = TRUE)
plot(loo_randomslopes)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_randomslopes, type="loo_pit_qq", ndraws=100) 
pp_check(mod1.1_randomslopes, type ="loo_pit_overlay", ndraws=100) 


## 4.6 Random slopes + quadratic delta-age -------------------------------------
#where population level slope is a quadratic function (but assumes all females share the same quadratic trajectory)
mod1.1_randomslopes_quadratic<-brm(total_lifespan|cens(1-event)~
                                     avg_age+ 
                                     poly(delta.age,2)+
                                     Temp + 
                                     (1+delta.age|PairID),
                                   family=weibull,
                                   data=data1,
                                   iter=5000,
                                   save_pars = save_pars(all = TRUE),
                                   control=list(adapt_delta=0.95),
                                   core=4)
summary(mod1.1_randomslopes_quadratic)
saveRDS(mod1.1_randomslopes_quadratic , 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_randomslopes_quadratic.rda")
mod1.1_randomslopes_quadratic<-readRDS("scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_randomslopes_quadratic.rda")

#diagnostics
loo_randomslopes_quadratic<-loo(mod1.1_randomslopes_quadratic, save_psis=TRUE, moment_match = TRUE)
plot(loo_randomslopes)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_randomslopes_quadratic, type="loo_pit_qq", ndraws=100) 
pp_check(mod1.1_randomslopes_quadratic, type ="loo_pit_overlay", ndraws=100) 


## 4.7 LOO comparison: random-effects and quadratic structure ------------------
LOO_ageeffects<-loo(default_weibull,mod1.1_randomslopes, 
                    mod1.1_randomslopes_quadratic, mod1.1_quadraticage,
                    moment_match = TRUE)
saveRDS(LOO_ageeffects, 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/LOO_ageeffects.rda")
LOO_ageeffects<-readRDS("scripts/model_outputs/Offspring Trait Models/total lifespan/LOO_ageeffects.rda")

loomod5<-loo(mod1.1_quadraticage)
loomod6<-loo(mod1.1_randomslopes)
loomod7<-loo(mod1.1_randomslopes_quadratic)

looageeffects2<-loo_compare(loomod6,
                            loomod1,
                            loomod5,
                            loomod7)
saveRDS(looageeffects2, file = "scripts/model_outputs/Offspring Trait Models/total lifespan/looageeffects2.rda")


## 4.8 Distributional models: evidence for shape(k) submodel ------------------
# Does parental age (and/or temperature) predict the shape of the hazard curve?-----
#4.9. Testing for a temperature effect on the shape parameter----
mod1.1_shape_temp<-brm(bf(total_lifespan|cens(1-event)~
                            avg_age+ 
                            delta.age +
                            Temp + 
                            (1+delta.age|PairID),
                          shape~
                            delta.age+
                            Temp),
                       family=weibull(),
                       data=data1,
                       iter=5000,
                       save_pars = save_pars(all = TRUE),
                       control=list(adapt_delta=0.98),
                       core=4) #Model struggles to converge
summary(mod1.1_shape_temp)
saveRDS(mod1.1_shape_temp , 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_shape_temp.rda")
mod1.1_shape_temp<-readRDS("scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_shape_temp.rda")

#diagnostics
loo_shapetemp<-loo(mod1.1_shape_temp, save_psis = TRUE, moment_match = TRUE)
plot(loo_shapetemp)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_shape_temp, type="loo_pit_qq", ndraws=100) 
pp_check(mod1.1_shape_temp, type ="loo_pit_overlay", ndraws=100) 


#4.10. The effect of just parental age on shape parameter--------------------------
mod1.1_shape<-brm(bf(total_lifespan|cens(1-event)~
                       avg_age+ 
                       delta.age +
                       Temp + 
                       (1+delta.age|PairID),
                     shape~
                       delta.age),
                  family=weibull(),
                  data=data1,
                  iter=5000,
                  save_pars = save_pars(all = TRUE),
                  control=list(adapt_delta=0.95),
                  core=4) #model is stuggling to converge, estimating shape may not be feasible in brms
summary(mod1.1_shape)
saveRDS(mod1.1_shape , file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_shape.rda")
mod1.1_shape<-readRDS("scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_shape.rda")

#diagnostics
loo_shape<-loo(mod1.1_shape, save_psis = TRUE)
plot(loo_shape)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_shape, type="loo_pit_qq", ndraws=100) 
pp_check(mod1.1_shape, type ="loo_pit_overlay", ndraws=100) 


#---------Both of the models with covariates on the shape parameter fail to converge-----
##------------Now assume a constant shape parameter across fixed effects--------------------
#This is unlikley to be the case, hence we also used our BaSTA framework

# ===========================================================================
# 5. PRIOR SPECIFICATION-----------------------------------------------------
# ===========================================================================
# Priors assume all continuous predictors are mean-centred and SD-scaled.
# Intercept prior informed by observed mean (18.3 weeks) and SD (3.97 weeks).

#---------------------Setting up the priors-------------------------------------
#what's the range of our data?
range(log(data1$total_lifespan)) 
mean(data1$total_lifespan)#18.30 weeks
sd(data1$total_lifespan)#3.965274 weeks

#SD on the empirical log scale
sd(log10(data1$total_lifespan))*2.302585 #=0.2585665

#Mean on the log scale
log(18.31) #2.907447 mean

#SD on theoretical log scale
sqrt(log(1 + (3.965274^2 / 18.30^2))) #0.214 on log scale

#No correction needed for the log(mean) in the weibull model

#---------------------SETTING MODEL PRIORS--------------------------------------

#Diffuse priors-----------------------------------------------------------------
#very similar to the flat brms priors, but the intercept is slightly more regularised
Diffusepriors <- c(
  prior(normal(2.907, 0.214), class = "Intercept"),             
  prior(normal(0,1), class = "b"),          
  prior(normal(0, 2.5), class ="shape", lb = 0),
  prior(student_t(3, 0, 2.5), class = "sd", lb = 0),           
  prior(lkj(2), class = "cor")                                 
)


#Weakly informative priors [SELECTED]-------------------------------------------
#assumes most effects are within 1SD of the response, and rarely beyond 1.96SDs
weakpriors <- c(
  prior(normal(2.907,0.214), class = "Intercept"),                
  prior(normal(0,0.214), class = "b"), 
  prior(normal(1, 2.5), class = "shape", lb=0), #brms model struggles to add parameters to the shape of the hazard         
  prior(exponential(10), class = "sd", lb = 0),                      
  prior(lkj(2), class = "cor")                                   
)


#Moderate prior-----------------------------------------------------------------
moderatepriors <- c(
  prior(normal(2.907, 0.214), class = "Intercept"),                 
  prior(normal(0,0.1), class = "b"),  
  prior(normal(1,1.5), class="shape", lb=0),  
  prior(exponential(10), class = "sd", lb = 0),                         
  prior(lkj(2), class = "cor")                                    
)

#Constraints prior --> assuming parental age reduces lifespan (i.e., an expected Lansing effect)
constraintpriors <- c(
  prior(normal(2.907, 0.214), class = "Intercept"),                 
  prior(normal(0,0.1), class = "b"),                               
  prior(normal(-0.186, 0.03), class = "b", coef="delta.age"), #estimate based on meta-analysis by Ed.Ivemy cook and Moorad (2023)
  prior(normal(1, 2.5), class ="shape"), 
  prior(exponential(10), class = "sd", lb = 0),                          
  prior(lkj(2), class = "cor")                                    
)


# ===========================================================================
# 6. PRIOR PREDICTIVE CHECKS---------------------------------------------------
# ===========================================================================

## 6.1 Prior-only model: diffuse ---------------------------------------------------
diffuse_prior_model<-brm(total_lifespan|cens(1-event)~
                           avg_age+ 
                           delta.age +
                           Temp + 
                           (1+delta.age|PairID),
                         family=weibull,
                         data=data1,
                         iter=5000,
                         #With specified priors added
                         prior = Diffusepriors,
                         sample_prior = "only",
                         cores = 4)

saveRDS(diffuse_prior_model , 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/diffusepriormodel.rda")
diffuse_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/diffusepriormodel.rda")

#setting the colour scheme
color_scheme_set("teal")

#Prior draws (I.e., the prior cumulative density function)
diffuseprior_cumvdis <- pp_check(diffuse_prior_model,
                                 ndraws = 100,
                                 type = "ecdf_overlay",
                                 discrete = FALSE) + scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 50)) + theme_classic()


## 6.2 Prior-only model: weak ----------------------------------------------------
weak_prior_model <-brm(bf(total_lifespan|cens(1-event)~
                            avg_age+ 
                            delta.age +
                            Temp + 
                            (1+delta.age|PairID)),
                       family=weibull,
                       data=data1,
                       iter=5000,
                       #With specified priors added
                       prior = weakpriors,
                       sample_prior = "only",
                       cores = 4)
summary(weak_prior_model)
saveRDS(weak_prior_model , file = "scripts/model_outputs/Offspring Trait Models/total lifespan/weakpriors.rda")
weak_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/weakpriors.rda")

#prior draws (for the weak priors)
weakprior_cumvdis <- pp_check(weak_prior_model,
                              ndraws = 100,
                              type = "ecdf_overlay",
                              discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 50)) +
  theme_classic()

##6.3 Prior-only model: moderate --------------------------------------------------
moderate_prior_model <-brm(total_lifespan|cens(1-event)~
                             avg_age+ 
                             delta.age +
                             Temp + 
                             (1+delta.age|PairID),
                           family=weibull,
                           data=data1,
                           iter=5000,
                           prior = moderatepriors,
                           sample_prior = "only",
                           cores = 4)
summary(moderate_prior_model)
saveRDS(moderate_prior_model , 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/moderatepriors.rda")
moderate_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/moderatepriors.rda")

#prior draws
moderateprior_cumvdis <- pp_check(moderate_prior_model,
                                  ndraws = 100,
                                  type = "ecdf_overlay",
                                  discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 50)) +
  theme_classic()

##6.4 Constraint prior draws-----------------------------------------------------
constraint_prior_model <-brm(total_lifespan|cens(1-event)~
                               avg_age+ 
                               delta.age +
                               Temp + 
                               F1_sex + 
                               (1+delta.age|PairID),
                             family=weibull,
                             data=data1,
                             iter=5000,
                             prior =constraintpriors,
                             sample_prior = "only",
                             cores = 4)
summary(constraint_prior_model)
saveRDS(constraint_prior_model , 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/constraintpriors.rda")
constraint_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/constraintpriors.rda")

#posterior draws
constraintprior_cumvdis<- pp_check(constraint_prior_model,
                                   ndraws = 100,
                                   type = "ecdf_overlay",
                                   discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 50)) +
  theme_classic()

# ===========================================================================
# 7. POSTERIOR MODELS---------------------------------------------------------
# ===========================================================================

## 7.1 Diffuse priors ---------------------------------------------------------
mod1.1_diffusepriors<-brm(total_lifespan|cens(1-event)~
                            avg_age+ 
                            delta.age +
                            Temp + 
                            (1+delta.age|PairID),
                          family=weibull,
                          data=data1,
                          prior = Diffusepriors,
                          iter=5000,
                          control=list(adapt_delta=0.98),
                          save_pars = save_pars(all = TRUE),
                          core=4)
summary(mod1.1_diffusepriors) 
saveRDS(mod1.1_diffusepriors, file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_diffusepriors.rda")
mod1.1_diffusepriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_diffusepriors.rda")

color_scheme_set("pink")

#posterior draws
mod1.1_diffusepriors_cumvdis<-
  pp_check(mod1.1_diffusepriors,
           ndraws = 100,
           type = "ecdf_overlay",
           discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 50)) +
  theme_classic()

## 7.2 Weak priors  [SELECTED MODEL] -------------------------------------------
mod1.1_weakpriors<-brm(total_lifespan|cens(1-event)~
                         avg_age+ 
                         delta.age +
                         Temp + 
                         (1+delta.age|PairID),
                       family=weibull,
                       data=data1,
                       prior = weakpriors,
                       iter=5000,
                       control=list(adapt_delta=0.98),
                       save_pars = save_pars(all = TRUE),
                       core=4)

summary(mod1.1_weakpriors)
saveRDS(mod1.1_weakpriors, 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_weakpriors.rda")
mod1.1_weakpriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_weakpriors.rda")

#posterior draws
mod1.1_moderatepriors_cumvdis<-
  pp_check(mod1.1_moderatepriors,
           ndraws = 100,
           type = "ecdf_overlay",
           discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 50)) +
  theme_classic()

## 7.3 Moderate priors ---------------------------------------------------------
mod1.1_moderatepriors<-brm(total_lifespan|cens(1-event)~
                             avg_age+ 
                             delta.age +
                             Temp + 
                             (1+delta.age|PairID),
                           family=weibull,
                           data=data1,
                           #With specified priors added
                           prior = moderatepriors,
                           iter=5000,
                           control=list(adapt_delta=0.98),
                           core=4)
summary(mod1.1_moderatepriors) #model has quite a large residual standard deviation (shape)--> need to identify why? This effects the cohens D estimates.
saveRDS(mod1.1_moderatepriors, file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_moderatepriors.rda")
mod1.1_moderatepriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_moderatepriors.rda")

#Posterior draws
color_scheme_set("pink")

mod1.1_weakpriors_cumvdis<-
  pp_check(mod1.1_weakpriors,
           ndraws = 100,
           type = "ecdf_overlay",
           discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 50)) +
  theme_classic()

### 7.4 default (flat) priors ---------------------------------------------------------
mod1.1_defaultprior<- brm(total_lifespan|cens(1-event)~
                            avg_age+ 
                            delta.age +
                            Temp + 
                            (1+delta.age|PairID),
                          family = weibull,
                          data = data1,
                          #No selected priors here, using the improper, default flat priors
                          iter=5000,
                          control=list(adapt_delta=0.98),
                          cores = 4)
summary(mod1.1_defaultprior)
saveRDS(mod1.1_defaultprior , 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_defaultprior.rda")
mod1.1_defaultprior<-readRDS("scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_defaultprior.rda")

# 7.5 Constraints priors -------------------------------------------------------
mod1.1_constraintpriors<-brm(total_lifespan|cens(1-event)~
                               avg_age+ 
                               delta.age +
                               Temp +
                               (1+delta.age|PairID),
                             family=weibull,
                             data=data1,
                             prior = constraintpriors,
                             iter=5000,
                             control=list(adapt_delta=0.98),
                             core=4)
summary(mod1.1_constraintpriors) 
saveRDS(mod1.1_constraintpriors, 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_constraintpriors.rda")
mod1.1_constraintpriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_constraintpriors.rda")

#posterior draws
constraintsprior_ecdf<-
  pp_check(mod1.1_weakpriors,
           ndraws = 100,
           type = "ecdf_overlay",
           discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 50)) +
  theme_classic()


## 7.6 Prior sensitivity: LOO comparison across prior specifications ---------
Loo_prior_performance<-LOO(mod1.1_weakpriors, mod1.1_diffusepriors, 
                           mod1.1_moderatepriors, mod1.1_constraintpriors, 
                           moment.match = TRUE)
saveRDS(Loo_prior_performance, 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/priorfitsummary.rda")
#Priors have almost no say over the posterior...
Loo_prior_performance<-readRDS("scripts/model_outputs/Offspring Trait Models/total lifespan/priorfitsummary.rda")

#selecting the model with weakly informative priors


# ===========================================================================
# 8. SELECTED MODEL DIAGNOSTICS  (mod1.1_weakpriors)--------------------------
# ===========================================================================

## LOO-CV --------------------------------------------------------------------
loo_weakpriors <- loo(mod1.1_weakpriors, save_psis = TRUE)
plot(loo_weakpriors)
psis_weakpriors <- loo_weakpriors$psis_object
psis_weights    <- weights(psis_weakpriors)

## Posterior predictive checks -----------------------------------------------
#posterior simulated and empirical mean
weakprior_mean     <- ppc_stat(yrep = posterior_predict(mod1.1_weakpriors),
                               y    = data1$total_lifespan,
                               stat = "mean") + theme_classic()
#posterior cumulative density
weakprior_ecdf     <- pp_check(mod1.1_weakpriors, ndraws = 100,
                               type = "ecdf_overlay") +
  xlim(0, 40) + theme_classic()

#posterior probability density
weakprior_pdf      <- pp_check(mod1.1_weakpriors, ndraws = 100) +
  xlim(0, 40) + theme_classic()

#posterior LOO Q-Q plot
weakprior_loo_qq   <- pp_check(mod1.1_weakpriors, type = "loo_pit_qq",
                               ndraws = 100) + theme_classic()
#LOO-PIT values fall well along the diagonal line (i.e. the expected uniform distribution of PIT integrals)

#LOO uniformity plot
weakprior_loo_unif <- pp_check(mod1.1_weakpriors, type = "loo_pit_overlay",
                               ndraws = 100) + theme_classic()
#LOO predictive interval plots
weakprior_intervals <- ppc_loo_intervals(
  y    = data1$total_lifespan,
  yrep = posterior_predict(mod1.1_weakpriors),
  psis_weakpriors
) + theme_classic()

## Combined fit plot ---------------------------------------------------------
fit_plots <- ggarrange(
  weakprior_mean, weakprior_ecdf, weakprior_pdf,
  weakprior_loo_qq, weakprior_loo_unif, weakprior_intervals,
  nrow = 2, ncol = 3,
  labels = c("A", "B", "C", "D", "E", "F")
)

ggsave("./bayesian_plots/model fit plots/total lifespan/selectedmodelfit.png",
       plot   = fit_plots,
       device = "png",
       width  = 500, height = 400, units = "mm")


## Bayes R² ------------------------------------------------------------------
bayes_R2(mod1.1_weakpriors, re.form = NA)   # fixed effects only (i.e., marginal)
bayes_R2(mod1.1_weakpriors, re.form = NULL)  # including random effects (i.e., conditional)


## MAP estimates and pd ------------------------------------------------------
MAP_totallifespan<- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  test        = "pd",
  ci_method   = "HDI",
  centrality  = "MAP",
  component   = "all",
  effects     = "full",
  ci          = 0.95
)
saveRDS(MAP_totallifespan,
        "scripts/model_outputs/Offspring Trait Models/total lifespan/MAP_totallifespan.rda")


## ROPE -----------------------------------------------------------------------
rope_totallifespan<- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  rope_range = c(-0.026 , 0.026), #(0.1*SD of the reponse)
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1 #estimating ROPE seperately as it requires the full posterior distribution
)
saveRDS(rope_totallifespan,
        "scripts/model_outputs/Offspring Trait Models/total lifespan/rope_totallifespan.rda")


# ===========================================================================
# 9. PRIOR vs. POSTERIOR SENSITIVITY PLOTS------------------------------------
# ===========================================================================

## Combined prior–posterior ECDF panel -----------------------------------------

#Saving the prior and posterior plots for default versus moderate priors
priorpostplots<-ggarrange(diffuseprior_cumvdis, weakprior_cumvdis, 
                          moderateprior_cumvdis, constraintprior_cumvdis,
                          mod1.1_diffusepriors_cumvdis, mod1.1_weakpriors_cumvdis,  
                          mod1.1_moderatepriors_cumvdis, constraintsprior_ecdf, nrow = 2, ncol = 4,
                          labels = c("A", "B", "c", "D", "E", "F", "G", "H"))

install.packages("ragg")
library(ragg)

ggsave(filename = "./bayesian_plots/model fit plots/total lifespan/priorvpost_ecdf.png",
       plot = priorpostplots, 
       device = png,
       dpi = 300,
       width = 490, 
       height = 340, 
       units = "mm")



# ===========================================================================
# 10. HYPOTHESIS TESTING: INTERACTION MODELS---------------------------------
# ===========================================================================
#Tests the core interaction asked by our paper

#10.1. Two-way interaction: delta age × temperature (mu only) ------------------------
#I.e., are parental age effects on offspring total lifespan conditional on the parents temperature treatment?
mod2.1<-brm(total_lifespan|cens(1-event)~
              avg_age+
              delta.age+
              Temp+
              delta.age:Temp+
              (1+delta.age|PairID),
            family=weibull(),
            data=data1,
            prior = weakpriors,
            iter=5000,
            control=list(adapt_delta=0.98),
            core=4)
summary(mod2.1)
saveRDS(mod2.1, file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod2.1.rda")
mod2.1<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod2.1.rda")

#diagnostics
loo_interaction<-loo(mod2.1, moment_match =TRUE, save_psis = TRUE)
plot(loo_interaction)

## MAP estimates for interaction model------------------------------
#For reporting in the main manuscript
MAP_totallifespan_mod2.1 <- bayestestR::describe_posterior(
  mod2.1,
  test       = "pd",
  ci_method  = "HDI",
  centrality = "MAP",
  effects    = "full",
  component = "all",
  ci         = 0.95
)
saveRDS(MAP_totallifespan_mod2.1,
        "scripts/model_outputs/Offspring Trait Models/total lifespan/MAP_totallifespan_mod2.1")

#ROPE estimates for interaction
rope_totallifespan_mod2.1 <- bayestestR::describe_posterior(
  mod2.1,
  rope_range = c(-0.026 , 0.026), #(0.1*SD of the reponse)
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1
)
saveRDS(rope_totallifespan_mod2.1,
        "scripts/model_outputs/Offspring Trait Models/total lifespan/rope_total_lifespan_mod2.1.rda")

# ===========================================================================
# 11. MODEL SELECTION----------------------------------------------------------
# ===========================================================================
hypothesis_fit<-loo(mod1.1_weakpriors,
                    mod2.1, moment_match = TRUE) 
saveRDS(hypothesis_fit, file = "scripts/model_outputs/Offspring Trait Models/total lifespan/hypothesis_fit.rda")
#Selecting model 1: No interactions --> Just single effect predictors

loocomparisons<-loo_compare(looweak,
                            loo_hypothesis2)
saveRDS(loocomparisons, file = "scripts/model_outputs/Offspring Trait Models/total lifespan/loocomparisons.rda")


## Selective disappearance test -----------------------------------------------
hypothesis(mod1.1_weakpriors, "avg_age - delta.age > 0")
hypothesis(mod1.1_weakpriors, "avg_age - delta.age < 0")

## Marginal means via emmeans ------------------------------------------------
#Effect of parental age
pairwise_estimates<- emmeans(mod1.1_weakpriors, "delta.age", type = "response", at = list(delta.age = c(-1.96, 1.96)))
pairs(pairwise_estimates) 
#for our post-hoc inference

# ===========================================================================
# 12. POSTERIOR DISTRIBUTION PLOTS  (mod1.1_weakpriors)-----------------------
# ===========================================================================
# Plots styled using stat_halfeye (ggdist) for consistency with other traits

post1 <- as_draws_df(mod1.1_weakpriors)
post2<-as_draws_df(mod3.1)

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
posterior_plot_totalife <- ggplot(
  posterior_df_mod1.1,
  aes(x = value, y = parameter, fill = parameter)
) +
  annotate("rect",
           xmin = -0.026, xmax = 0.026,
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
  xlim(c(-0.1, 0.1)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title      = element_text(size = 30),
    axis.text       = element_text(size = 30)
  ) +
  labs(x = "Fixed-effect posterior estimates\n(μ = mean total lifespan (weeks), log-scale)",
       y = NULL)

ggsave("./bayesian_plots/posterior plots/total lifespan/001_hypothesis1_halfeye.png",
       plot   = posterior_plot_totalife,
       device = "png",
       width  = 580, height = 400, units = "mm")


## Post-hoc standardised mean difference (SMD) --------------------------------
# Compares extreme ends of within-individual parental age; standardised by
# observed SD of F1 total lifespan. Descriptive only — not used for inference.

#Post-hoc standardised mean difference (on unit scale)
pairwise_estimates<- emmeans(mod1.1_weakpriors, "delta.age", type = "response", at = list(delta.age = c(-1.96, 1.96)))
#This is now converted to the raw scale where the extreme ends of parental age are compared
sd_life<-sd(data1$total_lifespan)

Total_lifespan_SMD<- summary(pairwise_age) %>%
  summarise(
    SMD       = (emmean[2]    - emmean[1])    / sd_life,
    SMD_lower = (lower.HPD[2] - upper.HPD[1]) / sd_life,
    SMD_upper = (upper.HPD[2] - lower.HPD[1]) / sd_life
  )

#Saving for use in combined forest plot of all traits
saveRDS(Total_lifespan_SMD, file = "scripts/model_outputs/Offspring Trait Models/Total_lifespan_SMD")

# ===========================================================================
# 13. MODEL PREDICTIONS AND RESULTS PLOTS-------------------------------------
# ===========================================================================

## 13.1 Prediction grids: within- and between-individual age effects ---------

# Within-individual: vary delta.age, hold avg_age constant
df_predict_within<-expand.grid(
  avg_age=mean(data1$avg_age),
  delta.age=unique(data1$delta.age),
  Temp=unique(data1$Temp),
  PairID=unique(data1$PairID[1]))

#calculating model predictions and extracting standard errors
pred_within <- fitted(mod1.1_weakpriors, df_predict_within, re_formula=NA) #NA generates estimates for "new" group levels (so not over your individuals, that is set by NULL)
pred_within<-as.data.frame(pred_within)
df_predict_within$prediction <- pred_within$Estimate
df_predict_within$lower <- pred_within$Q2.5
df_predict_within$upper<-pred_within$Q97.5

# Calculate mean predictions and SE for plotting
newdat <- df_predict_within %>%
  group_by(delta.age) %>%
  summarise(
    mean = mean(prediction),
    lower= mean(lower),
    upper=mean(upper)) %>% 
  ungroup() %>% 
  mutate(within_subject_age = delta.age*sd(data1$within_subject_age)) #converts standardised scale back to weeks (mean-centred)

# Between-individual: vary avg_age, hold delta.age constant
df_predict_between<-expand.grid(
  avg_age= unique(data1$avg.age),  
  delta.age=mean(data1$delta.age),
  Temp=unique(data1$Temp),
  PairID=unique(data1$PairID[1]))

#calculating model predictions and extracting standard errors
pred_between <- fitted(mod1.1_weakpriors, df_predict_between, re_formula=NA)
pred_between<-as.data.frame(pred_between)
df_predict_between$prediction <- pred_between$Estimate
df_predict_between$lower <- pred_between$Q2.5
df_predict_between$upper<-pred_between$Q97.5

# generating average age-slopes
avg.age_slopes<- df_predict_between %>%
  group_by(avg_age) %>% 
  summarise(prediction= mean(prediction),
            lower= mean(lower),
            upper=mean(upper)) %>% 
  mutate(
    avg.age2 = avg_age * sd(data1$avg.age)+mean(data1$avg.age),
    avg_timepoint_centered = avg.age2 - mean(avg.age2)
  )


## 13.2 Raw within-individual means (for overlaying on plot) ----------------
breaks <- seq(-4.5, 4.5, by = 1)  
labels <- -4:4  # label bins with week numbers

data1<- data1 %>%
  mutate(delta_age_bin = cut(within_subject_age, breaks = breaks, labels = labels))


#Raw delta age means
raw_deltaage <- data1%>%
  select(-within_subject_age) %>%        
  filter(!is.na(delta_age_bin)) %>% 
  rename(within_subject_age = delta_age_bin) %>% 
  mutate(within_subject_age = as.numeric(as.character(within_subject_age))) %>%  # convert factor to numeric
  group_by(within_subject_age) %>%
  summarise(
    n = sum(!is.na(total_lifespan)),                       # number of non-missing observations
    mean_life = mean(total_lifespan, na.rm = TRUE),        # group mean
    se_life = ifelse(n > 1, sd(total_lifespan, na.rm = TRUE)/sqrt(n), NA_real_), # SE only if n>1
    .groups = "drop"
  )

## 14.3 Main parental age effect plot ----------------------------------------
lifespan_plot <- ggplot(data = data1,
                        aes(x = within_subject_age, 
                            y = total_lifespan,
                            colour=Temp)) +
  geom_point(position = position_jitter(width = 0.2, height=0.1),
             shape=21,
             size=6,
             stroke=1.8,
             alpha = 0.7,
             colour="white",
             fill= "grey")+
  geom_line(data = newdat,
            aes(x =within_subject_age, 
                y = mean,
                linetype = "Within-individual",
                fill = "Within-individual",
                colour = "Within-individual"),
            linewidth=4) +
  geom_ribbon(data = newdat, 
              aes(y=NULL, 
                  ymin = lower, 
                  ymax = upper,
                  fill = "Within-individual",
                  colour="Within-individual"),
              alpha = 0.1,
              linewidth =1,
              show.legend = FALSE) +
  geom_line(data=avg.age_slopes,
            aes(x= avg_timepoint_centered,
                y=prediction,
                linetype = "Between-individual",
                fill = "Between-individual",
                colour = "Between-individual"),
            linewidth=4)+
  geom_ribbon(data = avg.age_slopes, 
              aes(x= avg_timepoint_centered,
                  y=NULL, 
                  ymin = lower, 
                  ymax = upper,
                  fill = "Between-individual",
                  colour="Between-individual"),
              alpha = 0.1,
              linewidth =1,
              show.legend = FALSE)+
  geom_linerange(data = raw_deltaage,
                 aes(y=mean_life,
                     ymin = mean_life-se_life,
                     ymax = mean_life+se_life,
                     colour="Within-individual"),
                 linewidth = 3,
                 show.legend = FALSE)+
  geom_point(data = raw_deltaage,
             aes(x = within_subject_age,
                 y = mean_life,
                 fill="Within-individual"),
             shape=21, 
             stroke=1.8,
             size= 18,
             alpha=1,
             colour="white",
             show.legend = FALSE)+
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
  )+
  scale_y_continuous(breaks=c(5, 10,15,20,25, 30))+
  scale_x_continuous(breaks=c(-4, -3, -2, -1, 0, 1,2,3, 4))+
  scale_fill_manual(name = "Parental age effect",
                    values = c("Within-individual" = "#4A6479", 
                               "Between-individual" = "#8C2F4B"))+
  scale_linetype_manual(name = "Parental age effect",
                        values = c("Within-individual" = "solid",
                                   "Between-individual" = "dashed"))+
  scale_colour_manual(name = "Parental age effect",
                      values = c("Within-individual" = "#4A6479", 
                                 "Between-individual" = "#8C2F4B"))+
  guides(linetype = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                                 override.aes = list(size = 15)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  theme(
    legend.position = "top", 
    legend.box = "horizontal"
  ) +
  labs(x = "Parents' adult age at reproduction (weeks; mean-centred)", y = "Offspring's total lifespan (weeks)")


# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 total lifespan/001_lifespan_plot.png",
       plot = lifespan_plot,
       bg="transparent",
       device = "png", 
       width = 500, 
       height = 700, 
       units = "mm")


## 13.4 Plotting effects with survival and hazard curves-----------------

#Code adapted from: https://github.com/Agasax/brms_survival_plot

#Functions to transform u to the scale parameter--------------------------------
median_weibull <- function(scale, shape) {
  (log(2)^(1/shape)) * scale
}


weibull_mu_to_scale <- function(mu, shape) {
  mu/gamma(1 + 1/shape)
}

weibull_survival <- function(scale, shape, time) {
  exp(-(time/scale)^shape)
}


weibull_hazard <- function(scale, shape, time) {
  (shape/scale) * (time/scale)^(shape - 1)
}


#Creating a dummy avriable for delta age that is on the SD scale

#Used to generate model predictions
newdat<-expand.grid(
  avg_age=mean(data1$avg_age),
  delta.age=c(-1.95, 0, 1.85),
  Temp=unique(data1$Temp),
  PairID=unique(data1$PairID))

# Creating a sequence over which to evaluate the survival and hazard
life_seq <- seq(min(0), 
                max(data1$total_lifespan, na.rm=TRUE), 
                length.out = 100)



#Estimating survival and hazard estimates over the dummy dataset
proportional_posterior <- mod1.1_weakpriors %>% 
  linpred_draws(
    newdat,
    value = "mu",
    allow_new_levels=TRUE,
    transform = TRUE,
    re_formula=NA,
    dpar = "shape",
    ndraws = 3000, 
    seed = 123
  ) 

#survival and hazard estimates
dempars<-proportional_posterior%>% 
  mutate(scale = weibull_mu_to_scale(mu, shape)) %>%
  ungroup() %>%
  crossing(time = life_seq) %>%
  mutate(
    S = weibull_survival(scale, shape, time),
    h = weibull_hazard(scale, shape, time)
  ) %>%
  group_by(time, .draw, delta.age) %>%
  summarise(S = mean(S), h = mean(h), .groups = "drop")

#Getting the median survival times
median_survival <-  mod1.1_weakpriors%>% 
  linpred_draws(
    newdat,
    allow_new_levels=TRUE,
    value = "mu",
    re_formula = NA,
    transform = TRUE,
    dpar = "shape", 
    ndraws = 3000 , 
    seed = 123
  ) %>% 
  mutate(scale = weibull_mu_to_scale(mu,shape)) %>% 
  group_by(.draw, delta.age) %>% 
  summarise(median_survival = mean(median_weibull(scale, shape))) %>% 
  group_by(delta.age) %>% 
  median_hdi(median_survival)

#summarising across draws for the survival and hazard curves
summary_survival <- dempars%>%
  group_by(time, delta.age) %>%
  median_hdi(S, .width = 0.95)

#Summarising to get the hazard
summary_hazard <- dempars %>%
  group_by(time, delta.age) %>%
  median_hdi(h, .width = 0.95)

#converting delta.age to a categorical variable for plotting
summary_survival<- summary_survival%>%
  mutate(timepoint= case_when(
    delta.age %in% -1.95 ~ "Early-Aged",
    delta.age %in% 0 ~ "Middle-Aged",
    delta.age %in% 1.85 ~ "Late-Aged"
  ))

summary_survival$timepoint<-factor(summary_survival$timepoint, 
                                   levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))
#converting delta.age
median_survival<- median_survival%>%
  mutate(timepoint= case_when(
    delta.age %in% -1.95 ~ "Early-Aged",
    delta.age %in% 0 ~ "Middle-Aged",
    delta.age %in% 1.85 ~ "Late-Aged"
  ))

#For the hazrad curves
#converting delta.age to a categorical variable for plotting
summary_hazard<- summary_hazard%>%
  mutate(timepoint= case_when(
    delta.age %in% -1.95 ~ "Early-Aged",
    delta.age %in% 0 ~ "Middle-Aged",
    delta.age %in% 1.85 ~ "Late-Aged"
  ))
summary_hazard$timepoint<-factor(summary_hazard$timepoint, 
                                 levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))



#Plotting the survival curve
survivalposterior <-ggplot(data=summary_survival,
                           aes(x = time, 
                               y = S, 
                               colour= timepoint)) + 
  geom_line(data = summary_survival, 
            aes(x = time, 
                y = S,
                colour = timepoint),
            linewidth=4, 
            alpha=1)+
  geom_ribbon(data = summary_survival, 
              aes(y=NULL, 
                  ymin = .lower, 
                  ymax = .upper, 
                  color= timepoint,
                  fill=timepoint,
                  alpha = timepoint),
              linetype="dashed",
              linewidth = 1)+
  geom_segment(data = median_survival,
               aes(x = 0, xend = median_survival, y = 0.5, yend = 0.5, color = timepoint),
               linetype = "dashed", linewidth = 2, alpha = 0.9) +
  geom_segment(data = median_survival,
               aes(x = median_survival, xend = median_survival, y = 0, yend = 0.5, color = timepoint),
               linetype = "dashed", linewidth = 2, alpha = 0.9)+
  theme_classic()+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), )+
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
  scale_color_manual(name = "Parents' adult age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' adult age at reproduction", values = c(0.6,0.1, 0.6))+
  scale_fill_manual(name = "Parents' adult age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  labs(x = "Offspring age (weeks)", y = "Cumulative survival probability, S(x)")

# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 total lifespan/001_totallifespan_survivalplot.png",
       plot = survivalposterior, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 250, 
       units = "mm")

#Plotting the hazard curves for each group --> this model assumes the shape is the same among groups...
hazardposterior <-ggplot(data=summary_hazard,
                         aes(x = time, 
                             y = h, 
                             colour= timepoint)) + 
  geom_line(data = summary_hazard, 
            aes(x = time, 
                y = h,
                colour = timepoint),
            linewidth=4, 
            alpha=1)+
  geom_ribbon(data = summary_hazard, 
              aes(y=NULL, 
                  ymin = .lower, 
                  ymax = .upper, 
                  color= timepoint,
                  fill=timepoint,
                  alpha = timepoint),
              linetype="dashed",
              linewidth = 1)+
  theme_classic()+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), )+
  ylim(c(0, 2.0))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
  scale_color_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' age at reproduction", values = c(0.6,0.1, 0.6))+
  scale_fill_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  labs(x = "Offspring age (weeks)", y = "Instantaneous hazard rate, μ (X)")

# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 total lifespan/001_totallifespan_hazard.png",
       plot = hazardposterior, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 250, 
       units = "mm")


#------------------COMBINATION PLOT FOR PUBLICATION------------------------------
new_total_mort_plot<- hazardposterior+
  labs( y = "Hazard rate, μ(x)")+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "none")

new_total_survival_plot<-survivalposterior +
  labs( y = "Survival, S(x)")+
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top")


right_panel <- new_total_survival_plot/ new_total_mort_plot

lifespan_inference<-(lifespan_plot | right_panel) +
  plot_layout(widths = c(0.8, 1)) +
  theme(plot.tag = element_text(size = 50, face = "bold"),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), )

ggsave(filename = "./bayesian_plots/offspring trait plots/F1 total lifespan/001_inference_plots2.png",
       plot = lifespan_inference, 
       device = "png", 
       width = 920, 
       height = 630, 
       units = "mm")

#Creating a table to export------------------------------------------------------------------------------

# ── helper ────────────────────────────────────────────────────────────────────
#Model draws from selected models
base_model_draws<-as_draws_df(mod1.1_weakpriors)#draws from model with no interaction
interaction_draws<-as_draws_df(mod2.1) #draws from the model with the interaction

#estimating the difference between average age and delta age
differences<-data.frame(
  selectivedis = base_model_draws$b_avg_age - base_model_draws$b_delta.age
)

#function for generating estimates
summarise_param <- function(draws, param, section,
                            rope = FALSE, rope_range = c(-0.026, 0.026), 
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
  summarise_param(base_model_draws$sd_PairID__delta.age, "σ slope Δage [High]",    "Random effects (Location)", pd = FALSE),
  summarise_param(base_model_draws$cor_PairID__Intercept__delta.age, "r intercept ~ slope Δage",    "Random effects (Location)", pd = TRUE)
)

#---4. Shape parameter──────────────────────────────────────────────────
shape<- bind_rows(
  summarise_param(base_model_draws$shape, "Shape parameter (K)", "Distributional parameters (scale component)", pd= FALSE)
)

#---5. Conditional and marginal Bayes R² ------------------------------------------------------------------
marginal<-bayes_R2(mod1.1_weakpriors, re.form = NA, summary = FALSE)   # fixed effects only (i.e., marginal)
conditional<-bayes_R2(mod1.1_weakpriors, re.form = NULL, summary = FALSE)  # including random effects (i.e., conditional)

bayes<-bind_rows(
  summarise_param(marginal, "Marginal R²", "Bayes R²", pd= FALSE),
  summarise_param(conditional, "Conditional R²", "Bayes R²", pd= FALSE)
)


# ── 5. Combine and render ─────────────────────────────────────────────────────
bind_rows(fe, fe_interaction, re_sd_explore, shape, bayes) %>%
  gt(groupname_col = "Section") %>%
  cols_label(Parameter = "Parameter", MAP = "MAP",
             `95% HDI` = "95% HDI", `% in ROPE` = "% in ROPE", pd = "pd") %>%
  tab_header(
    title    = "Total lifespan model summary: Weibull model"
  ) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) %>%
  tab_footnote("ROPE = [−0.026, 0.026] on log scale; applied to location fixed effects only.",
               locations = cells_column_labels(`% in ROPE`)) %>%
  tab_footnote("pd = probability of direction.",
               locations = cells_column_labels(pd)) %>%
  opt_stylize(style = 1)%>%
  gtsave("./tables/totalifespan_updatedtable.docx")


# ===========================================================================
# 14. SUBSET ANALYSIS: INCLUDING OFFSPRING SEX (REMOVING UNSEXED ANIMALS)
# ===========================================================================

#Filtering so only sexed adults are included
data2<-data1 %>%
  filter(!Include_in_total == "N") %>% 
  filter(!is.na(F1_sex))
length(unique(data2$F1_ID))
#Analysis includes 947 sexed offspring
#Assuming the weakly-informative priors fit the data the best...

#14.1: No interaction with offspring sex----------------------------------------
mod1.1_offspring_sex<- brm(total_lifespan|cens(1-event)~
                             avg_age+ 
                             delta.age +
                             Temp + 
                             F1_sex+
                             (1+delta.age|PairID),
                           family=weibull,
                           data=data2,
                           #With specified priors added
                           prior = weakpriors,
                           iter=5000,
                           control=list(adapt_delta=0.98),
                           save_pars = save_pars(all = TRUE),
                           core=4)

summary(mod1.1_offspring_sex)
saveRDS(mod1.1_offspring_sex, 
        file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_offspring_sex.rda")
mod1.1_offspring_sex<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod1.1_offspring_sex.rda")

bayes_R2(mod1.1_offspring_sex, re.form = NA)#model only explains 3% of variation in total lifespan...
bayes_R2(mod1.1_offspring_sex, re.form = NULL) #conditional R2

#Getting model parameters and pd values for this subset model
fixed_offspring<- parameters::parameters(mod1.1_offspring_sex, effects = "fixed")
parameters::parameters(mod1.1_offspring_sex)

#Testing for selective disappearance within this model------
hypothesis(mod1.1_offspring_sex, "avg_age - delta.age > 0")
hypothesis(mod1.1_offspring_sex, "avg_age - delta.age < 0")

#obtaining SD for the ROPE inference
mean_y<-mean(data2$total_lifespan)
sd_y<-sd(log10(data2$total_lifespan))*2.302585 #=0.2585665 SD on the log scale
rope_response <- c(-0.1 * sd_y, 0.1 * sd_y)

#ROPE 
ROPETOTALLIFESPAN2<-bayestestR::describe_posterior(mod1.1_offspring_sex, 
                                                   test = c("p_direction", "rope"),
                                                   rope_range = c(-0.02 , 0.02),
                                                   ci = 1)
saveRDS(ROPETOTALLIFESPAN2, file = "scripts/model_outputs/Offspring Trait Models/total lifespan/ROPETOTALLIFESPAN2.rda")
ROPETOTALLIFESPAN2<- readRDS("scripts/model_outputs/Offspring Trait Models/total lifespan/ROPETOTALLIFESPAN2.rda")


#Model 14.2: Three-way interaction between offspring sex, parental age, and temperature treatment----
mod2.1_offspring_sex<- brm(total_lifespan|cens(1-event)~
                             avg_age+ 
                             delta.age +
                             Temp + 
                             F1_sex+
                             F1_sex*delta.age*Temp+
                             (1+delta.age|PairID),
                           family=weibull,
                           data=data2,
                           #With specified priors added
                           prior = weakpriors,
                           iter=5000,
                           control=list(adapt_delta=0.98,
                                        max_treedepth = 10),#sampler struggles at lower tree depths
                           save_pars = save_pars(all = TRUE),
                           core=4)
summary(mod2.1_offspring_sex)
plot(mod2.1_offspring_sex)
saveRDS(mod2.1_offspring_sex, file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod2.1_offspring_sex.rda")
mod2.1_offspring_sex<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod2.1_offspring_sex.rda")

bayes_R2(mod2.1_offspring_sex, re.form = NA)#model only explains 3% of variation in total lifespan...
bayes_R2(mod2.1_offspring_sex, re.form = NULL) #conditiuonal R2


#Model 14.3: Two-way interaction  between offspring sex and parental age----
mod3.1_offspring_sex<- brm(total_lifespan|cens(1-event)~
                             avg_age+ 
                             delta.age +
                             Temp + 
                             F1_sex+
                             F1_sex:delta.age+
                             (1+delta.age|PairID),
                           family=weibull,
                           data=data2,
                           #With specified priors added
                           prior = weakpriors,
                           iter=5000,
                           control=list(adapt_delta=0.98),
                           save_pars = save_pars(all = TRUE),
                           core=4)
summary(mod3.1_offspring_sex)
plot(mod3.1_offspring_sex)
saveRDS(mod3.1_offspring_sex, file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod3.1_offspring_sex.rda")
mod3.1_offspring_sex<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/total lifespan/mod3.1_offspring_sex.rda")
bayes_R2(mod3.1_offspring_sex, re.form = NA)#model only explains 3% of variation in total lifespan...
bayes_R2(mod3.1_offspring_sex, re.form = NULL) #conditional R2

fixed_offspring_interaction<- parameters::parameters(mod3.1_offspring_sex, effects = "fixed")

#Getting the ROPE estimates for parental age effects
mean_y<-mean(data2$total_lifespan)
sd_y<-sd(log10(data2$total_lifespan))*2.302585 #=0.2585665 SD on the log scale
rope_response <- c(-0.1 * sd_y, 0.1 * sd_y)
ROPETOTALLIFESPAN3<-bayestestR::describe_posterior(mod3.1_offspring_sex, 
                                                   test = c("p_direction", "rope"),
                                                   rope_range = c(-0.02 , 0.02),
                                                   ci = 1)


#Comparing the fit of the models --> there's not really any evidence for any interaction
LOOoffspringsex<-LOO(mod1.1_offspring_sex, mod2.1_offspring_sex,mod3.1_offspring_sex)
loooffspringsex<-loo(mod1.1_offspring_sex, reloo=TRUE)
loooffspringsex_fullinteraction<-loo(mod2.1_offspring_sex, reloo=TRUE)
loooffspringsex_agesexinteraction<-loo(mod3.1_offspringsex, reloo=TRUE)



#Creating a table to export------------------------------------------------------------------------------
#For the subset of animals that we know became adults

# ── helper ────────────────────────────────────────────────────────────────────
#Model draws from selected models
base_model_draws2<-as_draws_df(mod1.1_offspring_sex)#draws from model with no interaction
interaction_draws2<-as_draws_df(mod2.1_offspring_sex) #draws from the model with the three-way interaction
interaction_draws3<-as_draws_df(mod3.1_offspring_sex) #draws from the model with the two-way interaction

#estimating the difference between average age and delta age
differences2<-data.frame(
  selectivedis = base_model_draws2$b_avg_age - base_model_draws2$b_delta.age
)

#function for generating estimates
summarise_param <- function(draws, param, section,
                            rope = FALSE, rope_range = c(-0.026, 0.026), 
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
fe2<- bind_rows(
  summarise_param(base_model_draws2$b_Intercept,        "Intercept",            "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws2$b_avg_age,   "Parents' average age (z-scaled)", "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws2$b_delta.age, "Parents' Δage (z-scaled)",        "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws2$b_Temp28, "Temperature: 28.0°C",        "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws2$b_Temp30.5, "Temperature: 30.5°C",        "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(base_model_draws2$b_F1_sexM, "Offspring sex (Male)", "Location submodel (Fixed effects)", rope = TRUE),
  summarise_param(differences2$selectivedis, "Average age - Δage", "Location submodel (Fixed effects)", rope = FALSE)
  )

# ── 1. dropped interaction (location) — ROPE + pd ───────────────────────────────────
fe_interaction2 <- bind_rows(summarise_param(interaction_draws2$`b_delta.age:Temp28:F1_sexM`, "Parents' Δage x  Temperature: 28.0°C x Offspring sex (Male)","Location submodel (Dropped interactions)", rope = TRUE),
                             summarise_param(interaction_draws2$`b_delta.age:Temp30.5:F1_sexM`, "Parents' Δage x  Temperature: 30.5°C x Offspring sex (Male)","Location submodel (Dropped interactions)", rope = TRUE),
                            summarise_param(interaction_draws2$`b_delta.age:Temp28`, "Parents' Δage x  Temperature: 28.0°C", "Location submodel (Dropped interactions)", rope = TRUE),
                            summarise_param(interaction_draws2$`b_delta.age:Temp30.5`, "Parents' Δage x  Temperature: 30.5°C", "Location submodel (Dropped interactions)", rope = TRUE),
                            summarise_param(interaction_draws3$`b_delta.age:F1_sexM`, "Parents' Δage x Offspring sex: Male","Location submodel (Dropped interactions)", rope = TRUE)
)

# ── 3.Random effects ──────────────────────────────────────────────────
# SDs 
re_sd_explore2<- bind_rows(
  summarise_param(base_model_draws2$sd_PairID__Intercept, "σ intercept",    "Random effects (Location)", pd = FALSE),
  summarise_param(base_model_draws2$sd_PairID__delta.age, "σ slope Δage [High]",    "Random effects (Location)", pd = FALSE),
  summarise_param(base_model_draws2$cor_PairID__Intercept__delta.age, "r intercept ~ slope Δage",    "Random effects (Location)", pd = TRUE)
)

#---4. Shape parameter──────────────────────────────────────────────────
shape2<- bind_rows(
  summarise_param(base_model_draws2$shape, "Shape parameter (K)", "Distributional parameters (scale component)", pd= FALSE)
)

#---5. Conditional and marginal Bayes R² ------------------------------------------------------------------
marginal2<-bayes_R2(mod1.1_offspring_sex, re.form = NA, summary = FALSE)   # fixed effects only (i.e., marginal)
conditional2<-bayes_R2(mod1.1_offspring_sex, re.form = NULL, summary = FALSE)  # including random effects (i.e., conditional)

bayes2<-bind_rows(
  summarise_param(marginal, "Marginal R²", "Bayes R²", pd= FALSE),
  summarise_param(conditional, "Conditional R²", "Bayes R²", pd= FALSE)
)


# ── 5. Combine and render ─────────────────────────────────────────────────────
bind_rows(fe2, fe_interaction2, re_sd_explore2, shape2, bayes2) %>%
  gt(groupname_col = "Section") %>%
  cols_label(Parameter = "Parameter", MAP = "MAP",
             `95% HDI` = "95% HDI", `% in ROPE` = "% in ROPE", pd = "pd") %>%
  tab_header(
    title    = "Total lifespan model summary: Weibull model"
  ) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) %>%
  tab_footnote("ROPE = [−0.026, 0.026] on log scale; applied to location fixed effects only.",
               locations = cells_column_labels(`% in ROPE`)) %>%
  tab_footnote("pd = probability of direction.",
               locations = cells_column_labels(pd)) %>%
  opt_stylize(style = 1)%>%
  gtsave("./tables/totalifespan_withsex.docx")


###############################################################################
##END OF SCRIPT----------------------------------------------------------------
###############################################################################
