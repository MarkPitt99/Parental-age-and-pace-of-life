
###############################################################################
##  Script: Effect of parental age on F1 adult lifespan
##  Note:   Analysis restricted to offspring that survived  to adulthood
#Script provides an insight into average adult lifespan only - see BaSTA scripts for how adult mortality
#parameters are affected by increasing parental age.
###############################################################################

# ===========================================================================
# 1. SETUP--------------------------------------------------------------------
# ===========================================================================

## Libraries ------------------------------------------------------------------
library(dplyr)  #v.1.1.4
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

# ===========================================================================
# 2. DATA ORGANISATION--------------------------------------------------------
# ===========================================================================

#read data for animals that are to be included in this analysis#
F1data<-readRDS("./raw data/F1_filtered_data06022025.RDS")
length(unique(F1data$F1_ID))#1317 offspring 

#Only retaining animals that survived to adulthood
data1<-F1data%>%
  filter(!is.na(adult_lifespan)) %>% 
  filter(Include_in_adult == "Y")

#Creating an adult lifespan variable: Calculating the difference between adult emergence and death dates
data1$adult_lifespan <- as.numeric(difftime(data1$F1_death_date, data1$F1_Adult_emergence, units = "weeks"))

#Creating an "event" column to signal when the offspring died
data1$event<-ifelse(is.na(data1$adult_lifespan), 0,1)

#Scaling continous predictors----------------------------
#Scaled to have mean = 0, and SD = 1 (i.e., z-transformed)
data1$avg_age<-as.numeric(scale(data1$avg.age, center = TRUE, scale = TRUE))

# Within-individual (delta) age: scaled but NOT mean-centred to preserve
# the within-subject structure (since this varibale is already mean centered); NAs replaced with 0
data1 <- data1 %>%
  group_by(PairID) %>%
  mutate(
    delta.age = as.numeric(scale(within_subject_age,
                                 center = FALSE,
                                 scale = TRUE)),
    delta.age = replace(delta.age, is.na(delta.age), 0)
  ) %>%
  ungroup()
sd(data1$within_subject_age)
sd(data1$delta.age)


2.062967 # A single SD increase in delta age reflects a 2.06 week increase in parental age

## Sample sizes ---------------------------------------------------------------
length(unique(data1$Mother_ID)) #77 parent pairs
length(unique(data1$F1_ID)) #939 offspring that reached adulthood


# ===========================================================================
# 3. SUMMARY STATISTICS----------------------------------------------------------
# ===========================================================================

## Distribution of F1 adult lifespan------------------------------------------
ggplot(data1, aes(x = adult_lifespan)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") + 
  facet_wrap(~Timepoint) +  # Create panels for each Temp level
  labs(
    x = "F1 adult lifespan (Weeks)",
    y = "Count"
  ) +
  theme_classic()

## Overall mean and SD (used to build intercept priors later) ----------------------
sum_dat<-data1 %>% 
  summarise(F1_lifespan = mean(adult_lifespan),
            sd_total= sd(adult_lifespan),
            n_animals=n(),
            n_parents =n_distinct(PairID))
#mean lifespan of 9.53 weeks, SD of 3.16 weeks, n =939, no.parents = 77

#mean adult lifespan in relation to each timepoint----------------------------
sum_dat1<-data1 %>% 
  group_by(Timepoint)%>%
  summarise(F1_lifespan = median (adult_lifespan),
            se_total= sd(adult_lifespan)/sqrt(n()),
            n_animals=n())



## Collinearity checks --------------------------------------------------------
X <- model.matrix(
  ~ avg_age + delta.age + cum_succesful_matings +
    Temp + F1_sex,
  data = data.frame(
    avg_age = data1$avg_age,
    delta.age = data1$delta.age,
    Temp = data1$Temp,
    cum_succesful_matings = data1$cum_successful_matings,
    F1_sex = data1$F1_sex
  )
)
# Correlation matrix
cor(X[, -1])  # correlation matrix (excluding intercept column)


# Variance inflation factors (VIFs)
vif(lm(adult_lifespan ~ avg_age + delta.age + cum_successful_matings+
         Temp + F1_sex, data = data1))

#removing cum_succesful matings (i.e., cumulative number of succesful matings)
vif(lm(adult_lifespan ~ avg_age + delta.age +
         Temp + F1_sex, data = data1))


# ===========================================================================
# 4. MODEL STRUCTURE SELECTION------------------------------------------------
# ===========================================================================
#Data is survival data (discrete time-to-event)
#Using brms, the following distributions can be used to models survival times
#Gamma, Weibull, Exponential, lognormal and cox


#4.1.Default Weibull model---------------------------------------------------
data1<-data1 %>% 
  filter(!adult_lifespan < 0.28) #brms Weibull model can't handle zero data (this removes the single individual where lifespan = 0)

#Weibull brms model --> requires responses > 0, 
default_weibull<-brm(adult_lifespan|cens(1-event)~
                       avg_age+
                       delta.age+
                       Temp+
                       F1_sex+
                       (1|PairID),
                     family=weibull(),
                     data=data1 ,
                     iter=5000,
                     save_pars = save_pars(all = TRUE),
                     cores = 4)
summary(default_weibull) 
saveRDS(default_weibull , file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/defaultweibull.rda")
default_weibull<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/defaultweibull.rda")

#Model diagnostics
Loo_Weibull<-LOO(default_weibull, save_psis = TRUE, moment_match=TRUE)
plot(Loo_Weibull)

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(default_weibull, type="loo_pit_qq", ndraws=100) #Very good fit
pp_check(default_weibull, type ="loo_pit_overlay", ndraws=100) #Within the expected uniform distribution

#Inspecting the fit with a KM curve
default_weibull_KM<- pp_check(default_weibull,
                              status_y=data1$event,
                              type="km_overlay", ndraws=100) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  ylab("Cumulative Survival Probability, S(x)") + coord_cartesian(xlim=c(0, 20))+
  xlab("Adult age (weeks)") + coord_cartesian(xlim=c(0, 20))


#4.2. Default lognormal model------------------------------------------------
default_lognormal<-brm(adult_lifespan|cens(1-event)~
                         avg.age+
                         delta.age+
                         F1_sex+
                         mothers_age_at_entry_scaled+ 
                         Temp+
                         (1|PairID),
                       family=lognormal(),
                       data=data1,
                       iter=5000,
                       save_pars = save_pars(all = TRUE),
                       cores = 4)
summary(default_lognormal) 
saveRDS(default_lognormal , file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/defaultlognormal.rda")
default_lognormal<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/defaultlognormal.rda")

#Model diagnostics
loo_lognormal<- loo(default_lognormal, save_psis = TRUE, moment_match = TRUE)
plot(loo_lognormal)

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(default_lognormal, type="loo_pit_qq", ndraws=100) 
pp_check(default_lognormal, type ="loo_pit_overlay", ndraws=100) 

#Lognormal kaplan meier curve
default_lognormal_KM<- pp_check(default_lognormal,
                                status_y=data1$event,
                                type="km_overlay", ndraws=100) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  ylab("Cumulative Survival Probability, S(x)") + coord_cartesian(xlim=c(0, 20))+
  xlab("Adult age (weeks)") + coord_cartesian(xlim=c(0, 20))


#4.3. Default exponential model-----------------------------------
default_exponential<-brm(adult_lifespan|cens(1-event)~
                           avg_age+
                           delta.age+
                           Temp+
                           F1_sex+
                           (1|PairID),
                         family=exponential(),
                         data=data1,
                         iter=5000,
                         save_pars = save_pars(all = TRUE),
                         cores = 4)
summary(default_exponential) 
saveRDS(default_exponential, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/defaultexponential.rda")
default_exponential<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/defaultexponential.rda")

#model diagnostics
loo_exponential<- loo(default_exponential, save_psis = TRUE, moment_match = TRUE)
plot(loo_exponential)

#Creating the LOO-PIT plots
#Loo probability integral transform (PIT)
pp_check(default_exponential, type="loo_pit_qq", ndraws=100) #hugely underdispersed
pp_check(default_exponential, type ="loo_pit_overlay", ndraws=100) 


#exponential kaplan meier curve
default_exponential_KM<- pp_check(default_exponential,
                                  status_y=data1$event,
                                  type="km_overlay", ndraws=100) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  ylab("Cumulative Survival Probability, S(x)") + coord_cartesian(xlim=c(0, 35))+
  xlab("Offspring age (weeks)") + coord_cartesian(xlim=c(0, 35))


#--------------------Comparison of the time to event models------------------------------------------
LOO_allsurvivalmodels<-loo(default_weibull, default_lognormal, default_exponential, moment_match = TRUE)
saveRDS(LOO_allsurvivalmodels, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/LOO_allsurvivalmodels.rda")
LOO_allsurvivalmodels<-readRDS(LOO_allsurvivalmodels, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/LOO_allsurvivalmodels.rda")
#Choosing the Weibull model as the best fitting

#Deciding on whether the effect of delta age should be a linear or a quadratic term-----
#4.4. Adding the quadratic age term into the model
mod1.1_quadraticage<-brm(adult_lifespan|cens(1-event)~
                           avg_age+
                           poly(delta.age,2)+ 
                           Temp+
                           F1_sex+
                           (1|PairID),
                         family=weibull,
                         data=data1,
                         iter=5000,
                         control=list(adapt_delta=0.98),
                         save_pars = save_pars(all = TRUE),
                         core=4)
summary(mod1.1_quadraticage) 
conditional_effects(mod1.1_quadraticage, effects = "delta.age") #Not much support for a quadratic age effect
saveRDS(mod1.1_quadraticage, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_quadraticage.rda")
mod1.1_quadraticage<-readRDS("scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_quadraticage.rda")

#MODEL DIAGNOSTICS
loo_quadratic<-loo(mod1.1_quadraticage,  save_psis = TRUE, moment_match = TRUE)
plot(loo_quadratic)

#Loo probability integral transform (PIT)
pp_check(mod1.1_quadraticage, type="loo_pit_qq", ndraws=100)
pp_check(mod1.1_quadraticage, type ="loo_pit_overlay", ndraws=100) 


#4.4. Quadratic age +Random slopes model------------------
#Model does not allow the shape of the slopes of individual parents to vary
mod1.1_quadraticslopes<-brm(adult_lifespan|cens(1-event)~
                              avg_age+ 
                              poly(delta.age, 2)+ 
                              Temp + 
                              F1_sex+
                              (1+delta.age|PairID),
                            family=weibull,
                            data=data1,
                            iter=5000,
                            save_pars = save_pars(all = TRUE),
                            control=list(adapt_delta=0.95),
                            core=4)
summary(mod1.1_quadraticslopes)
saveRDS(mod1.1_quadraticslopes, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_quadraticslopes.rda")
mod1.1_quadraticslopes<-readRDS("scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_quadraticslopes.rda")

#diagnostics
loo_quadratic_slopes<-loo(mod1.1_quadraticslopes, moment_match = TRUE, save_psis = TRUE)
plot(loo_quadratic_slopes)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_quadraticslopes, type="loo_pit_qq", ndraws=100)
pp_check(mod1.1_quadraticslopes, type ="loo_pit_overlay", ndraws=100) 

#4.5. Linear age + Random slopes model----------------------------
mod1.1_randomslopes<-brm(adult_lifespan|cens(1-event)~
                           avg_age+ 
                           delta.age + 
                           Temp + 
                           F1_sex+
                           (1+delta.age|PairID),
                         family=weibull,
                         data=data1,
                         iter=5000,
                         save_pars = save_pars(all = TRUE),
                         control=list(adapt_delta=0.95),
                         core=4)
summary(mod1.1_randomslopes)
saveRDS(mod1.1_randomslopes , file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_randomslopes.rda")
mod1.1_randomslopes<-readRDS("scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_randomslopes.rda")

#diagnostics
loo_randomslopes<-loo(mod1.1_randomslopes, save_psis=TRUE, moment_match = TRUE)
plot(loo_randomslopes)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_randomslopes, type="loo_pit_qq", ndraws=100)
pp_check(mod1.1_randomslopes, type ="loo_pit_overlay", ndraws=100) 


#--------------Comparing the fit of the random slopes and quadratic age terms---
LOO_ageeffects<-loo(mod1.1_randomslopes, mod1.1_quadraticage, default_weibull, mod1.1_quadraticslopes, moment_match = TRUE)
saveRDS(LOO_ageeffects, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/LOO_ageeffects.rda")
LOO_ageeffects<-readRDS("scripts/model_outputs/Offspring Trait Models/adult lifespan/LOO_ageeffects.rda")

## 4.6 Distributional models: evidence for shape(k) submodel ------------------
#Testing for a temperature effect + age effect on the shape parameter
mod1.1_shape_temp<-brm(bf(adult_lifespan|cens(1-event)~
                            avg_age+ 
                            delta.age +
                            Temp + 
                            F1_sex+
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
saveRDS(mod1.1_shape_temp , file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_shape_temp.rda")
mod1.1_shape_temp<-readRDS("scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_shape_temp.rda")

#diagnostics
loo_shapetemp<-loo(mod1.1_shape_temp, save_psis = TRUE, moment_match = TRUE)
plot(loo_shapetemp)


#4.7. effect of just parental age on the shape parameter--------------------------
mod1.1_shape<-brm(bf(adult_lifespan|cens(1-event)~
                       avg_age+ 
                       delta.age +
                       Temp + 
                       F1_sex+
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
saveRDS(mod1.1_shape , file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_shape.rda")
mod1.1_shape<-readRDS("scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_shape.rda")

#MODEL DIAGNOSTICS 
loo_shape<-loo(mod1.1_shape,  save_psis = TRUE, moment_match=TRUE)
plot(mod1.1_shape) 

#---------------------------Model setup-----------------------------------------
#---------Both of the models with covariates on the shape parameter fail to converge-----
##------------Now assume a constant shape parameter across fixed effects--------------------
#This is unlikely to be the case, hence we also used our BaSTA framework


# ===========================================================================
# 5. PRIOR SPECIFICATION-----------------------------------------------------
# ===========================================================================
# Priors assume all continuous predictors are mean-centered and SD-scaled.
# Intercept prior informed by observed mean (9.54 weeks) and SD (3.15 weeks).

#---------------------Setting up the priors-------------------------------------
#what's the range of our data?
range(log(data1$adult_lifespan)) 
mean(data1$adult_lifespan) #9.54 weeks
sd(data1$adult_lifespan) #3.145468 weeks

#SD on the empirical log scale
sd(log10(data1$adult_lifespan))*2.302585 #=0.42

#mean on log scale
log(9.54) #2.255493 weeks on the log scale

#SD on log scale
sqrt(log(1 + ( 3.145468^2 / 9.54^2))) #0.32 on log scale

#Setting a prior on shape --> maybe a shape centered on 1.5 (assumes an increase in hazard)
#this could also have a broad SD, permitting the hazard curve to be flat (shape = 1) or accelerating increases (shape >2)

#---------------------SETTING MODEL PRIORS--------------------------------------

#Diffuse priors---------------------------------
Diffusepriors <- c(
  prior(normal(2.26, 0.32), class = "Intercept"),             
  prior(normal(0,1), class = "b"),          
  prior(normal(0, 2.5), class ="shape", lb = 0),
  prior(student_t(3, 0, 2.5), class = "sd", lb = 0),           
  prior(lkj(2), class = "cor")                                 
)


#Weakly informative priors [SELECTED]------------
weakpriors <- c(
  prior(normal(2.26,0.32), class = "Intercept"),                
  prior(normal(0,0.32), class = "b"), 
  prior(normal(1, 2.5), class = "shape", lb=0), #brms model struggles to add parameters to the shape of the hazard         
  prior(exponential(10), class = "sd", lb = 0),                      
  prior(lkj(2), class = "cor")                                   
)


#Moderate prior----------------------------------
moderatepriors <- c(
  prior(normal(2.26, 0.32), class = "Intercept"),                 
  prior(normal(0,0.16), class = "b"),  
  prior(normal(1,1.5), class="shape", lb=0),  
  prior(exponential(10), class = "sd", lb = 0),                         
  prior(lkj(2), class = "cor")                                    
)

# ===========================================================================
# 6. PRIOR PREDICTIVE CHECKS---------------------------------------------------
# ===========================================================================

## 6.1 Prior-only model: diffuse ---------------------------------------------------
diffuse_prior_model<-brm(adult_lifespan|cens(1-event)~
                           avg_age+ 
                           delta.age +
                           Temp + 
                           F1_sex+
                           (1+delta.age|PairID),
                         family=weibull,
                         data=data1,
                         iter=5000,
                         #With specified priors added
                         prior = Diffusepriors,
                         sample_prior = "only",
                         cores = 4)
summary(diffuse_prior_model)
saveRDS(diffuse_prior_model , file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/diffusepriormodel.rda")
diffuse_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/diffusepriormodel.rda")

#setting the colour scheme
color_scheme_set("teal")

#Prior draws (I.e., the prior cumulative density function) 
diffuseprior_cumvdis <- pp_check(diffuse_prior_model,
                                 ndraws = 100,
                                 type = "ecdf_overlay",
                                 discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 40)) +
  theme_classic()


##6.2. weak prior draws--------------------------------------------------------------
weak_prior_model <-brm(bf(adult_lifespan|cens(1-event)~
                            avg_age+ 
                            delta.age +
                            Temp + 
                            F1_sex+
                            (1+delta.age|PairID)),
                       family=weibull,
                       data=data1,
                       iter=5000,
                       #With specified priors added
                       prior = weakpriors,
                       sample_prior = "only",
                       cores = 4)
summary(weak_prior_model)
saveRDS(weak_prior_model , file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/weakpriors.rda")
weak_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/weakpriors.rda")

#Prior draws (I.e., the prior cumulative density function) 
weakprior_cumvdis <- pp_check(weak_prior_model,
                              ndraws = 100,
                              type = "ecdf_overlay",
                              discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 40)) +
  theme_classic()

##6.3. moderate prior draws-----------------------------------------------------------------
moderate_prior_model <-brm(adult_lifespan|cens(1-event)~
                             avg_age+ 
                             delta.age +
                             Temp + 
                             F1_sex+
                             (1+delta.age|PairID),
                           family=weibull,
                           data=data1,
                           iter=5000,
                           #With specified priors added
                           prior = moderatepriors,
                           sample_prior = "only",
                           cores = 4)
summary(moderate_prior_model)
saveRDS(moderate_prior_model , file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/moderatepriors.rda")
moderate_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/moderatepriors.rda")

#Prior draws (I.e., the prior cumulative density function) 
moderateprior_cumvdis <- pp_check(moderate_prior_model,
                                  ndraws = 100,
                                  type = "ecdf_overlay",
                                  discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 40)) +
  theme_classic()


# ===========================================================================
# 7. POSTERIOR MODELS---------------------------------------------------------
# ===========================================================================

## 7.1 Diffuse priors ---------------------------------------------------------
mod1.1_diffusepriors<-brm(adult_lifespan|cens(1-event)~
                            avg_age+ 
                            delta.age +
                            Temp +
                            F1_sex+
                            (1+delta.age|PairID),
                          family=weibull,
                          data=data1,
                          #With specified priors added
                          prior = Diffusepriors,
                          iter=5000,
                          control=list(adapt_delta=0.98),
                          save_pars = save_pars(all = TRUE),
                          core=4)
summary(mod1.1_diffusepriors) 
saveRDS(mod1.1_diffusepriors, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_diffusepriors.rda")
mod1.1_diffusepriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_diffusepriors.rda")

color_scheme_set("pink")

#Posterior draws 
mod1.1_diffusepriors_cumvdis<-
  pp_check(mod1.1_diffusepriors,
           ndraws = 100,
           type = "ecdf_overlay",
           discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 40)) +
  theme_classic()

#7.2. Weak priors [SELECTED MODEL]---------------------------------------------------------------
mod1.1_weakpriors<-brm(adult_lifespan|cens(1-event)~
                         avg_age+ 
                         delta.age +
                         Temp + 
                         F1_sex+
                         (1+delta.age|PairID),
                       family=weibull,
                       data=data1,
                       #With specified priors added
                       prior = weakpriors,
                       iter=5000,
                       control=list(adapt_delta=0.98),
                       save_pars = save_pars(all = TRUE),
                       core=4)
summary(mod1.1_weakpriors)
saveRDS(mod1.1_weakpriors, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_weakpriors.rda")
mod1.1_weakpriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_weakpriors.rda")

#Posterior draws 
color_scheme_set("pink")
mod1.1_weakpriors_cumvdis<-
  pp_check(mod1.1_weakpriors,
           ndraws = 100,
           type = "ecdf_overlay",
           discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 40)) +
  theme_classic()

#7.3. Moderate priors---------------------------------------------------------------
mod1.1_moderatepriors<-brm(adult_lifespan|cens(1-event)~
                             avg_age+ 
                             delta.age +
                             Temp + 
                             F1_sex+
                             (1+delta.age|PairID),
                           family=weibull,
                           data=data1,
                           #With specified priors added
                           prior = moderatepriors,
                           iter=5000,
                           control=list(adapt_delta=0.98),
                           core=4)
summary(mod1.1_moderatepriors)
saveRDS(mod1.1_moderatepriors, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_moderatepriors.rda")
mod1.1_moderatepriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_moderatepriors.rda")

#Posterior draws
color_scheme_set("pink")
mod1.1_moderatepriors_cumvdis<-
  pp_check(mod1.1_moderatepriors,
           ndraws = 100,
           type = "ecdf_overlay",
           discrete = FALSE) +
  scale_y_reverse(breaks=0.5) +
  coord_cartesian(xlim = c(0, 40)) +
  theme_classic()

#7.4. default prior model-------------------------------
#default brms priors
mod1.1_defaultprior<- brm(adult_lifespan|cens(1-event)~
                            avg_age+ 
                            delta.age +
                            Temp + 
                            F1_sex+
                            (1+delta.age|PairID),
                          family = weibull,
                          data = data1,
                          #No selected priors here, using the improper, default flat priors
                          iter=5000,
                          control=list(adapt_delta=0.98),
                          cores = 4)
summary(mod1.1_defaultprior)
saveRDS(mod1.1_defaultprior , file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_defaultprior.rda")
mod1.1_defaultprior<-readRDS("scripts/model_outputs/Offspring Trait Models/adult lifespan/mod1.1_defaultprior.rda")


## 7.5 Prior sensitivity: LOO comparison across prior specifications ---------
Loo_prior_performance<-LOO(mod1.1_weakpriors, mod1.1_diffusepriors, 
                           mod1.1_moderatepriors, mod1.1_defaultprior,
                           moment.match = TRUE)
saveRDS(Loo_prior_performance, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/priorfitsummary.rda")
Loo_prior_performance<-readRDS("scripts/model_outputs/Offspring Trait Models/adult lifespan/priorfitsummary.rda")

#Selecting the weakly informative priors


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
                               y    = data1$adult_lifespan,
                               stat = "mean") + theme_classic()
#posterior cumulative density
weakprior_ecdf     <- pp_check(mod1.1_weakpriors, ndraws = 100,
                               type = "ecdf_overlay") +
  xlim(0, 20) + theme_classic()

#posterior probability density
weakprior_pdf      <- pp_check(mod1.1_weakpriors, ndraws = 100) +
  xlim(0, 20) + theme_classic()

#posterior LOO Q-Q plot
weakprior_loo_qq   <- pp_check(mod1.1_weakpriors, type = "loo_pit_qq",
                               ndraws = 100) + theme_classic()
#LOO-PIT values fall well along the diagonal line (i.e. the expected uniform distribution of PIT integrals)

#LOO uniformity plot
weakprior_loo_unif <- pp_check(mod1.1_weakpriors, type = "loo_pit_overlay",
                               ndraws = 100) + theme_classic()
#LOO predictive interval plots
weakprior_intervals <- ppc_loo_intervals(
  y    = data1$adult_lifespan,
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
#both fixed and random effects explain very little total variation in the offspring's adult lifespan

## MAP estimates and pd ------------------------------------------------------
MAP_adultlifespan<- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  test        = "pd",
  ci_method   = "HDI",
  centrality  = "MAP",
  component   = "all",
  effects     = "full",
  ci          = 0.95
)
saveRDS(MAP_totallifespan,
        "scripts/model_outputs/Offspring Trait Models/adult lifespan/MAP_adultlifespan.rda")

#calculating the rope range (10% of the response SD)
sd_y <- data1 %>%
  filter(adult_lifespan > 0) %>%
  summarise(sd_log10 = sd(log10(adult_lifespan), na.rm = TRUE)*2.302585) %>%
  pull(sd_log10)
rope_response <- c(-0.1 * sd_y, 0.1 * sd_y)


## ROPE -----------------------------------------------------------------------
rope_adultlifespan<- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  rope_range = c(-0.042 , 0.042), #(0.1*SD of the reponse)
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1 #estimating ROPE seperately as it requires the full posterior distribution
)
saveRDS(rope_adultlifespan,
        "scripts/model_outputs/Offspring Trait Models/adult lifespan/rope_adultlifespan.rda")

# ===========================================================================
# 9. PRIOR vs. POSTERIOR SENSITIVITY PLOTS------------------------------------
# ===========================================================================
## Combined prior–posterior ECDF panel -----------------------------------------

#Saving the prior and posterior plots for default versus moderate priors
priorpostplots<-ggarrange(diffuseprior_cumvdis, weakprior_cumvdis, moderateprior_cumvdis,
                          mod1.1_diffusepriors_cumvdis, mod1.1_weakpriors_cumvdis,  mod1.1_moderatepriors_cumvdis, nrow = 2, ncol = 4,
                          labels = c("A", "B", "c", "D", "E", "F"))
install.packages("ragg")
library(ragg)

ggsave(filename = "./bayesian_plots/model fit plots/adult lifespan/priorvpost_ecdf.png",
       plot = priorpostplots, 
       device = png,
       dpi = 300,
       width = 490, 
       height = 340, 
       units = "mm")


# ===========================================================================
# 10. HYPOTHESIS TESTING: INTERACTION MODELS---------------------------------
# ===========================================================================
#Tests the core interaction asked by our paper, does the parental temperature mediate the observed age effect?


#10.1. Three-way interaction: delta age × temperature x offspring sex (mu only) ------------------------
#I.e., are parental age effects on offspring adult lifespan conditional on both the parents temperature treatment and the offspring's sex?
mod2.1<-brm(adult_lifespan|cens(1-event)~
              avg_age+
              delta.age+
              Temp+
              F1_sex+
              delta.age*Temp*F1_sex+
              (1+delta.age|PairID),
            family=weibull(),
            data=data1,
            #With specified priors added
            prior = weakpriors,
            iter=5000,
            control=list(adapt_delta=0.98),
            core=4)
summary(mod2.1)
saveRDS(mod2.1, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod2.1.rda")
mod2.1<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod2.1.rda")

#diagnostics
loo_hypothesis2<-loo(mod2.1)
plot(loo_hypothesis2)

#Loo probability integral transform (PIT) plots
pp_check(mod2.1, type="loo_pit_qq", ndraws=100)
pp_check(mod2.1, type ="loo_pit_overlay", ndraws=100) 


#10.2: Testing the two-way interaction between Parental age: Temperature------------------------------------------------------------------------
mod3.1<-brm(adult_lifespan|cens(1-event)~
              avg_age+ 
              delta.age +
              Temp + 
              F1_sex + 
              delta.age:Temp+
              (1+delta.age|PairID),
            family=weibull,
            data=data1,
            #With specified priors added
            prior = weakpriors,
            iter=5000,
            control=list(adapt_delta=0.98),
            core=4)
summary(mod3.1)
saveRDS(mod3.1, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod3.1.rda")
mod3.1<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod3.1.rda")

## MAP estimates for interaction model--for reporting in the main manuscript
MAP_adultlifespan_mod3.1 <- bayestestR::describe_posterior(
  mod3.1,
  test       = "pd",
  ci_method  = "HDI",
  centrality = "MAP",
  effects    = "full",
  component = "all",
  ci         = 0.95
)
saveRDS(MAP_adultlifespan_mod3.1,
        "scripts/model_outputs/Offspring Trait Models/adult lifespan/MAP_adultlifespan_mod3.1")

#ROPE estimates for interaction
rope_adultlifespan_mod3.1 <- bayestestR::describe_posterior(
  mod3.1,
  rope_range = c(-0.042 , 0.042), #(0.1*SD of the reponse)
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1
)
saveRDS(rope_adultlifespan_mod3.1,
        "scripts/model_outputs/Offspring Trait Models/adult lifespan/rope_adultlifespan_mod3.1.rda")


#10.3: Testing the two-way interaction between parental age and offspring sex--------------------------------------------------------------
mod4.1<-brm(adult_lifespan|cens(1-event)~
              avg_age+ 
              delta.age +
              Temp + 
              F1_sex + 
              delta.age:F1_sex+
              (1+delta.age|PairID),
            family=weibull,
            data=data1,
            #With specified priors added
            prior = weakpriors,
            iter=5000,
            control=list(adapt_delta=0.98),
            core=4)
summary(mod4.1)
saveRDS(mod4.1, file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod4.1.rda")
mod4.1<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/adult lifespan/mod4.1.rda")


## MAP estimates for interaction model--for reporting in the main manuscript
MAP_adultlifespan_mod4.1 <- bayestestR::describe_posterior(
  mod3.1,
  test       = "pd",
  ci_method  = "HDI",
  centrality = "MAP",
  effects    = "full",
  component = "all",
  ci         = 0.95
)
saveRDS(MAP_adultlifespan_mod4.1,
        "scripts/model_outputs/Offspring Trait Models/adult lifespan/MAP_adultlifespan_mod4.1")

#ROPE estimates for interaction
rope_adultlifespan_mod4.1 <- bayestestR::describe_posterior(
  mod3.1,
  rope_range = c(-0.042 , 0.042), #(0.1*SD of the reponse)
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1
)
saveRDS(rope_adultlifespan_mod4.1,
        "scripts/model_outputs/Offspring Trait Models/adult lifespan/rope_adultlifespan_mod4.1.rda")


# ===========================================================================
# 11. MODEL SELECTION----------------------------------------------------------
# ===========================================================================
hypothesis_fit<-loo(mod1.1_weakpriors,
                    mod2.1,
                    mod3.1,
                    mod4.1, 
                    moment_match = TRUE) 
saveRDS(hypothesis_fit, file = "scripts/model_outputs/Offspring Trait Models/total lifespan/hypothesis_fit.rda")
#Selecting model 1: No interactions --> Just single effect predictors

## Selective disappearance test -----------------------------------------------
hypothesis(mod1.1_weakpriors, "avg_age - delta.age = 0")
hypothesis(mod1.1_weakpriors, "avg_age - delta.age > 0")
hypothesis(mod1.1_weakpriors, "avg_age - delta.age < 0")

## Marginal means via emmeans ------------------------------------------------
#Effect of parental age
pairwise_estimates<- emmeans(mod1.1_weakpriors, "delta.age", type = "response", at = list(delta.age = c(-1.96, 1.96)))
pairs(pairwise_estimates) 

# ===========================================================================
# 12. POSTERIOR DISTRIBUTION PLOTS  (mod1.1_weakpriors)-----------------------
# ===========================================================================
# Plots styled using stat_halfeye (ggdist) for consistency with other traits

post1 <- as_draws_df(mod1.1_weakpriors)
post2<-as_draws_df(mod3.1)
post3<-as_draws_df(mod4.1)

## Build combined posterior data frame ----------------------------------------
posterior_df_mod1.1 <- data.frame(
  "μ: Parents' Δage"                   = post1$b_delta.age,
  "μ: Parents' average age"            = post1$b_avg_age,
  "μ: Parent Temperature (28.0°C)"     = post1$b_Temp28,
  "μ: Parent Temperature (30.5°C)"     = post1$b_Temp30.5,
  "μ: Offspring sex (Male)"            = post1$b_F1_sexM,
  "μ: Δage × Temperature (28.0°C)"    = post2$`b_delta.age:Temp28`,
  "μ: Δage × Temperature (30.5°C)"    = post2$`b_delta.age:Temp30.5`,
  "μ: Δage × Offspring sex (male)"    = post3$`b_delta.age:F1_sexM`,
  check.names = FALSE
) %>%
  pivot_longer(everything(),
               names_to  = "parameter",
               values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = c(
    "μ: Parents' Δage",
    "μ: Parents' average age",
    "μ: Offspring sex (Male)",
    "μ: Parent Temperature (28.0°C)",
    "μ: Parent Temperature (30.5°C)",
    "μ: Δage × Temperature (28.0°C)",
    "μ: Δage × Temperature (30.5°C)",
    "μ: Δage × Offspring sex (male)"
  )))

y_levels <- rev(c(
  "μ: Parents' Δage",
  "μ: Parents' average age",
  "μ: Parent Temperature (28.0°C)",
  "μ: Parent Temperature (30.5°C)",
  "μ: Offspring sex (Male)",
  "μ: Δage × Temperature (28.0°C)",
  "μ: Δage × Temperature (30.5°C)",
  "μ: Δage × Offspring sex (male)"
))

## Halfeye posterior plot ------------------------------------------------------
posterior_plot_adultlife <- ggplot(
  posterior_df_mod1.1,
  aes(x = value, y = parameter, fill = parameter)
) +
  annotate("rect",
           xmin = -0.042, xmax = 0.042,
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
    "#e5a9f5", "#d2e9f5", "#d7c9a3", "#4c9279"
  )) +
  xlim(c(-0.1, 0.20)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title      = element_text(size = 30),
    axis.text       = element_text(size = 30)
  ) +
  labs(x = "Fixed-effect posterior estimates\n(on μ = mean adult lifespan (weeks), log-scale)",
       y = NULL)

ggsave("./bayesian_plots/posterior plots/adult lifespan/001_hypothesis1_halfeye.png",
       plot   = posterior_plot_adultlife,
       device = "png",
       width  = 580, height = 400, units = "mm")


## Post-hoc standardised mean difference (SMD) --------------------------------
# Compares extreme ends of within-individual parental age; standardised by
# observed SD of F1 adult lifespan. Descriptive only — not used for inference.

#Post-hoc standardised mean difference (on unit scale)
pairwise_estimates<- emmeans(mod1.1_weakpriors, "delta.age", type = "response", at = list(delta.age = c(-1.96, 1.96)))
#This is now converted to the raw scale where the extreme ends of parental age are compared
sd_life<-sd(data1$adult_lifespan)

#standardised emmeans --> estimating difference across the entire range of parental age
emm_life_scaled <- summary(pairwise_estimates) %>%
  mutate(emmean_sd = response / sd_life,
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
adult_lifespan_SMD<-data.frame(
  SMD = SMD,
  SMD_lower = SMD_lower,
  SMD_upper = SMD_upper
)

#Saving for use in combined forest plot of all traits
saveRDS(adult_lifespan_SMD, file = "scripts/model_outputs/Offspring Trait Models/adult_lifespan_SMD")


# ===========================================================================
# 13. MODEL PREDICTIONS AND RESULTS PLOTS-------------------------------------
# ===========================================================================

## 13.1 Prediction grids: within- and between-individual age effects ---------

#13.1.1 : Allowing predictions to vary for delta age. Holding average age constant.
df_predict_within<-expand.grid(
  avg_age=mean(data1$avg_age),
  delta.age=unique(data1$delta.age),
  Temp=unique(data1$Temp),
  F1_sex=unique(data1$F1_sex),
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
  mutate(within_subject_age= delta.age*sd(data1$within_subject_age)) #converts standardised scale back to weeks (mean-centred)

#13.2: Allowing predictions to vary for average age. Holding within parental age constant.
df_predict_between<-expand.grid(
  avg_age=unique(data1$avg_age), 
  delta.age=mean(data1$delta.age),
  Temp=unique(data1$Temp),
  F1_sex = unique(data1$F1_sex),
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
    avg.age2 = avg_age * sd(data1$avg.age) + mean(data1$avg.age),
    avg_timepoint_centered = avg.age2 - mean(avg.age2)
  )

## 13.2 Raw within-individual means (for overlaying on plot) ----------------
breaks <- seq(-4.5, 4.5, by = 1)  # creates bins: [-3.5, -2.5], [-2.5, -1.5], ..., [2.5, 3.5]
labels <- -4:4  # label bins with week numbers

data1<- data1 %>%
  mutate(delta_age_bin = cut(within_subject_age, breaks = breaks, labels = labels))


#Raw delta age means
raw_deltaage <- data1%>%
  select(-within_subject_age) %>%        
  filter(!is.na(delta_age_bin)) %>% 
  rename(within_subject_age = delta_age_bin) %>% 
  mutate(within_subject_age= as.numeric(as.character(within_subject_age))) %>%  # convert factor to numeric
  group_by(within_subject_age) %>%
  summarise(
    n = sum(!is.na(adult_lifespan)),                       # number of non-missing observations
    mean_life = mean(adult_lifespan, na.rm = TRUE),        # group mean
    se_life = ifelse(n > 1, sd(adult_lifespan, na.rm = TRUE)/sqrt(n), NA_real_), # SE only if n>1
    .groups = "drop"
  )

## 13.5 Main parental age effect plot on adult lifespan----------------------------------------
adult_lifespan_plot <- ggplot(data = data1,
                              aes(x = within_subject_age, 
                                  y = adult_lifespan,
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
  labs(x = "Parents' adult age at reproduction (weeks; mean-centred)", y = "Offspring's adult lifespan (weeks)")

# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 adult lifespan/001_lifespan_plot.png",
       plot = adult_lifespan_plot,
       bg="transparent",
       device = "png", 
       width = 500, 
       height = 700, 
       units = "mm")


#13.6. Plotting the effect of parental age on adult mortality and survival----------------------------
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
#part 1: Allowing predictions to vary for within_subject_timepoint. Holding average age constant.
newdat<-expand.grid(
  avg_age=mean(data1$avg_age),
  delta.age=c(-1.95, 0, 1.85),
  Temp=unique(data1$Temp),
  F1_sex = unique(data1$F1_sex),
  PairID=unique(data1$PairID))

# Create a sequence over which to evaluate the survival and hazard
life_seq <- seq(min(0), 
                max(data1$adult_lifespan, na.rm=TRUE), 
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
    ndraws = 1000, 
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
  scale_x_continuous(limits = c(0, 18), breaks = c(0, 3, 6, 9, 12, 15))+
  scale_color_manual(name = "Parents' adult age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' adult age at reproduction", values = c(0.6,0.1, 0.6))+
  scale_fill_manual(name = "Parents' adult age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  labs(x = "Offspring adult age (weeks)", y = "Cumulative survival probability, S(x)")

# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 adult lifespan/001_adultlifespan_survivalplot.png",
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
  ylim(c(0, 2.5))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  scale_x_continuous(limits = c(0, 18), breaks = c(0, 3, 6, 9, 12, 15))+
  scale_color_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' age at reproduction", values = c(0.6,0.1, 0.6))+
  scale_fill_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  labs(x = "Offspring adult age (weeks)", y = "Instantaneous hazard rate, μ (X)")

# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 adult lifespan/001_adultlifespan_hazard.png",
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

lifespan_inference<-(adult_lifespan_plot | right_panel) +
  plot_layout(widths = c(0.8, 1)) +
  theme(plot.tag = element_text(size = 50, face = "bold"),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), )

ggsave(filename = "./bayesian_plots/offspring trait plots/F1 adult lifespan/001_inference_plots2.png",
       plot = lifespan_inference, 
       device = "png", 
       width = 920, 
       height = 630, 
       units = "mm")

#Creating a table to export------------------------------------------------------------------------------

# ── helper ────────────────────────────────────────────────────────────────────
#Model draws from selected models
base_model_draws<-as_draws_df(mod1.1_weakpriors)
interaction_draws<-as_draws_df(mod2.1) 
interaction_draws2<-as_draws_df(mod3.1) 
interaction_draws3<-as_draws_df(mod4.1)

#estimating the difference between average age and delta age
differences<-data.frame(
  selectivedis = base_model_draws$b_avg_age - base_model_draws$b_delta.age
)

#function for generating estimates
summarise_param <- function(draws, param, section,
                            rope = FALSE, rope_range = c(-0.042, 0.042), 
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
bind_rows(fe, fe_interaction2, re_sd_explore, shape, bayes) %>%
  gt(groupname_col = "Section") %>%
  cols_label(Parameter = "Parameter", MAP = "MAP",
             `95% HDI` = "95% HDI", `% in ROPE` = "% in ROPE", pd = "pd") %>%
  tab_header(
    title    = "Adult lifespan model summary: Weibull model"
  ) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) %>%
  tab_footnote("ROPE = [−0.042, 0.042] on log scale; applied to location fixed effects only.",
               locations = cells_column_labels(`% in ROPE`)) %>%
  tab_footnote("pd = probability of direction.",
               locations = cells_column_labels(pd)) %>%
  opt_stylize(style = 1)%>%
  gtsave("./tables/adultlifespan_updated.docx")

###############################################################################
##END OF SCRIPT----------------------------------------------------------------
###############################################################################
