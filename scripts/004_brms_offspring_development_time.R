###############################################################################
##  Script: Effect of parental age on offspring development time
##  Note:   Analysis restricted to offspring that survived to adulthood
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
F1data<-readRDS("./raw data/F1_filtered_data06022025.RDS")
length(unique(F1data$F1_ID)) 

#filtering out animals that didn't reach adulthood/had uncertain adult emergence dates
data1<-F1data%>%
  filter(Include_in_adult=="Y") %>% 
  filter(adult_surv==1) 

#Creating a development time variable: Calculating the difference between birth dates and death dates
data1$development_time_weeks<- as.numeric(difftime(data1$F1_Adult_emergence, data1$F1_hatch, units = "weeks"))


#Creating an "event" column to signal when the offspring reached adulthood
#reuqired for time-to-event models (i.e., Lognormal, Weibull, Gamma, exp.)
data1$event<-ifelse(is.na(data1$development_time_weeks), 0,1)


## Scale continuous predictors ------------------------------------------------
# All predictors scaled to mean = 0, SD = 1; a 1-unit change reflects 1 SD
data1$avg_age<-as.numeric(scale(data1$avg.age, center = TRUE, scale = TRUE))


# Within-individual (delta) age: scaled but NOT mean-centred to preserve
# the within-subject structure; NAs replaced with 0
data1 <- data1 %>%
  group_by(PairID) %>%
  mutate(
    delta.age = as.numeric(scale(within_subject_age,
                                 center = FALSE,
                                 scale = TRUE)),
    delta.age = replace(delta.age, is.na(delta.age), 0)
  ) %>%
  ungroup()

#coding sex as a factor
data1$F1_sex <- factor(data1$F1_sex, levels = c("F", "M"))


## Sample sizes ---------------------------------------------------------------
length(unique(data1$Mother_ID))  # number of mothers (N = 77)
length(unique(data1$F1_ID))      # number of F1 offspring (n = 939)


# ===========================================================================
# 3. SUMMARY STATISTICS----------------------------------------------------------
# ===========================================================================

## Distribution of F1 development time -----------------------------------------
ggplot(data1, aes(x = development_time_weeks)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +  # Create histogram
  facet_wrap(~Timepoint) +  
  labs(
    x = "F1 developement time (weeks)",
    y = "Count"
  ) +
  theme_classic()


## Overall mean and SD (used to build intercept priors later) ----------------------
sum_overall <- data1 %>%
  summarise(
    F1_mass     = mean(development_time_weeks),
    sd_total    = sd(development_time_weeks),
    n_offspring = n()
  )
#mean development time of 9.19 weeks, SD of 1.48 weeks, n =939, no.parents = 77


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
vif(lm(development_time_weeks ~ avg_age + delta.age + cum_successful_matings+
         Temp + F1_sex, data = data1))

#removing cum_succesful matings
vif(lm(development_time_weeks ~ avg_age + delta.age +
         Temp + F1_sex, data = data1))

# ===========================================================================
# 4. MODEL STRUCTURE SELECTION------------------------------------------------
# ===========================================================================

#---------------Deciding on a family to use for the models----------------------
#Data is time-to-event (testing the Weibull, Lognormal, and exponential distributions here)

## 4.1: Weibull model-----------------------------------------------------------
default_weibull<-brm(development_time_weeks|cens(1-event)~
                       avg_age+
                       delta.age+
                       Temp+
                       F1_sex+
                       (1|PairID),
                     family=weibull,
                     save_pars = save_pars(all = TRUE),
                     data=data1,
                     iter=5000,
                     cores = 4)
summary(default_weibull) 
saveRDS(default_weibull , 
        file = "scripts/model_outputs/Offspring Trait Models/development time/defaultweibull.rda")
default_weibull<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/defaultweibull.rda")

#diagnostics
Loo_Weibull<-LOO(default_weibull, save_psis = TRUE, moment_match=TRUE)
plot(Loo_Weibull)

#Loo probability integral transform (PIT) plots
pp_check(default_weibull, type="loo_pit_qq", ndraws=100)
pp_check(default_weibull, type ="loo_pit_overlay", ndraws=100) 

#Inspecting the fit with a KM curve
default_weibull_KM<- pp_check(default_weibull,
                              status_y=data1$event,
                              type="km_overlay", ndraws=100) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  ylab("S(x)") + coord_cartesian(xlim=c(0, 20))+
  xlab("development time (weeks)") + coord_cartesian(xlim=c(0, 20))

## 4.2: Lognormal model---------------------------------------------------------
default_lognormal<-brm(development_time_weeks|cens(1-event)~
                         avg_age+
                         delta.age+
                         Temp+
                         F1_sex+
                         (1|PairID),
                       family=brms::lognormal,
                       data=data1,
                       iter=5000,
                       save_pars = save_pars(all = TRUE),
                       cores = 4)
summary(default_lognormal)
saveRDS(default_lognormal , file = "scripts/model_outputs/Offspring Trait Models/development time/defaultlognormal.rda")
default_lognormal<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/defaultlognormal.rda")

#MODEL DIAGNOSTICS
loo_lognormal<-loo(default_lognormal,  save_psis = TRUE, moment_match =TRUE)
plot(loo_lognormal)

#Loo probability integral transform (PIT) plots
pp_check(default_lognormal, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(default_lognormal, type ="loo_pit_overlay", ndraws=100) 


#Inspecting the fit with a KM curve
default_lognormal_KM<- pp_check(default_lognormal,
                              status_y=data1$event,
                              type="km_overlay", ndraws=100) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  ylab("S(x)") + coord_cartesian(xlim=c(0, 20))+
  xlab("Development time (weeks)") + coord_cartesian(xlim=c(0, 20))



#4.3: Default exponential-------------------------------------------------------
default_exponential<-brm(development_time_weeks|cens(1-event)~
                           avg_age+
                           delta.age+
                           Temp+
                           F1_sex+
                           (1|PairID),
                         save_pars = save_pars(all = TRUE),
                         family=exponential,
                         data=data1,
                         iter=5000,
                         cores = 4)
summary(default_exponential) #Robust = TRUE gives the posterior median
saveRDS(default_exponential, file = "scripts/model_outputs/Offspring Trait Models/development time/defaultexponential.rda")
default_exponential<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/defaultexponential.rda")

#model diagnostics
loo_exponential<-loo(default_exponential,  save_psis = TRUE, moment_match=TRUE)
plot(loo_exponential)

#Loo probability integral transform (PIT) plots
pp_check(loo_exponential, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(loo_exponential, type ="loo_pit_overlay", ndraws=100) 

#Inspecting distribution
default_exponential_KM<- pp_check(default_lognormal,
                              status_y=data1$event,
                              type="km_overlay", ndraws=100) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  ylab("S(x)") + coord_cartesian(xlim=c(0, 20))+
  xlab("Development time (weeks)") + coord_cartesian(xlim=c(0, 20))


#--------------------Comparison of time to event models------------------------------------------
LOO_allsurvivalmodels<-loo(default_weibull, default_lognormal, 
                           default_exponential, moment_match = TRUE)
saveRDS(LOO_allsurvivalmodels, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/LOO_allsurvivalmodels.rda")
LOO_allsurvivalmodels<-readRDS("scripts/model_outputs/Offspring Trait Models/development time/LOO_allsurvivalmodels.rda")

#Choosing the lognormal model as the best fitting


## 4.4 Quadratic delta-age term ------------------------------------------------
#Adding the quadratic delta age term into the model
mod1.1_quadraticage<-brm(development_time_weeks|cens(1-event)~
                           avg_age+
                           poly(delta.age,2)+
                           Temp+
                           F1_sex+
                           (1|PairID),
                         family=brms::lognormal,
                         data=data1,
                         iter=5000,
                         control=list(adapt_delta=0.98),
                         save_pars = save_pars(all = TRUE),
                         core=4)
summary(mod1.1_quadraticage) 
marginal_effects(mod1.1_quadraticage, effects = "delta.age") #No support for a quadratic age effect
saveRDS(mod1.1_quadraticage, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_quadraticage.rda")
mod1.1_quadraticage<-readRDS("scripts/model_outputs/Offspring Trait Models/development time/mod1.1_quadraticage.rda")

#MODEL DIAGNOSTICS
loo_quadratic<-loo(mod1.1_quadraticage, save_psis = TRUE, moment_match = TRUE)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_quadraticage, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(mod1.1_quadraticage, type ="loo_pit_overlay", ndraws=100) 


## 4.5 Random slopes (intercept + delta age) -----------------------------------
mod1.1_randomslopes<-brm(development_time_weeks|cens(1-event)~
                           avg_age+ 
                           delta.age +
                           Temp + 
                           F1_sex + 
                           (1+delta.age|PairID),
                         family=brms::lognormal,
                         data=data1,
                         iter=5000,
                         save_pars = save_pars(all = TRUE),
                         control=list(adapt_delta=0.98),
                         core=4)
summary(mod1.1_randomslopes)
saveRDS(mod1.1_randomslopes , 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_randomslopes.rda")
mod1.1_randomslopes<-readRDS("scripts/model_outputs/Offspring Trait Models/development time/mod1.1_randomslopes.rda")

#model diagnostics
loo_randomslopes<-loo(mod1.1_randomslopes, save_psis=TRUE, moment_match = TRUE)
plot(loo_randomslopes)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_randomslopes, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(mod1.1_randomslopes, type ="loo_pit_overlay", ndraws=100) 


## 4.6 Random slopes + quadratic delta-age -------------------------------------
#where population level slope is a quadratic function (but assumes all females share the same quadratic trajectory)
mod1.1_qudratic_randomslopes<-brm(development_time_weeks|cens(1-event)~
                                    avg_age+ 
                                    poly(delta.age,2)+
                                    Temp + 
                                    F1_sex + 
                                    (1+delta.age|PairID),
                                  family=brms::lognormal,
                                  data=data1,
                                  iter=5000,
                                  save_pars = save_pars(all = TRUE),
                                  control=list(adapt_delta=0.98),
                                  core=4)
summary(mod1.1_qudratic_randomslopes)
saveRDS(mod1.1_qudratic_randomslopes, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_qudratic_randomslopes.rda")
mod1.1_qudratic_randomslopes<-readRDS("scripts/model_outputs/Offspring Trait Models/development time/mod1.1_qudratic_randomslopes.rda")

#Model diagnostics
loo_quadratic_randomslopes<-loo(mod1.1_quadratic_randomslopes, save_psis = TRUE, moment_match = TRUE)
plot(loo_quadratic)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_qudratic_randomslopes, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(mod1.1_qudratic_randomslopes, type ="loo_pit_overlay", ndraws=100) 


## 4.7 LOO comparison: random-effects and quadratic delta age structure ------------------
LOO_ageeffects<-loo(default_lognormal, mod1.1_randomslopes, 
                    mod1.1_quadraticage,mod1.1_qudratic_randomslopes, 
                    moment_match = TRUE)
saveRDS(LOO_ageeffects, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/LOO_ageeffects.rda")
LOO_ageeffects<-readRDS("scripts/model_outputs/Offspring Trait Models/development time/LOO_ageeffects.rda")

#Selecting the linear delta age model with random slopes (removing quadratic effect of delta age)


## 4.8 Distributional models: evidence for a sigma submodel ---------------------------------------------------
# Does parental age (and/or temperature) predict residual variance?---------

#4.9. Effect of delta age and temp on sigma
mod1.1_sigma1<-brm(bf(development_time_weeks|cens(1-event)~
                        avg_age+ 
                        delta.age +
                        Temp + 
                        F1_sex + 
                        (1+delta.age|PairID),
                      sigma~delta.age+
                        Temp), #effect of temperature weakly identified
                   family=brms::lognormal,
                   data=data1,
                   iter=5000,
                   save_pars = save_pars(all = TRUE),
                   control=list(adapt_delta=0.98),
                   core=4)
summary(mod1.1_sigma1)
saveRDS(mod1.1_sigma1 , 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_sigma1.rda")
mod1.1_sigma1<-readRDS("scripts/model_outputs/Offspring Trait Models/development time/mod1.1_sigma1.rda")

#model diagnostics
loo_sigma<-loo(mod1.1_sigma1, moment_match = TRUE, save_psis = TRUE)
plot(loo_sigma)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_sigma1, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(mod1.1_sigma1, type ="loo_pit_overlay", ndraws=100) 


#4.10. Effect of just delta age on sigma----------------------------------------------
mod1.1_sigma<-brm(bf(development_time_weeks|cens(1-event)~
                       avg_age+ 
                       delta.age +
                       Temp + 
                       F1_sex + 
                       (1+delta.age|PairID),
                     sigma~delta.age),#Not enough data per timepoint per female to add random effects on the sigma submodel
                  family=brms::lognormal,
                  data=data1,
                  iter=5000,
                  save_pars = save_pars(all = TRUE),
                  control=list(adapt_delta=0.98),
                  core=4)
summary(mod1.1_sigma)
saveRDS(mod1.1_sigma , 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_sigma.rda")
mod1.1_sigma<-readRDS("scripts/model_outputs/Offspring Trait Models/development time/mod1.1_sigma.rda")

#model diagnostics
loo_sigma_age<-loo(mod1.1_sigma, moment_match = TRUE, save_psis = TRUE)
plot(loo_sigma_age)

#Loo probability integral transform (PIT) plots
pp_check(mod1.1_sigma, type="loo_pit_qq", ndraws=100) #very good fit
pp_check(mod1.1_sigma, type ="loo_pit_overlay", ndraws=100) 

#--------Comparison of model with distributional parameters------------------
modelsetup<-loo(mod1.1_sigma, mod1.1_sigma1,
                default_lognormal) 
saveRDS(modelsetup, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/modelsetup.rda")


# ===========================================================================
# 5. PRIOR SPECIFICATION------------------------------------------------------
# ===========================================================================
# Priors assume all continuous predictors are mean-centered and SD-scaled.

#---------------------Setting up the priors-------------------------------------
#what's the range of our data?
range(log(data1$development_time_weeks)) #Crickets range in log dev time from 1.76 to 2.87 weeks

mean(data1$development_time_weeks)#9.19 weeks
sd(data1$development_time_weeks)#1.48 weeks

#mean on log scale --> log is a non-linear transformation, this needs to be corrected later
log(9.19) # 2.22 

#SD on log scale
sqrt(log(1 + (1.48^2 /9.19^2))) #0.16 on log scale

#correcting mean on log scale (since it's a non-linear transformation)
2.22 - ((0.16^2)/2) #=2.21 
#intercept should be set at 2.21

#Setting a prior on sigma
#Priors for sigma
sigma_emp <- sd(log(data1$development_time_weeks))
#The empirical sigma is 0.154 on the log scale, a good starting value for sigma intercept


#---------------------SETTING MODEL PRIORS--------------------------------------
# Priors assume all continuous predictors are mean-centred and SD-scaled.
# Intercept prior informed by observed mean (19.9 weeks) and SD (1.48 weeks).

#--------------------setting up the priors--------------------------------------

#Diffuse priors-----------------------------------------------------------------
#similar to the brms defaults
Diffusepriors <- c(
  prior(normal(2.21, 0.16), class = "Intercept"),             
  prior(normal(0,1), class = "b"),          
  prior(normal(log(0.154), 0.5), class = "Intercept", dpar="sigma"), #sigma goes through two rounds of log transformations in the brms lognormal parameterisation...(hence why I'm taking the log again here)
  prior(normal(0,1), class = "b", dpar="sigma"),   
  prior(student_t(3, 0, 2.5), class = "sd", lb = 0),           
  prior(lkj(2), class = "cor")                                 
)


#Weakly informative priors [SELECTED]-------------------------------------------
weakpriors <- c(
  prior(normal(2.21, 0.16), class = "Intercept"),                
  prior(normal(0,0.16), class = "b"), 
  prior(normal(log(0.154), 0.5), class = "Intercept", dpar="sigma"), #model uses log link again
  prior(normal(0,0.5), class = "b", dpar="sigma"),               
  prior(exponential(10), class = "sd", lb = 0),                      
  prior(lkj(2), class = "cor")                                   
)


#Moderate prior-----------------------------------------------------------------
moderatepriors <- c(
  prior(normal(2.21, 0.16), class = "Intercept"),                 
  prior(normal(0,0.08), class = "b"),  
  prior(normal(log(0.15), 0.5), class = "Intercept", dpar="sigma"),
  prior(normal(0,0.5), class = "b", dpar="sigma"),  
  prior(exponential(10), class = "sd", lb = 0),                         
  prior(lkj(2), class = "cor")                                    
)

#Constraints prior --> assuming parental age extends development time-----------
constraintpriors <- c(
  prior(normal(2.21, 0.16), class = "Intercept"),                 
  prior(normal(0,0.08), class = "b"),                               
  prior(normal(0.186, 0.03), class = "b", coef="delta.age"), 
  prior(normal(log(0.15), 0.5), class = "Intercept", dpar="sigma"),
  prior(normal(0,0.5), class = "b", dpar="sigma"),  
  prior(exponential(10), class = "sd", lb = 0),                          
  prior(lkj(2), class = "cor")                                    
)

#Fast pace-of-life prior [assuming parental age accelerates pace-of-life]------
fastprior <- c(
  prior(normal(2.21, 0.16), class = "Intercept"),                
  prior(normal(0,0.08), class = "b"),                               
  prior(normal(-0.186, 0.03), class = "b", coef="delta.age"),  
  prior(normal(log(0.15), 0.5), class = "Intercept", dpar="sigma"),
  prior(normal(0,0.5), class = "b", dpar="sigma"), 
  prior(exponential(10), class = "sd", lb = 0),                          
  prior(lkj(2), class = "cor")                                   
)


# ===========================================================================
# 6. PRIOR PREDICTIVE CHECKS---------------------------------------------------
# ===========================================================================

## 6.1 Prior-only model: diffuse ---------------------------------------------------
diffuse_prior_model<-brm(bf(development_time_weeks|cens(1-event)~
                              avg_age+ 
                              delta.age +
                              Temp + 
                              F1_sex + 
                              (1+delta.age|PairID),
                            sigma~delta.age),
                         family=brms::lognormal,
                         data=data1,
                         iter=5000,
                         #With specified priors added
                         prior = Diffusepriors,
                         sample_prior = "only",
                         cores = 4)
summary(diffuse_prior_model)
saveRDS(diffuse_prior_model , 
        file = "scripts/model_outputs/Offspring Trait Models/development time/diffusepriormodel.rda")
diffuse_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/diffusepriormodel.rda")

#setting the colour scheme
color_scheme_set("teal")

#Prior draws (I.e., the prior cumulative density function)
diffuseprior_cumvdis<-pp_check(diffuse_prior_model,
                               ndraws=100, type = "ecdf_overlay",
                               discrete=FALSE)+ coord_cartesian(xlim=c(0, 20))+
  theme_classic()



## 6.2 Prior-only model: weak ----------------------------------------------------
weak_prior_model <-brm(bf(development_time_weeks|cens(1-event)~
                            avg_age+ 
                            delta.age +
                            Temp + 
                            F1_sex + 
                            (1+delta.age|PairID),
                          sigma~delta.age),
                       family=brms::lognormal,
                       data=data1,
                       iter=5000,
                       #With specified priors added
                       prior = weakpriors,
                       sample_prior = "only",
                       cores = 4)
summary(weak_prior_model)
saveRDS(weak_prior_model , 
        file = "scripts/model_outputs/Offspring Trait Models/development time/weakpriors.rda")
weak_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/weakpriors.rda")

#prior draws (for the weak priors)
weakprior_cumvdis<-pp_check(weak_prior_model, 
                            ndraws=100, type = "ecdf_overlay",
                            discrete=FALSE)+ coord_cartesian(xlim=c(0, 20))+
  theme_classic()

##6.3 Prior-only model: moderate --------------------------------------------------
moderate_prior_model <-brm(bf(development_time_weeks|cens(1-event)~
                                avg_age+ 
                                delta.age +
                                Temp + 
                                F1_sex + 
                                (1+delta.age|PairID),
                              sigma~delta.age),
                           family=brms::lognormal,
                           data=data1,
                           iter=5000,
                           #With specified priors added
                           prior = moderatepriors,
                           sample_prior = "only",
                           cores = 4)
summary(moderate_prior_model)
saveRDS(moderate_prior_model , 
        file = "scripts/model_outputs/Offspring Trait Models/development time/moderatepriors.rda")
moderate_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/moderatepriors.rda")

#prior draws (for the weak priors)
moderateprior_cumvdis<-pp_check(moderate_prior_model, 
                                ndraws=100, type = "ecdf_overlay",
                                discrete=FALSE)+coord_cartesian(xlim=c(0, 20))+
  theme_classic()


##6.4 Constraint prior draws-----------------------------------------------------
constraint_prior_model <-brm(bf(development_time_weeks|cens(1-event)~
                                  avg_age+ 
                                  delta.age +
                                  Temp + 
                                  F1_sex + 
                                  (1+delta.age|PairID),
                                sigma~delta.age),
                             family=brms::lognormal,
                             data=data1,
                             iter=5000,
                             #With specified priors added
                             prior =constraintpriors,
                             sample_prior = "only",
                             cores = 4)
summary(constraint_prior_model)
saveRDS(constraint_prior_model , 
        file = "scripts/model_outputs/Offspring Trait Models/development time/constraintpriors.rda")
constraint_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/constraintpriors.rda")

#Prior draws
constraintprior_cumvdis<-pp_check(constraint_prior_model, 
                                  ndraws=100, type = "ecdf_overlay", discrete=FALSE)+
  coord_cartesian(xlim=c(0, 20))+theme_classic()

##6.5 Fast-pace-of-life prior draws---------------------------------------------
fastlife_prior_model <-brm(bf(development_time_weeks|cens(1-event)~
                                avg_age+ 
                                delta.age +
                                Temp + 
                                F1_sex + 
                                (1+delta.age|PairID),
                              sigma~delta.age),
                           family=brms::lognormal,
                           data=data1,
                           iter=5000,
                           #With specified priors added
                           prior = fastprior,
                           sample_prior = "only",
                           cores = 4)
summary(fastlife_prior_model)
saveRDS(fastlife_prior_model ,
        file = "scripts/model_outputs/Offspring Trait Models/development time/fastpriors.rda")
fastlife_prior_model<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/fastpriors.rda")

#Prior draws
fastprior_cumvdis<-pp_check(fastlife_prior_model, ndraws=100, type = "ecdf_overlay", discrete=FALSE)+
  coord_cartesian(xlim=c(0, 20))+theme_classic()


# ===========================================================================
# 7. POSTERIOR MODELS---------------------------------------------------------
# ===========================================================================

## 7.1 Diffuse priors ---------------------------------------------------------
mod1.1_diffusepriors<-brm(bf(development_time_weeks|cens(1-event)~
                               avg_age+ 
                               delta.age +
                               Temp + 
                               F1_sex + 
                               (1+delta.age|PairID),
                             sigma~delta.age),
                          family=brms::lognormal,
                          data=data1,
                          #With specified priors added
                          prior = Diffusepriors,
                          iter=5000,
                          control=list(adapt_delta=0.98),
                          save_pars = save_pars(all = TRUE),
                          core=4)
summary(mod1.1_diffusepriors)
saveRDS(mod1.1_diffusepriors, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_diffusepriors.rda")
mod1.1_diffusepriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_diffusepriors.rda")

color_scheme_set("pink")

#Posterior draws
mod1.1_diffusepriors_cumvdis<-pp_check(mod1.1_diffusepriors, ndraws=100, type = "ecdf_overlay", discrete=FALSE)+
  coord_cartesian(xlim=c(0, 20))+theme_classic()


## 7.2 Weak priors  [SELECTED MODEL] -------------------------------------------
mod1.1_weakpriors<-brm(bf(development_time_weeks|cens(1-event)~
                            avg_age+ 
                            delta.age +
                            Temp + 
                            F1_sex + 
                            (1+delta.age|PairID),
                          sigma~delta.age),
                       family=brms::lognormal,
                       data=data1,
                       #With specified priors added
                       prior = weakpriors,
                       iter=5000,
                       control=list(adapt_delta=0.98),
                       save_pars = save_pars(all = TRUE),
                       core=4)
summary(mod1.1_weakpriors)
saveRDS(mod1.1_weakpriors, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_weakpriors.rda")
mod1.1_weakpriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_weakpriors.rda")

#posterior draws
weakprior_ecdf1<-pp_check(mod1.1_weakpriors, 
                          ndraws=100, type = "ecdf_overlay", discrete=FALSE)+
  coord_cartesian(xlim=c(0, 20))+theme_classic()


## 7.3 Moderate priors----------------------------------------------------------
mod1.1_moderatepriors<-brm(bf(development_time_weeks|cens(1-event)~
                                avg_age+ 
                                delta.age +
                                Temp + 
                                F1_sex + 
                                (1+delta.age|PairID),
                              sigma~delta.age),
                           family=brms::lognormal,
                           data=data1,
                           #With specified priors added
                           prior = moderatepriors,
                           iter=5000,
                           control=list(adapt_delta=0.98),
                           core=4)
summary(mod1.1_moderatepriors)
saveRDS(mod1.1_moderatepriors, file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_moderatepriors.rda")
mod1.1_moderatepriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_moderatepriors.rda")

#posterior draws
moderateprior_ecdf<-pp_check(mod1.1_moderatepriors,
                             ndraws=100, type = "ecdf_overlay")+
  xlim(0, 20)+theme_classic()  

## 7.4 default (flat) priors----------------------------------------------------
mod1.1_defaultprior<- brm(bf(development_time_weeks|cens(1-event)~
                               avg_age+ 
                               delta.age +
                               Temp + 
                               F1_sex + 
                               (1+delta.age|PairID),
                             sigma~delta.age),
                          family = brms::lognormal,
                          data = data1,
                          #No selected priors here, using the improper, default flat priors
                          iter=5000,
                          control=list(adapt_delta=0.98),
                          cores = 4)
summary(mod1.1_defaultprior)
saveRDS(mod1.1_defaultprior , 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_defaultprior.rda")
mod1.1_defaultprior<-readRDS("scripts/model_outputs/Offspring Trait Models/development time/mod1.1_defaultprior.rda")

## 7.5 Constraint priors (assume parental age extends dev time)-----------------
mod1.1_constraintpriors<-brm(bf(development_time_weeks|cens(1-event)~
                                  avg_age+ 
                                  delta.age +
                                  Temp + 
                                  F1_sex + 
                                  (1+delta.age|PairID),
                                sigma~delta.age),
                             family=brms::lognormal,
                             data=data1,
                             #With specified priors added
                             prior = constraintpriors,
                             iter=5000,
                             control=list(adapt_delta=0.98),
                             core=4)
summary(mod1.1_constraintpriors) #model has quite a large residual standard deviation (sigma)--> need to identify why? This effects the cohens D estimates.
saveRDS(mod1.1_constraintpriors, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_constraintpriors.rda")
mod1.1_constraintpriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_constraintpriors.rda")

#posterior draws
constraintsprior_ecdf<-pp_check(mod1.1_constraintpriors, ndraws=100, type = "ecdf_overlay")+xlim(0, 20)+theme_classic()  

## 7.6 Fast pace-of-life priors (assume parental age reduces devtime)-----------
mod1.1_fastlifepriors<-brm(bf(development_time_weeks|cens(1-event)~
                                avg_age+ 
                                delta.age +
                                Temp + 
                                F1_sex + 
                                (1+delta.age|PairID),
                              sigma~delta.age),
                           family=brms::lognormal,
                           data=data1,
                           #With specified priors added
                           prior = fastprior,
                           iter=5000,
                           control=list(adapt_delta=0.98),
                           core=4)
summary(mod1.1_fastlifepriors) #model has quite a large residual standard deviation (sigma)--> need to identify why? This effects the cohens D estimates.
saveRDS(mod1.1_fastlifepriors, file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_fastlifepriors.rda")
mod1.1_fastlifepriors<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/mod1.1_fastlifepriors.rda")

#Posterior draws
fastpaceoflifeprior_ecdf<-pp_check(mod1.1_fastlifepriors, ndraws=100, 
                                   type = "ecdf_overlay")+xlim(0, 20)+
  theme_classic()  


## 7.7 Prior sensitivity: LOO comparison across prior specifications -----------
Loo_prior_performance<-LOO(mod1.1_weakpriors, mod1.1_diffusepriors, 
                           mod1.1_moderatepriors, mod1.1_constraintpriors, 
                           mod1.1_fastlifepriors, moment_match = TRUE)
saveRDS(Loo_prior_performance, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/priorfitsummary.rda")
Loo_prior_performance<-readRDS("scripts/model_outputs/Offspring Trait Models/development time/priorfitsummary.rda")

#Priors have almost no say over the posterior...
#...choosing the weakly informative priors as a default


# ==============================================================================
# 8. SELECTED MODEL DIAGNOSTICS  (mod1.1_weakpriors)-----------------------------
# ==============================================================================

## LOO-CV --------------------------------------------------------------------
loo_weakpriors <- loo(mod1.1_weakpriors, save_psis = TRUE)
plot(loo_weakpriors)
psis_weakpriors <- loo_weakpriors$psis_object
psis_weights    <- weights(psis_weakpriors)


## Posterior predictive checks -----------------------------------------------

#posterior simulated and empirical mean
weakprior_mean     <- ppc_stat(yrep = posterior_predict(mod1.1_weakpriors),
                               y    = data1$development_time_weeks,
                               stat = "mean") + theme_classic()
#posterior cumulative density
weakprior_ecdf     <- pp_check(mod1.1_weakpriors, ndraws = 100,
                               type = "ecdf_overlay") +
  xlim(0, 25) + theme_classic()

#posterior probability density
weakprior_pdf      <- pp_check(mod1.1_weakpriors, ndraws = 100) +
  xlim(0, 25) + theme_classic()

#posterior LOO Q-Q plot
weakprior_loo_qq   <- pp_check(mod1.1_weakpriors, type = "loo_pit_qq",
                               ndraws = 100) + theme_classic()
#LOO-PIT values fall well along the diagonal line (i.e. the expected uniform distribution of PIT integrals)

#LOO uniformity plot
weakprior_loo_unif <- pp_check(mod1.1_weakpriors, type = "loo_pit_overlay",
                               ndraws = 100) + theme_classic()
#LOO predictive interval plots
weakprior_intervals <- ppc_loo_intervals(
  y    = data1$development_time_weeks,
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

ggsave("./bayesian_plots/model fit plots/development time/selectedmodelfit.png",
       plot   = fit_plots,
       device = "png",
       width  = 500, height = 400, units = "mm")


## Bayes R² ------------------------------------------------------------------
bayes_R2(mod1.1_weakpriors, re.form = NA)   # fixed effects only (i.e., marginal)
bayes_R2(mod1.1_weakpriors, re.form = NULL)  # including random effects (i.e., conditional)


## MAP estimates and pd ------------------------------------------------------
MAP_devtime<- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  test        = "pd",
  ci_method   = "HDI",
  centrality  = "MAP",
  component   = "all",
  effects     = "full",
  ci          = 0.95
)
saveRDS(MAP_devtime,
        "scripts/model_outputs/Offspring Trait Models/development time/MAP_devtime.rda")


## ROPE -----------------------------------------------------------------------
rope_devtime<- bayestestR::describe_posterior(
  mod1.1_weakpriors,
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1
)#estimating ROPE seperately as it requires the full posterior distribution

#For Gaussian and lognormal models, the ROPE range is automatically set as 10% of the SD

saveRDS(rope_devtime,
        "scripts/model_outputs/Offspring Trait Models/development time/rope_devtime.rda")

# ===========================================================================
# 9. PRIOR vs. POSTERIOR SENSITIVITY PLOTS------------------------------------
# ===========================================================================

## Combined prior–posterior ECDF panel -----------------------------------------

priorpostplots<-ggarrange(diffuseprior_cumvdis, weakprior_cumvdis, 
                          moderateprior_cumvdis, constraintprior_cumvdis, 
                          fastprior_cumvdis,
                          mod1.1_diffusepriors_cumvdis, weakprior_ecdf1,
                          moderateprior_ecdf, constraintsprior_ecdf,
                          fastpaceoflifeprior_ecdf, nrow = 2, ncol = 5,
                          labels = c("A", "B", "c", "D", "E", "F", "G", "H", "I", "J"))
install.packages("ragg")

library(ragg)
ggsave(filename = "./bayesian_plots/model fit plots/development time/priorvpost_ecdf.png",
       plot = priorpostplots, 
       device = png,
       dpi = 300,
       width = 490, 
       height = 340, 
       units = "mm")

# ===========================================================================
# 10. HYPOTHESIS TESTING: INTERACTION MODELS----------------------------------
# ===========================================================================

#10.1 Three-way interaction: delta age × temperature × offspring sex -------
mod2.1<-brm(bf(development_time_weeks|cens(1-event)~
                 avg_age+ 
                 delta.age +
                 Temp + 
                 F1_sex + 
                 F1_sex*delta.age*Temp+
                 (1+delta.age|PairID),
               sigma~delta.age),
            family=brms::lognormal,
            data=data1,
            #With specified priors added
            prior = weakpriors,
            iter=5000,
            control=list(adapt_delta=0.98),
            core=4)
summary(mod2.1)
saveRDS(mod2.1 , 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod2.1.rda")
mod2.1<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/mod2.1.rda")

## 10.2. Two-way interaction: delta age × temperature (mu only) -----------------
mod3.1<-brm(bf(development_time_weeks|cens(1-event)~
                 avg_age+ 
                 delta.age +
                 Temp + 
                 F1_sex + 
                 delta.age: Temp+
                 (1+delta.age|PairID),
               sigma~delta.age),
            family=brms::lognormal,
            data=data1,
            #With specified priors added
            prior = weakpriors,
            iter=5000,
            control=list(adapt_delta=0.98),
            core=4)
summary(mod3.1) 
saveRDS(mod3.1, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod3.1.rda")
mod3.1<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/mod3.1.rda")

## MAP estimates for interaction of interest model (mod3.1)
MAP_devtime_mod3.1 <- bayestestR::describe_posterior(
  mod3.1,
  test       = "pd",
  ci_method  = "HDI",
  centrality = "MAP",
  effects    = "full",
  ci         = 0.95
)
saveRDS(MAP_devtime_mod3.1,
        "scripts/model_outputs/Offspring Trait Models/development time/MAP_devtime_mod3.1.rda")

#ROPE estimates
devbayestest2<-bayestestR::describe_posterior(
  mod3.1,
  test       = "rope",
  ci_method  = "HDI",
  centrality = "MAP",
  ci         = 1
)
saveRDS(devbayestest2, file = "scripts/model_outputs/Offspring Trait Models/development time/devbayestest2.rda")
devbayestest2<-readRDS("scripts/model_outputs/Offspring Trait Models/development time/devbayestest2.rda")

#10.3 Two-way interaction: delta age × temperature (mu and sigma) ------------
#core interaction of the paper: does temperature mediate the observed parental age effects?
mod3.1_sigmainteraction<-brm(bf(development_time_weeks|cens(1-event)~
                                  avg_age+ 
                                  delta.age +
                                  Temp + 
                                  F1_sex + 
                                  delta.age: Temp+
                                  (1+delta.age|PairID),
                                sigma~delta.age+
                                  Temp+
                                  delta.age:Temp), #testing interaction on sigma as well
                             family=brms::lognormal,
                             data=data1,
                             #With specified priors added
                             prior = weakpriors,
                             iter=5000,
                             control=list(adapt_delta=0.98),
                             core=4)
summary(mod3.1_sigmainteraction) #No evidence for interaction bteween delta age and temperature on sigma
saveRDS(mod3.1_sigmainteraction, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod3.1_sigmainteraction.rda")
mod3.1_sigmainteraction<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/mod3.1_sigmainteraction.rda")

## 10.4 Two-way interaction: delta age × offspring sex (mu only) -------------
#weak evidence for the effect on sigma, restricting inference to mu only
mod4.1<-brm(bf(development_time_weeks|cens(1-event)~
                 avg_age+ 
                 delta.age +
                 Temp + 
                 F1_sex + 
                 delta.age:F1_sex+
                 (1+delta.age|PairID),
               sigma~delta.age),
            family=brms::lognormal,
            data=data1,
            #With specified priors added
            prior = weakpriors,
            iter=5000,
            control=list(adapt_delta=0.98),
            core=4)
summary(mod4.1)
saveRDS(mod4.1, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod4.1.rda")
mod4.1<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/mod4.1.rda")


## 10.5 Two-way interaction: delta age × offspring sex (mu and sigma) ----------
#interaction asks: do parental age effects differ between the offspring sexes?
mod4.1_sigmainteraction<-brm(bf(development_time_weeks|cens(1-event)~
                                  avg_age+ 
                                  delta.age +
                                  Temp + 
                                  F1_sex + 
                                  delta.age:F1_sex+
                                  (1+delta.age|PairID),
                                sigma~delta.age+
                                  F1_sex+
                                  delta.age:F1_sex),
                             family=brms::lognormal,
                             data=data1,
                             #With specified priors added
                             prior = weakpriors,
                             iter=5000,
                             control=list(adapt_delta=0.98),
                             core=4)
summary(mod4.1_sigmainteraction) 
saveRDS(mod4.1_sigmainteraction, 
        file = "scripts/model_outputs/Offspring Trait Models/development time/mod4.1_sigmainteraction.rda")
mod4.1_sigmainteraction<-readRDS(file = "scripts/model_outputs/Offspring Trait Models/development time/mod4.1_sigmainteraction.rda")



# ===========================================================================
# 11. MODEL SELECTION----------------------------------------------------------
# ===========================================================================

#selecting the model with the lowest ELPD value (e.g. the best fitting model)###
hypothesis_fit<-loo(mod1.1_weakpriors,mod2.1,
                    mod3.1,mod4.1, mod4.1_sigmainteraction,
                    mod3.1_sigmainteraction, moment.match = TRUE) 

saveRDS(hypothesis_fit, file = "scripts/model_outputs/Offspring Trait Models/development time/hypothesis_fit.rda")

#Selecting model 1: No interactions --> Just single effect predictors
#still report the ROPE and the MAP for the delta age and temperature interaction (Answers a core question of the paper!)

# ===========================================================================
# 12. POSTERIOR DISTRIBUTION PLOTS  (mod1.1_weakpriors)-----------------------
# ===========================================================================
# Plots styled using stat_halfeye (ggdist) for consistency with other traits
post1 <- as_draws_df(mod1.1_weakpriors)
post2<-as_draws_df(mod3.1)
post2.2<-as_draws_df(mod3.1_sigmainteraction) #for effects on sigma
post3<-as_draws_df(mod4.1)
post3.2<-as_draws_df(mod4.1_sigmainteraction) #for effects on sigma

## Extract posterior draws ---------------------------------------------------
posterior_df_mod1.1 <-
  data.frame(
    "μ: Parents' Δage"              = post1$b_delta.age,
    "μ: Parents' average age"       = post1$b_avg_age,
    "μ: Parent Temperature (28.0°C)"       = post1$b_Temp28,
    "μ: Parent Temperature (30.5°C)"       = post1$b_Temp30.5,
    "μ: Offspring sex (Male)"       = post1$b_F1_sexM,
    "μ: Δage × Temperature (28.0°C)"            = post2$`b_delta.age:Temp28`,
    "μ: Δage × Temperature (30.5°C)"            = post2$`b_delta.age:Temp30.5`,
    "μ: Δage × Offspring sex (Male)"            = post3$`b_delta.age:F1_sexM`,
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
  "μ: Δage × Temperature (28.0°C)",
  "μ: Δage × Temperature (30.5°C)",
  "μ: Δage × Offspring sex (Male)",
  "σ: Parents' Δage"
))



## Halfeye posterior plot ----------------------------------------------------
posterior_plot_mod1.1 <- ggplot(
  posterior_df_mod1.1,
  aes(x = value, y = parameter, fill = parameter)
) +
  annotate("rect",
           xmin = -0.015, xmax = 0.015,
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
  scale_fill_manual(values = c("#4A6479", "#A9BAC2", "#5f9289", "#e653a1",
                               "#f16122", "#d4a5a5", "#e5a9f5", "#d2e9f5",
                               "#8C2F4B", "#9b7fc7"))+
  xlim(c(-0.2,0.2))+
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 30)
  ) +
  labs(x = "Fixed-effect posterior estimates\n(on μ = mean development time (log-weeks); on σ = residual SD, log scale)",
       y = NULL)

ggsave("./bayesian_plots/posterior plots/development time/001_hypothesis1_halfeye.png",
       plot   = posterior_plot_mod1.1,
       device = "png",
       width  = 580, height = 400, units = "mm")


## Selective disappearance test -----------------------------------------------
hypothesis(mod1.1_weakpriors, "avg_age - delta.age > 0")
hypothesis(mod1.1_weakpriors, "avg_age - delta.age < 0")

## Marginal means via emmeans ------------------------------------------------
pairwise_estimates<- emmeans(mod1.1_weakpriors, "delta.age", type = "response", at = list(delta.age = c(-1.96, 1.96)))
pairs(pairwise_estimates) 


## Post-hoc standardised mean difference (SMD) --------------------------------
# Compares extreme ends of within-individual parental age; standardised by
# observed SD of F1 development time. Descriptive only — not used for inference.

pairwise_estimates<- emmeans(mod1.1_weakpriors, "delta.age", type = "response", at = list(delta.age = c(-1.96, 1.96)))
#This is now converted to the raw scale where the extreme ends of parental age are compared
sd_life<-sd(data1$development_time)

#standardised emmeans --> estimating difference across the entire range of parental age
#First getting rid of the link transformation
emm_life_exp <- summary(pairwise_estimates) %>%
  mutate(emmean = exp(emmean),
         lower.HPD = exp(lower.HPD),
         upper.HPD= exp(upper.HPD))

#Then standardising by SD
emm_life_scaled <- emm_life_exp %>%
  mutate(emmean_sd = emmean/ sd_life,
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
dev_time_SMD<-data.frame(
  SMD = SMD,
  SMD_lower = SMD_lower,
  SMD_upper = SMD_upper
)

#Saving for use in combined forest plot of all traits
saveRDS(dev_time_SMD, 
        file = "scripts/model_outputs/Offspring Trait Models/dev_time_SMD")
dev_time_SMD<-readRDS("scripts/model_outputs/Offspring Trait Models/dev_time_SMD")



# ===========================================================================
# 13. MODEL PREDICTIONS AND RESULTS PLOTS-------------------------------------
# ===========================================================================

## 13.1 Prediction grids: within- and between-individual age effects ---------

# Within-individual: vary delta.age, hold avg_age constant
df_predict_within<-expand.grid(
  avg_age=mean(data1$avg_age),
  delta.age=unique(data1$delta.age),
  Temp=unique(data1$Temp),
  F1_sex=unique(data1$F1_sex),
  PairID=unique(data1$PairID[1]))


#calculating model predictions and extracting standard errors
pred_within <- fitted(mod1.1_weakpriors, df_predict_within, re_formula=NA)
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


# 13.2. Between-individual: vary avg_age, hold delta.age constant
df_predict_between<-expand.grid(
  avg_age=unique(data1$avg_age),
  delta.age=mean(data1$delta.age, na.rm =TRUE),
  Temp=unique(data1$Temp),
  F1_sex=unique(data1$F1_sex),
  PairID=unique(data1$PairID[1]))

df_predict_between$F1_sex <- factor(df_predict_between$F1_sex, levels = c("F", "M"))

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

## 13.3 Raw within-individual means (for overlaying on plot) ----------------
breaks <- seq(-4.5, 4.5, by = 1)  # creates bins: [-3.5, -2.5], [-2.5, -1.5], ..., [2.5, 3.5]
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
    n = sum(!is.na(development_time_weeks)),                       # number of non-missing observations
    mean_dev = mean(development_time_weeks, na.rm = TRUE),        # group mean
    se_dev = ifelse(n > 1, sd(development_time_weeks, na.rm = TRUE)/sqrt(n), NA_real_), # SE only if n>1
    .groups = "drop"
  )

## 13.4 Main parental age effect plot ----------------------------------------
development_plot <- ggplot(data = data1,
                           aes(x = within_subject_age, 
                               y = development_time_weeks,
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
                 aes(y=mean_dev,
                     ymin = mean_dev-se_dev,
                     ymax = mean_dev+se_dev,
                     colour="Within-individual"),
                 linewidth = 3,
                 show.legend = FALSE)+
  geom_point(data = raw_deltaage,
             aes(x = within_subject_age,
                 y = mean_dev,
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
  )
  scale_x_continuous(breaks=c(-4, -3, -2, -1, 0, 1,2,3, 4)) +
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
  labs(x = "Parents' adult age at reproduction (weeks; mean-centred)", y = "Offspring's development time (weeks)")


# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 development time/001_development_plot.png",
       plot = development_plot,
       bg="transparent",
       device = "png", 
       width = 480, 
       height = 700, 
       units = "mm")



## 13.4 Cumulative incidence function for F1 development time -----------------

lognormal_incidence <- function(mu, sigma, time) {
  plnorm(time, meanlog=mu, sdlog=sigma)
}

lognormal_pdf <- function(mu, sigma, time) {
  dlnorm(time, meanlog=mu, sdlog=sigma)
} #Calculates the probability density function based on the models mu and sigma parameters

lognormal_hazard <- function(mu, sigma, time) {
  lognormal_pdf(mu, sigma, time) / lognormal_survival(mu, sigma, time)
} #Instantaneous hazard rate = f(t)/S(t)

#Proxy data
newdat <- expand.grid(
  avg_age = mean(data1$avg_age),                         
  delta.age = c(-1.96, 0, 1.96), 
  Temp = unique(data1$Temp),                          
  F1_sex = unique(data1$F1_sex),
  PairID = unique(data1$PairID)[1]
)

# Create a sequence over which to evaluate
dev_seq <- seq(min(4), 
               max(data1$development_time_weeks, na.rm=TRUE), 
               length.out = 250)

#Creating the incidendce curve for F1 devtime
proportional_posterior <- mod1.1_weakpriors %>% 
  linpred_draws(
    newdat,
    value = "mu",
    allow_new_levels = TRUE, 
    transform = TRUE,
    re_formula = NA,  # population-level predictions
    dpar = "sigma",
    ndraws = 3000, 
    seed = 123
  ) 

#calculating S across posterior draws
posterior_draws<-proportional_posterior%>% 
  ungroup() %>%
  crossing(x = dev_seq) %>%   
  mutate(Fx = lognormal_incidence(mu, sigma, x)) %>% 
  group_by(x, .draw, delta.age) %>% 
  summarise(Fx=median(Fx),.groups = "drop")



#Incidence function---------------------------------------------------
summary_incidence<- posterior_draws%>%
  group_by(x, delta.age) %>%
  summarise(
    median = median(Fx),  
    lower  = quantile(Fx, 0.025),
    upper  = quantile(Fx, 0.975),
    .groups = "drop"
  )

#converting delta.age to a categorical variable for plotting
summary_incidence<- summary_incidence%>%
  mutate(timepoint= case_when(
    delta.age %in% -1.96 ~ "Early-Aged",
    delta.age %in% 0 ~ "Middle-Aged",
    delta.age %in% 1.96 ~ "Late-Aged"
  ))

summary_incidence$timepoint<-factor(summary_incidence$timepoint, 
                                    levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))


#Getting the median development times
median_development_draws<- mod1.1_weakpriors %>% 
  linpred_draws(
    newdat,
    value = "mu",
    allow_new_levels=TRUE,
    transform = TRUE,
    re_formula=NA,#population level predictions
    dpar = "sigma",
    ndraws = 5000, 
    seed = 123
  ) %>% 
  ungroup() %>%
  group_by(.draw, delta.age) %>% 
  summarise(median_dev= exp(median(mu))) %>% 
  group_by(delta.age) %>% 
  median_hdi(median_dev)

#converting delta.age
median_development<- median_development_draws%>%
  mutate(timepoint= case_when(
    delta.age %in% -1.96 ~ "Early-Aged",
    delta.age %in% 0 ~ "Middle-Aged",
    delta.age %in% 1.96 ~ "Late-Aged"
  ))


#Cumulative incidence curve-----------------------------------------------------

#Plotting the survival curves for each group --> entire range of development time
lognormaldevelopmenttime<-ggplot(data=summary_incidence,
                                 aes(x = x, 
                                     y = median, 
                                     colour= timepoint)) +
  
  geom_line(data = summary_incidence, 
            aes(x = x, 
                y = median,
                colour = timepoint),
            linewidth=4, 
            alpha=1)+
  geom_ribbon(data = summary_incidence, 
              aes(y=NULL, 
                  ymin = lower, 
                  ymax = upper, 
                  color= timepoint,
                  fill=timepoint,
                  alpha = timepoint),
              linetype="dashed",
              linewidth =1)+
  geom_segment(data = median_development,
               aes(x = 4, xend = median_dev, y = 0.5, yend = 0.5, color = timepoint),
               linetype = "dashed", linewidth = 3, alpha = 0.9) +
  geom_segment(data = median_development,
               aes(x = median_dev, xend = median_dev, y = 0, yend = 0.5, color = timepoint),
               linetype = "dashed", linewidth = 3, alpha = 0.9)+
  theme_classic()+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50, face="bold"),
        panel.background = element_rect(fill = "transparent", color = NA),  # Make panel background transparent
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), # Make plot background transparent
  ) +
  scale_x_continuous(limits = c(4, 15), breaks = c(4, 6, 8, 10, 12, 14))+
  scale_color_manual(name = "Parents' adult age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' adult age at reproduction", values = c(0.6,0.1, 0.6))+
  scale_fill_manual(name = "Parents' adult age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  theme(
    legend.position = "top", 
    legend.box = "horizontal"
  ) +
  labs(x = "Offspring age (weeks)", y = "Cumulative incidence of adult emergence, F(x)")

# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/F1 development time/002_incidenceofadultemergence.png",
       plot = lognormaldevelopmenttime, 
       bg="transparent",
       device = "png", 
       width = 500, 
       height = 400, 
       units = "mm")


#Combining the mu and scale plot into one (for the results section)-------------

#Saving the prior and poosterior plots for default versus moderate priors
inference4<-ggarrange(
  development_plot,
  lognormaldevelopmenttime,
  label.x = c(0.05, 0.05),
  ncol = 2,
  nrow = 1,
  labels = c("A", "B"),
  font.label = list(size = 50, face = "bold"),
  widths = c(0.7, 1) 
)

ggsave(filename = "./bayesian_plots/offspring trait plots/F1 development time/inference_plots.png",
       plot = inference4, 
       device = "png", 
       width = 990, 
       height = 570, 
       units = "mm")


#Creating a table to export------------------------------------------------------------------------------
#For the subset of animals that we know became adults

# ── helper ────────────────────────────────────────────────────────────────────
#Model draws from selected models
base_model_draws<-as_draws_df(mod1.1_weakpriors) #just single effects, and the final selected model
sigma_dropped_terms<-as_draws_df(mod1.1_sigma1)
interaction_draws<-as_draws_df(mod2.1) 
interaction_draws2<-as_draws_df(mod3.1) 
interaction_draws2.1<-as_draws_df(mod3.1_sigmainteraction) 
interaction_draws3<-as_draws_df(mod4.1)
interaction_draws3.1<-as_draws_df(mod4.1_sigmainteraction)


#estimating the difference between average age and delta age
differences<-data.frame(
  selectivedis = base_model_draws$b_avg_age - base_model_draws$b_delta.age
)

SDy<-sd(log(data1$development_time_weeks))*0.1 #(ROPE range should be 0.015)

#function for generating estimates
summarise_param <- function(draws, param, section,
                            rope = TRUE, #using default ROPE range automatically calculated from bayestest
                            rope_range = c(-0.015, +0.015),
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
  summarise_param(base_model_draws$sd_PairID__delta.age, "σ slope Δage",    "Random effects (Location)", pd = FALSE),
  summarise_param(base_model_draws$cor_PairID__Intercept__delta.age, "r intercept ~ slope Δage",    "Random effects (Location)", pd = TRUE)
)

#---4. Sigma parameter──────────────────────────────────────────────────
sigma<- bind_rows(
  summarise_param(base_model_draws$b_sigma_Intercept, "Sigma intercept (σ)", "Distributional parameters (scale component)",  rope = FALSE),
  summarise_param(base_model_draws$b_sigma_delta.age, "Parents' Δage (z-scaled)", "Distributional parameters (scale component)", rope= FALSE)
)

#---5. Sigma dropped terms----------------------------------------------
sigma_dropped<-bind_rows(
  summarise_param(sigma_dropped_terms$b_sigma_Temp28, "Temperature: 28.0°C", "Distributional parameters (Dropped terms)", rope = FALSE),
  summarise_param(sigma_dropped_terms$b_sigma_Temp30.5, "Temperature: 30.5°C", "Distributional parameters (Dropped terms)", rope = FALSE),
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
bind_rows(fe, fe_interaction2, re_sd_explore, sigma, sigma_dropped, bayes) %>%
  gt(groupname_col = "Section") %>%
  cols_label(Parameter = "Parameter", MAP = "MAP",
             `95% HDI` = "95% HDI", `% in ROPE` = "% in ROPE", pd = "pd") %>%
  tab_header(
    title    = "Development time model summary: Lognormal model"
  ) %>%
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) %>%
  tab_footnote("ROPE = [−0.015, 0.015] on log scale; applied to location fixed effects only.",
               locations = cells_column_labels(`% in ROPE`)) %>%
  tab_footnote("pd = probability of direction.",
               locations = cells_column_labels(pd)) %>%
  opt_stylize(style = 1)%>%
  gtsave("./tables/developmenttime_updated.docx")


###############################################################################
##  END OF SCRIPT--------------------------------------------------------------------
###############################################################################
