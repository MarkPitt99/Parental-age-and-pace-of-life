###############################################################################
##  Script: BaSTA analysis looking at the total lifespan of the offspring, 
#Analysis aims to assess how parental age affects different mortality parameters of the offspring (i.e., baseline and age-specific mortality)
#This script removes early-life mortality (as this period is not well-captured by the bathtub model)
###############################################################################

# ===========================================================================
# 1. SETUP--------------------------------------------------------------------
# ===========================================================================

## Libraries ------------------------------------------------------------------
library(dplyr)  #v.1.1.4
library(tidyr) #v.1.3.1
library(BaSTA) #v.2.0.2
library(snowfall) #v.1.84.6.3
library(ggplot2) #v.4.0.0
library(survival) #v.3.8.3
library(survminer) #v.0.5.1
library(ggpubr) #v.0.6.1
library(stringr) #v.1.5.1
library(cowplot) #v.1.1.3

# ===========================================================================
# 2. DATA ORGANISATION--------------------------------------------------------
# ===========================================================================

#2.1. Importing data------------
F1data<-readRDS("./raw data/BaSTA data/total_basta_data.RDS")
F1data<-as.data.frame(F1data)
length(unique(F1data$ID))
#987 offspring from 77 parent pairs that survived past the first four weeks of life


#2.2. grouping the parental ages into categorical groups (rather than treating it as a continuous variable)----
F1data <- F1data%>%
  mutate(Timepoint_binned = case_when(
    Timepoint %in% 1:2 ~ "early",
    Timepoint %in% 3:5 ~ "middle",
    Timepoint %in% 6:8 ~ "late"
  ))

#2.3. creating dataset with the required covariates-------------
totcovdata<- F1data %>% 
  select(ID, Birth.Date, Min.Birth.Date, Max.Birth.Date, Entry.Date, Depart.Date, Depart.Type, Timepoint, Timepoint_binned, Temp,F1_sex)

#DATA CHECK: re-checking whether the filtered data passes BaSTA's built-in data check function
checkedDataCens <- DataCheck(object = totcovdata, dataType = "census",
                             silent = FALSE)#No inconsistencies between dates
print(checkedDataCens)#data seems to all be coded correctly

#How many offspring do we have in each age category at each temperature?
table_temp_lifespan <- table(F1data$Temp, F1data$Timepoint_binned)
print(table_temp_lifespan)#this table gives counts for each level of age and temperature
#       early late middle
#25.5    88  125    147
#28      78   91    121
#30.5    99   88    150

#How many parents produced offspring in each age category? Irrespective of temperature
table_lifespan_pairID<-F1data %>% 
  group_by(Timepoint_binned) %>% 
  summarise(n_parents = n_distinct(PairID))
print(table_lifespan_pairID)

# Timepoint_binned n_parents
#early                   75
#late                    60
#middle                  72

# ===========================================================================
# 3. SUMMARY STATISTICS-------------------------------------------------------
# ===========================================================================
#creating a total lifespan variable
totcovdata$Lifespan_weeks <- as.numeric(difftime(totcovdata$Depart.Date, totcovdata$Birth.Date, units = "weeks"))

#Raw observed mean lifespan
sum_dat1<-totcovdata %>% 
  summarise(F1_lifespan = median (Lifespan_weeks),
            sd_total= sd(Lifespan_weeks),
            n_animals=n())
#median lifespan    sd_total    n_animals
#18.57143         3.965274       987

#Median lifespan of the offspring from each parental age class
sum_dat<-totcovdata %>%
  group_by(Timepoint_binned)%>%
  summarise(F1_lifespan = median (Lifespan_weeks),
            sd_total= sd(Lifespan_weeks),
            n_animals=n())
#Timepoint_binned F1_lifespan sd_total n_animals

#early                   18.6     3.96       265
#late                    18.7     4.10       304
#middle                  18.1     3.85       418

#median lifespan per temperature treatment
sum_dat2<-totcovdata %>%
  group_by(Temp)%>%
  summarise(F1_lifespan = median (Lifespan_weeks),
            sd_total= sd(Lifespan_weeks),
            n_animals=n(),)

#Temp  F1_lifespan sd_total n_animals
#25.5         18.6     3.97       360
#28           18.6     3.91       290
#30.5         18.4     4.01       337

# ===========================================================================
# 4.TESTING MORTALITY DISTRIBUTIONS----------------------------------------
# ===========================================================================
#Using a series of null models with no covariates
#Deciding on which parametric survival models fit our survival distribution the best

#4.1. Gompertz Model------------------------------
gomptzmod <- basta(object = totcovdata, 
                        dataType = "census",
                        niter=60000,
                        burnin=1001,
                        thinning=50,
                        nsim = 4, parallel = TRUE, ncpus = 4)
summary(gomptzmod)

#Saving the model
saveRDS(gomptzmod, file = "scripts/model_outputs/BaSTA/total lifespan/lifespan_gompertz.rda")
gomptzmod<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/lifespan_gompertz.rda")

#Assessing model fit and plotting the model outputs
plot(gomptzmod)
plot(gomptzmod , plot.type = "demorates")
plot(gomptzmod, plot.type = "gof")
plot(gomptzmod, densities=TRUE)


#4.2 Gompertz model with a Makeham shape ----
#adds a Makeham constant term (accounts for age-independent background mortality)
gomptzmod.2 <- basta(object = totcovdata, 
                          dataType = "census",
                          shape="Makeham",
                          niter=60000,
                          burnin=1001,
                          thinning=50,
                          nsim = 4, parallel = TRUE, ncpus = 4)
summary(gomptzmod.2)

#Saving the model
saveRDS(gomptzmod.2, file = "scripts/model_outputs/BaSTA/total lifespan/gompertz_makeham.rda")
gomptzmod.2<-readRDS("scripts/model_outputs/BaSTA/total lifespan/gompertz_makeham.rda")

#assessing model fit
plot(gomptzmod.2, plot.type = "gof")
plot(gomptzmod.2, plot.type = "demorates")
plot(gomptzmod.2, densities=TRUE)

#4.3. Gompertz Model with a bathtub function-------------------------
#adds two terms to capture the rate of early-life mortality--
gomp.3 <- basta(object = totcovdata, 
                     dataType = "census",
                     shape="bathtub",
                     niter=60000,#trying to double the burnin and the number of iterations
                     burnin=1001,
                     thinning=50,
                     nsim = 4, parallel = TRUE, ncpus = 4)
summary(gomp.3)

#Saving the model
saveRDS(gomp.3, file = "scripts/model_outputs/BaSTA/total lifespan/gompertz_bathtub.rda")
gomp.3 <-readRDS("scripts/model_outputs/BaSTA/total lifespan/gompertz_bathtub.rda")

#Assessing model fit and plotting the model outputs
plot(gomp.3, plot.type = "gof")
plot(gomp.3, plot.type = "demorates")
plot(gomp.3, densities=TRUE)


#4.4. Weibull model----------------------------------------
bastaweibull <- basta(object = totcovdata, 
                           dataType = "census", 
                           model="WE",
                           niter=60000,
                           burnin=1001,
                           thinning=50,
                           nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastaweibull)

#Saving the model
saveRDS(bastaweibull, file = "scripts/model_outputs/BaSTA/total lifespan/lifespan_weibull.rda")
bastaweibull<-readRDS("scripts/model_outputs/BaSTA/total lifespan/lifespan_weibull.rda")

#Assessing model fit and plotting the model outputs
plot(bastaweibull, plot.type = "gof")
plot(bastaweibull, plot.type = "demorates")
plot(bastaweibull, densities=TRUE)

#4.5. Weibull model with a Makeham term --------------------------
bastaweibull.2 <- basta(object = totcovdata, 
                             dataType = "census",
                             model="WE",
                             shape="Makeham",
                             niter=60000,#trying to double the burnin and the number of iterations
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastaweibull.2)

#Saving the model
saveRDS(bastaweibull.2, file = "scripts/model_outputs/BaSTA/total lifespan/weibull_makeham.rda")
bastaweibull.2<-readRDS("scripts/model_outputs/BaSTA/total lifespan/weibull_makeham.rda")

#Assessing model fit and plotting the model outputs
plot(bastaweibull.2, plot.type = "gof")
plot(bastaweibull.2, plot.type = "demorates")
plot(bastaweibull.2, densities=TRUE)


#4.6. Weibull model with a bathtub term ----------------------
bastaweibull.3 <- basta(object = totcovdata, 
                             dataType = "census",
                             model="WE",
                             shape="bathtub",
                             niter=60000,
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastaweibull.3)

#saving the model
saveRDS(bastaweibull.3, file = "scripts/model_outputs/BaSTA/total lifespan/weibull_bathtub.rda")
bastaweibull.3<-readRDS("scripts/model_outputs/BaSTA/total lifespan/weibull_bathtub.rda")

#Assessing model fit and plotting the model outputs
plot(bastaweibull.3, plot.type = "gof")
plot(bastaweibull.3, plot.type = "demorates")
plot(bastaweibull.3, densities=TRUE)


#4.7.exponential model-------------------------
bastaexp <- basta(object = totcovdata, 
                       dataType = "census",
                       model="EX",
                       niter=60000,
                       burnin=1001,
                       thinning=50,
                       nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastaexp)
plot(bastaexp, plot.type = "gof")

#saving the model
saveRDS(bastaexp, file = "scripts/model_outputs/BaSTA/total lifespan/lifespan_exp.rda")
bastaexp<-readRDS("scripts/model_outputs/BaSTA/total lifespan/lifespan_exp.rda")

#Assessing model fit and plotting the model outputs
plot(bastaexp, plot.type = "gof")
plot(bastaexp, plot.type = "demorates")
plot(bastaexp, densities=TRUE)

#4.8.logistic model-----------------------------
bastalog <- basta(object = totcovdata, 
                       dataType = "census",
                       model="LO",
                       niter=60000,
                       burnin=1001,
                       thinning=50,
                       nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastalog)

#saving the model
saveRDS(bastalog, file = "scripts/model_outputs/BaSTA/total lifespan/lifespan_log.rda")
bastalog<-readRDS("scripts/model_outputs/BaSTA/total lifespan/lifespan_log.rda")

#Assessing model fit and plotting the model outputs
plot(bastalog, plot.type = "gof")
plot(bastalog, plot.type = "demorates")
plot(bastalog, densities=TRUE)

#4.9. Logistic model with a makeham shape------------------------------
bastalog.2 <- basta(object = totcovdata, 
                         dataType = "census", 
                         model="LO",
                         shape="Makeham",
                         niter=60000,
                         burnin=1001,
                         thinning=50,
                         nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastalog.2)

#saving the model
saveRDS(bastalog.2, file = "scripts/model_outputs/BaSTA/total lifespan/lifespan_log_makeham.rda")
bastalog.2<-readRDS("scripts/model_outputs/BaSTA/total lifespan/lifespan_log_makeham.rda")

#Model fit and convergence plots
plot(bastalog.2, plot.type = "gof")
plot(bastalog.2, plot.type = "demorates")
plot(bastalog.2, densities=TRUE)


#4.10. Logistic model with a bathtub shape--------------------------------
bastalog.3 <- basta(object = totcovdata, 
                         dataType = "census",
                         model="LO",
                         shape="bathtub",
                         niter=60000,
                         burnin=1001,
                         thinning=50,
                         nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastalog.3)

#saving the model
saveRDS(bastalog.3, file = "scripts/model_outputs/BaSTA/total lifespan/lifespan_log_bathtub.rda")
bastalog.3<-readRDS("scripts/model_outputs/BaSTA/total lifespan/lifespan_log_bathtub.rda")

#Model fit and convergence plots
plot(bastalog.3, plot.type = "gof")
plot(bastalog.3, plot.type = "demorates")
plot(bastalog.3, densities=TRUE)
#The Weibull Makeham model has the best fit, selecting this moving forward

# ===========================================================================
# 5. PRIOR SPECIFICATION FOR THE NULL MODEL----------------------------------
# ===========================================================================
#Reminder that the null model is fitted with no covariates

#-----------NULL MODEL---------------------------------------------
mean(totcovdata$Lifespan_weeks) #18.30598 weeks 
sd(totcovdata$Lifespan_weeks) #3.965274 weeks

#Converting mean years into the characteristic lifespan (b1) 
#BaSTA estimates b1 in years, and inverts this estimate
#ensuring we dop the same when setting the prior on this parameter
mu_years<-18.30598/52.1429
prior_scale_mu<-gamma(1+1/1.5)/mu_years
#Gives an empirical prior value for b1 of 2.571387, using this to set the prior SD

#----------------------SETTING MODEL PRIORS-------------------
#Weakly informative priors-----------------
weakMean2 <- matrix(c(
  0, 1.5, prior_scale_mu
), nrow = 1, byrow = TRUE)
#mean for makeham (c), shape (b0), and scale (b1) parameter, respectively

weakSd2 <- matrix(c(
  1.0, 1.0, 0.5
), nrow = 1, byrow = TRUE)
#sd for makeham (c), shape (b0), and scale (b1) parameter, respectively

weakLower2<-matrix(c(
  0, 0, 0
), nrow = 1, byrow = TRUE)
#setting lower bounds of each prior to zero (to respect Weibull parameterisation)

#Moderate null priors------------------------
moderateMean2 <- matrix(c(
  0, 1.5, prior_scale_mu
), nrow = 1, byrow = TRUE)
#mean for makeham (c), shape (b0), and scale (b1) parameter, respectively

moderateSd2 <- matrix(c(
  1.0, 0.5, 0.25
), nrow = 1, byrow = TRUE)
#sd for makeham (c), shape (b0), and scale (b1) parameter, respectively

moderateLower2<-matrix(c(
  0, 0, 0
), nrow = 1, byrow = TRUE)
#setting lower bounds of each prior to zero (to respect Weibull parameterisation)


# ===========================================================================
# 6. PRIOR CHECKS---------------------------------------------------
# ===========================================================================

#6.1. Model with weak priors----------------------
Nullbastaweibull <- basta(object = totcovdata, 
                          dataType = "census",
                          model="WE",
                          shape="Makeham",
                          thetaPriorMean = weakMean2, #adding the weakpriors
                          thetaPriorSd = weakSd2,
                          thetaPriorLower = weakLower2,
                          niter=80000,
                          burnin=1001,
                          thinning=50,
                          nsim = 4, parallel = TRUE, ncpus = 4)
summary(Nullbastaweibull)
saveRDS(Nullbastaweibull, file = "scripts/model_outputs/BaSTA/total lifespan/Null_model_weibull.rda")
Nullbastaweibull<-readRDS("scripts/model_outputs/BaSTA/total lifespan/Null_model_weibull.rda")

#assessing how the priors affect model fit
plot(Nullbastaweibull, plot.type = "gof")

#6.2. Null model with moderate priors-------------------------------
Nullbastaweibull.moderate <- basta(object = totcovdata, 
                                   dataType = "census",
                                   model="WE",
                                   shape="Makeham",
                                   thetaPriorMean = moderateMean2, #adding the moderate priors
                                   thetaPriorSd = moderateSd2,
                                   thetaPriorLower = moderateLower2,
                                   niter=80000,
                                   burnin=1001,
                                   thinning=50,
                                   nsim = 4, parallel = TRUE, ncpus = 4)
summary(Nullbastaweibull.moderate)
saveRDS(Nullbastaweibull.moderate, file = "scripts/model_outputs/BaSTA/total lifespan/Nullbastaweibull.moderate.rda")
Nullbastaweibull.moderate<-readRDS("scripts/model_outputs/BaSTA/total lifespan/Nullbastaweibull.moderate.rda")

#assesing model fit
plot(Nullbastaweibull.moderate)
plot(Nullbastaweibull.moderate, plot.type = "gof")


# ===========================================================================
# 7. ADDING COVARIATES: ISOLATED EFFECT OF PARENTAL AGE----------------------
# ===========================================================================

#7.1. Plotting a Kaplein-Maier curve of offspring observed survival in each parental age category--------------

#Ensuring variables are correctly ordered
totcovdata$Timepoint_binned<-factor(totcovdata$Timepoint_binned, 
                                    levels = c("early", "middle", "late"))
#censoring variable for survival package
totcovdata$event<-ifelse(totcovdata$Depart.Type =="D", 1,0)

#Renaming levels for the plot
totcovdata$Timepoint_binnedplot<-factor(totcovdata$Timepoint_binned, 
                                        levels = c("early", "middle", "late"),
                                        labels=c("Early-Aged", "Middle-Aged", "Late-Aged"))

#1.Kaplain-Maier curve for parental age
#survival estimates
km_fit_parentalage <- survfit(Surv(Lifespan_weeks,event) ~ Timepoint_binnedplot, data = totcovdata)
#Kaplan-Meier curve
Km_survival_parentalage <- ggsurvplot(km_fit_parentalage, 
                                      data = totcovdata, 
                                      pval = FALSE,  
                                      conf.int = TRUE,  
                                      risk.table = FALSE,
                                      break.time.by = 5,
                                      surv.median.line = "hv",  
                                      xlab = "Offspring age (weeks)", 
                                      ylab = "Cumulative survival probability, S(x)",
                                      font.x = 35,
                                      font.y = 35,
                                      font.tickslab = c(35, "grey25"),
                                      font.legend = 25,
                                      risk.table.fontsize = 6,
                                      legend.title = "Parents' age at reproduction",
                                      legend.labs = c("Early-Aged", "Middle-Aged", "Late-Aged"),
                                      palette = c("#f96161", "#66b2b2", "#066594"))


#7.2. PARAMETRIC SURVIVAL MODELS FOR PARENTAL AGE EFFECTS------------------------------------------

#7.2.1. Weibull model------------------------------------------------------
totalbastaweibull <- basta(object = totcovdata, 
                           dataType = "census",
                           formulaMort=~Timepoint_binned-1, 
                           model="WE",
                           niter=60000,
                           burnin=1001,
                           thinning=50,
                           nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibull)
saveRDS(totalbastaweibull, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_weibull.rda")
totalbastaweibull<-readRDS("scripts/model_outputs/BaSTA/total lifespan/total_lifespan_weibull.rda")

#Assessing model fit and plotting the model outputs
plot(totalbastaweibull, plot.type = "gof")
plot(totalbastaweibull, plot.type = "demorates")
plot(totalbastaweibull, densities=TRUE)


#7.2.2. Weibull model with a Makeham term --------------------------------
totalbastaweibull.2 <- basta(object = totcovdata, 
                        dataType = "census",
                        formulaMort=~Timepoint_binned-1, 
                        model="WE",
                        shape="Makeham",
                        niter=60000,
                        burnin=1001,
                        thinning=50,
                        nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibull.2)
plot(totalbastaweibull.2)
saveRDS(totalbastaweibull.2, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_weibull_makeham.rda")
totalbastaweibull.2<-readRDS("scripts/model_outputs/BaSTA/total lifespan/total_lifespan_weibull_makeham.rda")

#Assessing model fit and plotting the model outputs
plot(totalbastaweibull.2, plot.type = "gof")
plot(totalbastaweibull.2, plot.type = "demorates")
plot(totalbastaweibull.2, densities=TRUE)


#7.2.3. Weibull model with a bathtub term --------------------------------------------------
totalbastaweibull.3 <- basta(object = totcovdata, 
                             dataType = "census",
                             formulaMort=~Timepoint_binned-1, 
                             model="WE",
                             shape="bathtub",
                             niter=60000,#trying to double the burnin and the number of iterations
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibull.3)
saveRDS(totalbastaweibull.3, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_weibull_bathtub.rda")
totalbastaweibull.3<-readRDS("scripts/model_outputs/BaSTA/total lifespan/total_lifespan_weibull_bathtub.rda")

#Assessing model fit and plotting the model outputs
plot(totalbastaweibull.3, plot.type = "gof")
plot(totalbastaweibull.3, plot.type = "demorates")
plot(totalbastaweibull.3, densities=TRUE)

#Selecting the Weibull Makeham model

# ===========================================================================
# 8. PRIOR SENSITIVITY ANALYSIS FOR THE PARENTAL AGE MODEL-------------------
# ===========================================================================

#-------------------Setting regularising priors-----------------------------

#Weak priors-------------------------------------------------

#Weakly informative priors (mean)
weakMean <- matrix(c(
  0, 1.5, prior_scale_mu
), nrow = 3, ncol =3, byrow = TRUE)
#Mean for the c, b0, and b1 terms, respectively

#weakly infromative priors (SD)
weakSd <- matrix(c(
  1.0, 1.0, 0.5
), nrow = 3, ncol = 3, byrow = TRUE)
#SD for the c, b0, and b1 terms, respectively

weakLower<-matrix(c(
  0, 0, 0
), nrow = 3, ncol=3, byrow = TRUE)
#lower bounds for the c, b0, and b1 terms, respectively


#---------------Moderately informative priors-----------------------------------

#moderate priors (mean)
moderateMean <- matrix(c(
  0, 1.5, prior_scale_mu
), nrow = 3, ncol=3, byrow = TRUE)
#Mean for the c, b0, and b1 terms, respectively

#moderate priors (SD)
moderateSd <- matrix(c(
  0.5, 0.5, 0.25
), nrow = 3, ncol=3, byrow = TRUE)
#SD for the c, b0, and b1 terms, respectively

#moderate priors
moderateLower<-matrix(c(
  0, 0, 0
), nrow = 3, ncol =3, byrow = TRUE)
#lower bounds for the c, b0, and b1 terms, respectively


#-----------------8.1. Fitting priors to the parental age model-----------------------------------

#8.1.1. Weibull model with a Makeham term---------------------
totalbastaweibull.weak <- basta(object = totcovdata, 
                                dataType = "census",
                                formulaMort=~Timepoint_binned-1, 
                                model="WE",
                                shape="Makeham",
                                thetaPriorMean = weakMean,
                                thetaPriorSd = weakSd,
                                thetaPriorLower = weakLower,
                                niter=80000,
                                burnin=1001,
                                thinning=50,
                                nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibull.weak)
saveRDS(totalbastaweibull.weak, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_weibull_weak.rda")
totalbastaweibull.weak<-readRDS("scripts/model_outputs/BaSTA/total lifespan/total_lifespan_weibull_weak.rda")

#Assessing model fit and plotting the model outputs
plot(totalbastaweibull.weak, plot.type = "gof")
plot(totalbastaweibull.weak, plot.type = "demorates")
plot(totalbastaweibull.weak, densities=TRUE)
plot(totalbastaweibull.weak, type="fancy")


#-------------Moderate priors---------------------------------------------------

#8.1.2. Weibull model with a Makeham term ----------------------------
totalbastaweibull.moderate<- basta(object = totcovdata, 
                                dataType = "census",
                                formulaMort=~Timepoint_binned-1, 
                                model="WE",
                                shape="Makeham",
                                thetaPriorMean = moderateMean,
                                thetaPriorSd = moderateSd,
                                thetaPriorLower = moderateLower,
                                niter=80000,
                                burnin=1001,
                                thinning=50,
                                nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibull.moderate)
plot(totalbastaweibull.moderate)
saveRDS(totalbastaweibull.moderate, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_weibull_moderate.rda")
totalbastaweibull.moderate<-readRDS("scripts/model_outputs/BaSTA/total lifespan/total_lifespan_weibull_moderate.rda")

#Assessing model fit and plotting the model outputs
plot(totalbastaweibull.moderate, plot.type = "gof")
plot(totalbastaweibull.moderate, plot.type = "demorates")
plot(totalbastaweibull.moderate, densities=TRUE)


#Selected the weakly informative prior for all models moving forward

# ===========================================================================
# 9. ISOLATED EFFECT OF PARENTAL AGE PLOTS-----------------------------------
# ===========================================================================
#Plotting estimates from weibull model with weakly-informative priors

#9.1. MORALITY PARAMETER POSTERIOR DISTRIBUTIONS (I.E., C, B1, AND B0)----

#Extracting the model parameters
total_theta_params <- totalbastaweibull.weak$params
total_theta_params<-as.data.frame(total_theta_params)

# Reshaping the data
total_posterior_df <- total_theta_params %>%
  pivot_longer(cols = everything(), 
               names_to = c("Type", "Parental_Age"),
               names_sep = "\\.", # Splitting at "."
               values_to = "Value") %>%
  mutate(
    Parental_Age = dplyr::recode(Parental_Age, 
                          Timepoint_binnedearly = "Early-Aged",
                          Timepoint_binnedmiddle = "Middle-Aged",
                          Timepoint_binnedlate = "Late-Aged")
  )

#making sure the levels are in  the right order for the plot:
total_posterior_df$Parental_Age<- factor(total_posterior_df$Parental_Age, 
                                   levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))


# Create separate datasets for "c", `b0` and `b1`
totb0_df <- total_posterior_df %>% filter(Type == "b0")
totb1_df <- total_posterior_df %>% filter(Type == "b1") 
totalcdf<-total_posterior_df %>% filter(Type == "c")


#9.2. B0 POSTERIOR---------------------------------------
tot_b0_plot <- ggplot(totb0_df, aes(x = Value, 
                                    fill = Parental_Age, 
                                    color = Parental_Age, 
                                    alpha=Parental_Age)) +
  geom_density(linewidth= 2) +
  theme_classic() +
  scale_fill_manual(values = c("Early-Aged" = "#f96161", 
                               "Middle-Aged" = "#66b2b2", 
                               "Late-Aged" = "#066594")) +
  scale_color_manual(values = c("Early-Aged" = "#f96161", 
                                "Middle-Aged" = "#66b2b2", 
                                "Late-Aged" = "#066594")) +
  scale_alpha_manual(values = c("Early-Aged" = 0.6, 
                                "Middle-Aged" = 0.1, 
                                "Late-Aged" = 0.6))+
  labs(
    x = "b0 parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' age at reproduction",
    color = "Parents' age at reproduction",
    alpha="Parents' age at reproduction")+
  scale_y_continuous(breaks= c(0.0, 0.6, 1.2, 1.80))+
  scale_x_continuous(limits = c(4.5, 6.6), breaks = c(4.6, 5.2, 5.8, 6.4))+
  theme(legend.position = "top",
        axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.title=element_text(size=45),
        legend.text = element_text(size = 45),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), 
  )+
  guides(alpha=guide_legend(order = 1, title.position = "top",title.hjust = 0.5,
                            ,override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) 

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/001_b0_parentalage_totallifespan.png",
       plot = tot_b0_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")


#9.3. B1 POSTERIOR plot----------------------------
tot_b1_plot <- ggplot(totb1_df, aes(x = Value, 
                                    fill = Parental_Age, 
                                    color = Parental_Age, 
                                    alpha=Parental_Age)) +
  geom_density(linewidth= 2) +
  theme_classic() +
  scale_fill_manual(values = c("Early-Aged" = "#f96161", 
                               "Middle-Aged" = "#66b2b2", 
                               "Late-Aged" = "#066594")) +
  scale_color_manual(values = c("Early-Aged" = "#f96161", 
                                "Middle-Aged" = "#66b2b2", 
                                "Late-Aged" = "#066594")) +
  
  scale_alpha_manual(values = c("Early-Aged" = 0.6, 
                                "Middle-Aged" = 0.1, 
                                "Late-Aged" = 0.6))+
  scale_x_continuous(limits = c(2.44, 2.77), breaks = c(2.45, 2.55, 2.65, 2.75))+
  labs(
    x = "b1 parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' age at reproduction",
    color = "Parents' age at reproduction",
    alpha="Parents' age at reproduction") +
  theme(legend.position = "top",
        axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.title=element_text(size=45),
        legend.text = element_text(size = 45),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), 
  )+
  guides(alpha=guide_legend(order = 1, title.position = "top",title.hjust = 0.5,
                            ,override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) 


ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/001_b1_parentalage_totallifespan.png",
       plot = tot_b1_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")


#9.4. C PARAMETER-------------------------
tot_c_plot <- ggplot(totalcdf, aes(x = Value, 
                                   fill = Parental_Age, 
                                   color = Parental_Age, 
                                   alpha=Parental_Age)) +
  geom_density(linewidth =2) +
  theme_classic() +
  scale_fill_manual(values = c("Early-Aged" = "#f96161", 
                               "Middle-Aged" = "#66b2b2", 
                               "Late-Aged" = "#066594")) +
  scale_color_manual(values = c("Early-Aged" = "#f96161", 
                                "Middle-Aged" = "#66b2b2", 
                                "Late-Aged" = "#066594")) +
  
  scale_alpha_manual(values = c("Early-Aged" = 0.6, 
                                "Middle-Aged" = 0.1, 
                                "Late-Aged" = 0.6))+
  scale_x_continuous(limits = c(0.0, 0.35), breaks = c(0.0, 0.1, 0.2, 0.3))+
  labs(
    x = "c parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' age at reproduction",
    color = "Parents' age at reproduction",
    alpha="Parents' age at reproduction") +
  
  theme(legend.position = "top",
        axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.title=element_text(size=45),
        legend.text = element_text(size = 45),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
  )+
  guides(alpha=guide_legend(order = 1, title.position = "top",title.hjust = 0.5,
                            ,override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) 

ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/001_c_parentalage_totallifespan.png",
       plot = tot_c_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")


#9.5--------CREATING THE SURVIVAL AND MORTALITY CURVES-----------------------------------

#Function for extracting mortality and survival estimates from the BaSTA model
plot_data_total <- list()
for (demv in c("mort", "surv", "dens")) {
  for (icat in seq_along(totalbastaweibull.weak[[demv]])) {
    cuts <-totalbastaweibull.weak$cuts[[icat]]
    minAge <- as.numeric(totalbastaweibull.weak$modelSpecs["min. age"])
    xx <- totalbastaweibull.weak$x[cuts] + minAge
    yy <- totalbastaweibull.weak[[demv]][[icat]][, cuts]
    
    df <- data.frame(
      Age = xx,
      Rate = yy[1, ], # Median
      LowerCI = yy[2, ], # Lower Credible interval Bound
      UpperCI = yy[3, ], # Upper Credible interval Bound
      Category = names(totalbastaweibull.weak[[demv]])[icat],
      Type = demv # Mortality or Survival (demographic parameters)
    )
    plot_data_total[[length(plot_data_total) + 1]] <- df
  }
}

# Combine all into one data frame
plot_data_total<- do.call(rbind, plot_data_total)

#renaming the values of Timepoint
plot_data_total$Category<- factor(plot_data_total$Category, 
                                  levels = c("Timepoint_binnedearly", "Timepoint_binnedmiddle", "Timepoint_binnedlate"), 
                                  labels = c("Early-Aged", "Middle-Aged", "Late-Aged"))


# Set the scaling factors (e.g., minimum and maximum age in days)
#Converts hazard so it is estimated in days rather than years
max_days <- 365.25 
min_days<-0

# Adjust Age column back to days
plot_data_total<- plot_data_total %>%
  mutate(Age_days = Age * (max_days - min_days) + min_days)

#ensuring the age categories are added as a factor
plot_data_total$Category<-as.factor(plot_data_total$Category)


#9.5.1. Splitting the data into the Mortality and Survival datasets------

#9.5.2. Hazard rate data frame--------------------------------
mort_df_total <- plot_data_total %>% filter(Type == "mort") %>% 
   mutate(Rate_weeks = Rate / 52.1429,
          weeks_LowerCI = LowerCI/52.1429,
          weeks_UpperCI = UpperCI/52.1429)
#Converting the hazard rate to be per unit week

#9.5.3. Survival data frame
surv_df_total<- plot_data_total%>% filter(Type == "surv")
dens_df_total<-plot_data_total %>% filter(Type=="dens")
  

#9.5.5. Median survival probability----------------------------------
surv_df_median<-surv_df_total %>%
  group_by(Category) %>%  
  arrange(Age_days) %>%
  mutate(CI = Rate) %>%  
  summarise(median_line = approx(CI, Age_days, xout = 0.5)$y)

#in weeks
surv_df_median_line<-surv_df_median %>%
  mutate(median_line_weeks = median_line/7)


#9.6. PLOTTING THE SURVIVAL AND MORTALITY CURVES-----------------------------------

#9.6.1. Cumulative survival probability plot---------------------------------------
total_survival_plot<- ggplot(data = surv_df_total,
                             aes(x = Age_days/7, #ensuring that age is in weeks
                                 y = Rate,
                                 colour = Category))+
  geom_line(data = surv_df_total, 
            aes(x = Age_days/7, 
                y = Rate,
                colour = Category),
            linewidth=4, 
            alpha=1)+
  geom_ribbon(data = surv_df_total, 
              aes(y=NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI, 
                  color= Category,
                  fill=Category, 
                  alpha = Category), 
              linetype="dashed",
              linewidth =1)+
  geom_segment(data = surv_df_median,
               aes(x = 0, xend = median_line/7, y = 0.5, yend = 0.5, color = Category),
               linetype = "dashed", linewidth = 2, alpha = 0.9) +
  geom_segment(data = surv_df_median,
               aes(x = median_line/7, xend = median_line/7, y = 0, yend = 0.5, color = Category),
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
        legend.background = element_rect(fill = "transparent", color = NA), ) +
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
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
  ) +labs(x = "Offspring age (weeks)", y = "Cumulative survival probability, S(x)")


# to save plot
  ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/001_survival_parentalage_totallifespan.png.png",
         plot = total_survival_plot, 
         bg="transparent",
         device = "png", 
         width = 520, 
         height = 320, 
         units = "mm")

  
#9.6.2. mortality risk plot---------------------------------------------------------------
total_mort_plot<- ggplot(data = mort_df_total,
                         aes(x = Age_days/7, 
                             y = Rate_weeks,
                             colour = Category))+
  geom_line(data = mort_df_total, 
            aes(x = Age_days/7, #ensuring that age is in weeks
                y = Rate_weeks,
                color = Category),
            linewidth =4, 
            alpha=1) +
  
  geom_ribbon(data = mort_df_total, 
              aes(y=NULL, 
                  ymin = weeks_LowerCI, 
                  ymax = weeks_UpperCI, 
                  color= Category,
                  fill=Category, 
                  alpha = Category),
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
        legend.background = element_rect(fill = "transparent", color = NA), ) +
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
  scale_y_continuous(breaks= c(0.00, 0.25, 0.50, 0.75, 1.00))+
  scale_color_manual(name = "Parents' adult age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' adult age at reproduction", values = c(0.6,0.1, 0.6))+
  scale_fill_manual(name = "Parents' adult age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  theme(
    legend.position = "top", 
    legend.box = "horizontal"
  ) +
  labs(x = "Offspring age (weeks)", y = "Instantaneous hazard rate, μ(x)")

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/001_mortality_parentalage_totallifespan.png",
       plot = total_mort_plot, 
       bg="transparent",
       device = "png", 
       width = 520, 
       height = 320, 
       units = "mm")


#------------------COMBINATION PLOT FOR PUBLICATION------------------------------
new_total_mort_plot<- total_mort_plot+
  labs( y = "Hazard rate, μ(x)")+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "none")
  
new_total_survival_plot<-total_survival_plot+
  labs( y = "Survival, S(x)")+
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top")
  

#timescale
new_tot_b1_plot<-tot_b1_plot+
  theme(axis.title.y = element_blank(),
    axis.title.x = element_blank(), 
        axis.title=element_text(size=50),
        axis.text=element_text(size=50),
    axis.text.y = element_blank(),
    axis.line.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))
  

#shape
new_tot_b0_plot<-tot_b0_plot+
  theme(axis.title.x = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50),
        axis.text.y = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))

#makeham
new_tot_c_plot<-tot_c_plot+
  labs( x = "Parameter value")+
  theme(axis.title.y = element_blank(),
                                 axis.title=element_text(size=50),
                                 axis.text=element_text(size=50),
        axis.line.y = element_blank(),
        axis.text.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))


#Combining all plots together as in "fancy basta" layout------------------------
lifespan_inference<-ggarrange(
  ggarrange(new_tot_b1_plot,new_tot_b0_plot, new_tot_c_plot, nrow = 3, labels = c("b1", "b0", "c"),
            font.label = list(size = 60, face = "bold"),
            label.x = c(0.03, 0.03, 0.05),
            heights = c(0.94, 0.94,  1),align = "v"),
  ggarrange(new_total_survival_plot, new_total_mort_plot, nrow = 2, labels = c("",""),
            font.label = list(size = 50, face = "bold"),
            heights = c(0.94,  1), align = "v"),
  ncol = 2,
  nrow=1,
  align = "h",
  labels = c("", ""),
  widths = c(0.9, 1.5))


ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/001_inference_plots.png",
       plot = lifespan_inference, 
       device = "png", 
       width = 900, 
       height = 630, 
       units = "mm")


# ===========================================================================
# 10. MODEL FIT PLOTS--------------------------------------------------------
# ===========================================================================

#Extracting the data from the Kaplan Meier curve
Km_survival_parentalage$data.survplot<-
  Km_survival_parentalage$data.survplot %>% 
  rename("Category" = "Timepoint_binnedplot")

# Create the zero-point for each group (so that the parametric and empirical curve start from the same point)
zero_point_001_age <- Km_survival_parentalage$data.survplot %>%
  group_by(Category) %>%
  summarise(time = 0, surv = 1, lower = 1, upper = 1)

# Combine zero-point and the actual KM data
Km_survival_parentalage$data.survplot<- Km_survival_parentalage$data.survplot %>%
  bind_rows(zero_point_001_age) %>%
  arrange(Category, time)  

#Combining the parametric and non-parametric estimates into one plot
Km_survival_combined_plot_001 <- ggplot() +
  #The Kaplan-Meier curve
  geom_step(data = Km_survival_parentalage$data.survplot, 
            aes(x = time, 
                y = surv), 
                colour = "#066594",
            linewidth = 2.5, alpha = 1)+
  geom_step(data = Km_survival_parentalage$data.survplot, 
            aes(x = time, 
                y = lower),
            colour = "#066594", 
            linetype="dashed",
            alpha = 0.5,
            size=1.2)+ 
  geom_step(data = Km_survival_parentalage$data.survplot, 
            aes(x = time, 
                y = upper),
            colour = "#066594", 
            linetype="dashed",
            alpha = 0.5,
            linewidth=1.2)+
  #Parametric survival curve
  geom_line(data = surv_df_total, 
            aes(x = Age_days / 7, #ensuring age is in weeks 
                y = Rate),
                colour = "#f96161", 
            linewidth = 2.5, alpha = 1) +
  geom_ribbon(data = surv_df_total, 
              aes(x = Age_days / 7, 
                  y = NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI),
                  fill = "#f96161", 
              alpha = 0.3)+
  facet_wrap(~ Category, ncol = 1) +  
  theme_classic()+
theme(axis.title=element_text(size=50),
      axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        strip.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), ) +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 5, 10, 15, 20, 25, 30))+
  labs(x = "Offspring age (weeks)", y = "Cumulative survival probability, S(x)")

#A good overall fit 

ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/001_modelfit_parentalage_totallifespan.png",
       plot = Km_survival_combined_plot_001, 
       bg="transparent",
       device = "png", 
       width = 620, 
       height = 480, 
       units = "mm")


# ===========================================================================
# 11. ADDING COVARIATES: ISOLATED EFFECT OF PARENTAL TEMPERATURE TREATMENT---
# ===========================================================================

#Plotting a Kaplein-Maier curve for parental temperature only
totcovdata$Temp<-factor(totcovdata$Temp, 
                                    levels = c("25.5", "28", "30.5"))
totcovdata$event<-ifelse(totcovdata$Depart.Type =="D", 1,0)

#11.1. Kaplain-Maier curve for parental temperature
km_fit_Temp <- survfit(Surv(Lifespan_weeks,event) ~ Temp, data = totcovdata)
Km_survival_Temp <- ggsurvplot(km_fit_Temp, 
                                      data = totcovdata, 
                                      pval = FALSE,  
                                      conf.int = TRUE,  
                                      risk.table = "abs_pct",
                                      break.time.by = 5,
                                      surv.median.line = "hv",  
                                      xlab = "Offspring age (weeks)", 
                                      ylab = "Cumulative survival probability, S(x)",
                                      font.x = 35,
                                      font.y = 35,
                                      font.tickslab = c(35, "grey25"),
                                      font.legend = 25,
                                      risk.table.fontsize = 6,
                                      legend.title = "Parents' Temperature Treatment",
                                      legend.labs = c("25.5°C", "28.0°C", "30.5°C"),
                                      palette = c("#2f4b7c","orange2", "#a64b61"))

#CombinING the survival plot and the risk table
Km_temp_combined_plot <- cowplot::plot_grid(Km_survival_Temp $plot, Km_survival_Temp$table, ncol = 1, rel_heights = c(3, 1))

# Save the combined plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/002_KM_temperature_totallifespan.png",
       plot = Km_temp_combined_plot,
       bg = "transparent",
       device = "png",
       width = 420,
       height = 360,
       units = "mm")

#11.2. PARAMETRIC SURVIVAL MODELS FOR PARENTAL TEMPERATURE EFFECTS------------------------------------------

#-----------------------Using weibull distributions only-------------------------

#11.2.1. Weibull model with Makeham term------------------
totalbastaweibulltemp.2 <- basta(object = totcovdata, 
                               dataType = "census",
                               formulaMort=~Temp-1, 
                               model="WE",
                               shape="Makeham",
                               niter=60000,
                               burnin=1001,
                               thinning=50,
                               nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibulltemp.2)
saveRDS(totalbastaweibulltemp.2, file = "scripts/model_outputs/BaSTA/total lifespan/temp_total_lifespan_weibull_Makeham.rda")
totalbastaweibulltemp.2<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/temp_total_lifespan_weibull_Makeham.rda")

#Assessing model fit and plotting the model outputs
plot(totalbastaweibulltemp, plot.type = "gof")
plot(totalbastaweibulltemp, plot.type = "demorates")
plot(totalbastaweibulltemp, densities=TRUE)


#11.2.2. weibull model----------------------------------------------
totalbastaweibulltemp <- basta(object = totcovdata, 
                               dataType = "census",
                               formulaMort=~Temp-1, 
                               model="WE",
                               niter=60000,#trying to double the burnin and the number of iterations
                               burnin=1001,
                               thinning=50,
                               nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibulltemp)
saveRDS(totalbastaweibulltemp, file = "scripts/model_outputs/BaSTA/total lifespan/temp_total_lifespan_weibull.rda")
totalbastaweibulltemp<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/temp_total_lifespan_weibull.rda")

#Assessing model fit and plotting the model outputs
plot(totalbastaweibulltemp, plot.type = "gof")
plot(totalbastaweibulltemp, plot.type = "demorates")
plot(totalbastaweibulltemp, densities=TRUE)


#11.2.3 weibull model with a bathtub shape--------------------------
totalbastaweibulltemp.3 <- basta(object = totcovdata, 
                                 dataType = "census",
                                 formulaMort=~Temp-1, 
                                 model="WE",
                                 shape="bathtub",
                                 niter=60000,
                                 burnin=1001,
                                 thinning=50,
                                 nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibulltemp.3)
saveRDS(totalbastaweibulltemp.3, file = "scripts/model_outputs/BaSTA/total lifespan/temp_total_lifespan_weibull.3.rda")
totalbastaweibulltemp.3<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/temp_total_lifespan_weibull.3.rda")

#Assessing model fit and plotting the model outputs
plot(totalbastaweibulltemp.3, plot.type = "gof")
plot(totalbastaweibulltemp.3, plot.type = "demorates")
plot(totalbastaweibulltemp.3, densities=TRUE)

#Selecting the Weibull Makeham model moving forwards

# ===========================================================================
# 12. PRIOR SENSITIVITY ANALYSIS FOR THE PARENTAL TEMPERATURE MODEL-----------
# ===========================================================================
#Using the same priors as specified for the parental age analysis

#12.1. Weibull Makeham model with weak priors---------------------------------
totalbastaweibulltemp.weak <- basta(object = totcovdata, 
                                 dataType = "census",
                                 formulaMort=~Temp-1, 
                                 model="WE",
                                 shape="Makeham",
                                 thetaPriorMean = weakMean,
                                 thetaPriorSd = weakSd,
                                 thetaPriorLower = weakLower,
                                 niter=80000,
                                 burnin=1001,
                                 thinning=50,
                                 nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibulltemp.weak)
saveRDS(totalbastaweibulltemp.weak, file = "scripts/model_outputs/BaSTA/total lifespan/temp_total_lifespan_weibull_Makeham_weakprior.rda")
totalbastaweibulltemp.weak<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/temp_total_lifespan_weibull_Makeham_weakprior.rda")#this is with the makeham shape to be confusing 

#Assessing model fit and plotting the model outputs
plot(totalbastaweibulltemp.weak, plot.type = "gof")
plot(totalbastaweibulltemp.weak, plot.type = "demorates")
plot(totalbastaweibulltemp.weak, densities=TRUE)

#12.2. With moderate priors---------------------------------------
totalbastaweibulltemp.moderate <- basta(object = totcovdata, 
                                    dataType = "census",
                                    formulaMort=~Temp-1, 
                                    model="WE",
                                    shape="Makeham",
                                    thetaPriorMean = moderateMean,
                                    thetaPriorSd = moderateSd,
                                    thetaPriorLower = moderateLower,
                                    niter=80000,
                                    burnin=1001,
                                    thinning=50,
                                    nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibulltemp.moderate)
saveRDS(totalbastaweibulltemp.moderate, file = "scripts/model_outputs/BaSTA/total lifespan/temp_total_lifespan_weibull_Makeham_moderate.rda")
totalbastaweibulltemp.moderate<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/temp_total_lifespan_weibull_Makeham_moderate.rda") 

# ===========================================================================
# 13. ISOLATED EFFECT OF PARENTAL TEMPERATURE PLOTS-----------------------------------
# ===========================================================================
#Plotting estimates from weibull model with weakly-informative priors

#13.1. Plotting the posterior distributions for the mortality parameters (c, b0, and b1)

#extracting parameters
temptotal_theta_params <- totalbastaweibulltemp.weak$params

#converting to a data frame
temptotal_theta_params<-as.data.frame(temptotal_theta_params)

#Reshaping the data
temptotal_posterior_df <- temptotal_theta_params %>%
  pivot_longer(
    cols = everything(),
    names_to = "Combined",
    values_to = "Value"
  ) %>%
  mutate(
    Type = str_extract(Combined, "^[^.]+"),  #Extracts everything before the first period
    Temp = str_extract(Combined, "(?<=Temp)\\d+\\.?\\d*"),  # Extract numbers after "Temp"
    Temp = case_when(
      Temp == "25.5" ~ "25.5°C",
      Temp == "28"   ~ "28.0°C",
      Temp == "30.5" ~ "30.5°C",
      TRUE ~ Temp
    )
  ) %>%
  select(Type, Temp, Value) 

#making sure the levels are in  the right order for the plot:
temptotal_posterior_df$Temp<- factor(temptotal_posterior_df$Temp, 
                                         levels = c("25.5°C", "28.0°C", "30.5°C"))


# 13.2. Creating separate datasets for `b0` and `b1`, and c mortality parameters:
temtotb0_df <- temptotal_posterior_df %>% filter(Type == "b0")
temtotb1_df <- temptotal_posterior_df %>% filter(Type == "b1")
temtotalcdf<-temptotal_posterior_df %>% filter(Type == "c")

#13.3 Plotting the b0 parameter----------------
temtot_b0_plot <- ggplot(temtotb0_df, aes(x = Value, 
                                          fill = Temp, 
                                          color = Temp, 
                                          alpha=Temp)) +
  geom_density(linewidth= 2) +
  theme_classic() +
  scale_fill_manual(values = c("25.5°C" = "#2f4b7c", 
                               "28.0°C"= "orange2", 
                               "30.5°C" = "#a64b61")) +
  scale_color_manual(values = c("25.5°C"= "#2f4b7c", 
                                "28.0°C" = "orange2", 
                                "30.5°C" = "#a64b61"))+
  scale_alpha_manual(values = c("25.5°C" = 0.7, 
                              "28.0°C" = 0.1, 
                              "30.5°C" = 0.4))+
  scale_x_continuous(limits = c(3.9, 7.1), breaks = c(4, 5, 6, 7))+
  labs(
    x = "b0 parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' temperature treatment",
    color = "Parents' temperature treatment",
    alpha="Parents' temperature treatment") +
  theme(legend.position = "top",
        axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.title=element_text(size=45),
        legend.text = element_text(size = 45),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), 
  )+
  guides(alpha=guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) 

#saving plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/002_b0_parentaltemp_totallifespan.png",
       plot = temtot_b0_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")

#13.2. Creating the posterior plot for the b1 posterior parameter values------------------
temtot_b1_plot <- ggplot(temtotb1_df, aes(x = Value, 
                                          fill = Temp,
                                          color = Temp,
                                          alpha=Temp)) +
  geom_density(linewidth = 2) +
  theme_classic() +
  scale_fill_manual(values = c("25.5°C" = "#2f4b7c", 
                               "28.0°C"= "orange2", 
                               "30.5°C" = "#a64b61")) +
  scale_color_manual(values = c("25.5°C"= "#2f4b7c", 
                                "28.0°C" = "orange2", 
                                "30.5°C" = "#a64b61"))+
scale_alpha_manual(values = c("25.5°C" = 0.7, 
                              "28.0°C" = 0.1, 
                              "30.5°C" = 0.4)) +
  scale_x_continuous(limits = c(2.49, 2.75))+
  labs(
    x = "b1 parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' temperature treatment",
    color = "Parents' temperature treatment",
    alpha="Parents' temperature treatment") +
  theme(legend.position = "top",
        axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.title=element_text(size=45),
        legend.text = element_text(size = 45),
        panel.background = element_rect(fill = "transparent", color = NA),  # Make panel background transparent
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), # Make plot background transparent
  )+
  guides(alpha=guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) 

#saving plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/002_b1_parentaltemp_totallifespan.png",
       plot = temtot_b1_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")

#13.3. Creating the posterior plot for the c parameter posterior values---------------
temtot_c_plot <- ggplot(temtotalcdf, aes(x = Value, fill = Temp, color = Temp, alpha=Temp)) +
  geom_density(linewidth = 2) +
  theme_classic() +
  scale_fill_manual(values = c("25.5°C" = "#2f4b7c", 
                               "28.0°C"= "orange2", 
                               "30.5°C" = "#a64b61")) +
  scale_color_manual(values = c("25.5°C"= "#2f4b7c", 
                                "28.0°C" = "orange2", 
                                "30.5°C" = "#a64b61"))+
  scale_alpha_manual(values = c("25.5°C" = 0.7, 
                                "28.0°C" = 0.1, 
                                "30.5°C" = 0.4))+
  scale_x_continuous(limits = c(0.0, 0.3), breaks = c(0.0, 0.1, 0.2, 0.3))+
  labs(
    x = "c parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' temperature treatment",
    color = "Parents' temperature treatment",
    alpha="Parents' temperature treatment") +
  theme(legend.position = "top",
        axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.title=element_text(size=45),
        legend.text = element_text(size = 45),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), 
  )+ 
    guides(alpha=guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) 

#saving plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/002_c_parentaltemp_totallifespan.png",
       plot = temtot_c_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")


#13.4. Creating the survival and Instantaneous hazard rate curves------

#13.4.1. function to extract necessary values from the BaSTA model-----
tempplot_data_total <- list()
for (demv in c("mort", "surv")) {
  for (icat in seq_along(totalbastaweibulltemp.weak[[demv]])) {
    cuts <-totalbastaweibulltemp.weak$cuts[[icat]]
    minAge <- as.numeric(totalbastaweibulltemp.weak$modelSpecs["min. age"])
    xx <- totalbastaweibulltemp.weak$x[cuts] + minAge
    yy <- totalbastaweibulltemp.weak[[demv]][[icat]][, cuts]
    #Converting to long format
    df <- data.frame(
      Age = xx,
      Rate = yy[1, ], # Median
      LowerCI = yy[2, ], # Lower Credible interval
      UpperCI = yy[3, ], # Upper Credible interval
      Category = names(totalbastaweibulltemp.weak[[demv]])[icat],
      Type = demv # Mortality or Survival (i.e., the demographic parameters)
    )
    tempplot_data_total[[length(tempplot_data_total) + 1]] <- df
  }
}

#Combining all into one data frame
tempplot_data_total<- do.call(rbind, tempplot_data_total)

#renaming the values of Temperature
tempplot_data_total$Category<- factor(tempplot_data_total$Category, 
                                  levels = c("Temp25.5", "Temp28", "Temp30.5"), 
                                  labels = c("25.5°C", "28.0°C", "30.5°C"))

#Setting the scaling factors (e.g., minimum and maximum age in days)
#ensures the axis is estimated in days not years
max_days <- 365.25 
min_days<-0

# Adjust the Age column back to days (from years)
tempplot_data_total<- tempplot_data_total %>%
  mutate(Age_days = Age * (max_days - min_days) + min_days)

tempplot_data_total$Category<-as.factor(tempplot_data_total$Category)


#13.4.2. Splitting the data into Mortality and Survival datasets (for plotting)

#mortality data frame------------
mort_df_temp <- tempplot_data_total %>% filter(Type == "mort") %>% 
  mutate(Rate_weeks = Rate / 52.1429,
         weeks_LowerCI = LowerCI/52.1429, #converting the estimate to be per week
         weeks_UpperCI = UpperCI/52.1429)

#survival data frame-------------
surv_df_temp<- tempplot_data_total%>% filter(Type == "surv")

#median survival probability-----
tempsurv_df_median<-surv_df_temp %>%
  group_by(Category) %>%  
  arrange(Age_days) %>%
  mutate(CI = 1 - Rate) %>% 
  summarise(median_line = approx(CI, Age_days, xout = 0.5)$y)

#13.4.3. cumulative survival plot - effect of parental temperature-----------
temptotal_survival_plot<- ggplot(data = surv_df_temp,
                             aes(x = Age_days/7, 
                                 y = Rate,
                                 colour = Category))+
  geom_line(data = surv_df_temp, 
            aes(x = Age_days/7, 
                y = Rate,
                colour = Category),
            linewidth = 4, 
            alpha=1)+
  geom_ribbon(data = surv_df_temp, 
              aes(y=NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI, 
                  color= Category,
                  fill=Category,
                  alpha=Category),
                 linewidth = 1,
                 linetype="dashed") +
  geom_segment(data = tempsurv_df_median,
               aes(x = 0, xend = median_line/7, y = 0.5, yend = 0.5, color = Category),
               linetype = "dashed", linewidth = 2, alpha = 0.9) +
  geom_segment(data = tempsurv_df_median,
               aes(x = median_line/7, xend = median_line/7, y = 0, yend = 0.5, color = Category),
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
        legend.background = element_rect(fill = "transparent", color = NA), ) +
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
  scale_color_manual(name = "Parents' temperature treatment", values = c("#26408B","orange2", "#a64b61"))+
  scale_fill_manual(name = "Parents' temperature treatment", values = c("#26408B","orange2", "#a64b61"))+
  scale_alpha_manual(name = "Parents' temperature treatment", values = c(0.7,0.1,0.4))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
         )+
  labs(x = "Offspring age (weeks)", y = "Cumulative survival probability, S(x)")

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/002_survival_parentaltemp_totallifespan.png",
       plot = temptotal_survival_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 220, 
       units = "mm")


#13.4.4. Instantaneous mortality risk plot--------------------
temptotal_mort_plot<- ggplot(data = mort_df_temp,
                         aes(x = Age_days/7,
                             y = Rate_weeks,
                             colour = Category))+
  geom_line(data = mort_df_temp, 
            aes(x = Age_days/7, 
                y = Rate_weeks,
                color = Category),
            linewidth = 4, alpha=1) +
  geom_ribbon(data = mort_df_temp, 
              aes(y=NULL, 
                  ymin = weeks_LowerCI, 
                  ymax = weeks_UpperCI, 
                  color= Category,
                  fill=Category,
                  alpha=Category),
              linewidth = 1,
              linetype="dashed")+
  theme_classic()+ 
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), ) +
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
  scale_y_continuous(limits = c(0, 1.22), breaks = c(0, 0.3, 0.6, 0.9, 1.2))+
  scale_color_manual(name = "Parents' temperature treatment", values = c("#2f4b7c","orange2", "#a64b61"))+
  scale_fill_manual(name = "Parents' temperature treatment", values = c("#2f4b7c","orange2", "#a64b61"))+
  scale_alpha_manual(name = "Parents' temperature treatment", values = c(0.7,0.1,0.4))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  labs(x = "Offspring age (weeks)", y = "Instantaneous hazard rate, μ(x)")

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/002_mortality_parentaltemp_totallifespan.png",
       plot = temptotal_mort_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 220, 
       units = "mm")


#-----------------------COMBINED PLOT FOR PAPER---------------------------------
new_temp_mort_plot<- temptotal_mort_plot+
  labs( y = "Hazard rate, μ(x)")+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "none")

new_temp_survival_plot<-temptotal_survival_plot+
  labs( y = "Survival, S(x)")+
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50))


#timescale
new_temp_b1_plot<-temtot_b1_plot+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.title=element_text(size=50),
        axis.text=element_text(size=50),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))+
  scale_y_continuous(limits=c(0,17), breaks = c(0,5,10,15))+
  scale_x_continuous(limits=c(2.49,2.75), breaks = c(2.5, 2.58, 2.66, 2.74))


#shape
new_temp_b0_plot<-temtot_b0_plot+
  theme(axis.title.x = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))+
  scale_y_continuous(limits=c(0,2), breaks = c(0, 0.6, 1.2, 1.8))

#makeham
new_temp_c_plot<-temtot_c_plot+
  labs( x = "Parameter value")+
  theme(axis.title.y = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))+
  scale_y_continuous(limits=c(0,29), breaks= c(0, 9, 18, 27))


#Combining all plots together as in "fancy basta" layout------------------------

lifespan_temp_inference<-ggarrange(
  ggarrange(new_temp_b1_plot,new_temp_b0_plot, new_temp_c_plot, nrow = 3, labels = c("b1", "b0", "c"),
            font.label = list(size = 60, face = "bold"),
            label.x = c(0.03, 0.03, 0.05),
            heights = c(0.94, 0.94,  1),align = "v"),
  ggarrange(new_temp_survival_plot, new_temp_mort_plot, nrow = 2, labels = c("",""),
            font.label = list(size = 50, face = "bold"),
            heights = c(0.94,  1),align = "v"),
  ncol = 2,
  nrow=1,
  align = "h",
  labels = c("", ""),
  widths = c(0.9, 1.5))


ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/002_inference_temperature.png",
       plot = lifespan_temp_inference, 
       device = "png", 
       width = 900, 
       height = 630, 
       units = "mm")

# ===========================================================================
# 14. MODEL FIT PLOTS FOR PARENTAL TEMPERATURE-------------------------------
# ===========================================================================

#Extracting the data from the Kaplan Meier curve
Km_survival_Temp$data.survplot<-
  Km_survival_Temp$data.survplot %>%
  rename(Category = Temp) %>% 
  mutate(Category = case_when(  
    Category == "25.5" ~ "25.5°C",
    Category == "28"   ~ "28.0°C",
    Category == "30.5" ~ "30.5°C",
    TRUE ~ Category 
  ))

# Create the zero-point for each group
zero_point_002_temp <- Km_survival_Temp$data.survplot %>%
  group_by(Category) %>%
  summarise(time = 0, surv = 1, lower = 1, upper = 1)

# Combine zero-point and the actual KM data
Km_survival_Temp$data.survplot<-Km_survival_Temp$data.survplot %>%
  bind_rows(zero_point_002_temp) %>%
  arrange(Category, time)  


#14.1. Plot Combining the parametric and non-parametric estimates into one---------
Km_survival_combined_plot_002 <- ggplot() +#
  #The Kaplan-Meier curve
  geom_step(data = Km_survival_Temp$data.survplot, 
            aes(x = time, 
                y = surv), 
            colour = "#066594",
            linewidth = 2.5, alpha = 1)+
  geom_step(data = Km_survival_Temp$data.survplot, 
            aes(x = time, 
                y = lower),
            colour = "#066594", 
            linetype="dashed",
            alpha = 0.5,
            size=1.2)+ 
  geom_step(data = Km_survival_Temp$data.survplot, 
            aes(x = time, 
                y = upper),
            colour = "#066594", 
            linetype="dashed",
            alpha = 0.5,
            linewidth=1.2)+
  #Parametric survival curve
  geom_line(data = surv_df_temp, 
            aes(x = Age_days / 7, 
                y = Rate),
            colour = "#f96161", 
            linewidth = 2.5, alpha = 1) +
  geom_ribbon(data = surv_df_temp, 
              aes(x = Age_days / 7, #ensuring that age is in weeks not days
                  y = NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI),
              fill = "#f96161", 
              alpha = 0.3)+
  facet_wrap(~ Category, ncol = 1) + 
  theme_classic()+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        strip.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), ) +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 5, 10, 15, 20, 25, 30))+
  labs(x = "Offspring age (weeks)", y = "Cumulative survival probability, S(x)")

ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/002_modelfit_Temp_totallifespan.png",
       plot = Km_survival_combined_plot_002 , 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 380, 
       units = "mm")

#14.2. Creating a combined model fit plot of parental temperature and parental age---------------------------
new_Km_survival_combined_plot_002<-Km_survival_combined_plot_002+
  theme(axis.title.y = element_blank())

#combining into one plot using ggarange
tempagefit<-ggarrange(
  Km_survival_combined_plot_001,
  new_Km_survival_combined_plot_002,
  ncol = 2,
  nrow=1,
  align = "h",
  labels = c("", ""),
  widths = c(0.6,0.5))


ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/001_modelfits.png",
       plot = tempagefit,
       device = "png", 
       width = 900, 
       height = 630, 
       units = "mm")


# ===========================================================================
# 15. ADDING COVARIATTES: INTERACTIVE EFFECT OF PARENTAL AGE & TEMPERATURE----
# ===========================================================================

#15.1 Creating the Kaplan Meier curve for parental age*temperature
KMdataagetemp<- totcovdata %>%
  mutate(temp_ParentalAge = interaction(Timepoint_binnedplot, Temp)) #using interaction() to create a new interactive term for age and temperature
#adding celsius to the label
KMdataagetemp$Temp_label <- paste0(KMdataagetemp$Temp, "°C")

#15.2.Kaplain-Maier curve for parental age and temperature interaction
#survival model
km_fit_parentalagetemp <- survfit(Surv(Lifespan_weeks,event) ~ temp_ParentalAge, data = KMdataagetemp)
#KM curve
Km_survival_parentalagetempplot <- ggsurvplot(
  km_fit_parentalagetemp,
  data = KMdataagetemp,
  pval = FALSE,
  conf.int = TRUE,
  xlab = "Offspring age (weeks)",
  ylab = "Cumulative survival probability, S(x)",
  ggtheme = theme_classic() +
    theme(
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      strip.text = element_text(size = 15),
      legend.position = "none",
    ),
  palette = c("#f96161", "#66b2b2", "#066594", "#f96162", "#66b2b3", "#066593", "#f96164", "#66b2b5", "#066597"),
 facet.by = "Temp_label",
 ncol=1,
 short.panel.labs =TRUE,
 legend.labs=NULL)

#adding my own manual theme to this plot
Km_survival_parentalagetempplot<- Km_survival_parentalagetempplot + theme(legend.position = "none")


# Save the combined plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/004_KM_parentalagetemp_totallifespan.png",
       plot = Km_survival_parentalagetempplot,
       bg = "transparent",
       device = "png",
       width = 380,
       height = 340,
       units = "mm")


#15.2. PARAMETRIC SURVIVAL MODELS--------------------------------------

#15.2.1. Weibull model with a Makeham term--------------
tempagebastaweibull.2<- basta(object = totcovdata, 
                             dataType = "census",
                             formulaMort=~Temp:Timepoint_binned-1,
                             model="WE",
                             shape="Makeham",
                             niter=60000,
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(tempagebastaweibull.2)
saveRDS(tempagebastaweibull.2, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_age_temp.rda")
tempagebastaweibull.2<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_age_temp.rda")

#Assessing model fit and plotting the model outputs
plot(tempagebastaweibull.2, plot.type = "gof")
plot(tempagebastaweibull, plot.type = "demorates")
plot(tempagebastaweibull, densities=TRUE)

#15.2.2.Weibull model---------------------------
tempagebastaweibull<- basta(object = totcovdata, 
                           dataType = "census",
                           formulaMort= ~Temp:Timepoint_binned-1, 
                           model="WE",
                           niter=60000,#trying to double the burnin and the number of iterations
                           burnin=1001,
                           thinning=50,
                           nsim = 4, parallel = TRUE, ncpus = 4)
summary(tempagebastaweibull)
saveRDS(tempagebastaweibull, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_tempage_weibull.rda")
tempagebastaweibull<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_tempage_weibull.rda")


#15.2.3.Weibull model with a bathtub term-------------------
tempagebastaweibull.3<- basta(object = totcovdata, 
                             dataType = "census",
                             formulaMort= ~Temp:Timepoint_binned-1, 
                             model="WE",
                             shape="bathtub",
                             niter=60000,#trying to double the burnin and the number of iterations
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(tempagebastaweibull.3)
saveRDS(tempagebastaweibull.3, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_tempage_weibull_bathtub.rda")
tempagebastaweibull.3<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_tempage_weibull_bathtub.rda")

#Selecting the Weibull Makeham model for the prior specification process


# ===========================================================================
# 16. PRIOR SENSITIVITY ANALYSIS FOR THE PARENTAL AGE X TEMPERATURE MODEL-----------
# ===========================================================================
#Using the same priors as specified for the parental age analysis

#16.1. Prior specification---------------------------

#Weakly informative priors (mean)
weakMean3 <- matrix(c(0, 1.5, prior_scale_mu), 
                    nrow = 9, ncol = 3, byrow = TRUE)
#Priors for the mean of the c, b0 and b1 terms

#weakly infromative priors (SD)
weakSd3 <- matrix(c(
  1.0, 1.0, 0.5
), nrow = 9, ncol = 3, byrow = TRUE)
#Priors for the sd of the c, b0 and b1 terms

weakLower3<-matrix(c(
  0, 0, 0
), nrow = 9, ncol = 3,byrow = TRUE)
#Priors for the lower bound of the c, b0 and b1 terms

#16.2. Fitting the priors to the model (just testing the weak priors here)------------
#Fitting weak priors to the weibull makeham model
tempagebastaweibull.weak<- basta(object = totcovdata, 
                              dataType = "census",
                              formulaMort=~Temp:Timepoint_binned-1, 
                              thetaPriorMean = weakMean3,
                              thetaPriorSd = weakSd3,
                              thetaPriorLower = weakLower3,
                              model="WE",
                              shape="Makeham",
                              niter=80000,
                              burnin=1001,
                              thinning=50,
                              nsim = 4, parallel = TRUE, ncpus = 4)
summary(tempagebastaweibull.weak)
saveRDS(tempagebastaweibull.weak, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_age_tempweak.rda")
tempagebastaweibull.weak<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_age_tempweak.rda")

#Assessing model fit and plotting the model outputs
plot(tempagebastaweibull.weak, plot.type = "gof")
plot(tempagebastaweibull.weak, plot.type = "demorates")
plot(tempagebastaweibull.weak, densities=TRUE)


# =======================================================================================================
# 17. PLOTTINGT THE MORTALITY PARAMETERS AND DEMOGRAPHIC CURVES FOR THE PARENTAL AGE X TEMP INTERACTION-----------
# ========================================================================================================

#17.1.PLOTTING THE POSTERIOR DISTRIBUTIONS FOR B0, B1 AND C------------------------------

#Extracting the model parameters
tempage_theta_params <- tempagebastaweibull.weak$params

#saving as a data frame
tempage_theta_params<-as.data.frame(tempage_theta_params)

#Reshaping the data
tempage_posterior_df <- tempage_theta_params %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Type", "Temp", "Timepoint_binned"),
    names_pattern = "(b[0-9]|c)\\.Temp(\\d{2}(?:\\.\\d)?):Timepoint_binned(early|middle|late)",
    values_to = "Value"
  ) %>%
  mutate(
    `Parental_Age + Temperature` = paste(Temp, Timepoint_binned, sep = "."),
    `Parental_Age + Temperature` = dplyr::recode(`Parental_Age + Temperature`, 
                                          "25.5.early" = "Early-Aged (25.5°C)",
                                          "25.5.middle" = "Middle-Aged (25.5°C)",
                                          "25.5.late" = "Late-Aged (25.5°C)",
                                          "28.early" = "Early-Aged (28°C)",
                                          "28.middle" = "Middle-Aged (28°C)",
                                          "28.late" = "Late-Aged (28°C)",
                                          "30.5.early" = "Early-Aged (30.5°C)",
                                          "30.5.middle" = "Middle-Aged (30.5°C)",
                                          "30.5.late" = "Late-Aged (30.5°C)")
  ) %>%
  separate(`Parental_Age + Temperature`, 
           into = c("Parental_Age", "Temperature"), 
           sep = " \\(", 
           extra = "merge") %>%
  mutate(
    Temperature = gsub("\\)$", "", Temperature) 
  ) %>%
  select(Type, Parental_Age, Temperature, Value)



#making sure the levels are in  the right order for the plot:
tempage_posterior_df$Parental_Age<- factor(tempage_posterior_df$Parental_Age, 
                                         levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))



#17.2. Create two separate datasets for `b0` and `b1`, and c------------
tottempb0_df <- tempage_posterior_df  %>% filter(Type == "b0")
tottempb1_df <- tempage_posterior_df  %>% filter(Type == "b1")
tottempcdf<-tempage_posterior_df  %>% filter(Type == "c")

#17.3. B0 PARAMETER PLOT--------------
tottemp_b0_plot <- ggplot(tottempb0_df, aes(x = Value, fill = Parental_Age, color = Parental_Age, alpha=Parental_Age)) +
  geom_density(size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("Early-Aged" = "#f96161", 
                               "Middle-Aged" = "#66b2b2", 
                               "Late-Aged" = "#066594")) +
  scale_color_manual(values = c("Early-Aged" = "#f96161", 
                                "Middle-Aged" = "#66b2b2", 
                                "Late-Aged" = "#066594")) +
  scale_alpha_manual(values = c("Early-Aged" = 0.6, 
                                "Middle-Aged" = 0.1, 
                                "Late-Aged" = 0.6))+
  scale_x_continuous(limits=c(3.5,7.5))+
  facet_wrap(~Temperature, ncol = 1)+
  labs(
    x = "b0 parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' age at reproduction",
    color = "Parents' age at reproduction",
    alpha="Parents' age at reproduction") +
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_text(size = 35, face = "bold"))+
guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                            override.aes = list(size = 10)),
       fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
       colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
)

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/004_b0_tempandage_totallifespan.png",
       plot = tottemp_b0_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 460, 
       units = "mm")


#17.4.-B1 PARAMETER PLOT-------------
tottemp_b1_plot <- ggplot(tottempb1_df, aes(x = Value, 
                                            fill = Parental_Age, 
                                            color = Parental_Age, alpha=Parental_Age)) +
  geom_density(size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("Early-Aged" = "#f96161", 
                               "Middle-Aged" = "#66b2b2", 
                               "Late-Aged" = "#066594")) +
  scale_color_manual(values = c("Early-Aged" = "#f96161", 
                                "Middle-Aged" = "#66b2b2", 
                                "Late-Aged" = "#066594")) +
  
  scale_alpha_manual(values = c("Early-Aged" = 0.6, 
                                "Middle-Aged" = 0.1, 
                                "Late-Aged" = 0.6))+
  facet_wrap(~Temperature, ncol=1)+
  scale_x_continuous(limits = c(2.39, 2.81), breaks = c(2.4, 2.5, 2.6, 2.7, 2.8))+
  labs(
    x = "b1 parameter value",
    y = "Probability Density, f(x)",
    fill = "Parents' age at reproduction",
    color = "Parents' age at reproduction",
    alpha="Parents' age at reproduction") +
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_text(size = 35, face = "bold"))+
guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                            override.aes = list(size = 10)),
       fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
       colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
)

ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/004_b1_tempandage_totallifespan.png",
       plot = tottemp_b1_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 460, 
       units = "mm")



#17.5. C PARAMETER PLOT----------------------
tottemp_c_plot <- ggplot(tottempcdf, aes(x = Value, fill = Parental_Age, color = Parental_Age, alpha=Parental_Age)) +
  geom_density(size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("Early-Aged" = "#f96161", 
                               "Middle-Aged" = "#66b2b2", 
                               "Late-Aged" = "#066594")) +
  scale_color_manual(values = c("Early-Aged" = "#f96161", 
                                "Middle-Aged" = "#66b2b2", 
                                "Late-Aged" = "#066594")) +
  
  scale_alpha_manual(values = c("Early-Aged" = 0.6, 
                                "Middle-Aged" = 0.1, 
                                "Late-Aged" = 0.6))+
  facet_wrap(~Temperature, ncol=1)+
  scale_x_continuous(limits = c(0, 0.3))+
  labs(
    x = "b1 parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' age at reproduction",
    color = "Parents' age at reproduction",
    alpha="Parents' age at reproduction") +
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_text(size = 35, face = "bold"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
  )

ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/004_c_tempandage_totallifespan.png",
       plot = tottemp_c_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 460, 
       units = "mm")


#17.5. PLOTTING THE SURVIVAL AND MORTALITY CURVES---------------------------

#17.5.1 Function for extracting demorgaphic parameters from the model----
plot_data_agetemp <- list()
for (demv in c("mort", "surv")) {
  for (icat in seq_along(tempagebastaweibull.weak[[demv]])) {
    cuts <-tempagebastaweibull.weak$cuts[[icat]]
    minAge <- as.numeric(tempagebastaweibull.weak$modelSpecs["min. age"])
    xx <- tempagebastaweibull.weak$x[cuts] + minAge
    yy <- tempagebastaweibull.weak[[demv]][[icat]][, cuts]
    
    # Convert to long format
    df <- data.frame(
      Age = xx,
      Rate = yy[1, ], # Median
      LowerCI = yy[2, ], # Lower Credible interval
      UpperCI = yy[3, ], # Upper Credible interval
      Category = names(tempagebastaweibull.weak[[demv]])[icat],
      Type = demv # Mortality or Survival
    )
    plot_data_agetemp[[length(plot_data_agetemp) + 1]] <- df
  }
}

#Combine all into one data frame
plot_data_agetemp<- do.call(rbind, plot_data_agetemp)

#splitting parental age and temperature into two seperate columns
plot_data_agetemp <- plot_data_agetemp %>%
  mutate(
    Temp = str_extract(Category, "\\d{2}\\.\\d|\\d{2}"),  # Extract temperature part
    Parental_Age = str_extract(Category, "early|middle|late")  # Extract parental age
  ) %>%
  mutate(
    Temp = factor(Temp, levels = c("25.5", "28", "30.5"), 
                  labels = c("25.5°C", "28°C", "30.5°C")),
    Parental_Age = factor(Parental_Age, 
                          levels = c("early", "middle", "late"),
                          labels = c("Early-Aged", "Middle-Aged", "Late-Aged"))
  )

#making sure the levels are inthe right order for the plot:
plot_data_agetemp$Parental_Age<- factor(plot_data_agetemp$Parental_Age, 
                                           levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))


# Set the scaling factors (e.g., minimum and maximum age in days)
max_days <- 365.25 
min_days<-0

# Adjust Age column back to days
plot_data_agetemp<- plot_data_agetemp %>%
  mutate(Age_days = Age * (max_days - min_days) + min_days)
plot_data_agetemp$Parental_Age<-as.factor(plot_data_agetemp$Parental_Age)


#17.2. Splitting the data into Mortality and Survival datasets----------

#mortality data set
mort_df_agetemp <- plot_data_agetemp %>% filter(Type == "mort")

#survival data set
surv_df_agetemp<- plot_data_agetemp%>% filter(Type == "surv")


#17.3. Median survival probability-------
surv_df_median<-surv_df_agetemp %>%
  group_by(Temp, Parental_Age) %>%  
  arrange(Age_days) %>%
  mutate(CI = 1 - Rate) %>% 
  summarise(median_line = approx(CI, Age_days, xout = 0.5)$y)

#17.3.1. cumulative survival plot-----------
agetemp_survival_plot<- ggplot(data = surv_df_agetemp,
                             aes(x = Age_days/7, 
                                 y = Rate,
                                 colour = Parental_Age))+
  geom_line(data = surv_df_agetemp, 
            aes(x = Age_days/7, 
                y = Rate,
                colour = Parental_Age),
            linewidth=4, 
            alpha=1)+
  geom_ribbon(data = surv_df_agetemp, 
              aes(y=NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI, 
                  color= Parental_Age,
                  fill=Parental_Age, 
                  alpha = Parental_Age), linetype="dashed") +
  geom_segment(data = surv_df_median,
               aes(x = 0, xend = median_line/7, y = 0.5, yend = 0.5, color = Parental_Age),
               linetype = "dashed", linewidth = 2, alpha = 0.9) +
  geom_segment(data = surv_df_median,
               aes(x = median_line/7, xend = median_line/7, y = 0, yend = 0.5, color = Parental_Age),
               linetype = "dashed", linewidth =2, alpha = 0.9)+
  facet_wrap(~Temp, ncol =1)+
  theme_classic()+ 
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_text(size = 35, face = "bold"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
  )+
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
  scale_color_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' age at reproduction", values = c(0.6,0.1, 0.6))+
  scale_fill_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  labs(x = "Offspring age (weeks)", y = "Cumulative survival probability, S(x)")


#saving plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/004_survival_tempandage_totallifespan.png",
       plot = agetemp_survival_plot, 
       bg="transparent",
       device = "png", 
       width = 370, 
       height = 440, 
       units = "mm")


#17.3.2. Mortality risk plot---------------------
agetemp_mort_plot<- ggplot(data = mort_df_agetemp,
                         aes(x = Age_days/7,
                             y = Rate/52.1429,
                             colour = Parental_Age))+
  geom_line(data = mort_df_agetemp, 
            aes(x = Age_days/7, 
                y = Rate/52.1429,
                color = Parental_Age),
            linewidth=4, alpha=1)+
  geom_ribbon(data = mort_df_agetemp, 
              aes(y=NULL, 
                  ymin = LowerCI/52.1429, 
                  ymax = UpperCI/52.1429, 
                  color= Parental_Age,
                  fill= Parental_Age, 
                  alpha = Parental_Age),
              linewidth = 1,
              linetype="dashed") +
  facet_wrap(~Temp, ncol=1)+
  theme_classic()+ 
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_text(size = 35, face = "bold"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
  )+
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
  scale_color_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' age at reproduction", values = c(0.6,0.1, 0.6))+
  scale_fill_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  labs(x = "Offspring age (weeks)", y = "Instantaneous hazard rate, μ(x)")

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/004_mortality_tempandage_totallifespan.png",
       plot = agetemp_mort_plot, 
       bg="transparent",
       device = "png", 
       width = 370, 
       height = 440, 
       units = "mm")

#--------------COMBINED PLOT FOR APPENDIX---------------------------------------

agetemp_b1<-tottemp_b1_plot+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        legend.text = element_text(size=70))

agetemp_b0<-tottemp_b0_plot+
  labs( x= "Parameter value")+
  theme(axis.title.y = element_blank())
  
agetemp_c<- tottemp_c_plot+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


agetemp_surv<- agetemp_survival_plot+
  labs( y = "Survival, S(x)")+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "none")

agetemp_mort<- agetemp_mort_plot+
  labs( y = "Hazard rate, μ(x)")+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "none")


#-------Combination plot for appendix------------------
tempage_inference<-ggarrange(
  ggarrange(agetemp_b1,agetemp_b0, agetemp_c, ncol= 3, labels = c("b1", "b0", "c"),
            font.label = list(size = 70, face = "bold"),
            label.y = 1.01, 
            label.x = .05,
            heights = c(0.8, 1, 0.8),
            widths=c(1, 0.9, 0.9),
            align = "h"),
  ggarrange(agetemp_surv, agetemp_mort, ncol= 2, labels = c("",""),
            font.label = list(size = 50, face = "bold"),
            heights = c(0.94,  1),align = "h"),
  nrow=2,
  align = "v",
  labels = c("", ""),
  heights = c(0.9, 1))


ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/003_tempageinference.png",
       plot = tempage_inference, 
       device = "png", 
       width = 1000, 
       height = 1000, 
       units = "mm")

# ==============================================================================
# 18. ADDING COVARIATTES: INTERACTIVE EFFECT OF PARENTAL AGE & OFFSPRING SEX----
# ==============================================================================
#Uses a different subset of the data (filtering out animals that weren't sexed, i.e those that died before adulthood)

#18.1.----------------------------
sexdata<-totcovdata %>% 
  filter(F1_sex %in% c("F", "M"))
length(unique(sexdata$ID))
#should be 947 animals for this subset analysis

#Manually setting the interaction term for the KM plot
KMdataagesex<- totcovdata %>%
  mutate(sex_ParentalAge = interaction(Timepoint_binnedplot, F1_sex))

#18.2. Kaplain-Maier curve for parental age------------

#KM survival model
km_fit_parentalagesex<- survfit(Surv(Lifespan_weeks,event) ~ sex_ParentalAge, data = KMdataagesex)
#KM plot
Km_survival_parentalagesex <- ggsurvplot(
  km_fit_parentalagesex,
  data = KMdataagesex,
  pval = FALSE,
  conf.int = TRUE,
  xlab = "Offspring age (weeks)",
  ylab = "Cumulative survival probability, S(x)",
  ggtheme = theme_classic() +
    theme(
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      axis.text.x = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      strip.text = element_text(size = 15),
      legend.position = "none",
    ),
  palette = c("#f96161", "#66b2b2", "#066594", "#f96162", "#66b2b3", "#066593"),
  facet.by = "F1_sex",
  ncol=1,
  short.panel.labs =TRUE,
  legend.labs=NULL)

Km_survival_parentalagesex<- Km_survival_parentalagesex + theme(legend.position = "none")

# Save the combined plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/005_KM_parentalagesex_totallifespan.png",
       plot = Km_survival_parentalagesex,
       bg = "transparent",
       device = "png",
       width = 380,
       height = 340,
       units = "mm")


#This analysis removes 40 animals that were included in the first two (since these did not have a sex included)

#18.3. PARAMETRIC SURVIVAL MODELS FOR PARENTAL AGE X OFFSPRING SEX-----------------------------------------

#18.3.1. A weibull model with just a makeham term------
sexagebastaweibull.2<- basta(object = sexdata, 
                            dataType = "census",
                            formulaMort= ~F1_sex:Timepoint_binned-1, 
                            model="WE",
                            shape="Makeham",
                            niter=60000,#trying to double the burnin and the number of iterations
                            burnin=1001,
                            thinning=50,
                            nsim = 4, parallel = TRUE, ncpus = 4)
summary(sexagebastaweibull.2)
saveRDS(sexagebastaweibull.2, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_sexage_weibullmakeham.rda")
sexagebastaweibull.2<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_sexage_weibullmakeham.rda")

#Assessing model fit and plotting the model outputs
plot(sexagebastaweibull.2, plot.type = "gof")
plot(sexagebastaweibull.2, plot.type = "demorates")
plot(sexagebastaweibull.2, densities=TRUE)

#18.3.2.Weibull model--------------------------------
sexagebastaweibull<- basta(object = sexdata, 
                             dataType = "census",
                             formulaMort= ~F1_sex:Timepoint_binned-1, 
                             model="WE",
                             niter=60000,#trying to double the burnin and the number of iterations
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(sexagebastaweibull)
saveRDS(sexagebastaweibull, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_sexage_weibull.rda")
sexagebastaweibull<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_sexage_weibull.rda")

#18.3.4.Weibull model with a bathtub family-----------------
sexagebastaweibull.3<- basta(object = sexdata, 
                           dataType = "census",
                           formulaMort= ~F1_sex:Timepoint_binned-1, 
                           model="WE",
                           shape="bathtub",
                           niter=60000,#trying to double the burnin and the number of iterations
                           burnin=1001,
                           thinning=50,
                           nsim = 4, parallel = TRUE, ncpus = 4)
summary(sexagebastaweibull.3)
saveRDS(sexagebastaweibull.3, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_sexage_weibull_bathtub.rda")
sexagebastaweibull.3<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_sexage_weibull_bathtub.rda")

#Selecting the Weibull Makeham model


# ==============================================================================
# 19. ADDING WEAKLY INFORMATIVE PRIORS TO THE PARENTAL AGE X OFFSPRING SEX MODEL----
# ==============================================================================

#19.1. Adding weak priors into the model------------------------------

#Weakly informative priors (mean)
weakMean4 <- matrix(c(0, 1.5, prior_scale_mu), 
                    nrow = 6, ncol = 3, byrow = TRUE)

#weakly informative priors (SD)
weakSd4 <- matrix(c(
  1.0, 1.0, 0.5
), nrow = 6, ncol = 3, byrow = TRUE)

weakLower4<-matrix(c(
  0, 0, 0
), nrow = 6, ncol = 3,byrow = TRUE)


#19.2. Fitting weak priors to the weibull makeham model-----------
sexagebastaweibull.weak<- basta(object = sexdata, 
                                 dataType = "census",
                                 formulaMort=~F1_sex:Timepoint_binned-1, 
                                 thetaPriorMean = weakMean4,
                                 thetaPriorSd = weakSd4,
                                 thetaPriorLower = weakLower4,
                                 model="WE",
                                 shape="Makeham",
                                 niter=80000,
                                 burnin=1001,
                                 thinning=50,
                                 nsim = 4, parallel = TRUE, ncpus = 4)
summary(sexagebastaweibull.weak)
saveRDS(sexagebastaweibull.weak, file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_age_sexweak.rda")
sexagebastaweibull.weak<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/total_lifespan_age_sexweak.rda")

#Assessing model fit and plotting the model outputs
plot(sexagebastaweibull.weak, plot.type = "gof")
plot(sexagebastaweibull.weak, plot.type = "demorates")
plot(sexagebastaweibull.weak, densities=TRUE)

# ============================================================================================
# 20. PLOTTING THE POSTERIOR DISTRIBUTIONS FOR MORTALITY PARAMETERS AND DEMOGRAPHIC CURVES----
# ============================================================================================

#20.1. POSTERIOR DISTRIBUTIONS FOR B0, B1 AND C----------------------------------------------------

#Extracting the model parameters
sexage_theta_params <- sexagebastaweibull.weak$params

#Converting to a data frame
sexage_theta_params<-as.data.frame(sexage_theta_params)

# Reshaping the data
sexageposterior_df_long <- sexage_theta_params %>%
  pivot_longer(
    cols = everything(),
    names_to = "Category",
    values_to = "Value"
  ) %>%
  mutate(
    # Extracting F1_sex (Female or Male)
    F1_sex = str_extract(Category, "F1_sexF|F1_sexM"),
    
    # Extracting Parental_Age (early, middle, or late)
    Parental_Age = str_extract(Category, "early|middle|late"),
    
    # Extracting Type (everything before 'sexage')
    Type = str_extract(Category, "^.*(?=.F1_sex)")
  ) %>% 
  mutate(
    Parental_Age = case_when(
      Parental_Age == "early" ~ "Early-Aged",
      Parental_Age== "middle"   ~ "Middle-Aged",
      Parental_Age == "late" ~ "Late-Aged",
      TRUE ~ Parental_Age
    )) %>%
  mutate(F1_sex = case_when(
    F1_sex == "F1_sexF" ~ "Female",
    F1_sex== "F1_sexM"   ~ "Male",
    TRUE ~ F1_sex
  )) %>% 
  select(Type, F1_sex, Parental_Age, Value)


#making sure the levels are in  the right order for the plot:
sexageposterior_df_long$Parental_Age<- factor(sexageposterior_df_long$Parental_Age, 
                                           levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))



#20.2. Creating seperate data sets for c, `b0` and `b1`
sexageb0_df <- sexageposterior_df_long  %>% filter(Type == "b0")
sexageb1_df <- sexageposterior_df_long  %>% filter(Type == "b1")
sexagecdf<- sexageposterior_df_long %>% filter(Type =="c")


#20.3. Posterior plots for the b0 parameter (the shape of ageing)-------------------
sexage_b0_plot <- ggplot(sexageb0_df, aes(x = Value, fill = Parental_Age, color = Parental_Age, alpha=Parental_Age)) +
  geom_density(size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("Early-Aged" = "#f96161", 
                               "Middle-Aged" = "#66b2b2", 
                               "Late-Aged" = "#066594")) +
  scale_color_manual(values = c("Early-Aged" = "#f96161", 
                                "Middle-Aged" = "#66b2b2", 
                                "Late-Aged" = "#066594")) +
  scale_alpha_manual(values = c("Early-Aged" = 0.6, 
                                "Middle-Aged" = 0.1, 
                                "Late-Aged" = 0.6))+
  scale_x_continuous(limits = c(4.4, 7.6), breaks = c(4.5, 5.5, 6.5, 7.5))+
  facet_wrap(~F1_sex, ncol=1)+
  labs(
    x = "b0 parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' age at reproduction",
    color = "Parents' age at reproduction",
    alpha="Parents' age at reproduction") +
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_text(size = 35, face = "bold"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
  )

#saving plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/005_b0_sexandage_totallifespan.png",
       plot = sexage_b0_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 460, 
       units = "mm")


#20.4. Posterior plot for the b1 parameter values (the rate of ageing)---------------
sexage_b1_plot <- ggplot(sexageb1_df , aes(x = Value, fill = Parental_Age, color = Parental_Age, alpha=Parental_Age)) +
  geom_density(size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("Early-Aged" = "#f96161", 
                               "Middle-Aged" = "#66b2b2", 
                               "Late-Aged" = "#066594")) +
  scale_color_manual(values = c("Early-Aged" = "#f96161", 
                                "Middle-Aged" = "#66b2b2", 
                                "Late-Aged" = "#066594")) +
  
  scale_alpha_manual(values = c("Early-Aged" = 0.6, 
                                "Middle-Aged" = 0.1, 
                                "Late-Aged" = 0.6))+
  scale_x_continuous(limits = c(2.39, 2.81), breaks = c(2.4, 2.5, 2.6, 2.7, 2.8))+
  facet_wrap(~F1_sex, ncol=1)+
  labs(
    x = "b1 parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' age at reproduction",
    color = "Parents' age at reproduction",
    alpha="Parents' age at reproduction") +
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_text(size = 35, face = "bold"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
  )

ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/005_b1_sexandage_totallifespan.png",
       plot = sexage_b1_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 460, 
       units = "mm")



#20.5. creating a posterior plot for the c parameter value------
sexage_c_plot <- ggplot(sexagecdf, aes(x = Value, fill = Parental_Age, color = Parental_Age, alpha=Parental_Age)) +
  geom_density(size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("Early-Aged" = "#f96161", 
                               "Middle-Aged" = "#66b2b2", 
                               "Late-Aged" = "#066594")) +
  scale_color_manual(values = c("Early-Aged" = "#f96161", 
                                "Middle-Aged" = "#66b2b2", 
                                "Late-Aged" = "#066594")) +
  
  scale_alpha_manual(values = c("Early-Aged" = 0.6, 
                                "Middle-Aged" = 0.1, 
                                "Late-Aged" = 0.6))+
  facet_wrap(~F1_sex, ncol=1)+
  labs(
    x = "c parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' age at reproduction",
    color = "Parents' age at reproduction",
    alpha="Parents' age at reproduction") +
  scale_x_continuous(limits = c(0, 0.21), breaks = c(0.00, 0.05, 0.10, 0.15, 0.20))+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_text(size = 35, face = "bold"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
  )


ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/004_posteriors_c_sexandage.png",
       plot = sexage_c_plot, 
       bg="transparent",
       device = "png", 
       width = 310, 
       height = 420, 
       units = "mm")




#20.6. PARENTAL AGE + OFFSPRING SEX: Plotting the survival and mortality curves from the parental age + offspring sex analysis--------

#20.6.1. Function to extract mortality parameters-------------
plot_data_agesex <- list()
for (demv in c("mort", "surv")) {
  for (icat in seq_along(sexagebastaweibull.weak[[demv]])) {
    cuts <-sexagebastaweibull.weak$cuts[[icat]]
    minAge <- as.numeric(sexagebastaweibull.weak$modelSpecs["min. age"])
    xx <- sexagebastaweibull.weak$x[cuts] + minAge
    yy <- sexagebastaweibull.weak[[demv]][[icat]][, cuts]
    
    # Convert to long format
    df <- data.frame(
      Age = xx,
      Rate = yy[1, ], # Median
      LowerCI = yy[2, ], # Lower Credible interval Bound
      UpperCI = yy[3, ], # Upper Credible interval Bound
      Category = names(sexagebastaweibull.weak[[demv]])[icat],
      Type = demv # Mortality or Survival
    )
    plot_data_agesex[[length(plot_data_agesex) + 1]] <- df
  }
}

# Combine all into one data frame
plot_data_agesex<- do.call(rbind, plot_data_agesex)



#splitting parental age and temperature into two seperate columns
plot_data_agesex <- plot_data_agesex %>%
  mutate(
    F1_sex = str_extract(Category, "F1_sexF|F1_sexM"),  # Extract temperature part
    Parental_Age = str_extract(Category, "early|middle|late")  # Extract parental age
  ) %>%
  mutate(
    F1_sex = factor(F1_sex, levels = c("F1_sexF", "F1_sexM"), 
                    labels = c("Female", "Male")),
    Parental_Age = factor(Parental_Age, 
                          levels = c("early", "middle", "late"),
                          labels = c("Early-Aged", "Middle-Aged", "Late-Aged"))
  )


#making sure the levels are in  the right order for the plot:
plot_data_agesex$Parental_Age<- factor(plot_data_agesex$Parental_Age, 
                                       levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))


# Set the scaling factors (e.g., minimum and maximum age in days)
max_days <- 365.25 
min_days<-0

# Adjust Age column back to days
plot_data_agesex<- plot_data_agesex %>%
  mutate(Age_days = Age * (max_days - min_days) + min_days)
plot_data_agesex$Parental_Age<-as.factor(plot_data_agesex$Parental_Age)


#20.6.2. Split the data into Mortality and Survival datasets--------
mort_df_agesex <- plot_data_agesex %>% filter(Type == "mort") %>% 
  mutate(Rate_weeks = Rate/52.1429,
         LowerCI_weeks = LowerCI/52.1429,
         UpperCI_weeks = UpperCI/52.1429)
surv_df_agesex<- plot_data_agesex%>% filter(Type == "surv")

#20.6.3. Median survival probability----------------
surv_df_median<-surv_df_agesex %>%
  group_by(F1_sex, Parental_Age) %>%  
  arrange(Age_days) %>%
  mutate(CI = 1 - Rate) %>%  
  summarise(median_line = approx(CI, Age_days, xout = 0.5)$y)


#20.6.4. cumulative survival plot---------------------
sexage_survival_plot<- ggplot(data = surv_df_agesex,
                               aes(x = Age_days/7, 
                                   y = Rate,
                                   colour = Parental_Age))+
  geom_line(data = surv_df_agesex, 
            aes(x = Age_days/7, 
                y = Rate,
                colour = Parental_Age),
            linewidth=4, alpha=1)+
  geom_ribbon(data = surv_df_agesex, 
              aes(y=NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI, 
                  color= Parental_Age,
                  fill=Parental_Age, 
                  alpha = Parental_Age), 
              linewidth = 1, 
              linetype="dashed") +
  geom_segment(data = surv_df_median,
               aes(x = 0, xend = median_line/7, y = 0.5, yend = 0.5, color = Parental_Age),
               linetype = "dashed", size = 1.75, alpha = 0.9) +
  geom_segment(data = surv_df_median,
               aes(x = median_line/7, xend = median_line/7, y = 0, yend = 0.5, color = Parental_Age),
               linetype = "dashed", size = 1.75, alpha = 0.9)+
  facet_wrap(~F1_sex, ncol=1)+
  theme_classic()+ 
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_text(size = 35, face = "bold"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
  )+
  scale_color_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' age at reproduction", values = c(0.5,0.1, 0.5))+
  scale_fill_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  labs(x = "Offspring age (weeks)", y = "Cumulative survival probability, S(x)")


#save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/005_survival_sexandage_totallifespan.png",
       plot = sexage_survival_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 460, 
       units = "mm")



#20.6.5. Mortality risk plot-----------------------
agesex_mort_plot<- ggplot(data = mort_df_agesex,
                          aes(x = Age_days/7, 
                              y = Rate_weeks,
                              colour = Parental_Age))+
  geom_line(data = mort_df_agesex, 
            aes(x = Age_days/7, 
                y = Rate_weeks,
                color = Parental_Age),
            linewidth=4, alpha=1)+
  geom_ribbon(data = mort_df_agesex, 
              aes(y=NULL, 
                  ymin = LowerCI_weeks, 
                  ymax = UpperCI_weeks, 
                  color= Parental_Age,
                  fill= Parental_Age, 
                  alpha = Parental_Age),
              linewidth= 1,
              linetype="dashed") +
  facet_wrap(~F1_sex, ncol=1)+
  theme_classic()+ 
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_text(size = 35, face = "bold"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5),
  )+
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
  scale_color_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' age at reproduction", values = c(0.5,0.1, 0.5))+
  scale_fill_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  labs(x = "Offspring age (weeks)", y = "Instantaneous hazard rate, μ(x)")

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/005_mortality_sexandage_totallifespan.png",
       plot = agesex_mort_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 460, 
       units = "mm")


#--------------COMBINED PLOT FOR APPENDIX---------------------------------------


agesex_b1<-sexage_b1_plot+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        legend.text = element_text(size=70))

agesex_b0<-sexage_b0_plot+
  labs( x= "Parameter value")+
  theme(axis.title.y = element_blank())

agesex_c<- sexage_c_plot+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


agesex_surv<- sexage_survival_plot+
  labs( y = "Survival, S(x)")+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "none")
theme(legend.position = "none")

agesex_mort<- agesex_mort_plot+
  labs( y = "Hazard rate, μ(x)")+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "none")


#-------Combination plot for appendix------------------
sexage_inference<-ggarrange(
  ggarrange(agesex_b1,agesex_b0,agesex_c, ncol= 3, labels = c("b1", "b0", "c"),
            font.label = list(size = 70, face = "bold"),
            label.y = 1.01, 
            label.x = .01,
            heights = c(0.8, 1, 0.8),
            widths=c(1, 0.9, 0.9),
            align = "h"),
  ggarrange(agesex_surv, agesex_mort, ncol= 2, labels = c("",""),
            font.label = list(size = 50, face = "bold"),
            heights = c(0.94,  1),align = "h"),
  nrow=2,
  align = "v",
  labels = c("", ""),
  heights = c(0.9, 1))


ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/004_sexageinference.png",
       plot = sexage_inference, 
       device = "png", 
       width = 900, 
       height = 800, 
       units = "mm")

###############################################################################
##END OF MAIN SCRIPT----------------------------------------------------------------
###############################################################################






#BONUS ANALYSIS

#6. Part 3: The independent effect of offspring sex on offspring longevity####

#6.1. creating a new variable that groups individuals by parental age and parental temperature (creates a new grouped variable for the interaction)
sexdata<-totcovdata %>% 
  filter(F1_sex %in% c("F", "M"))
length(unique(sexdata$ID))#947 animals from 77 parent pairs-->this filters out the 40 animals that died before reaching adulthood



#1.Kaplain-Maier curve for offspring sex
km_fit_sex<- survfit(Surv(Lifespan_weeks,event) ~ F1_sex, data = sexdata)#fitting the curve for Timepoint
Km_survival_sex<- ggsurvplot(km_fit_sex, 
                             data = totcovdata, 
                             pval = FALSE,  
                             conf.int = TRUE,  
                             risk.table = "abs_pct",
                             break.time.by = 5,
                             surv.median.line = "hv",  
                             xlab = "Offspring age (weeks)", 
                             ylab = "Cumulative survival probability, S(x)",
                             theme = theme_classic(),
                             font.x = 35,
                             font.y = 35,
                             font.tickslab = c(35, "grey25"),
                             font.legend = 25,
                             risk.table.fontsize = 6,
                             tables.theme = theme_classic() +
                               theme(
                                 axis.title.x = element_text(size = 25), 
                                 axis.title.y = element_blank(),  
                                 axis.text.x = element_text(size = 25),  
                                 axis.text.y = element_text(size = 25), 
                                 legend.position = "none",
                                 legend.title = element_blank(),  
                                 legend.text = element_blank(),   
                                 strip.text = element_text(size = 25)
                               ),
                             legend.title = "Offspring sex",
                             legend.labs = c("Female", "Male"),
                             palette = c("#A8C2A1","#6A7FDB"))

# Combine the survival plot and the risk table
Km_sex_combined_plot <- cowplot::plot_grid(Km_survival_sex$plot, Km_survival_sex$table, ncol = 1, rel_heights = c(3, 1))

# Save the combined plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/003_KM_offspringsex_totallifespan.png",
       plot = Km_sex_combined_plot,
       bg = "transparent",
       device = "png",
       width = 420,
       height = 360,
       units = "mm")


############Statistical models###########

#6.2. Weibull model with Makeham term
totalbastaweibullsex.2 <- basta(object = sexdata, 
                                dataType = "census",
                                formulaMort=~F1_sex-1, 
                                model="WE",
                                shape="Makeham",
                                niter=60000,#trying to double the burnin and the number of iterations
                                burnin=1001,
                                thinning=50,
                                nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibullsex.2)
saveRDS(totalbastaweibullsex.2, file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_weibull_Makeham.rda")
totalbastaweibullsex.2<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_weibull_Makeham.rda")#this is with the makeham shape to be confusing 


#Assessing model fit and plotting the model outputs
plot(totalbastaweibullsex.2, plot.type = "gof")
plot(totalbastaweibullsex.2, plot.type = "demorates")
plot(totalbastaweibullsex.2, densities=TRUE)


#6.3 weibull model
totalbastaweibullsex <- basta(object = sexdata, 
                              dataType = "census",
                              formulaMort=~F1_sex-1, 
                              model="WE",
                              niter=60000,#trying to double the burnin and the number of iterations
                              burnin=1001,
                              thinning=50,
                              nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibullsex)
saveRDS(totalbastaweibullsex, file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_weibull.rda")
totalbastaweibullsex<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_weibull.rda")


#Assessing model fit and plotting the model outputs
plot(totalbastaweibullsex, plot.type = "gof")
plot(totalbastaweibullsex, plot.type = "demorates")
plot(totalbastaweibullsex, densities=TRUE)


#6.4 weibull model with a bathtub shape
totalbastaweibullsex.3 <- basta(object = sexdata, 
                                dataType = "census",
                                formulaMort=~F1_sex-1, 
                                model="WE",
                                shape="bathtub",
                                niter=60000,#trying to double the burnin and the number of iterations
                                burnin=1001,
                                thinning=50,
                                nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibullsex.3)
saveRDS(totalbastaweibullsex.3, file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_weibull.3.rda")
totalbastaweibullsex.3<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_weibull.3.rda")



#6.5 weibull model
totalbastagompsex<- basta(object = sexdata, 
                          dataType = "census",
                          formulaMort=~F1_sex-1, 
                          niter=60000,#trying to double the burnin and the number of iterations
                          burnin=1001,
                          thinning=50,
                          nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastagompsex)
saveRDS(totalbastagompsex, file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_gompertz.rda")
totalbastagompsex<-readRDS("scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_gompertz.rda")

#Assessing model fit and plotting the model outputs
plot(totalbastagompsex, plot.type = "gof")
plot(totalbastagompsex, plot.type = "demorates")
plot(totalbastagompsex, densities=TRUE)



#6.6. BaSTA gompertz model with makeham shape
totalbastagompsex.2<- basta(object = sexdata, 
                            dataType = "census",
                            formulaMort=~F1_sex-1, 
                            shape="Makeham",
                            niter=60000,#trying to double the burnin and the number of iterations
                            burnin=1001,
                            thinning=50,
                            nsim = 4, parallel = TRUE, ncpus = 4)
saveRDS(totalbastagompsex.2, file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_gompertz_makeham.rda")
totalbastagompsex.2<-readRDS("scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_gompertz_makeham.rda")

#Assessing model fit and plotting the model outputs
plot(sexbastagomptemp.2, plot.type = "gof")
plot(totalbastagomptemp.2, plot.type = "demorates")
plot(totalbastagomptemp.2, densities=TRUE)



#6.7 gompertz model with a bathtub shape
totalbastagompsex.3 <- basta(object = sexdata, 
                             dataType = "census",
                             formulaMort=~F1_sex-1, 
                             model="GO",
                             shape="bathtub",
                             niter=60000,#trying to double the burnin and the number of iterations
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastagompsex.3)
saveRDS(totalbastagompsex.3, file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_gompertz.3.rda")
totalbastagompsex.3 <-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_gompertz.3.rda")



#6.7. Logarithmic model
totalbastalogsex <- basta(object = sexdata, 
                          dataType = "census",
                          formulaMort=~F1_sex-1, 
                          model="LO",
                          niter=60000,#trying to double the burnin and the number of iterations
                          burnin=1001,
                          thinning=50,
                          nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastalogsex)
saveRDS(totalbastalogsex, file = "scripts/model_outputs/BaSTA/total lifespan/Sex_total_lifespan_log.rda")
totalbastalogsex<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/Sex_total_lifespan_log.rda")#this is with the makeham shape to be confusing 


#Assessing model fit and plotting the model outputs
plot(totalbastalogsex, plot.type = "gof")
plot(totalbastalogsex, plot.type = "demorates")
plot(totalbastalogsex, densities=TRUE)


#6.8.Logarithmic model with a makeham shape
totalbastalogsex.2 <- basta(object = sexdata, 
                            dataType = "census",
                            formulaMort=~F1_sex-1, 
                            model="LO",
                            shape="Makeham",
                            niter=60000,#trying to double the burnin and the number of iterations
                            burnin=1001,
                            thinning=50,
                            nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastalogsex.2)
saveRDS(totalbastalogsex.2, file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_log_makeham.rda")
totalbastalogsex.2<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_log_makeham.rda")#this is with the makeham shape to be confusing 



#6.9. logarithmic model with a bathtub shape
totalbastalogsex.3 <- basta(object = sexdata, 
                            dataType = "census",
                            formulaMort=~F1_sex-1, 
                            model="LO",
                            shape="bathtub",
                            niter=60000,#trying to double the burnin and the number of iterations
                            burnin=1001,
                            thinning=50,
                            nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastalogsex.3)
saveRDS(totalbastalogsex.3, file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_log_bathtub.rda")
totalbastalogsex.3<-readRDS(file = "scripts/model_outputs/BaSTA/total lifespan/sex_total_lifespan_log_bathtub.rda")#this is with the makeham shape to be confusing 




#6.9.1 DIC values
#Comparing fits of all models

# Extract DIC for multiple models
sex_dic_weibull <- totalbastaweibullsex$DIC["DIC"]
sex_dic_weibull_2 <- totalbastaweibullsex.2$DIC["DIC"]
sex_dic_weibull_3 <- totalbastaweibullsex.3$DIC["DIC"]
sex_dic_gompertz <- totalbastagompsex$DIC["DIC"]
sex_dic_gompertz_2 <- totalbastagompsex.2$DIC["DIC"]
sex_dic_total_gompertz_3 <- totalbastagompsex.3$DIC["DIC"]
sex_dic_total_log<-totalbastalogsex$DIC["DIC"]
sex_dic_total_log_2<-totalbastalogsex.2$DIC["DIC"]
sex_dic_total_log_3<-totalbastalogsex.3$DIC["DIC"]


# Comparing the DIC values
totalsex_dic_values <- data.frame(
  Model = c("totalbastalog", "totalbastalog.2", "totalbastalog.3",  "totalbastaweibull", "totalbastaweibull.2", "totalbastaweibull.3", 
            "totalgomptzmod", "totalgomptzmod.2", "totalgomp.3"),
  DIC = c(sex_dic_total_log, sex_dic_total_log_2, sex_dic_total_log_3, sex_dic_weibull, sex_dic_weibull_2 , sex_dic_weibull_3, 
          sex_dic_gompertz, sex_dic_gompertz_2, sex_dic_total_gompertz_3 )
)

print(totalsex_dic_values)

#6.9.2. Offspring sex: Plotting the Weibull with a Makeham shape#########

#Plotting the posterior distributions from the weibull model with a makeham shape
sextotal_theta_params <- totalbastaweibullsex$params
sextotal_theta_params<-as.data.frame(sextotal_theta_params)

# Reshaping the data
sextotal_posterior_df <- sextotal_theta_params %>%
  pivot_longer(
    cols = everything(),
    names_to = "Category",
    values_to = "Value"
  ) %>%
  mutate(
    # Extracting F1_sex (Female or Male)
    F1_sex = str_extract(Category, "F1_sexF|F1_sexM"),
    
    # Extracting Type (everything before the offspring's sex)
    Type = str_extract(Category, "^.*(?=F1_sex)")
  ) %>% 
  mutate(
    F1_sex = case_when(
      F1_sex == "F1_sexF" ~ "Female",
      F1_sex== "F1_sexM"   ~ "Male",
      TRUE ~ F1_sex
    )) %>%
  select(Type, F1_sex, Value)

# Create two separate datasets for `b0` and `b1`, and c:
sextotb0_df <- sextotal_posterior_df %>% filter(Type == "b0.")
sextotb1_df <- sextotal_posterior_df %>% filter(Type == "b1.")


#creating posterior plot for b0 values:
sextot_b0_plot <- ggplot(sextotb0_df, aes(x = Value, fill = F1_sex, color = F1_sex, alpha=F1_sex)) +
  geom_density(size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("Male" = "#6A7FDB", 
                               "Female"= "#A8C2A1")) +
  scale_color_manual(values = c("Male" = "#6A7FDB", 
                                "Female"= "#A8C2A1"))+
  scale_alpha_manual(values = c("Male" = 0.6, 
                                "Female" = 0.6))+
  scale_x_continuous(limits = c(3.5, 8), breaks = c(4, 5, 6, 7, 8))+
  labs(
    x = "b0 parameter value",
    y = "Density",
    fill = "Offspring sex",
    color = "Offspring sex",
    alpha="Offspring sex") +
  theme(legend.position = "top",
        axis.title = element_text(size = 35),
        axis.text = element_text(size = 35),
        legend.title=element_text(size=35),
        legend.text = element_text(size = 35))

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/003_b0_offspringsex_totallifespan.png",
       plot = sextot_b0_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")

#creating a posterior plot for the b1 parameter values
sextot_b1_plot <- ggplot(sextotb1_df, aes(x = Value, fill = F1_sex, color = F1_sex, alpha=F1_sex)) +
  geom_density(size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("Male" = "#6A7FDB", 
                               "Female"= "#A8C2A1")) +
  scale_color_manual(values = c("Male" = "#6A7FDB", 
                                "Female"= "#A8C2A1"))+
  scale_alpha_manual(values = c("Male" = 0.6, 
                                "Female" = 0.6))+
  scale_x_continuous(limits = c(2.4, 2.8), breaks = c(2.4, 2.5, 2.6, 2.7, 2.8))+
  labs(
    x = "b1 parameter value",
    y = "Density",
    fill = "Offspring sex",
    color = "Offspring sex",
    alpha="Offspring sex") +
  theme(legend.position = "top",
        axis.title = element_text(size = 35),
        axis.text = element_text(size = 35),
        legend.title=element_text(size=25),
        legend.text = element_text(size = 25))

ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/003_b1_Offspringsex_totallifespan.png",
       plot = sextot_b1_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")



#6.9.3. POffspring sex: Plotting the survival and mortality curves for The weibull model with makeham shape has the lowest DIC value#######
plot_data_total <- list()
for (demv in c("mort", "surv")) {
  for (icat in seq_along(totalbastaweibullsex[[demv]])) {
    cuts <-totalbastaweibullsex$cuts[[icat]]
    minAge <- as.numeric(totalbastaweibullsex$modelSpecs["min. age"])
    xx <- totalbastaweibullsex$x[cuts] + minAge
    yy <- totalbastaweibullsex[[demv]][[icat]][, cuts]
    
    # Convert to long format
    df <- data.frame(
      Age = xx,
      Rate = yy[1, ], # Median
      LowerCI = yy[2, ], # Lower Confidence Bound
      UpperCI = yy[3, ], # Upper Confidence Bound
      Category = names(totalbastaweibullsex[[demv]])[icat],
      Type = demv # Mortality or Survival
    )
    plot_data_total[[length(plot_data_total) + 1]] <- df
  }
}

# Combine all into one data frame
plot_data_total<- do.call(rbind, plot_data_total)

# Set the scaling factors (e.g., minimum and maximum age in days)
max_days <- 365.25 
min_days<-0

# Adjust Age column back to days
plot_data_total<- plot_data_total %>%
  mutate(Age_days = Age * (max_days - min_days) + min_days) %>% 
  mutate(
    Category = case_when(
      Category == "F1_sexF" ~ "Female",
      Category== "F1_sexM"   ~ "Male",
      TRUE ~ Category
    ))

plot_data_total$Category<-as.factor(plot_data_total$Category)


# Split the data into Mortality and Survival datasets
mort_df_total <- plot_data_total %>% filter(Type == "mort")
surv_df_total<- plot_data_total%>% filter(Type == "surv")

#3.4 Survival probability###
surv_df_median<-surv_df_total %>%
  group_by(Category) %>%  
  arrange(Age_days) %>%
  mutate(CI = 1 - Rate) %>%  # This assumes Rate is a survival function; adjust if it's already cumulative incidence
  summarise(median_line = approx(CI, Age_days, xout = 0.5)$y)

#5.2. cumulative survival plot
sextotal_survival_plot<- ggplot(data = surv_df_total,
                                aes(x = Age_days/7, #this suggests that animals are developing into totals earlier than expected?
                                    y = Rate,
                                    colour = Category))+
  geom_line(data = surv_df_total, 
            aes(x = Age_days/7, 
                y = Rate,
                colour = Category),
            size=1.5, alpha=1)+
  geom_ribbon(data = surv_df_total, 
              aes(y=NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI, 
                  color= Category,
                  fill=Category,
                  alpha=Category),
              linetype="dashed") +
  theme_classic()+ 
  geom_segment(data = surv_df_median,
               aes(x = 0, xend = median_line/7, y = 0.5, yend = 0.5, color = Category),
               linetype = "dashed", size = 1.75, alpha = 0.9) +
  geom_segment(data = surv_df_median,
               aes(x = median_line/7, xend = median_line/7, y = 0, yend = 0.5, color = Category),
               linetype = "dashed", size = 1.75, alpha = 0.9)+
  theme(axis.title=element_text(size=35),
        axis.text=element_text(size=35))+
  theme(legend.position = "top",
        legend.title=element_text(size=35),
        legend.text=element_text(size=35),
        strip.text = element_text(size=35)) +
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
  scale_color_manual(name = "Offspring sex", values = c("#A8C2A1","#6A7FDB"))+
  scale_fill_manual(name = "Offspring sex", values = c("#A8C2A1","#6A7FDB"))+
  scale_alpha_manual(name = "Offspring sex", values = c(0.5,0.5))+
  labs(x = "Offspring age (weeks)", y = "Cumulative survival probability, s(X)")


# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/003_survival_offspringsex_totallifespan.png",
       plot = sextotal_survival_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 220, 
       units = "mm")


#5.2: mortality risk plot

#Logged mortality risk
sexlog_total_mort_plot<- ggplot(data = mort_df_total,
                                aes(x = Age_days/7, #this suggests that animals are developing into adults earlier than expected?
                                    y = log(Rate+ 1e-6),
                                    colour = Category))+
  geom_line(data = mort_df_total, 
            aes(x = Age_days/7, 
                y = log(Rate+ 1e-6),
                color = Category),
            size=1.5, alpha=1)+
  geom_ribbon(data = mort_df_total, 
              aes(y=NULL, 
                  ymin = log(LowerCI+ 1e-6), 
                  ymax = log(UpperCI+ 1e-6), 
                  color= Category,
                  fill=Category), 
              alpha = .05, linetype="dashed") +
  theme_classic()+ 
  theme(axis.title=element_text(size=35),
        axis.text=element_text(size=35))+
  theme(legend.position = "top",
        legend.title=element_text(size=35),
        legend.text=element_text(size=35),
        strip.text = element_text(size=35)) +
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
  scale_color_manual(name = "Offspring sex", values = c("#A8C2A1","#6A7FDB"))+
  scale_fill_manual(name = "Offspring sex", values = c("#A8C2A1","#6A7FDB"))+
  scale_alpha_manual(name = "Offspring sex", values = c(0.5,0.5))+
  labs(x = "Offspring age (weeks)", y = "log (mortality risk), log (μ (X))")

#to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/003_logmortality_offspringsex_totallifespan.png",
       plot = sexlog_total_mort_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 220, 
       units = "mm")



#5.3. mortality risk (not logged)
sextotal_mort_plot<- ggplot(data = mort_df_total,
                            aes(x = Age_days/7, 
                                y = Rate,
                                colour = Category))+
  geom_line(data = mort_df_total, 
            aes(x = Age_days/7, 
                y = Rate,
                color = Category),
            size=1.5, alpha=1)+
  geom_ribbon(data = mort_df_total, 
              aes(y=NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI, 
                  color= Category,
                  fill=Category,
                  alpha=Category),
              linetype="dashed") +
  theme_classic()+ 
  theme(axis.title=element_text(size=35),
        axis.text=element_text(size=35))+
  theme(legend.position = "top",
        legend.title=element_text(size=35),
        legend.text=element_text(size=35),
        strip.text = element_text(size=35)) +
  scale_x_continuous(limits = c(0, 28), breaks = c(0, 5, 10, 15, 20, 25))+
  scale_color_manual(name = "Offspring sex", values = c("#A8C2A1","#6A7FDB"))+
  scale_fill_manual(name = "Offspring sex", values = c("#A8C2A1","#6A7FDB"))+
  scale_alpha_manual(name = "Offspring sex", values = c(0.5,0.5))+
  labs(x = "Offspring age (weeks)", y = "Instantaneous hazard rate, μ (X)")

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/003_mortality_offspringsex_totallifespan.png",
       plot = sextotal_mort_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 220, 
       units = "mm")



#Assessing how well the survival curves overlap with the Kaplan-Meier curve estimates####

#Extracting the data from the Kaplan Meier curve
Km_survival_sex$data.survplot<-
  Km_survival_sex$data.survplot %>%
  rename(Category=F1_sex) %>% 
  mutate(Category = case_when(  
    Category == "F" ~ "Female",
    Category == "M"   ~ "Male",
    TRUE ~ Category 
  ))

# Create the zero-point for each group
zero_point_003_sex <- Km_survival_sex$data.survplot %>%
  group_by(Category) %>%
  summarise(time = 0, surv = 1, lower = 1, upper = 1)

# Combine zero-point and the actual KM data
Km_survival_sex$data.survplot<-Km_survival_sex$data.survplot %>%
  bind_rows(zero_point_003_sex) %>%
  arrange(Category, time)  # Make sure points are sorted for step plot


#Combining the parametric and non-parametric estimates into one plot
Km_survival_combined_plot_003<- ggplot() +#
  #The Kaplan-Meier curve
  geom_step(data = Km_survival_sex$data.survplot, 
            aes(x = time, 
                y = surv), 
            colour = "#066594",
            size = 1.2, alpha = 1)+
  #Kaplan-Meir CI
  geom_ribbon(data = Km_survival_sex$data.survplot, 
              aes(x = time, 
                  y = NULL, 
                  ymin = lower, 
                  ymax = upper), 
              fill = "#066594", 
              alpha = 0.3)+ 
  #Parametric survival curve
  geom_line(data = surv_df_total, 
            aes(x = Age_days / 7, 
                y = Rate),
            colour = "#f96161", 
            size = 1.2, alpha = 1) +
  #Parametric survival cI
  geom_ribbon(data = surv_df_total, 
              aes(x = Age_days / 7, 
                  y = NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI),
              fill = "#f96161", 
              alpha = 0.3)+
  facet_wrap(~ Category, ncol = 1) +  # Ensure one column of facets
  theme_classic() +
  theme(axis.title=element_text(size=35),
        axis.text=element_text(size=35),
        legend.position = "",
        strip.text = element_text(size=25))+
  scale_x_continuous(limits = c(0, 32), breaks = c(0, 5, 10, 15, 20, 25, 30))+
  labs(x = "Offspring age (weeks)", y = "Cumulative survival probability, s(X)")

ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/003_modelfit_sex_totallifespan.png",
       plot = Km_survival_combined_plot_003 , 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 380, 
       units = "mm")
