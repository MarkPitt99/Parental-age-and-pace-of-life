###############################################################################
##  Script: BaSTA analysis looking at the adult lifespan of the offspring, 
#Only including animals which had reliably estimated adult emergence dates
#Analysis also answers a core question of our manuscript
#Does increasing parental age deterimentally affect offspring adult baseline or age-specific mortality?
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
F1data<-readRDS("./raw data/BaSTA data/adult_basta_data.RDS")
F1data<-as.data.frame(F1data)
length(unique(F1data$ID))
#939 offspring from 77 parent pairs that had reliably estimated adult emergence dates

#2.2. grouping the parental ages into categorical groups (rather than treating it as a continuous variable)----
F1data <- F1data%>%
  mutate(Timepoint_binned = case_when(
    Timepoint %in% 1:2 ~ "early",
    Timepoint %in% 3:5 ~ "middle",
    Timepoint %in% 6:8 ~ "late"
  ))

#2.3. creating data with the required covariates-------
adcovdata<- F1data %>% 
  select(ID, Birth.Date, Min.Birth.Date, Max.Birth.Date, Entry.Date, Depart.Date, Depart.Type, Timepoint, Timepoint_binned, Temp,F1_sex)

#DATA CHECK: re-checking whether the filtered data passes BaSTA's built-in data check function
checkedDataCens <- DataCheck(object = adcovdata, dataType = "census",
                             silent = FALSE)#No inconsistencies between dates
print(checkedDataCens)#data seems to all be coded correctly, and I've includes all necessary covariates

#Calculating adult lifespan (in weeks)
adcovdata$adult_lifespan <- as.numeric(difftime(adcovdata$Depart.Date, adcovdata$Birth.Date, units = "weeks"))

#How many offspring do we have in each age category at each temperature?
table_temp_lifespan <- table(F1data$Temp, F1data$Timepoint_binned)
print(table_temp_lifespan)
#       early late middle
#25.5    78  120    139
#28      75   89    115
#30.5    95   83    145

#How many offspring per paremt pair at each level?
table_lifespan_pairID<-table(F1data$PairID, F1data$Timepoint_binned)

#how many parents produced offspring in each age category?
table_lifespan_pairID<-F1data %>% 
  group_by(Timepoint_binned) %>% 
  summarise(n_parents = n_distinct(PairID))
  
#  early                   74
#  late                    60
#  middle                  72

# ===========================================================================
# 3. SUMMARY STATISTICS-------------------------------------------------------
# ===========================================================================

#3.1. Median adult lifespans------
sum_dat<-adcovdata %>%
  group_by(Timepoint_binned) %>% 
  summarise(F1_lifespan = median (adult_lifespan),
            sd_adult= sd(adult_lifespan),
            n_animals=n())

#Timepoint_binned F1_lifespan sd_adult   n_animals
#early                    9       2.87       248
#late                    10.4     3.37       292
#middle                   9       3.10       399


#Median adult lifespans per temperature treatment
sum_dat2<-totcovdata%>%
  group_by(Temp)%>%
  summarise(F1_lifespan = median (adult_lifespan),
            sd_adult= sd(adult_lifespan),
            n_animals=n())

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
gomptzmod <- basta(object = adcovdata, 
                   dataType = "census",
                   niter=60000,
                   burnin=1001,
                   thinning=50,
                   nsim = 4, parallel = TRUE, ncpus = 4)
summary(gomptzmod)

#saving the model
saveRDS(gomptzmod, file = "scripts/model_outputs/BaSTA/adult lifespan/lifespan_gompertz.rda")
gomptzmod<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/lifespan_gompertz.rda")

#Assessing model fit and plotting the model outputs
plot(gomptzmod)
plot(gomptzmod , plot.type = "demorates")
plot(gomptzmod, plot.type = "gof")
plot(gomptzmod, densities=TRUE)


#4.2 Gompertz model with a Makeham shape ----
#adds a Makeham constant term (accounts for age-independent background mortality)
gomptzmod.2 <- basta(object = adcovdata, 
                     dataType = "census",
                     shape="Makeham",
                     niter=60000,
                     burnin=1001,
                     thinning=50,
                     nsim = 4, parallel = TRUE, ncpus = 4)
summary(gomptzmod.2)
#Saving the model
saveRDS(gomptzmod.2, file = "scripts/model_outputs/BaSTA/adult lifespan/gompertz_makeham.rda")
gomptzmod.2<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/gompertz_makeham.rda")

#assessing model fit
plot(gomptzmod.2, plot.type = "gof")
plot(gomptzmod.2, plot.type = "demorates")
plot(gomptzmod.2, densities=TRUE)


#4.3. Gompertz Model with a bathtub function-------------------------
#adds two terms to capture the rate of early-life mortality--
gomp.3 <- basta(object = adcovdata, 
                dataType = "census",
                shape="bathtub",
                niter=60000,
                burnin=1001,
                thinning=50,
                nsim = 4, parallel = TRUE, ncpus = 4)
summary(gomp.3)

#Saving the model
saveRDS(gomp.3, file = "scripts/model_outputs/BaSTA/adult lifespan/gompertz_bathtub.rda")
gomp.3 <-readRDS("scripts/model_outputs/BaSTA/adult lifespan/gompertz_bathtub.rda")

#Assessing model fit and plotting the model outputs
plot(gomp.3, plot.type = "gof")
plot(gomp.3, plot.type = "demorates")
plot(gomp.3, densities=TRUE)

#4.4. Weibull model----------------------------------------
bastaweibull <- basta(object = adcovdata, 
                      dataType = "census", 
                      model="WE",
                      niter=60000,
                      burnin=1001,
                      thinning=50,
                      nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastaweibull)
#saving the model
saveRDS(bastaweibull, file = "scripts/model_outputs/BaSTA/adult lifespan/lifespan_weibull.rda")
bastaweibull<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/lifespan_weibull.rda")

#Assessing model fit and plotting the model outputs
plot(bastaweibull, plot.type = "gof")
plot(bastaweibull, plot.type = "demorates")
plot(bastaweibull, densities=TRUE)

#4.5. Weibull model with a Makeham term --------------------------
bastaweibull.2 <- basta(object = adcovdata, 
                        dataType = "census",
                        model="WE",
                        shape="Makeham",
                        niter=60000,
                        burnin=1001,
                        thinning=50,
                        nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastaweibull.2)
saveRDS(bastaweibull.2, file = "scripts/model_outputs/BaSTA/adult lifespan/weibull_makeham.rda")
bastaweibull.2<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/weibull_makeham.rda")

#Assessing model fit and plotting the model outputs
plot(bastaweibull.2, plot.type = "gof")
plot(bastaweibull.2, plot.type = "demorates")
plot(bastaweibull.2, densities=TRUE)

#4.6. Weibull model with a bathtub term ----------------------
bastaweibull.3 <- basta(object = adcovdata, 
                        dataType = "census",
                        model="WE",
                        shape="bathtub",
                        niter=60000,
                        burnin=1001,
                        thinning=50,
                        nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastaweibull.3)
saveRDS(bastaweibull.3, file = "scripts/model_outputs/BaSTA/adult lifespan/weibull_bathtub.rda")
bastaweibull.3<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/weibull_bathtub.rda")

#Assessing model fit and plotting the model outputs
plot(bastaweibull.3, plot.type = "gof")
plot(bastaweibull.3, plot.type = "demorates")
plot(bastaweibull.3, densities=TRUE)


#4.7.exponential model-------------------------
bastaexp <- basta(object = adcovdata, 
                  dataType = "census",
                  model="EX",
                  niter=60000,#trying to double the burnin and the number of iterations
                  burnin=1001,
                  thinning=50,
                  nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastaexp)
saveRDS(bastaexp, file = "scripts/model_outputs/BaSTA/adult lifespan/lifespan_exp.rda")
bastaexp<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/lifespan_exp.rda")

#Assessing model fit and plotting the model outputs
plot(bastaexp, plot.type = "gof")
plot(bastaexp, plot.type = "demorates")
plot(bastaexp, densities=TRUE)


#4.8.logistic model-----------------------------
bastalog <- basta(object = adcovdata, 
                  dataType = "census",
                  model="LO",
                  niter=60000,
                  burnin=1001,
                  thinning=50,
                  nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastalog)
saveRDS(bastalog, file = "scripts/model_outputs/BaSTA/adult lifespan/lifespan_log.rda")
bastalog<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/lifespan_log.rda")

#Assessing model fit and plotting the model outputs
plot(bastalog, plot.type = "gof")
plot(bastalog, plot.type = "demorates")
plot(bastalog, densities=TRUE)


#4.9. Logistic model with a makeham shape------------------------------
bastalog.2 <- basta(object = adcovdata, 
                    dataType = "census", 
                    model="LO",
                    shape="Makeham",
                    niter=60000,
                    burnin=1001,
                    thinning=50,
                    nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastalog.2)
saveRDS(bastalog.2, file = "scripts/model_outputs/BaSTA/adult lifespan/lifespan_log_makeham.rda")
bastalog.2<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/lifespan_log_makeham.rda")

#Model fit and convergence plots
plot(bastalog.2, plot.type = "gof")
plot(bastalog.2, plot.type = "demorates")
plot(bastalog.2, densities=TRUE)


#4.10. Logistic model with a bathtub shape--------------------------------
bastalog.3 <- basta(object = adcovdata, 
                    dataType = "census",
                    model="LO",
                    shape="bathtub",
                    niter=60000,
                    burnin=1001,
                    thinning=50,
                    nsim = 4, parallel = TRUE, ncpus = 4)
summary(bastalog.3)
saveRDS(bastalog.3, file = "scripts/model_outputs/BaSTA/adult lifespan/lifespan_log_bathtub.rda")
bastalog.3<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/lifespan_log_bathtub.rda")

#Model fit and convergence plots
plot(bastalog.3, plot.type = "gof")
plot(bastalog.3, plot.type = "demorates")
plot(bastalog.3, densities=TRUE)

#The Weibull Makeham model has the best fit, selecting this moving forward

# ===========================================================================
# 5. PRIOR SPECIFICATION FOR THE NULL MODEL----------------------------------
# ===========================================================================
#Reminder that the null model is fitted with no covariates

mean(adcovdata$adult_lifespan) #9.534763 weeks
sd(adcovdata$adult_lifespan) # 3.159184 weeks

#Converting mean years into the characteristic lifespan (b1)
#BaSTA estimates b1 in years, and inverts this estimate
#ensuring we dop the same when setting the prior on this parameter
mu_years<-9.534763/52.1429
prior_scale_mu<-gamma(1+1/1.5)/mu_years
#setting the mean for b1 to 4.936857
hist(rnorm(1000, 4.94, 1))


#----------------------SETTING MODEL PRIORS-------------------

#Weakly informative priors-----------------
weakMean2 <- matrix(c(
  0, 1.5, prior_scale_mu
), nrow = 1, byrow = TRUE)
#mean for makeham (c), shape (b0), and scale (b1) parameter, respectively

weakSd2 <- matrix(c(
  1.0, 1.0, 1.0
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
  1.0, 0.5, 0.5
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
Nullbastaweibull <- basta(object = adcovdata, 
                          dataType = "census",
                          model="WE",
                          shape="Makeham",
                          thetaPriorMean = weakMean2, #adding the weak priors
                          thetaPriorSd = weakSd2,
                          thetaPriorLower = weakLower2,
                          niter=80000,
                          burnin=1001,
                          thinning=50,
                          nsim = 4, parallel = TRUE, ncpus = 4)
summary(Nullbastaweibull)
saveRDS(Nullbastaweibull, file = "scripts/model_outputs/BaSTA/adult lifespan/Null_model_weibull.rda")
Nullbastaweibull<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/Null_model_weibull.rda")

#assessing how the priors affect model fit and convergence
plot(Nullbastaweibull)
plot(Nullbastaweibull, plot.type = "gof")

#6.2. Null model with moderate priors---------------
Nullbastaweibull.moderate <- basta(object = adcovdata, 
                                   dataType = "census",
                                   model="WE",
                                   shape="Makeham",
                                   thetaPriorMean = moderateMean2,
                                   thetaPriorSd = moderateSd2,
                                   thetaPriorLower = moderateLower2,
                                   niter=80000,#trying to double the burnin and the number of iterations
                                   burnin=1001,
                                   thinning=50,
                                   nsim = 4, parallel = TRUE, ncpus = 4)
summary(Nullbastaweibull.moderate)
saveRDS(Nullbastaweibull.moderate, file = "scripts/model_outputs/BaSTA/adult lifespan/Nullbastaweibull.moderate.rda")
Nullbastaweibull.moderate<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/Nullbastaweibull.moderate.rda")

#assessing how the priors affect model fit and convergence
plot(Nullbastaweibull.moderate)
plot(Nullbastaweibull.moderate, plot.type = "gof")

# ===========================================================================
# 7. ADDING COVARIATES: ISOLATED EFFECT OF PARENTAL AGE----------------------
# ===========================================================================

#7.1. Plotting a Kaplein-Maier curve of offspring observed survival in each parental age category--------------

#Ensuring variables are correctly ordered
adcovdata$Timepoint_binned<-factor(adcovdata$Timepoint_binned, 
                                levels = c("early", "middle", "late"))
#censoring variable for survival package
adcovdata$event<-ifelse(adcovdata$Depart.Type =="D", 1,0)

#Renaming levels for the plot
adcovdata$Timepoint_binnedplot<-factor(adcovdata$Timepoint_binned, 
                                       levels = c("early", "middle", "late"),
                                       labels=c("Early-Aged", "Middle-Aged", "Late-Aged"))



#1.Kaplain-Maier curve for parental age
#survival estimates
km_fit_parentalage<- survfit(Surv(adult_lifespan,event) ~ Timepoint_binnedplot, data = adcovdata)#fitting the curve for Timepoint
#Kaplan-Meier curve
Km_adult_parentalage <- ggsurvplot(km_fit_parentalage, 
                                         data = adcovdata, 
                                         pval = FALSE,  
                                         conf.int = TRUE, 
                                         risk.table = "abs_pct",
                                         break.time.by = 3,
                                         surv.median.line = "hv",  
                                         xlab = "Offspring adult age (weeks)", 
                                         ylab = "Cumulative adult survival probability, S(x)",
                                         font.x = 35,
                                         font.y = 35,
                                         font.tickslab = c(35, "grey25"),
                                         font.legend = 25,
                                         risk.table.fontsize = 6,
                                         legend.title = "Parents' age at reproduction",
                                         legend.labs = c("Early-Aged", "Middle-Aged", "Late-Aged"),
                                         palette = c("#f96161", "#66b2b2", "#066594"))

# Combine the survival plot and the risk table
Km_adult_combined_plot <- cowplot::plot_grid(Km_adult_parentalage$plot, Km_adult_parentalage$table, ncol = 1, rel_heights = c(3, 1))

# Save the combined plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/001_KM_parentalage_adultlifespan.png",
       plot = Km_adult_combined_plot,
       bg = "transparent",
       device = "png",
       width = 460,
       height = 380,
       units = "mm")


#7.2. PARAMETRIC SURVIVAL MODELS FOR PARENTAL AGE EFFECTS------------------------------------------


#7.2.1. Weibull model------------------------------------------------------
adultbastaweibull <- basta(object = adcovdata, 
                           dataType = "census",
                           formulaMort=~Timepoint_binned-1, 
                           model="WE",
                           niter=60000,
                           burnin=1001,
                           thinning=50,
                           nsim = 4, parallel = TRUE, ncpus = 4)
summary(adultbastaweibull)
saveRDS(adultbastaweibull, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_weibull.rda")
adultbastaweibull<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_weibull.rda")

#Assessing model fit and plotting the model outputs
plot(adultbastaweibull, plot.type = "gof")
plot(adultbastaweibull, plot.type = "demorates")
plot(adultbastaweibull, densities=TRUE)


#7.2.2. Weibull model with a Makeham term --------------------------------
adultbastaweibull.2 <- basta(object = adcovdata, 
                             dataType = "census",
                             formulaMort=~Timepoint_binned-1, 
                             model="WE",
                             shape="Makeham",
                             niter=60000,
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(adultbastaweibull.2)
saveRDS(adultbastaweibull.2, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_weibull_makeham.rda")
adultbastaweibull.2<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_weibull_makeham.rda")

#Assessing model fit and plotting the model outputs
plot(adultbastaweibull.2, plot.type = "gof")
plot(adultbastaweibull.2, plot.type = "demorates")
plot(adultbastaweibull.2, densities=TRUE)


#7.2.3. Weibull model with a bathtub term --------------------------------------------------
adultbastaweibull.3 <- basta(object = adcovdata, 
                             dataType = "census",
                             formulaMort=~Timepoint_binned-1, 
                             model="WE",
                             shape="bathtub",
                             niter=60000,#trying to double the burnin and the number of iterations
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(adultbastaweibull.3)
saveRDS(adultbastaweibull.3, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_weibull_bathtub.rda")
adultbastaweibull.3<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_weibull_bathtub.rda")

#Assessing model fit and plotting the model outputs
plot(adultbastaweibull.3, plot.type = "gof")
plot(adultbastaweibull.3, plot.type = "demorates")
plot(adultbastaweibull.3, densities=TRUE)

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
  1.0, 1.0, 1.0
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
  0.5, 0.5, 0.5
), nrow = 3, ncol=3, byrow = TRUE)
#SD for the c, b0, and b1 terms, respectively

#moderate priors
moderateLower<-matrix(c(
  0, 0, 0
), nrow = 3, ncol =3, byrow = TRUE)
#lower bounds for the c, b0, and b1 terms, respectively


#-----------------8.1. Fitting priors to the parental age model-----------------------------------

#---------------Weak priors-----------------------------------------------------

#8.1.1. Weibull model with a Makeham term---------------------
adultbastaweibull.weak <- basta(object = adcovdata, 
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
summary(adultbastaweibull.weak)

saveRDS(adultbastaweibull.weak, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_weibull_weak.rda")
adultbastaweibull.weak<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_weibull_weak.rda")

#Assessing model fit and plotting the model outputs
plot(adultbastaweibull.weak, plot.type = "gof")
plot(adultbastaweibull.weak, plot.type = "demorates")
plot(adultbastaweibull.weak, densities=TRUE)
plot(adultbastaweibull.weak, type="fancy")


#-------------Moderate priors---------------------------------------------------

#8.1.2. Weibull model with a Makeham term ----------------------------
adultbastaweibull.moderate<- basta(object = adcovdata, 
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
summary(adultbastaweibull.moderate)
plot(adultbastaweibull.moderate)
saveRDS(adultbastaweibull.moderate, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_weibull_moderate.rda")
adultbastaweibull.moderate<-readRDS("scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_weibull_moderate.rda")

#Assessing model fit and plotting the model outputs
plot(adultbastaweibull.moderate, plot.type = "gof")
plot(adultbastaweibull.moderate, plot.type = "demorates")
plot(adultbastaweibull.moderate, densities=TRUE)


#Selected the weakly informative prior for all models moving forward

# ===========================================================================
# 9. ISOLATED EFFECT OF PARENTAL AGE PLOTS-----------------------------------
# ===========================================================================
#Plotting estimates from weibull model with weakly-informative priors

#9.1. MORALITY PARAMETER POSTERIOR DISTRIBUTIONS (I.E., C, B1, AND B0)----

#Extracting the model parameters
adult_theta_params <- adultbastaweibull.weak$params

#saving as a data frame
adult_theta_params<-as.data.frame(adult_theta_params)

# Reshaping the data
adult_posterior_df <- adult_theta_params %>%
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
adult_posterior_df$Parental_Age<- factor(adult_posterior_df$Parental_Age, 
                                         levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))

# Create separate datasets for "c", `b0` and `b1`
adb0_df <- adult_posterior_df %>% filter(Type == "b0")
adb1_df <- adult_posterior_df %>% filter(Type == "b1")
adcdf<-adult_posterior_df %>% filter(Type == "c")

#9.2. B0 POSTERIOR---------------------------------------
ad_b0_plot <- ggplot(adb0_df, aes(x = Value, fill = Parental_Age, color = Parental_Age, alpha=Parental_Age)) +
  geom_density(linewidth = 2) +
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
  scale_y_continuous(limits = c(0, 2.2), breaks = c(0, 0.7, 1.4, 2.1))+
  scale_x_continuous(limits = c(2.9, 5.1), breaks = c(3, 3.6, 4.2, 4.8))+
  labs(
    x = "b0 parameter value",
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

#save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/001_b0_parentalage_adultlifespan.png",
       plot = ad_b0_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")


#9.3. B1 POSTERIOR plot----------------------------
ad_b1_plot <- ggplot(adb1_df, aes(x = Value, 
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
  
  scale_y_continuous(limits = c(0, 5.2), breaks = c(0, 1.7, 3.4, 5.1 ))+
  scale_x_continuous(limits = c(4.09, 5.3), breaks = c(4.2, 4.53, 4.86, 5.20))+
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


ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/001_b1_parentalage_adultlifespan.png",
       plot = ad_b1_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")


#9.4. C PARAMETER-------------------------
ad_c_plot <- ggplot(adcdf, aes(x = Value, fill = Parental_Age, 
                               color = Parental_Age, 
                               alpha=Parental_Age)) +
  geom_density(linewidth = 2) +
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
  scale_x_continuous(limits = c(0, 1.22), breaks = c(0.0, 0.4, 0.8, 1.2))+
  labs(
    x = "c parameter value",
    y = "Density",
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

#save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/001_c_parentalage_adultlifespan.png",
       plot = ad_c_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")



#9.5--------CREATING THE SURVIVAL AND MORTALITY CURVES-----------------------------------

#Function for extracting mortality and survival estimates from the BaSTA model
plot_data_adult <- list()
for (demv in c("mort", "surv")) {
  for (icat in seq_along(adultbastaweibull.weak[[demv]])) {
    cuts <-adultbastaweibull.weak$cuts[[icat]]
    minAge <- as.numeric(adultbastaweibull.weak$modelSpecs["min. age"])
    xx <- adultbastaweibull.weak$x[cuts] + minAge
    yy <- adultbastaweibull.weak[[demv]][[icat]][, cuts]
    
    # Convert to long format
    df <- data.frame(
      Age = xx,
      Rate = yy[1, ], # Median
      LowerCI = yy[2, ], # Lower Confidence Bound
      UpperCI = yy[3, ], # Upper Confidence Bound
      Category = names(adultbastaweibull.weak[[demv]])[icat],
      Type = demv # Mortality or Survival
    )
    plot_data_adult[[length(plot_data_adult) + 1]] <- df
  }
}

#Combine all into one data frame
plot_data_adult<- do.call(rbind, plot_data_adult)

#renaming the values of Timepoint
plot_data_adult$Category<- factor(plot_data_adult$Category, 
                                  levels = c("Timepoint_binnedearly", "Timepoint_binnedmiddle", "Timepoint_binnedlate"), 
                                  labels = c("Early-Aged", "Middle-Aged", "Late-Aged"))


# Set the scaling factors (e.g., minimum and maximum age in days)
max_days <- 365.25 
min_days<-0

# Adjust Age column back to days
plot_data_adult<- plot_data_adult %>%
  mutate(Age_days = Age * (max_days - min_days) + min_days)

#ensuring the age categories are added as a factor
plot_data_adult$Category<-as.factor(plot_data_adult$Category)



#9.5.1. Splitting the data into the Mortality and Survival datasets------

#mortality dataset
mort_df_adult <- plot_data_adult %>% filter(Type == "mort") %>% 
  mutate(Rate_weeks = Rate / 52.1429,
         weeks_LowerCI = LowerCI/52.1429,
         weeks_UpperCI = UpperCI/52.1429) #Converting the hazard rate to be per unit week

#survival dataset
surv_df_adult<- plot_data_adult%>% filter(Type == "surv")

#median Survival probability
  surv_df_median<-surv_df_adult %>%
    group_by(Category) %>%  
    arrange(Age_days) %>%
    mutate(CI = Rate) %>%  
    summarise(median_line = approx(Rate, Age_days, xout = 0.50)$y)


#9.6. PLOTTING THE SURVIVAL AND MORTALITY CURVES-----------------------------------
  
#9.6.1. Cumulative survival probability plot---------------------------------------
adult_survival_plot<- ggplot(data = surv_df_adult,
                             aes(x = Age_days/7, 
                                 y = Rate,
                                 colour = Category))+
  geom_line(data = surv_df_adult, 
            aes(x = Age_days/7, 
                y = Rate,
                colour = Category),
           linewidth=4, alpha=1)+
  geom_ribbon(data = surv_df_adult, 
              aes(y=NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI, 
                  color= Category,
                  fill=Category, 
                  alpha = Category), 
              linewidth = 1,
              linetype="dashed") +
  geom_segment(data = surv_df_median,
               aes(x = 0, xend = median_line/7, y = 0.5, yend = 0.5, color = Category),
               linetype = "dashed", linewidth = 2, alpha = 0.9) +
  geom_segment(data = surv_df_median,
               aes(x = median_line/7, xend = median_line/7, y = 0, yend = 0.5, color = Category),
               linetype = "dashed", linewidth= 2, alpha = 0.9)+
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
  scale_color_manual(name = "Parents' adult age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' adult age at reproduction", values = c(0.6,0.1, 0.6))+
  scale_fill_manual(name = "Parents' adult age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15))+
  theme(
    legend.position = "top", 
    legend.box = "horizontal"
  )+
  labs(x = "Offspring adult age (weeks)", y = "Cumulative adult survival probability, S(x)")


#save plot
  ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/001_survival_parentalage_adultlifespan.png",
         plot = adult_survival_plot, 
         bg="transparent",
         device = "png", 
         width = 440, 
         height = 240, 
         units = "mm")

  
  
  #9.6.2. mortality risk plot---------------------------------------------------------------
adult_mort_plot<- ggplot(data = mort_df_adult,
                         aes(x = Age_days/7, 
                             y = Rate_weeks,
                             colour = Category))+
  geom_line(data = mort_df_adult, 
            aes(x = Age_days/7, 
                y = Rate_weeks,
                color = Category),
            linewidth=4, 
            alpha=1)+
  geom_ribbon(data = mort_df_adult, 
              aes(y=NULL, 
                  ymin = weeks_LowerCI, 
                  ymax = weeks_UpperCI, 
                  color= Category,
                  fill=Category, 
                  alpha = Category),
              linewidth = 1,
              linetype="dashed") +
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
  ) + scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15))+
  scale_y_continuous(breaks=c(0, 0.33, 0.66, 0.99, 1.3))+
  labs(x = "Offspring adult age (weeks)", y = "Instantaneous hazard rate, μ(x)")

ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/001_rate_parentalage_adultlifespan.png",
       plot = adult_mort_plot, 
       bg="transparent",
       device = "png", 
       width = 440, 
       height = 240, 
       units = "mm")

#---------------Creating combination plot for manuscript-----------------------------
new_adult_mort_plot<- adult_mort_plot+
  labs( y = "Hazard rate, μ(x)")+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "none")

new_adult_survival_plot<-adult_survival_plot+
  labs( y = "Survival, S(x)")+
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top")


#timescale
new_ad_b1_plot<-ad_b1_plot+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.title=element_text(size=50),
        axis.text=element_text(size=50),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))


#shape
new_ad_b0_plot<-ad_b0_plot+
  theme(axis.title.x = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))

#makeham
new_ad_c_plot<-ad_c_plot+
  labs( x = "Parameter value")+
  theme(axis.title.y = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))


#Combining all plots together as in "fancy basta" layout------------------------

adult_lifespan_inference<-ggarrange(
  ggarrange(new_ad_b1_plot,new_ad_b0_plot, new_ad_c_plot, nrow = 3, labels = c("b1", "b0", "c"),
            font.label = list(size = 60, face = "bold"),
            label.x = c(0.03, 0.03, 0.05),
            heights = c(0.94, 0.94,  1),align = "v"),
  ggarrange(new_adult_survival_plot, new_adult_mort_plot, nrow = 2, labels = c("",""),
            font.label = list(size = 50, face = "bold"),
            heights = c(0.94,  1),align = "v"),
  ncol = 2,
  nrow=1,
  align = "h",
  labels = c("", ""),
  widths = c(0.9, 1.5))


ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/001_inference_plots.png",
       plot = adult_lifespan_inference, 
       device = "png", 
       width = 900, 
       height = 630, 
       units = "mm")

# ===========================================================================
# 10. MODEL FIT PLOTS--------------------------------------------------------
# ===========================================================================

#Extracting the data from the Kaplan Meier curve
Km_adult_parentalage$data.survplot<-
  Km_adult_parentalage$data.survplot %>%
  rename(Category = Timepoint_binnedplot)

# Create the zero-point for each group
zero_point_001 <- Km_adult_parentalage$data.survplot %>%
  group_by(Category) %>%
  summarise(time = 0, surv = 1, lower = 1, upper = 1)

# Combine zero-point and the actual KM data
Km_adult_parentalage$data.survplot<- Km_adult_parentalage$data.survplot %>%
  bind_rows(zero_point_001) %>%
  arrange(Category, time)  


#Combining the parametric and non-parametric estimates into one plot
Km_adult_combined_plot_001 <- ggplot() +
  #The Kaplan-Meier curve
  geom_step(data = Km_adult_parentalage$data.survplot, 
            aes(x = time, 
                y = surv), 
            colour = "#066594",
            linewidth = 2.5, alpha = 1)+
  geom_step(data = Km_adult_parentalage$data.survplot, 
            aes(x = time, 
                y = lower),
            colour = "#066594", 
            linetype="dashed",
            alpha = 0.5,
            linewidth = 1.2)+ 
  geom_step(data = Km_adult_parentalage$data.survplot, 
            aes(x = time, 
                y = upper),
            colour = "#066594", 
            linetype="dashed",
            alpha = 0.5,
            linewidth=1.2)+
  #Parametric survival curve
  geom_line(data = surv_df_adult, 
            aes(x = Age_days / 7, 
                y = Rate),
            colour = "#f96161", 
            linewidth = 2.5, alpha = 1) +
  geom_ribbon(data = surv_df_adult, 
              aes(x = Age_days / 7, 
                  y = NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI),
              fill = "#f96161", 
              alpha = 0.3)+
  facet_wrap(~ Category, ncol = 1) +  
  theme_classic() +
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        strip.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), )+
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20))+
  labs(x = "Offspring adult age (weeks)", y = "Cumulative adult survival probability, S(x)")

ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/001_modelfit_parentalage_adult.png",
       plot = Km_adult_combined_plot_001, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 380, 
       units = "mm")



# ===========================================================================
# 11. ADDING COVARIATES: ISOLATED EFFECT OF PARENTAL TEMPERATURE TREATMENT---
# ===========================================================================

#Plotting a Kaplein-Maier curve for parental temperature
adcovdata$Temp<-factor(adcovdata$Temp, 
                                   levels = c("25.5", "28", "30.5"))

#1.Kaplain-Maier curve for parental age
km_fit_Temp <- survfit(Surv(adult_lifespan,event) ~ Temp, data =adcovdata)#fitting the curve for Timepoint
Km_ad_Temp <- ggsurvplot(km_fit_Temp, 
                          data = adcovdata,
                          pval = FALSE,  
                          conf.int = TRUE,  
                          risk.table = "abs_pct",
                          break.time.by = 3,
                          surv.median.line = "hv",  
                          xlab = "Offspring age (weeks)", 
                          ylab = "Cumulative adult survival probability, S(x)",
                          font.x = 35,
                          font.y = 35,
                          font.tickslab = c(35, "grey25"),
                          font.legend = 25,
                          risk.table.fontsize = 6,
                          legend.title = "Parents' Temperature Treatment",
                          legend.labs = c("25.5°C", "28.0°C", "30.5°C"),
                          palette = c("#2f4b7c","orange2", "#a64b61"))

# Combine the survival plot and the risk table
Km_ad_temp_combined_plot <- cowplot::plot_grid(Km_ad_Temp$plot, Km_ad_Temp$table, ncol = 1, rel_heights = c(3, 1))

# Save the combined plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/002_KM_temperature_adulthood.png",
       plot = Km_ad_temp_combined_plot,
       bg = "transparent",
       device = "png",
       width = 420,
       height = 360,
       units = "mm")

#11.2. PARAMETRIC SURVIVAL MODELS FOR PARENTAL TEMPERATURE EFFECTS------------------------------------------

#-----------------------Using weibull distributions only-------------------------

#11.2.1. Weibull model with Makeham term------------------
adultbastaweibulltemp.2 <- basta(object = adcovdata, 
                                 dataType = "census",
                                 formulaMort=~Temp-1, 
                                 model="WE",
                                 shape="Makeham",
                                 niter=60000,
                                 burnin=1001,
                                 thinning=50,
                                 nsim = 4, parallel = TRUE, ncpus = 4)
summary(adultbastaweibulltemp.2)
saveRDS(adultbastaweibulltemp.2, file = "scripts/model_outputs/BaSTA/adult lifespan/temp_adult_lifespan_weibull_Makeham.rda")
adultbastaweibulltemp.2<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/temp_adult_lifespan_weibull_Makeham.rda")#this is with the makeham shape to be confusing 


#Assessing model fit and plotting the model outputs
plot(adultbastaweibulltemp.2, plot.type = "gof")
plot(adultbastaweibulltemp.2, plot.type = "demorates")
plot(adultbastaweibulltemp.2, densities=TRUE)



#11.2.2. weibull model----------------------------------------------
adultbastaweibulltemp <- basta(object = adcovdata, 
                               dataType = "census",
                               formulaMort=~Temp-1, 
                               model="WE",
                               niter=60000,
                               burnin=1001,
                               thinning=50,
                               nsim = 4, parallel = TRUE, ncpus = 4)
summary(adultbastaweibulltemp)
saveRDS(adultbastaweibulltemp, file = "scripts/model_outputs/BaSTA/adult lifespan/temp_adult_lifespan_weibull.rda")
adultbastaweibulltemp<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/temp_adult_lifespan_weibull.rda")


#Assessing model fit and plotting the model outputs
plot(adultbastaweibulltemp, plot.type = "gof")
plot(adultbastaweibulltemp, plot.type = "demorates")
plot(adultbastaweibulltemp, densities=TRUE)


#11.2.3 weibull model with a bathtub shape--------------------------
adultbastaweibulltemp.3 <- basta(object = adcovdata, 
                                 dataType = "census",
                                 formulaMort=~Temp-1, 
                                 model="WE",
                                 shape="bathtub",
                                 niter=60000,
                                 burnin=1001,
                                 thinning=50,
                                 nsim = 4, parallel = TRUE, ncpus = 4)
summary(adultbastaweibulltemp.3)
saveRDS(adultbastaweibulltemp.3, file = "scripts/model_outputs/BaSTA/adult lifespan/temp_adult_lifespan_weibull.3.rda")
adultbastaweibulltemp.3<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/temp_adult_lifespan_weibull.3.rda")

#Selecting the Weibull Makeham model moving forwards

# ===========================================================================
# 12. PRIOR SENSITIVITY ANALYSIS FOR THE PARENTAL TEMPERATURE MODEL-----------
# ===========================================================================
#Using the same priors as specified for the parental age analysis

#12.1. Weibull Makeham model with weak priors---------------------------------
adultbastaweibulltemp.weak <- basta(object = adcovdata, 
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
summary(adultbastaweibulltemp.weak)
saveRDS(adultbastaweibulltemp.weak, file = "scripts/model_outputs/BaSTA/adult lifespan/temp_adult_lifespan_weibull_Makeham_weakprior.rda")
adultbastaweibulltemp.weak<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/temp_adult_lifespan_weibull_Makeham_weakprior.rda")#this is with the makeham shape to be confusing 


#Assessing model fit and plotting the model outputs
plot(adultbastaweibulltemp.weak, plot.type = "gof")
plot(adultbastaweibulltemp.weak, plot.type = "demorates")
plot(adultbastaweibulltemp.weak, densities=TRUE)

#12.2. With moderate priors---------------------------------------
adultbastaweibulltemp.moderate <- basta(object = adcovdata, 
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
summary(adultbastaweibulltemp.moderate)
saveRDS(adultbastaweibulltemp.moderate, file = "scripts/model_outputs/BaSTA/adult lifespan/temp_adult_lifespan_weibull_Makeham_moderate.rda")
adultbastaweibulltemp.moderate<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/temp_adult_lifespan_weibull_Makeham_moderate.rda") 

#Assessing model fit and plotting the model outputs
plot(adultbastaweibulltemp.moderate, plot.type = "gof")
plot(adultbastaweibulltemp.moderate, plot.type = "demorates")
plot(adultbastaweibulltemp.moderate, densities=TRUE)


# ===========================================================================
# 13. ISOLATED EFFECT OF PARENTAL TEMPERATURE PLOTS-----------------------------------
# ===========================================================================
#Plotting estimates from weibull model with weakly-informative priors

#13.1. Plotting the posterior distributions for the mortality parameters (c, b0, and b1)

#extracting parameters
tempadult_theta_params <- adultbastaweibulltemp.weak$params

#converting to a data frame
tempadult_theta_params<-as.data.frame(tempadult_theta_params)

#Reshaping the data
tempadult_posterior_df <- tempadult_theta_params %>%
  pivot_longer(
    cols = everything(),
    names_to = "Combined", 
    values_to = "Value"
  ) %>%
  mutate(
    Type = str_extract(Combined, "^[^.]+"),  # Extract everything before the first period
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
tempadult_posterior_df$Temp<- factor(tempadult_posterior_df$Temp, 
                                     levels = c("25.5°C", "28.0°C", "30.5°C"))


# 13.2. Creating separate datasets for `b0` and `b1`, and c mortality parameters:
temadb0_df <- tempadult_posterior_df %>% filter(Type == "b0")
temadb1_df <- tempadult_posterior_df %>% filter(Type == "b1")
temadcdf<-tempadult_posterior_df %>% filter(Type == "c")

#13.3 Plotting the b0 parameter----------------
tempad_b0_plot <- ggplot(temadb0_df, aes(x = Value, fill = Temp, color = Temp, alpha=Temp)) +
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
  scale_x_continuous(limits = c(2.8, 4.7), breaks = c(3.0, 3.5, 4.0, 4.5))+
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
  guides(alpha=guide_legend(order = 1, title.position = "top",title.hjust = 0.5,
                            ,override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) 

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/002_b0_temperature_adultlifespan.png",
       plot = tempad_b0_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")

#13.2. Creating the posterior plot for the b1 posterior parameter values------------------
tempad_b1_plot <- ggplot(temadb1_df, aes(x = Value, fill = Temp, color = Temp, alpha=Temp)) +
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
  scale_x_continuous(limits = c(4.3, 5.1), breaks = c(4.4, 4.6, 4.8, 5.0))+
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
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), 
  )+
  guides(alpha=guide_legend(order = 1, title.position = "top",title.hjust = 0.5,
                            ,override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) 

ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/002_b1_temperature_adultlifespan.png",
       plot = tempad_b1_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")

#13.3. Creating the posterior plot for the c parameter posterior values---------------
tempad_c_plot <- ggplot(temadcdf, aes(x = Value, fill = Temp, color = Temp, alpha=Temp)) +
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
  scale_x_continuous(limits = c(0.0, 1.0), breaks = c(0.00, 0.33, 0.66, 1.00))+
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
  guides(alpha=guide_legend(order = 1, title.position = "top",title.hjust = 0.5,
                            ,override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) 


ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/002_c_temperature_adultlifespan.png",
       plot = tempad_c_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
       units = "mm")



#---------------------SURVIVAL AND MORTALITY CURVES-----------------------------
#13.4. Creating the survival and Instantaneous hazard rate curves------

#13.4.1. function to extract necessary values from the BaSTA model-----
temp_plot_data_adult <- list()
for (demv in c("mort", "surv")) {
  for (icat in seq_along(adultbastaweibulltemp.weak[[demv]])) {
    cuts <-adultbastaweibulltemp.weak$cuts[[icat]]
    minAge <- as.numeric(adultbastaweibulltemp.weak$modelSpecs["min. age"])
    xx <- adultbastaweibulltemp.weak$x[cuts] + minAge
    yy <- adultbastaweibulltemp.weak[[demv]][[icat]][, cuts]
    # Convert to long format
    df <- data.frame(
      Age = xx,
      Rate = yy[1, ], #median rate
      LowerCI = yy[2, ],#lower CI
      UpperCI = yy[3, ], #upper CI
      Category = names(adultbastaweibulltemp.weak[[demv]])[icat],
      Type = demv # Mortality or Survival
    )
    temp_plot_data_adult[[length(temp_plot_data_adult) + 1]] <- df
  }
}

# Combine all into one data frame
temp_plot_data_adult<- do.call(rbind, temp_plot_data_adult)

#renaming the values of Temperature
temp_plot_data_adult$Category<- factor(temp_plot_data_adult$Category, 
                                  levels = c("Temp25.5", "Temp28", "Temp30.5"), 
                                  labels = c("25.5°C", "28.0°C", "30.5°C"))


# Set the scaling factors (e.g., minimum and maximum age in days)
max_days <- 365.25 
min_days<-0

# Adjust Age column back to days
temp_plot_data_adult<- temp_plot_data_adult %>%
  mutate(Age_days = Age * (max_days - min_days) + min_days)

temp_plot_data_adult$Category<-as.factor(temp_plot_data_adult$Category)


#13.4.2. Splitting the data into Mortality and Survival datasets (for plotting)

#mortality data frame------------
temp_mort_df_adult <- temp_plot_data_adult %>% filter(Type == "mort") %>% 
  mutate(Rate_weeks = Rate / 52.1429,
         weeks_LowerCI = LowerCI/52.1429,
         weeks_UpperCI = UpperCI/52.1429)

#survival data frame-------------
temp_surv_df_adult<- temp_plot_data_adult%>% filter(Type == "surv")

#median survival probability-----
temp_surv_df_median<-temp_surv_df_adult %>%
  group_by(Category) %>%  
  arrange(Age_days) %>%  
  summarise(median_line = approx(Rate, Age_days, xout = 0.5)$y)


#13.4.3. cumulative survival plot - effect of parental temperature-----------
tempadult_survival_plot<- ggplot(data = temp_surv_df_adult,
                                 aes(x = Age_days/7,
                                     y = Rate,
                                     colour = Category))+
  geom_line(data = temp_surv_df_adult, 
            aes(x = Age_days/7, 
                y = Rate,
                colour = Category),
            linewidth=4, alpha=1)+
  geom_ribbon(data = temp_surv_df_adult, 
              aes(y=NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI, 
                  color= Category,
                  fill=Category,
                  alpha=Category),
              linewidth =1,
              linetype="dashed") +
  geom_segment(data = temp_surv_df_median,
               aes(x = 0, xend = median_line/7, y = 0.5, yend = 0.5, color = Category),
               linetype = "dashed", linewidth = 2, alpha = 0.9) +
  geom_segment(data = temp_surv_df_median,
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
  scale_x_continuous(limits = c(0, 18), breaks = c(0, 3, 6, 9, 12, 15, 18))+
  scale_color_manual(name = "Parents' Temperature Treatment", values = c("#2f4b7c","orange2", "#a64b61"))+
  scale_fill_manual(name = "Parents' Temperature Treatment", values = c("#2f4b7c","orange2", "#a64b61"))+
  scale_alpha_manual(name = "Parents' Temperature Treatment", values = c(0.7,0.1,0.4))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  theme(
    legend.position = "top", 
    legend.box = "horizontal"
  )+labs(x = "Offspring adult age (weeks)", y = "Cumulative adult survival probability, S(x)")


# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/002_survival_temperature_adultlifespan.png",
       plot = tempadult_survival_plot, 
       bg="transparent",
       device = "png", 
       width = 440, 
       height = 240, 
       units = "mm")


#13.4.4. Instantaneous mortality risk plot--------------------
tempadult_mort_plot<- ggplot(data = temp_mort_df_adult,
                             aes(x = Age_days/7, 
                                 y = Rate_weeks,
                                 colour = Category))+
  geom_line(data = temp_mort_df_adult, 
            aes(x = Age_days/7, 
                y = Rate_weeks,
                color = Category),
            linewidth=4, alpha=1)+
  geom_ribbon(data = temp_mort_df_adult, 
              aes(y=NULL, 
                  ymin = weeks_LowerCI, 
                  ymax = weeks_UpperCI, 
                  color= Category,
                  fill=Category,
                  alpha=Category),
              linewidth = 1,
              linetype="dashed") +
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
  scale_y_continuous(limits = c(0, 1.32), breaks = c(0, 0.33, 0.66, 0.99, 1.30))+
  scale_x_continuous(limits = c(0, 18), breaks = c(0, 3, 6, 9, 12, 15, 18))+
  scale_color_manual(name = "Parents' Temperature Treatment", values = c("#2f4b7c","orange2", "#a64b61"))+
  scale_fill_manual(name = "Parents' Temperature Treatment", values = c("#2f4b7c","orange2", "#a64b61"))+
  scale_alpha_manual(name = "Parents' Temperature Treatment", values = c(0.7,0.1,0.4))+
  guides(alpha = guide_legend(order = 1,title.position = "top",title.hjust = 0.5,
                              override.aes = list(size = 10)),
         fill = guide_legend(order = 1, title.position = "top",title.hjust = 0.5),
         colour = guide_legend(order = 1,title.position = "top",title.hjust = 0.5)) +
  theme(
    legend.position = "top", 
    legend.box = "horizontal"
  )+labs(x = "Offspring adult age (weeks)", y = "Instantaneous hazard rate, μ(x)")

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/002_rate_temperature_adultlifespan.png.png",
       plot = tempadult_mort_plot, 
       bg="transparent",
       device = "png", 
       width = 440, 
       height = 240, 
       units = "mm")


#-----------------------COMBINED PLOT FOR PAPER---------------------------------
new_temp_mort_plot<- tempadult_mort_plot+
  labs( y = "Hazard rate, μ(x)")+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "none")

new_temp_survival_plot<-tempadult_survival_plot+
  labs( y = "Survival, S(x)")+
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50))


#timescale
new_temp_b1_plot<-tempad_b1_plot+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.title=element_text(size=50),
        axis.text=element_text(size=50),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))

#shape
new_temp_b0_plot<-tempad_b0_plot+
  theme(axis.title.x = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))

#makeham
new_temp_c_plot<-tempad_c_plot+
  labs( x = "Parameter value")+
  theme(axis.title.y = element_blank(),
        axis.title=element_text(size=50),
        axis.text=element_text(size=50),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank())+
  theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 20, 10))

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


ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/002_inference_temperature.png",
       plot = lifespan_temp_inference, 
       device = "png", 
       width = 900, 
       height = 630, 
       units = "mm")


# ===========================================================================
# 14. MODEL FIT PLOTS FOR PARENTAL TEMPERATURE-------------------------------
# ===========================================================================

#Extracting the data from the Kaplan Meier curve
Km_ad_Temp$data.survplot<-
  Km_ad_Temp$data.survplot %>%
  rename(Category = Temp) %>% 
  mutate(Category = case_when(  
    Category == "25.5" ~ "25.5°C",
    Category == "28"   ~ "28.0°C",
    Category == "30.5" ~ "30.5°C",
    TRUE ~ Category 
  ))


# Create the zero-point for each group
zero_point_002 <- Km_ad_Temp$data.survplot %>%
  group_by(Category) %>%
  summarise(time = 0, surv = 1, lower = 1, upper = 1)

# Combine zero-point and the actual KM data
Km_ad_Temp$data.survplot<- Km_ad_Temp$data.survplot %>%
  bind_rows(zero_point_002) %>%
  arrange(Category, time)  


#14.1. Plot Combining the parametric and non-parametric estimates into one---------
Km_ad_combined_plot_002 <- ggplot() +
  #The Kaplan-Meier curve
  geom_step(data = Km_ad_Temp$data.survplot, 
            aes(x = time, 
                y = surv), 
            colour = "#066594",
            linewidth = 2.5, alpha = 1)+
  geom_step(data = Km_ad_Temp$data.survplot, 
            aes(x = time, 
                y = lower),
            colour = "#066594", 
            linetype="dashed",
            alpha = 0.5,
            linewidth = 1.2)+ 
  geom_step(data = Km_ad_Temp$data.survplot, 
            aes(x = time, 
                y = upper),
            colour = "#066594", 
            linetype="dashed",
            alpha = 0.5,
            linewidth =1.2)+
  #Parametric survival curve
  geom_line(data = temp_surv_df_adult, 
            aes(x = Age_days / 7, 
                y = Rate),
            colour = "#f96161", 
            linewidth = 2.5, alpha = 1) +
  geom_ribbon(data = temp_surv_df_adult, 
              aes(x = Age_days / 7, 
                  y = NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI),
              fill = "#f96161", 
              alpha = 0.3)+
  facet_wrap(~ Category, ncol = 1) +  
  theme_classic() +
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "top",
        legend.title=element_text(size=50),
        legend.text=element_text(size=50),
        strip.text = element_text(size=50),
        strip.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),  
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA), )+
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20))+
  labs(x = "Offspring adult age (weeks)", y = "Cumulative adult survival probability, S(x)")

#saving plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/002_modelfit_temp_adultlife.png",
       plot = Km_ad_combined_plot_002, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 380, 
       units = "mm")

#14.2. Creating a combined model fit plot of parental temperature and parental age---------------------------
new_Km_survival_combined_plot_002<-Km_ad_combined_plot_002+
  theme(axis.title.y = element_blank())

tempagefit<-ggarrange(
  Km_adult_combined_plot_001,
  new_Km_survival_combined_plot_002,
  ncol = 2,
  nrow=1,
  align = "h",
  labels = c("", ""),
  widths = c(0.6,0.5))


ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/001_modelfits.png",
       plot = tempagefit,
       device = "png", 
       width = 900, 
       height = 630, 
       units = "mm")


# ===========================================================================
# 15. ADDING COVARIATTES: INTERACTIVE EFFECT OF PARENTAL AGE & TEMPERATURE----
# ===========================================================================

#15.1 Creating the Kaplan Meier curve for parental age*temperature
KMdataagetemp<- adcovdata %>%
  mutate(temp_ParentalAge = interaction(Timepoint_binnedplot, Temp)) #using interaction() to create a new interactive term for age and temperature
#adding celsius to the label
KMdataagetemp$Temp_label <- paste0(KMdataagetemp$Temp, "°C")

#15.2.Kaplain-Maier curve for parental age and temperature interaction
#survival model
km_fit_parentalagetemp <- survfit(Surv(adult_lifespan,event) ~ temp_ParentalAge, data = KMdataagetemp)
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
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/004_KM_parentalagetemp_adultlifespan.png",
       plot = Km_survival_parentalagetempplot,
       bg = "transparent",
       device = "png",
       width = 380,
       height = 340,
       units = "mm")

#15.2. PARAMETRIC SURVIVAL MODELS--------------------------------------

#15.2.1. Weibull model with a Makeham term--------------
adtempagebastaweibull.2<- basta(object = adcovdata, 
                              dataType = "census",
                              formulaMort=~Temp:Timepoint_binned-1, 
                              model="WE",
                              shape="Makeham",
                              niter=60000,
                              burnin=1001,
                              thinning=50,
                              nsim = 4, parallel = TRUE, ncpus = 4)
summary(adtempagebastaweibull.2$fullpar)
saveRDS(adtempagebastaweibull.2, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_age_temp_wemakeham.rda")
adtempagebastaweibull.2<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_age_temp_wemakeham.rda")

#Assessing model fit and plotting the model outputs
plot(adtempagebastaweibull.2, plot.type = "gof")
plot(adtempagebastaweibull.2, plot.type = "demorates")
plot(adtempagebastaweibull.2, densities=TRUE)


#15.2.2.Weibull model---------------------------
adtempagebastaweibull<- basta(object = adcovdata, 
                            dataType = "census",
                            formulaMort= ~Temp:Timepoint_binned-1, 
                            model="WE",
                            niter=60000,
                            burnin=1001,
                            thinning=50,
                            nsim = 4, parallel = TRUE, ncpus = 4)
summary(adtempagebastaweibull)
saveRDS(adtempagebastaweibull, file = "scripts/model_outputs/BaSTA/adult lifespan/ad_lifespan_tempage_weibull.rda")
adtempagebastaweibull<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/ad_lifespan_tempage_weibull.rda")

#Assessing model fit and plotting the model outputs
plot(adtempagebastaweibull, plot.type = "gof")
plot(adtempagebastaweibull, plot.type = "demorates")
plot(adtempagebastaweibull, densities=TRUE)

#15.2.3.Weibull model with a bathtub term-------------------
adtempagebastaweibull.3<- basta(object = adcovdata, 
                              dataType = "census",
                              formulaMort= ~Temp:Timepoint_binned-1, 
                              model="WE",
                              shape="bathtub",
                              niter=60000,
                              burnin=1001,
                              thinning=50,
                              nsim = 4, parallel = TRUE, ncpus = 4)
summary(adtempagebastaweibull.3)
saveRDS(adtempagebastaweibull.3, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_tempage_weibull_bathtub.rda")
adtempagebastaweibull.3<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_tempage_weibull_bathtub.rda")


#Selecting the Weibull Makeham model for the prior specification process


#------------FITTING THE MODEL WITH WEAK PRIORS---------------------------------

#Weakly informative priors (mean)
weakMean3 <- matrix(c(0, 1.5, prior_scale_mu), 
                    nrow = 9, ncol = 3, byrow = TRUE)
#Priors for the mean of the c, b0 and b1 terms

#weakly infromative priors (SD)
weakSd3 <- matrix(c(
  1.0, 1.0, 1.0
), nrow = 9, ncol = 3, byrow = TRUE)
#Priors for the sd of the c, b0 and b1 terms

weakLower3<-matrix(c(
  0, 0, 0
), nrow = 9, ncol = 3,byrow = TRUE)
#Priors for the lower bound of the c, b0 and b1 terms


#16.2. Fitting the priors to the model (just testing the weak priors here)------------
#Fitting weak priors to the weibull makeham model
adtempagebastaweibull.weak<- basta(object = adcovdata, 
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
summary(adtempagebastaweibull.weak)
saveRDS(adtempagebastaweibull.weak, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_age_tempweak.rda")
adtempagebastaweibull.weak<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_age_tempweak.rda")

#Assessing model fit and plotting the model outputs
plot(adtempagebastaweibull.weak, plot.type = "gof")
plot(adtempagebastaweibull.weak, plot.type = "demorates")
plot(adtempagebastaweibull.weak, densities=TRUE)


# =======================================================================================================
# 17. PLOTTINGT THE MORTALITY PARAMETERS AND DEMOGRAPHIC CURVES FOR THE PARENTAL AGE X TEMP INTERACTION-----------
# ========================================================================================================

#17.1.PLOTTING THE POSTERIOR DISTRIBUTIONS FOR B0, B1 AND C------------------------------

#Extracting the model parameters
adtempage_theta_params <- adtempagebastaweibull.weak$params

#saving as a data frame
adtempage_theta_params<-as.data.frame(adtempage_theta_params)

# Reshaping the data
adtempage_posterior_df <- adtempage_theta_params %>%
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
adtempage_posterior_df$Parental_Age<- factor(adtempage_posterior_df$Parental_Age, 
                                           levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))

#17.2. Create two separate datasets for `b0` and `b1`, and c------------
adagetempb0_df <- adtempage_posterior_df  %>% filter(Type == "b0")
adagetempb1_df <- adtempage_posterior_df  %>% filter(Type == "b1")
adagetempcdf<-adtempage_posterior_df  %>% filter(Type == "c")

#17.3. B0 PARAMETER PLOT--------------
adagetemp_b0_plot <- ggplot(adagetempb0_df, aes(x = Value, fill = Parental_Age, color = Parental_Age, alpha=Parental_Age)) +
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
  
  scale_x_continuous(limits = c(2.5, 5.23), breaks = c(2.5, 3.4, 4.3, 5.2))+
  
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
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/004_b0_agetemp_adultlife.png",
       plot = adagetemp_b0_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 460, 
       units = "mm")


#17.4.-B1 PARAMETER PLOT-------------
adagetemp_b1_plot <- ggplot(adagetempb1_df, aes(x = Value, fill = Parental_Age, color = Parental_Age, alpha=Parental_Age)) +
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
  scale_x_continuous(limits = c(3.79, 5.45), breaks = c(3.80, 4.34, 4.88, 5.42))+
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

#saving plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/004_b1_agetemp_adultlife.png",
       plot = adagetemp_b1_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 460, 
       units = "mm")


#17.5. C PARAMETER PLOT----------------------
adagetempc_plot<- ggplot(adagetempcdf, aes(x = Value, fill = Parental_Age, color = Parental_Age, alpha=Parental_Age)) +
  geom_density(size = 1) +
  theme_classic() +
  scale_fill_manual(values = c("Early-Aged" = "#f96161", 
                               "Middle-Aged" = "#66b2b2", 
                               "Late-Aged" = "#066594")) +
  scale_color_manual(values = c("Early-Aged" = "#f96161", 
                                "Middle-Aged" = "#66b2b2", 
                                "Late-Aged" = "#066594")) +
  scale_x_continuous(limits = c(0.0, 1.5), breaks = c(0.0, 0.5, 1.0, 1.5))+
  
  scale_alpha_manual(values = c("Early-Aged" = 0.6, 
                                "Middle-Aged" = 0.1, 
                                "Late-Aged" = 0.6))+
  facet_wrap(~Temperature, ncol=1)+
  labs(
    x = "c parameter value",
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


ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/004_c_agetemp_adultlife.png",
       plot = adagetempc_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 460, 
       units = "mm")



#17.5. PLOTTING THE SURVIVAL AND MORTALITY CURVES---------------------------

#17.5.1 Function for extracting demorgaphic parameters from the model----
plot_data_agetemp <- list()
for (demv in c("mort", "surv")) {
  for (icat in seq_along(adtempagebastaweibull.weak[[demv]])) {
    cuts <-adtempagebastaweibull.weak$cuts[[icat]]
    minAge <- as.numeric(adtempagebastaweibull.weak$modelSpecs["min. age"])
    xx <- adtempagebastaweibull.weak$x[cuts] + minAge
    yy <- adtempagebastaweibull.weak[[demv]][[icat]][, cuts]
    
    # Convert to long format
    df <- data.frame(
      Age = xx,
      Rate = yy[1, ], # Median
      LowerCI = yy[2, ], # Lower Confidence Bound
      UpperCI = yy[3, ], # Upper Confidence Bound
      Category = names(adtempagebastaweibull.weak[[demv]])[icat],
      Type = demv # Mortality or Survival
    )
    plot_data_agetemp[[length(plot_data_agetemp) + 1]] <- df
  }
}

# Combine all into one data frame
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

#making sure the levels are in  the right order for the plot:
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

#median Survival times
surv_df_median<-surv_df_agetemp %>%
  group_by(Temp, Parental_Age) %>%  
  arrange(Age_days) %>%
  mutate(CI = 1 - Rate) %>% 
  summarise(median_line = approx(CI, Age_days, xout = 0.5)$y)

#17.3.1. cumulative survival plot-----------
adagetemp_survival_plot<- ggplot(data = surv_df_agetemp,
                               aes(x = Age_days/7, 
                                   y = Rate,
                                   colour = Parental_Age))+
  geom_line(data = surv_df_agetemp, 
            aes(x = Age_days/7, 
                y = Rate,
                colour = Parental_Age),
            linewidth=4, alpha=1)+
  geom_ribbon(data = surv_df_agetemp, 
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
  scale_x_continuous(limits = c(0, 20), breaks = c(0, 3, 6, 9, 12, 15, 18))+
  scale_color_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' age at reproduction", values = c(0.6,0.1, 0.6))+
  scale_fill_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  labs(x = "Offspring adult age (weeks)", y = "Cumulative adult survival probability, S(x)")


# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/004_survival_tempandage_adultlifespan.png",
       plot = adagetemp_survival_plot, 
       bg="transparent",
       device = "png", 
       width = 370, 
       height = 440, 
       units = "mm")


#17.3.2. Mortality risk plot---------------------
adagetemp_mort_plot<- ggplot(data = mort_df_agetemp,
                           aes(x = Age_days/7, 
                               y = Rate/52.1429,
                               colour = Parental_Age))+
  geom_line(data = mort_df_agetemp, 
            aes(x = Age_days/7, 
                y = Rate/52.1429,
                color = Parental_Age),
            linewidth = 4, alpha=1)+
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
  scale_x_continuous(limits = c(0, 20), breaks = c(0, 3, 6, 9, 12, 15, 18))+
  scale_color_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' age at reproduction", values = c(0.5,0.1, 0.5))+
  scale_fill_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  labs(x = "Offspring adult age (weeks)", y = "Instantaneous hazard rate, μ(x)")

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/004_rate_tempandage_adultlidfe.png",
       plot = adagetemp_mort_plot, 
       bg="transparent",
       device = "png", 
       width = 370, 
       height = 440, 
       units = "mm")


#Combination plot for the appendix----------------------------------------------


agetemp_b1<-adagetemp_b0_plot+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        legend.text = element_text(size=70))

agetemp_b0<-adagetemp_b1_plot+
  labs( x= "Parameter value")+
  theme(axis.title.y = element_blank())

agetemp_c<- adagetempc_plot+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


agetemp_surv<- adagetemp_survival_plot+
  labs( y = "Survival, S(x)")+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "none")

agetemp_mort<- adagetemp_mort_plot+
  labs( y = "Hazard rate, μ(x)")+
  theme(axis.title=element_text(size=50),
        axis.text=element_text(size=50))+
  theme(legend.position = "none")


#-------Combination plot for appendix------------------
tempage_inference<-ggarrange(
  ggarrange(agetemp_b1,agetemp_b0, agetemp_c, ncol= 3, labels = c("b1", "b0", "c"),
            font.label = list(size = 70, face = "bold"),
            label.y = 1.01, 
            label.x = .01,
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

#saving
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/003_tempageinference.png",
       plot = tempage_inference, 
       device = "png", 
       width = 1000, 
       height = 1000, 
       units = "mm")


# ==============================================================================
# 18. ADDING COVARIATTES: INTERACTIVE EFFECT OF PARENTAL AGE & OFFSPRING SEX----
# ==============================================================================

#Manually setting the interaction term for the KM plot
KMdataagesex<- adcovdata %>%
  mutate(sex_ParentalAge = interaction(Timepoint_binnedplot, F1_sex))

#18.2. Kaplain-Maier curve for parental age------------

#KM survival model
km_fit_parentalagesex<- survfit(Surv(adult_lifespan,event) ~ sex_ParentalAge, data = KMdataagesex)
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
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/005_KM_parentalagesex_adultlifespan.png",
       plot = Km_survival_parentalagesex,
       bg = "transparent",
       device = "png",
       width = 380,
       height = 340,
       units = "mm")


#18.3. PARAMETRIC SURVIVAL MODELS FOR PARENTAL AGE X OFFSPRING SEX-----------------------------------------

#18.3.1. A weibull model with just a makeham term------
adsexagebastaweibull.2<- basta(object = adcovdata, 
                             dataType = "census",
                             formulaMort= ~F1_sex:Timepoint_binned-1, 
                             model="WE",
                             shape="Makeham",
                             niter=60000,
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(adsexagebastaweibull.2)
saveRDS(adsexagebastaweibull.2, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_sexage_weibullmakeham.rda")
adsexagebastaweibull.2<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_sexage_weibullmakeham.rda")

#Assessing model fit and plotting the model outputs
plot(adsexagebastaweibull.2, plot.type = "gof")
plot(adsexagebastaweibull.2, plot.type = "demorates")
plot(adsexagebastaweibull.2, densities=TRUE)

#18.3.2.Weibull model--------------------------------
adsexagebastaweibull<- basta(object = adcovdata, 
                           dataType = "census",
                           formulaMort= ~F1_sex:Timepoint_binned-1, 
                           model="WE",
                           niter=60000,
                           burnin=1001,
                           thinning=50,
                           nsim = 4, parallel = TRUE, ncpus = 4)
summary(adsexagebastaweibull)
saveRDS(adsexagebastaweibull, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_adsexage_weibull.rda")
adsexagebastaweibull<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_adsexage_weibull.rda")


#18.3.4.Weibull model with a bathtub family-----------------
adsexagebastaweibull.3<- basta(object = adcovdata, 
                             dataType = "census",
                             formulaMort= ~F1_sex:Timepoint_binned-1, 
                             model="WE",
                             shape="bathtub",
                             niter=60000,
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(adsexagebastaweibull.3)
saveRDS(adsexagebastaweibull.3, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_adsexage_weibull_bathtub.rda")
adsexagebastaweibull.3<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_adsexage_weibull_bathtub.rda")

# ==============================================================================
# 19. ADDING WEAKLY INFORMATIVE PRIORS TO THE PARENTAL AGE X OFFSPRING SEX MODEL----
# ==============================================================================

#19.1. Adding weak priors into the model------------------------------

#Weakly informative priors (mean)
weakMean4 <- matrix(c(0, 1.5, prior_scale_mu), 
                    nrow = 6, ncol = 3, byrow = TRUE)

#weakly informative priors (SD)
weakSd4 <- matrix(c(
  1.0, 1.0, 1.0
), nrow = 6, ncol = 3, byrow = TRUE)

weakLower4<-matrix(c(
  0, 0, 0
), nrow = 6, ncol = 3,byrow = TRUE)
#lower bound

#19.2. Fitting weak priors to the weibull makeham model-----------
adsexagebastaweibull.weak<- basta(object = adcovdata, 
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
summary(adsexagebastaweibull.weak)
saveRDS(adsexagebastaweibull.weak, file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_age_sexweak.rda")
adsexagebastaweibull.weak<-readRDS(file = "scripts/model_outputs/BaSTA/adult lifespan/adult_lifespan_age_sexweak.rda")

#Assessing model fit and plotting the model outputs
plot(adsexagebastaweibull.weak, plot.type = "gof")
plot(adsexagebastaweibull.weak, plot.type = "demorates")
plot(adsexagebastaweibull.weak, densities=TRUE)


# ============================================================================================
# 20. PLOTTING THE POSTERIOR DISTRIBUTIONS FOR MORTALITY PARAMETERS AND DEMOGRAPHIC CURVES----
# ============================================================================================

#20.1. POSTERIOR DISTRIBUTIONS FOR B0, B1 AND C----------------------------------------------------

#extracting  mortality parameters
adsexage_theta_params <- adsexagebastaweibull.weak$params

#Converting to a data frame
adsexage_theta_params<-as.data.frame(adsexage_theta_params)

#Reshaping the data
adsexageposterior_df_long <- adsexage_theta_params %>%
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
adsexageposterior_df_long$Parental_Age<- factor(adsexageposterior_df_long$Parental_Age, 
                                               levels = c("Early-Aged", "Middle-Aged", "Late-Aged"))


#20.2. Creating seperate data sets for c, `b0` and `b1`
sexageb0_df <- adsexageposterior_df_long %>% filter(Type == "b0")
sexageb1_df <- adsexageposterior_df_long %>% filter(Type == "b1")
sexagecdf<- adsexageposterior_df_long%>% filter(Type == "c")

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
  facet_wrap(~F1_sex, ncol=1)+
  labs(
    x = "b0 parameter value",
    y = "Probability density, f(x)",
    fill = "Parents' age at reproduction",
    color = "Parents' age at reproduction",
    alpha="Parents' age at reproduction") +
  scale_x_continuous(limits = c(2.8, 6.2), breaks = c(3, 4, 5, 6))+
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
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/005-b0_tempandage_adultlifespan.png",
       plot = sexage_b0_plot, 
       bg="transparent",
       device = "png", 
       width = 420, 
       height = 200, 
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
  facet_wrap(~F1_sex, ncol=1)+
  scale_x_continuous(limits = c(3.8, 5.9), breaks = c(4, 4.6, 5.2, 5.8))+
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

ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/005_adult_posteriors_b1_sexandage.png",
       plot = sexage_b1_plot, 
       bg="transparent",
       device = "png", 
       width = 310, 
       height = 420, 
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

ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/005_posteriors_c_sexandage.png",
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
  for (icat in seq_along(adsexagebastaweibull.weak[[demv]])) {
    cuts <-adsexagebastaweibull.weak$cuts[[icat]]
    minAge <- as.numeric(adsexagebastaweibull.weak$modelSpecs["min. age"])
    xx <- adsexagebastaweibull.weak$x[cuts] + minAge
    yy <- adsexagebastaweibull.weak[[demv]][[icat]][, cuts]
    
    # Convert to long format
    df <- data.frame(
      Age = xx,
      Rate = yy[1, ], # Median
      LowerCI = yy[2, ], # Lower Confidence Bound
      UpperCI = yy[3, ], # Upper Confidence Bound
      Category = names(adsexagebastaweibull.weak[[demv]])[icat],
      Type = demv # Mortality or Survival
    )
    plot_data_agesex[[length(plot_data_agesex) + 1]] <- df
  }
}

#Combine all into one data frame
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
mort_df_agesex <- plot_data_agesex %>% filter(Type == "mort")
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
            linewidth = 4, alpha=1)+
  geom_ribbon(data = surv_df_agesex, 
              aes(y=NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI, 
                  color= Parental_Age,
                  fill=Parental_Age, 
                  alpha = Parental_Age), 
              linewidth =1,
              linetype="dashed") +
  geom_segment(data = surv_df_median,
               aes(x = 0, xend = median_line/7, y = 0.5, yend = 0.5, color = Parental_Age),
               linetype = "dashed", linewidth = 2, alpha = 0.9) +
  geom_segment(data = surv_df_median,
               aes(x = median_line/7, xend = median_line/7, y = 0, yend = 0.5, color = Parental_Age),
               linetype = "dashed", linewidth = 2, alpha = 0.9)+
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
  )+scale_x_continuous(limits = c(0, 18), breaks = c(0, 3, 6, 9, 12, 15, 18))+
  scale_color_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' age at reproduction", values = c(0.5,0.1, 0.5))+
  scale_fill_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  labs(x = "Offspring adult age (weeks)", y = "Cumulative survival probability, S(x)")


# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/005_adult_survival_ageandsex.png",
       plot = sexage_survival_plot, 
       bg="transparent",
       device = "png", 
       width = 320, 
       height = 430, 
       units = "mm")

#20.6.5. Mortality risk plot-----------------------
agesex_mort_plot<- ggplot(data = mort_df_agesex,
                          aes(x = Age_days/7, 
                              y = Rate/52.1429,
                              colour = Parental_Age))+
  geom_line(data = mort_df_agesex, 
            aes(x = Age_days/7, 
                y = Rate/52.1429,
                color = Parental_Age),
            linewidth=4, alpha=1)+
  geom_ribbon(data = mort_df_agesex, 
              aes(y=NULL, 
                  ymin = LowerCI/52.1429, 
                  ymax = UpperCI/52.1429, 
                  color= Parental_Age,
                  fill= Parental_Age, 
                  alpha = Parental_Age),
              linewidth =1,
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
  scale_x_continuous(limits = c(0, 18), breaks = c(0, 3, 6, 9, 12, 15, 18))+
  scale_color_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  scale_alpha_manual(name = "Parents' age at reproduction", values = c(0.5,0.1, 0.5))+
  scale_fill_manual(name = "Parents' age at reproduction", values = c("#f96161","#66b2b2", "#066594"))+
  labs(x = "Offspring adult age (weeks)", y = "Instantaneous hazard rate, μ (X)")

# to save plot
ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/005_adult_mort_plot_ageandsex.png",
       plot = agesex_mort_plot, 
       bg="transparent",
       device = "png", 
       width = 320, 
       height = 430, 
       units = "mm")


#-----------------Combined plot for supplementary results-----------------------

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


ggsave(filename = "./bayesian_plots/BaSTA plots/adult lifespan/004_sexageinference.png",
       plot = sexage_inference, 
       device = "png", 
       width = 900, 
       height = 800, 
       units = "mm")


###############################################################################
##END OF MAIN SCRIPT----------------------------------------------------------------
###############################################################################

