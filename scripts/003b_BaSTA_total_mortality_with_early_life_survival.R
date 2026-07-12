
###############################################################################
##  Script: BaSTA analysis looking at the total lifespan of the offspring, 
#including all the offspring who died in their first month
#Aim to capture steep early-life mortality using the BaSTA bathtub Weibull model
###############################################################################

# ===========================================================================
# 1. SETUP--------------------------------------------------------------------
# ===========================================================================

## Libraries ------------------------------------------------------------------
library(dplyr)  #v.1.1.4
library(BaSTA) #v.2.0.2
library(snowfall) #v.1.84.6.3
library(ggplot2) #v.4.0.0
library(survival) #v.3.8.3
library(survminer) #v.0.5.1

# ===========================================================================
# 2. DATA ORGANISATION--------------------------------------------------------
# ===========================================================================

#2.1. Importing data------------
F1data<-readRDS("./raw data/BaSTA data/total_basta_data.2.RDS")
F1data<-as.data.frame(F1data)
length(unique(F1data$ID))#1317 offspring from 78 parent pairs


#2.2. grouping the parental ages into categorical groups (rather than treating it as a continuous variable)
F1data <- F1data%>%
  mutate(Timepoint_binned = case_when(
    Timepoint %in% 1:2 ~ "early",
    Timepoint %in% 3:5 ~ "middle",
    Timepoint %in% 6:8 ~ "late"
  ))

#How many offspring do we have in each age category at each temperature?
table_temp_lifespan <- table(F1data$Temp, F1data$Timepoint_binned)
print(table_temp_lifespan)#this table gives counts for each level
        #early late middle
#25.5   131  175    189
#28     111  121    162
#30.5   133  108    187

#How many offspring per group? (averaged over temperature)
table_lifespan <- table(F1data$Timepoint_binned)
print(table_lifespan)
# early late middle 
#375    404    538 

#How many offspring in each age group per parent pair?
table_lifespan_pairID<-table(F1data$PairID, F1data$Timepoint_binned)
print(table_lifespan_pairID)


#2.3. creating dataset with the required covariates-------------
totcovdata<- F1data %>% 
  select(ID, Birth.Date, Min.Birth.Date, Max.Birth.Date, Entry.Date, Depart.Date, Depart.Type, Timepoint, Timepoint_binned, Temp,F1_sex)

#DATA CHECK: re-checking whether the filtered data passes BaSTA's built-in data check function
checkedDataCens <- DataCheck(object = totcovdata, dataType = "census",
                             silent = FALSE)
print(checkedDataCens) #Passed the BaSTA data check


# ===========================================================================
# 3. MODEL STRUCTURE SELECTION------------------------------------------------
# ===========================================================================
#Data is survival data (discrete time-to-event). Aim to capture early-life mortality so using the bathtub function
#Using , the following distributions Gompertz bathtub, Weibull bathtub, and logistic bathtub

#3.1--Gompertz Model with a bathtub function-------------------------
totalgomp.3 <- basta(object = totcovdata, 
                     dataType = "census",
                     shape="bathtub",
                     niter=60000,
                     burnin=1001,
                     thinning=50,
                     nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalgomp.3)
saveRDS(totalgomp.3, file = "scripts/model_outputs/BaSTA/total lifespan.2/total_lifespan_gompertz_bathtub.rda")
totalgomp.3 <-readRDS("scripts/model_outputs/BaSTA/total lifespan.2/total_lifespan_gompertz_bathtub.rda")


#Assessing model fit and plotting the model outputs
plot(totalgomp.3, plot.type = "gof")
plot(totalgomp.3, plot.type = "demorates")
plot(totalgomp.3, densities=TRUE)


#3.2. Weibull model with a bathtub term --------------------------
totalbastaweibull.3 <- basta(object = totcovdata, 
                             dataType = "census",
                             model="WE",
                             shape="bathtub",
                             niter=60000,#trying to double the burnin and the number of iterations
                             burnin=1001,
                             thinning=50,
                             nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastaweibull.3)
saveRDS(totalbastaweibull.3, file = "scripts/model_outputs/BaSTA/total lifespan.2/total_lifespan_weibull_bathtub.rda")
totalbastaweibull.3<-readRDS("scripts/model_outputs/BaSTA/total lifespan.2/total_lifespan_weibull_bathtub.rda")

#Assessing model fit and plotting the model outputs
plot(totalbastaweibull.3, plot.type = "gof") #Model does not capture early-life survival well
plot(totalbastaweibull.3, plot.type = "demorates")
plot(totalbastaweibull.3, densities=TRUE)

#3.3. Logistic model with a bathtub shape----------------------------------------------------------
totalbastalog.3 <- basta(object = totcovdata, 
                         dataType = "census",
                         model="LO",
                         shape="bathtub",
                         niter=60000,
                         burnin=1001,
                         thinning=50,
                         nsim = 4, parallel = TRUE, ncpus = 4)
summary(totalbastalog.3)

saveRDS(totalbastalog.3, file = "scripts/model_outputs/BaSTA/total lifespan.2/total_lifespan_log_bathtub.rda")
totalbastalog.3<-readRDS("scripts/model_outputs/BaSTA/total lifespan.2/total_lifespan_log_bathtub.rda")

# ===========================================================================
# 4. Model diagnostics--Kaplan Meir curve of observed survival compared to parametric survival curve-----------------------------------
# ===========================================================================
#Inspecting fit of the Weibull bathtub model (this seemed to cpature the distribution most effectively)

#Ensuring variables are correctly ordered
totcovdata$Timepoint_binned<-factor(totcovdata$Timepoint_binned, 
                                    levels = c("early", "middle", "late"))
#censoring variable for survival package
totcovdata$event<-ifelse(totcovdata$Depart.Type =="D", 1,0) 

#Renaming levels for the plot
totcovdata$Timepoint_binnedplot<-factor(totcovdata$Timepoint_binned, 
                                        levels = c("early", "middle", "late"),
                                        labels=c("Early-Aged", "Middle-Aged", "Late-Aged"))

#Calculating lifespan (in weeks)
totcovdata$Lifespan_weeks <- as.numeric(difftime(totcovdata$Depart.Date, totcovdata$Birth.Date, units = "weeks"))

#4.1.Kaplan-Maier curve of observed survival

#survival model
km_fit_parentalage <- survfit(Surv(Lifespan_weeks,event) ~ 1, data = totcovdata)

#KM plot
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
                                      risk.table.fontsize = 6)
#Gives the raw observed survival curve

#4.2. Function for pulling the survival estimates from the weibull bathtub model------
plot_data_total <- list()
for (demv in c("mort", "surv", "dens")) {
  for (icat in seq_along(totalbastaweibull.3[[demv]])) {
    cuts <-totalbastaweibull.3$cuts[[icat]]
    minAge <- as.numeric(totalbastaweibull.3$modelSpecs["min. age"])
    xx <- totalbastaweibull.3$x[cuts] + minAge
    yy <- totalbastaweibull.3[[demv]][[icat]][, cuts]
    
    # Convert to long format
    df <- data.frame(
      Age = xx,
      Rate = yy[1, ], # Median
      LowerCI = yy[2, ], # Lower Confidence Bound
      UpperCI = yy[3, ], # Upper Confidence Bound
      Category = names(totalbastaweibull.3[[demv]])[icat],
      Type = demv # Mortality or Survival
    )
    plot_data_total[[length(plot_data_total) + 1]] <- df
  }
}

# Combine all into one data frame
plot_data_total<- do.call(rbind, plot_data_total)

#renaming the values of Timepoint (early aged, middle aged, and late-aged)
plot_data_total$Category<- factor(plot_data_total$Category, 
                                    levels = c("Timepoint_binnedearly", "Timepoint_binnedmiddle", "Timepoint_binnedlate"), 
                                    labels = c("Early-Aged", "Middle-Aged", "Late-Aged"))


# Set the scaling factors (e.g., minimum and maximum age in days)
max_days <- 365.25 #for converting the timescale to days
min_days<-0

# Adjust Age column back to days
plot_data_total<- plot_data_total %>%
  mutate(Age_days = Age * (max_days - min_days) + min_days)

#ensuring age is a factor for plotting
plot_data_total$Category<-as.factor(plot_data_total$Category)


# Splitting the data into Mortality and Survival datasets
#For this script, just plotting the survival curves--
surv_df_total<- plot_data_total%>% filter(Type == "surv")


#Extracting the data from the Kaplan Meier curve
Km_survival_parentalage$data.survplot<-
  Km_survival_parentalage$data.survplot

# Create the zero-point for each group (ensuring both the parametric and survival curves start from the same position)
zero_point_001_age <- Km_survival_parentalage$data.survplot %>%
  summarise(time = 0, surv = 1, lower = 1, upper = 1)

# Combine zero-point and the actual KM data
Km_survival_parentalage$data.survplot<- Km_survival_parentalage$data.survplot %>%
  bind_rows(zero_point_001_age) %>%
  arrange(time)  # Make sure points are sorted for step plot

#4.3. Combining the parametric and non-parametric estimates into one plot---------------
Km_survival_combined_plot_002 <- ggplot() +
#Kaplan-Meier curve
  geom_step(data = Km_survival_parentalage$data.survplot, 
            aes(x = time, 
                y = surv), 
            colour = "#066594",
            linewidth = 3, alpha = 1)+
  geom_step(data = Km_survival_parentalage$data.survplot, 
            aes(x = time, 
                y = lower),
            colour = "#066594", 
            linetype="dashed",
            alpha = 0.5,
            linewidth=1.2)+ 
  geom_step(data = Km_survival_parentalage$data.survplot, 
            aes(x = time, 
                y = upper),
            colour = "#066594", 
            linetype="dashed",
            alpha = 0.5,
            linewidth=1.2)+
  #Parametric survival curve
  geom_line(data = surv_df_total, 
            aes(x = Age_days / 7, #converting the timescale from days to weeks
                y = Rate),
            colour = "#f96161", 
            size = 2.5, alpha = 1) +
  geom_ribbon(data = surv_df_total, 
              aes(x = Age_days / 7, 
                  y = NULL, 
                  ymin = LowerCI, 
                  ymax = UpperCI),
              fill = "#f96161", 
              alpha = 0.3)+
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

ggsave(filename = "./bayesian_plots/BaSTA plots/total lifespan/000_total_survival.png",
       plot = Km_survival_combined_plot_002, 
       bg="transparent",
       device = "png", 
       width = 620, 
       height = 480, 
       units = "mm")

###############################################################################
##END OF SCRIPT----------------------------------------------------------------
###############################################################################
