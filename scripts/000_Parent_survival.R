
###############################################################################
##  Script: Longevity and removal times of parents' in the three temperature treatments
##  Note:   Analysis to investigate parent longevity and removal times
#small sample sizes, and the complexities with the setup (i.e., removing parents to climate room),
#this functions more as a visual, exploratory confirmation that temperature affects removal/survival times of the parents
###############################################################################

# ===========================================================================
# 1. SETUP--------------------------------------------------------------------
# ===========================================================================

## Libraries ------------------------------------------------------------------
library(dplyr) #v.1.1.4
library(ggplot2) #v.4.0.0
library(survival) #v.3.8.3
library(survminer) #v.0.5.1
library(ggpubr) #v.0.6.1
library(cowplot) #v.1.1.3

# ===========================================================================
# 2. DATA ORGANISATION--------------------------------------------------------
# ===========================================================================
#read data for animals that are to be included in this analysis#
F0data<-readRDS('./raw data/F0_final_data_08012026.RDS')
length(unique(F0data$Cricket_ID))#156 parents (i.e., 78 parent pairs)

#Formatting data
F0data$Cricket_ID<-as.factor(F0data$Cricket_ID)
F0data$event<-as.numeric(F0data$event)
F0data$lifespan<-as.numeric(F0data$lifespan)
F0data$Temperature_treatment<-as.factor(F0data$Temperature_treatment)

# ===========================================================================
# 3. SUMMARY STATISTICS--------------------------------------------------------
# ===========================================================================

#Looking at the lifespan of parents that died within their experimental temperature period
#I.e., these animals actually died in the incubators, and weren't removed due to reproductive senescence
hist_died<-ggplot(F0data %>% 
               filter(event == 0), aes(x = time_in_study)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") + 
  facet_wrap(~Temperature_treatment) +  # Create panels for each Temp level
  labs(
    x = "Lifespan (weeks)",
    y = "Number died"
  ) +
  theme_classic() #A greater number of individuals died in the experimental period at 30.5 celsius


#Looking at the time spent in study for removed animals
#This tests reproductive senescence: i.e., what ages did individuals stop producing viable eggs at?
hist_removed<-ggplot(F0data %>% 
                    filter(event == 1) %>% 
                    filter(Sex == "F"), aes(x = time_in_study)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") + 
  facet_wrap(~Temperature_treatment) +  # Create panels for each Temp level
  labs(
    x = "Time in study (weeks)",
    y = "Number removed from incubators"
  ) +
  theme_classic()
#more individuals removed at week 10 (i.e., the end point of the study)
#individuals were removed from their temperature treatments at week 10 to free up the incubators


#1. Overall Mean F0 time in study and standard deviation for animals that didn't die at their set temperature treatments-----
sum_removal<-F0data %>% 
  filter(event ==1) %>%  #time in study for animals that didn't die in the incubators
  summarise(F0_time= mean(time_in_study),
            sd_time= sd(time_in_study),
            n_animals=n(),
            n_parents =n_distinct(Mother_ID))
#parents, on average, spent 8.22 weeks in the study
#F0_time     sd_time     n_animals n_parents
#1 8.226891 1.703805       102        12

#2. Overall mean F0 lifespan for animals that did die at their set temperature treatments----
sum_lifespan<-F0data %>% 
  filter(event ==0) %>% #lifespan of animals that did die in the incubators
  summarise(F0_time= mean(lifespan),
            sd_time= sd(lifespan),
            n_animals=n(),
            n_parents =n_distinct(Mother_ID))

#Mean lifespan for animals that died:
#F0_time     sd_time      n_animals n_parents
# 7.174603  2.263651        54        10

###1. Median time in study for censored animals at each temperature treatment----
  sum_dat1<-F0data %>%
  filter(event == 1) %>% 
  group_by(Temperature_treatment)%>%
  summarise(duration_in_study= mean(time_in_study),
            se_total= sd(time_in_study),
            n_animals=n(),
              n_parents = n_distinct(Mother_ID))

# Temperature_treatment duration_in_study se_total n_animals n_parents
#25.5                               8.89     1.46        46        10
#28                                 7.85     1.65        31         9
#30.5                               7.48     1.79        25        11

#Mean lifespan for animals that died at their temperature treatments-----
sum_dat2<-F0data %>% 
  filter(event == 0) %>% 
  group_by(Temperature_treatment)%>%
  summarise(mean_lifespan= mean(time_in_study),
            se_total= sd(time_in_study),
            n_animals=n(),
            n_parents = n_distinct(Mother_ID))

#Temperature_treatment mean_lifespan se_total n_animals n_parents
#25.5                           7.77     2.42         8         5
#28                             6.97     2.79        15         7
#30.5                           7.12     1.98        31         9


#When deaths and experimental exits are viewed as the same process
sum_dat2<-F0data %>% 
  group_by(Temperature_treatment)%>%
  summarise(mean_lifespan= mean(time_in_study),
            se_total= sd(time_in_study),
            n_animals=n(),
            n_parents = n_distinct(Mother_ID))

#Temperature_treatment mean_lifespan se_total n_animals n_parents
#25.5                           8.72     1.66        54        10
#28                             7.57     2.10        46         9
#30.5                           7.28     1.89        56        11



# ===========================================================================
# 4. Survival Curves for parent animals--------------------------------------
# ===========================================================================

#4.1: KM curve with censoring of removed animals-------
#Here animals that were removed from the study before death was observed were right cesnored in the survival curve

#Plotting a Kaplein-Maier curve for parental temperature 
F0data2<- F0data %>%
  mutate(Temperature_treatment = factor(case_when(
                            Temperature_treatment == "25.5" ~ "25.5°C",
                            Temperature_treatment == "28" ~ "28.0°C",
                            Temperature_treatment == "30.5" ~ "30.5°C"
                            ),
                            levels = c("25.5°C", "28.0°C", "30.5°C")))

#setting all events to observed deaths (animals that didn't die will be right censored here)
F0data3<-F0data %>% 
  mutate(event = 1-as.numeric(event)) #switching the censoring around, so those that were observed dead = 1


#Calculating lifespan (in weeks), with right censoring for animals that were removed to the climate room
#1.Kaplain-Maier curve for parental temperature
km_fit_Temp <- survfit(Surv(time_in_study,event) ~ Temperature_treatment, data = F0data3) #fitting the curve for Timepoint
Km_survival_Temp <- ggsurvplot(km_fit_Temp, 
                               data = F0data3, 
                               pval = FALSE,  
                               conf.int = TRUE,  
                               risk.table = "nrisk_cumcensor",
                               break.time.by = 3, 
                               xlim = c(0, 12),
                               xlab = "Age (weeks)", 
                               ylab = "Cumulative survival probability, S(x)",
                               font.x = 35,
                               font.y = 35,
                               font.tickslab = c(35, "grey25"),
                               font.legend = 25,
                               risk.table.fontsize = 10,
                               tables.theme = theme_survminer(
                                 font.main      = c(28, "bold"),  
                                 font.x         = c(26),           
                                 font.y         = c(26),           
                                 font.tickslab  = c(24)            
                               ),
                               legend.title = "Temperature Treatment",
                               legend.labs = c("25.5°C", "28.0°C", "30.5°C"),
                               palette = c("#2f4b7c","orange2", "#a64b61"))

# Combine the survival plot and the risk table
Km_temp_combined_plot <- cowplot::plot_grid(Km_survival_Temp$plot, Km_survival_Temp$table, ncol = 1, rel_heights = c(3, 1))

# Save the combined plot
ggsave(filename = "./bayesian_plots/Offspring trait plots/parent lifespan/002_KM_parentlifespan.png",
       plot = Km_temp_combined_plot,
       bg = "transparent",
       device = "png",
       width = 420,
       height = 360,
       units = "mm")

#4.2. Kaplan Meier curve when removal and death are treated as the same event (i.e., this includes natural detahs and removals due to reproductive senescence)
F0data4<-F0data %>% 
  mutate(event =1) #All animals included (no right censoring)

#Kaplain-Maier curve
km_fit_all <- survfit(Surv(time_in_study,event) ~ Temperature_treatment, data = F0data4) #fitting the curve for Timepoint
Km_survival_all <- ggsurvplot(km_fit_all, 
                               data = F0data4, 
                               pval = FALSE,  
                               conf.int = TRUE,  
                               risk.table = "absolute",
                               break.time.by = 3, 
                               xlim = c(0, 12),
                               xlab = "Time spent in study (weeks)", 
                               ylab = "Cumulative probability of experimental exit, S(x)",
                               font.x = 35,
                               font.y = 35,
                               font.tickslab = c(35, "grey25"),
                               font.legend = 25,
                               risk.table.fontsize = 10,
                               tables.theme = theme_survminer(
                                 font.main      = c(28, "bold"),  
                                 font.x         = c(26),           
                                 font.y         = c(26),           
                                 font.tickslab  = c(24)            
                               ),
                               legend.title = "Temperature Treatment",
                               legend.labs = c("25.5°C", "28.0°C", "30.5°C"),
                               palette = c("#2f4b7c","orange2", "#a64b61"))

# Combine the survival plot and the risk table
Km_survival_all <- cowplot::plot_grid(Km_survival_all$plot, Km_survival_all$table, ncol = 1, rel_heights = c(3, 1))

# Save the combined plot
ggsave(filename = "./bayesian_plots/Offspring trait plots/parent lifespan/002_KM_parentexit.png",
       plot = Km_survival_all,
       bg = "transparent",
       device = "png",
       width = 420,
       height = 360,
       units = "mm")


#Creating a combined plot
temp_inference<-ggarrange(
  ggarrange(Km_survival_all, 
            Km_temp_combined_plot, 
            nrow = 1,
            font.label = list(size = 30, face = "bold"),
  ncol = 2,
  labels = c("A", "B"),
  align = "h",
  widths = c(1,1)))


ggsave(filename = "./bayesian_plots/Offspring trait plots/parent lifespan/001_inference_temperature.png",
       plot = temp_inference, 
       device = "png", 
       width = 700, 
       height = 430, 
       units = "mm")

###############################################################################
##END OF SCRIPT----------------------------------------------------------------
###############################################################################


