###############################################################################
##  Script: Looking at correlations between variables (paternal and maternal age, number of mating attempts and delta age)
###############################################################################

# ===========================================================================
# 1. SETUP-------------------------------------------------------------------
# ===========================================================================

## Libraries ------------------------------------------------------------------

library(plotly) #v.4.10.4
library(ggstatsplot) #v.0.13.0
library(ggplot2) #v.4.0.0
library(dplyr) #v.1.1.4
library(tidyr)  #v.1.3.1
library(openxlsx) #v.4.2.5.2

## Load data ------------------------------------------------------------------
F1data <- readRDS('./raw data/F1_filtered_data06022025.RDS')
f1stats<-read.xlsx("./raw data/F1_sizetest_28012025.xlsx") #for subset of animals with pronotum measurements 

#ensuring f1stats are correctly formatted
f1stats$Pronotum_width<-as.numeric(f1stats$Pronotum_width)
f1stats$Pronotum_length<-as.numeric(f1stats$Pronotum_length)
f1stats$Nymph_ID<-as.factor(f1stats$Nymph_ID)
f1stats$Adult.mass<-as.numeric(f1stats$Adult.mass)

#Need to strip this back so there is only one offspring per parent per mating attempt
prepdata2<- F1data %>%
  group_by (Mother_ID)%>%
  arrange(Timepoint)%>%
  distinct(Timepoint, .keep_all=TRUE)#removes duplicates from each timepoint per mother

#Correlations between age measures and number of remaining mating attempts
cor.test(prepdata2$within_subject_age, prepdata2$cum_successful_matings) #0.9367288

#Plotting the correlation between within-subject age and the number of mating attempts remaining 
mating_v_age<-ggscatterstats(data = prepdata2, 
                             x = within_subject_age, 
                             y = cum_successful_matings,
                             point.height.jitter = 0.2,
                             point.width.jitter = 0.1,
                             point.label.args = list(alpha = 0.7, size = 4, color = "grey50"),
                             xsidehistogram.args = list(fill = "#CC79A7"), ## fill for marginals on the x-axis
                             ysidehistogram.args = list(fill = "#009E73"),
                             type ="bayes")+ ## fill for marginals on the y-axis
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8))+
  labs(x = "Parents' delta age (weeks; mean-centred)",
       y = "Cumulative number of mating attempts")

ggsave(filename = "./bayesian_plots/Offspring trait plots/mating_versus_age.png",
       plot = mating_v_age, 
       device = "png", 
       bg="transparent",
       width = 200, 
       height = 210, 
       units = "mm")


#Checking correlations between pronotum width, length, and size:
sizedata<-f1stats %>% 
  select(Nymph_ID, Pronotum_width, Pronotum_length, Adult.mass)

#How correlated are the F1 body size measurements
cor.test(f1stats$Adult.mass, f1stats$Pronotum_width,method = "pearson")#correlation test to see if egg weight and egg volume are correlated
cor.test(f1stats$Adult.mass, f1stats$Pronotum_length, method="pearson")#correlation test to see if egg weight and egg volume are correlated
cor.test(f1stats$Pronotum_width, f1stats$Pronotum_length, method="pearson")#correlation test to see if egg weight and egg volume are correlated

#correlation between pronotum length and adult mass
mass_pronolength<-ggscatterstats(data = f1stats, 
                                 x = Adult.mass, 
                                 y = Pronotum_length)+
  theme_classic()+labs(x = "F1 Adult Mass (g)", y = "F1 Pronotum Length (mm)")

ggsave(filename = "./bayesian_plots/Offspring trait plots/F1_mass_and_pronotum_length.png",
       plot = mass_pronolength, 
       device = "png", 
       bg="transparent",
       width = 200, 
       height = 210, 
       units = "mm")

#correlation between pronotum width and adult mass
mass_pronowidth<-ggscatterstats(data = f1stats, 
                                x = Adult.mass, 
                                y = Pronotum_width)+
  theme_classic()+labs(x = "F1 Adult Mass (g)", y = "F1 Pronotum width (mm)")

ggsave(filename = "./bayesian_plots/Offspring trait plots/F1_mass_and_pronotum_width.png",
       plot = mass_pronowidth, 
       device = "png", 
       bg="transparent",
       width = 200, 
       height = 210, 
       units = "mm")

#correlation between prnotum width and pronotum length
pronolength_pronowidth<-ggscatterstats(data = f1stats, 
                                       x = Pronotum_length, 
                                       y = Pronotum_width)+
  theme_classic()+labs(x = "F1 Pronotum Length (mm)", y = "F1 Pronotum width (mm)")

ggsave(filename = "./bayesian_plots/Offspring trait plots/F1_pronotum_length_and_pronotum_width.png",
       plot = pronolength_pronowidth, 
       device = "png", 
       bg="transparent",
       width = 200, 
       height = 210, 
       units = "mm")

###############################################################################
##  END OF SCRIPT----------------------------------------------------------------
###############################################################################

