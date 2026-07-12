
#####Script for creating a forest plot of standardised mean differences for every trait

#Libraries
library(ggplot2) 
library(dplyr) 
library(tidyverse)
library(tidyr)
library(brms) 
library(tidybayes)
library(lme4) 
library(bayesplot)
library(tidybayes)
library(projpred)


#----------Loading necessary SMD estimates----------------------
hatching_SMD <-readRDS("scripts/model_outputs/Offspring Trait Models/hatching_SMD")
fecundity_SMD<-readRDS("scripts/model_outputs/Offspring Trait Models/fecundity_SMD")
dev_time_SMD<-readRDS("scripts/model_outputs/Offspring Trait Models/dev_time_SMD")
dev_time_SMD<-readRDS("scripts/model_outputs/Offspring Trait Models/dev_time_SMD")
juvsurv_SMD<-readRDS("scripts/model_outputs/Offspring Trait Models/juvsurv_SMD")
early_life_surv_SMD<-readRDS("scripts/model_outputs/Offspring Trait Models/early_life_surv_SMD")
Total_lifespan_SMD<-readRDS("scripts/model_outputs/Offspring Trait Models/Total_lifespan_SMD")
mass_SMD<-readRDS("scripts/model_outputs/Offspring Trait Models/mass_SMD")
adult_lifespan_SMD<-readRDS("scripts/model_outputs/Offspring Trait Models/adult_lifespan_SMD")

#--------Combining all traits into one---------------------------
all_SMD <- bind_rows(
  early_life_surv_SMD %>% mutate(Trait = "Early life survival"),
  juvsurv_SMD         %>% mutate(Trait = "Juvenile survival"),
  Total_lifespan_SMD  %>% mutate(Trait = "Total lifespan"),
  dev_time_SMD        %>% mutate(Trait = "Development time"),
  mass_SMD            %>% mutate(Trait = "Adult mass"),
  fecundity_SMD       %>% mutate(Trait = "Fecundity"),
  hatching_SMD        %>% mutate(Trait = "Hatching success"),
  adult_lifespan_SMD  %>% mutate(Trait = "Adult lifespan")
)%>%
  mutate(
    Trait = factor(
      Trait,
      levels = c(
        "Early life survival",
        "Juvenile survival",
        "Total lifespan",
        "Development time",
        "Adult mass",
        "Fecundity",
        "Hatching success",
        "Adult lifespan"
      )
    )
  )


#Creating a ggplot with estimates

forestplot<-ggplot(data = all_SMD,
       aes(x = SMD, 
           y = Trait,
           colour = Trait)) +
  coord_flip()+
  # --- Shaded background for positive/negative SMD ---
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, 
           fill = "#FBBABA", alpha = 0.3) +   # light red for decrease
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "#BFFFC1", alpha = 0.3) +   # light green for increase
    geom_linerange(data = all_SMD,
                 aes(x=SMD,
                     xmin = SMD_lower,
                     xmax = SMD_upper,
                     colour=Trait),
                 linewidth = 3)+
  geom_point(data = all_SMD,
             aes(x = SMD,
                 y = Trait,
                 fill=Trait),
             shape=21, 
             stroke=1.8,
             size= 18,
             alpha=1,
             colour="white",
             show.legend = FALSE)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1)+
  theme_minimal() + 
  theme(axis.title=element_text(size=30),
        axis.text=element_text(size=28))+
  theme(legend.position = "none") +
  scale_color_manual(name = "Trait", values = c(  "#E89BD0", 
                                                  "#D4712E",
                                                  "#4F4F4F", 
                                                  "#5FA3B3",
                                                  "#4F6D7A", 
                                                  "#A8D5BA", 
                                                  "#3B7A57", 
                                                  "#F96161"))+
  scale_fill_manual(name = "Trait", values = c(  "#E89BD0", 
                                                 "#D4712E",
                                                 "#4F4F4F", 
                                                 "#5FA3B3",
                                                 "#4F6D7A", 
                                                 "#A8D5BA", 
                                                 "#3B7A57", 
                                                 "#F96161"))+
  labs(x = "Standardised mean difference (SMD)", y = "Offspring Trait")



# to save plot
ggsave(filename = "./bayesian_plots/offspring trait plots/forest_plot.png",
       plot = forestplot, 
       bg="transparent",
       device = "png", 
       width = 680, 
       height = 200, 
       units = "mm")
 
