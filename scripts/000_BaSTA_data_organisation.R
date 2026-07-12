
###############################################################################
##  Script: Tidying data and formatting it correctly for each BaSTA analysis
##  Note:   BaSTA requires data to be formatted correctly,
#ensuring our data frame has values in the right place to be fitted in BaSTA
###############################################################################

# ===========================================================================
# 1. SETUP--------------------------------------------------------------------
# ===========================================================================
library(dplyr) #v.1.1.4
library(tidyr) #v.1.3.1
library(lubridate) #v.1.9.3
library(BaSTA) #v.2.0.2
library(snowfall) #v.1.84.6.3

##Clearing the workspace
rm(list=ls())

# ===========================================================================
# 2. DATA ORGANISATION--------------------------------------------------------
# ===========================================================================

#reading in the F1 data, this will be used to create all BaSTA datasets
F1datas<-readRDS("./raw data/F1_filtered_data06022025.RDS")
F1datas<-as.data.frame(F1datas)


#2.1. Adult lifespan data set

#Here, the "birth date" is the adult emergence date, and "death dates" are the actual death dates
#this dataset excludes animals which a) died before adulthood, and b) have uncertain adult emergence dates

#Filtering out the high early life mortality-------
adultF1data<-F1datas %>% 
  rename(ID=F1_ID) %>% 
  filter(Include_in_adult=="Y") %>% 
  filter(adult_surv==1) %>% #only including adults here
  rename(Birth.Date=F1_Adult_emergence) %>% #Isolates the analysis to just look at the adult lifespan (difference between the adult emergence and death dates)
  rename(Depart.Date=F1_death_date)
length(unique(adultF1data$ID))#Only looking at 939 animals that reached adulthood and had complete life history data


#Adult mortality dataset ------------
adultBastaData<-adultF1data %>% 
  select(ID, Birth.Date, Depart.Date, F1_sex, Timepoint, Temp, PairID) %>% 
  filter(!is.na(Birth.Date))#this isolates the data to only look at adult mortality trajectories


#Animals that are still alive need a C in depart. type, animals that died need a D
adultBastaData$Depart.Type<-ifelse(!is.na(adultBastaData$Depart.Date), "D", "C")#assigned a value of D if the individual died in the study, C if not (didn't occur, we recorded all deaths)


#creating variables for maximum and minimum adult emergence date----------
adultBastaData<-adultBastaData %>% 
  mutate(Min.Birth.Date= (adultBastaData$Birth.Date-6)) %>% #accounts for the fact that they could have emerged as adult anywhere within that week
  mutate(Max.Birth.Date=adultBastaData$Birth.Date) %>% 
  mutate(Entry.Date=adultBastaData$Birth.Date)

#organising the columns so they are in the right order for basta to function
adultBastaData<-adultBastaData %>% 
  select(ID, Birth.Date, Min.Birth.Date, Max.Birth.Date,Entry.Date, Depart.Date, Depart.Type, F1_sex, Timepoint, Temp, PairID)


checkedDataCens <- DataCheck(object = adultBastaData, dataType = "census",
                             silent = FALSE)#No inconsistencies between dates
print(checkedDataCens)#data seems to all be coded correctly, and I've included all the required covariates for the model to work

###Saving data for final analysis
saveRDS(object=adultBastaData, file='./raw data/BaSTA data/adult_basta_data.RDS')


#2.2. BaSTA data for total lifespan (excluding early-life mortality)----------------------------------------------------------

#Here, the "birth date" is the hatching date (i.e., the actual birth date), and "death dates" are the actual death dates
#this dataset excludes animals which died within the first four weeks of life (but includes both adults and nymphs that survived past this period)

#Filtering out the high early life mortality
totalF1data<-F1datas %>% 
  rename(ID=F1_ID) %>% 
  filter(!total_lifespan<=4) %>% #filtering out the really high early life mortality
  filter(Include_in_total=="Y") %>% 
  rename(Birth.Date=F1_hatch) %>%
  rename(Depart.Date=F1_death_date)

length(unique(totalF1data$ID))#Only looking at 987 animals that survived past week 4 of life (40 of these died before adulthood)


#total mortality dataset --> animals that died before week 4 of life are not included in this analysis
TotalBastaData<-totalF1data %>% 
  select(ID, Birth.Date, Depart.Date, F1_sex, Timepoint, Temp, PairID) %>% 
  filter(!is.na(Birth.Date))


#Animals that died are given a value of D (should be all animals)
TotalBastaData$Depart.Type<-ifelse(!is.na(TotalBastaData$Depart.Date), "D", "C")#assigned a value of D if the individual died in the study


#creating variables for maximum and minimum birth date:
TotalBastaData<-TotalBastaData %>% 
  mutate(Min.Birth.Date= (TotalBastaData$Birth.Date)) %>% 
  mutate(Max.Birth.Date=TotalBastaData$Birth.Date) %>%  
  mutate(Entry.Date=TotalBastaData$Birth.Date)


#organising the columns so they are in the right order for BaSTA to function
TotalBastaData<-TotalBastaData %>% 
  select(ID, Min.Birth.Date, Max.Birth.Date, Birth.Date,Entry.Date,Depart.Date, Depart.Type, F1_sex, Timepoint, Temp, PairID)


checkedDataCens <- DataCheck(object = TotalBastaData, dataType = "census",
                             silent = FALSE)#No inconsistencies between dates
print(checkedDataCens)#data seems to all be coded correctly, and I've included all the required covariates for the model to work

###Saving data for final analysis
saveRDS(object=TotalBastaData, file='./raw data/BaSTA data/total_basta_data.RDS')


#2.3. BaSTA data for total lifespan (Including early-life mortality)####

#Here, the "birth date" is the hatching date, and "death dates" are the actual death dates
#this dataset includes all animals (including those that died in their first four weeks of life)


#Including all F1 offspring here
totalF1data2<-F1datas %>% 
  rename(ID=F1_ID) %>% 
  filter(Include_in_total=="Y") %>% 
  rename(Birth.Date=F1_hatch) %>%
  rename(Depart.Date=F1_death_date)

length(unique(totalF1data2$ID))#Looking at the 1317 offspring


#total mortality dataset --> Includes all animals in this analysis
TotalBastaData2<-totalF1data2 %>% 
  select(ID, Birth.Date, Depart.Date, F1_sex, Timepoint, Temp, PairID) %>% 
  filter(!is.na(Birth.Date))


#Animals that died are given a value of D (should be all animals)
TotalBastaData2$Depart.Type<-ifelse(!is.na(TotalBastaData2$Depart.Date), "D", "C")#assigned a value of D if the individual died in the study


#creating variables for maximum and minimum birth date:
TotalBastaData2<-TotalBastaData2 %>% 
  mutate(Min.Birth.Date= (TotalBastaData2$Birth.Date)) %>% 
  mutate(Max.Birth.Date=TotalBastaData2$Birth.Date) %>%  
  mutate(Entry.Date=TotalBastaData2$Birth.Date)


#organising the columns so they are in the right order for BaSTA to function
TotalBastaData2<-TotalBastaData2 %>% 
  select(ID, Min.Birth.Date, Max.Birth.Date, Birth.Date,Entry.Date,Depart.Date, Depart.Type, F1_sex, Timepoint, Temp, PairID)

checkedDataCens <- DataCheck(object = TotalBastaData2, dataType = "census",
                             silent = FALSE)#No inconsistencies between dates
print(checkedDataCens)#data seems to all be coded correctly, and I've included all the required covariates for the model to work

###Saving data for final analysis
saveRDS(object=TotalBastaData2, file='./raw data/BaSTA data/total_basta_data.2.RDS')

###############################################################################
##END OF SCRIPT----------------------------------------------------------------
###############################################################################

