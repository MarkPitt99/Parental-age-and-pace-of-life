---
# **Datasets for the paper:** "*Parental age accelerates offspring pace-of-life*"

This folder contains the data that was collected and analysed as part of a project investigating parental age effects 
on offspring life-history traits.

------------------------------------------------------------------------

Mark D. Pitt, Brendan O’Connor, Timothy D. Sheen, Davide M. Dominoni, Tom Tregenza, Jelle J. Boonekamp **Parental age at reproduction accelerates offspring pace of life**. 
*bioRxiv*. [DOI:10.64898/2026.04.14.718189](https://www.biorxiv.org/content/10.64898/2026.04.14.718189v2.full)

------------------------------------------------------------------------

For any further information, please contact: Mark Pitt, email: [markdavidpitt\@gmail.com](mailto:markdavidpitt@gmail.com)

## Description of data:

The data from this project is contained in three seperate RDS files.

------------------------------------------------------------------------

### RDS FILE 1: F1_filtered_data06022025.RDS

This sheet contains all the necessary information required to model the following offspring life-history traits:
early-life survival (i.e., first month survival), juvenile survival (i.e., survival from the first month to adulthood), 
total lifespan (i.e., survival from birth until death), development time (time between birth and adult emergence), adult lifespan (i.e., survival from adult emergence until death), 
and adult mass (in grams). This data was also used to create the data sheets required for the BaSTA analysis on the offspring's mortality parameters(see: 000_ BaSTA_data_organisation.R in scripts). This data sheet also includes additional information on the parent's identity and age at the time of mating, alongisde flags indicating which trait analyses individual offspring could be included within.

#### Key:

-   Column A: **F1_ID** - ID- The unique identifier assigned to each individual offspring.

-   Column B: **Include_in_total** - Y/N - Flag dictates whether the offspring had complete lifespan data recorded.
-   All offspring included had accurately recorded death dates.
  
-   Column C: **Include_in_adult** - Y/N - Dictates whether an individual offspring had accurately recorded adult emergence dates.
-   If **N**, then the offspring could not be included in the analysis of any adult traits (i.e., adult mass, adult lifespan, development time, reproductive success)

-   Column D: **Mother_ID** - ID - The unique identifier assigned to each individual mother. 

-   Column E: **Father_ID** - ID - The unique identifier given to each individual father.

-   Column F: **PairID** - ID - The combined unique identifier given to a pair of parents.
-   All parents were mated with the same partner for the full duration of the study.

-   Column G: **Temp** - Category - The temperature treatment that parent pairs were moved to following their first mating attempt.
-   Parents' were assigned to be maintained under one of the three following temperature treatments for the full duration of the study: 25.5°C, 28.0°C, or 30.5°C.

-   Column H: **Timepoint** - numeric - The numeric identifer given to each mating attempt.
-   Parents had up to 8 mating attempts from which we collected and housed offspring to be included in the study.

-   Column I: **age_at_mating** - numeric (weeks) - The exact adult age of each individual offspring's mother at their time of conception (i.e., laying) (recorded in weeks). 

-   Column J: **avg.age** - numeric (weeks) - The average age of the mother across all mating attempts resulting in succesfully hatched eggs (in weeks).
-   A fixed measure of the mother's average age of reproduction that is shared and stable across all offspring in a family. 

-   Column K: **within_subject_age** - numeric (weeks) - the mother's Δage at reproduction when the offspring was conceived, reflecting the parents actual age at mating (relative to their mean age).
-   This value was calculated by taking the following: **age_at_mating - avg.age**.
-   This value varies across offspring within individual families depending on the age of the mother at the time of laying.

-   Column L: **cum_successful_matings** - numeric - The number of succesful mating attempts (i.e., those matings that resulted in successfully hatched eggs) that the parent had undertaken when the individual offspring was conceived. This variable was included in an attempt to disentangle accumulated reproductive effort from Δage.

-   Column M: **F1_hatch** - date (yyyy-mm-dd) - The date that the individual offspring hatched.

-   Column N: **F1_Adult_emergence** - date(yyyy-mm-dd) - the date that the individual offspring reached adulthood.
-   Given a value of **NA** if the offspring failed to reach adulthood.
-   Eight offspring reached adulthood but had highly uncertain adult emergence dates, so were given an N in the Include_in_adult flag column, and excluded from the analysis of adult traits. 

-   Column O: **F1_sex** - M/F - The sex of the individual offspring. Either Male (M), Female (F), or NA if the offspring died before adulthood.

-   Column P: **F1_death_date** - date(yyyy-mm-dd)- The date when the individual offspring was observed dead. 

-   Column Q: **F1_adultmass** - grams - The weight of the offspring (in grams), measured at least 24 hours after we first observed adult emegrence. 
-   Given a value of NA if the offspring died before reaching adulthood.

-   Column R: **adult_surv** - 1/0 - The flag indicating if an individual offspring survived to adulthood.
-   1 = the offspring succesfully reached adulthood. 0 = the offspring died before adulthood.

-   Column S: **total_lifespan** - numeric (in weeks) - The total lifespan of the offspring (i.e., the difference between the death date and hatch date).

-   Column T: **development_time** - numeric (in weeks) - The time taken for the offspring to reach adulthood (i.e., the difference between the hatch date and adult emergence date).
-   Given a value of NA if the offspring did not reach adulthood.

-   Column U: **adult_lifespan** - numeric (in weeks) - The total time spent as an adult/total time spent having reached reproductive maturity (i.e, the difference between the adult emergence date and the death date).
-   Given a value of NA if the offspring did not reach adulthood. 

-   Column V: **Mother_bodymass** - numeric (in grams) - The adult weight of the offspring's mother (measured 24 hours after the mothers own adult emergence date).
-   This variable was included in the analysis of offspring adult mass.

-   Column W: **Father_bodymass** - numeric (in grams) - The adult weight of the offspring's father (as measured 24 hours after the fathers own adult emergence date).

-   Column X: **Maternal_grandfather_ID** - ID - The unique identifier given to the offsprings maternal grandfather.

-   Column Y: **Maternal_grandmother_ID** - ID - The unique identifier given to the offsprings maternal grandmother.

-   Column Z: **Paternal_grandmother_ID** - ID - The unique identifier given to the offsprings paternal grandfather.

-   Column BA: **Paternal_grandfather_ID** - ID - The unique identifier given to the offsprings paternal grandmother.

-   Column BB: **mothers_age_at_entry** - numeric  (in weeks) - The exact age of the mother at the **first** mating attempt (i.e., timepoint 1).

-   Column BC: **succesful_mate_count** - numeric- The total number of matings that the offspring's parents succesfully undertook over their full life course (i.e., the total number of matings resulting in succesfully hatched eggs).
