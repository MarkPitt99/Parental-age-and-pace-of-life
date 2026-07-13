---
# **R scripts for the paper:** "*Parental age accelerates offspring pace-of-life*"

This folder contains the R scripts used to analyse the data collected as part of a project investigating parental age effects 
on offspring life-history traits.

------------------------------------------------------------------------

Mark D. Pitt, Brendan O’Connor, Timothy D. Sheen, Davide M. Dominoni, Tom Tregenza, Jelle J. Boonekamp **Parental age at reproduction accelerates offspring pace of life**. 
*bioRxiv*. [DOI:10.64898/2026.04.14.718189](https://www.biorxiv.org/content/10.64898/2026.04.14.718189v2.full)

------------------------------------------------------------------------

For any further information, please contact: Mark Pitt, email: [markdavidpitt\@gmail.com](mailto:markdavidpitt@gmail.com)

## Description of scripts:

Scripts are ordered numerically by their order of appearance in the main manuscript. We used two different Bayesian modelling approaches throughout this project. First, we assessed how increasing parental age impacted the mean and variance of the offspring trait values using *brms* models (*brms v.2.21.0*), following the location-scale approach to trait modelling outlined by: 

Nakagawa, S., Ortega, S., Gazzea, E., Lagisz, M., Lenz, A., Lundgren, E. and Mizuno, A., 2026. **Location–scale models in ecology and evolution: Heteroscedasticity in continuous, count and proportion data**. Methods in Ecology and Evolution, 17(2), pp.554-566. Available at: [https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210x.70203]

Meanwhile we assesed how the mortality parameters (i.e., the Makeham term (*c*), shape parameter (*b0*), and scale parameter (*b1*)) were affected by parental age, across both the offspring's total and adult lifespan, using Bayesian Survival Trajectory Analysis (*BaSTA*) package (v.2.0.2).
