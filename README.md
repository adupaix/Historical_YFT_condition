# Drifting Fish Aggregating Devices (DFADs) as ecological traps: no evidence displayed by a long-term analysis of yellowfin tuna condition.

[![License](https://img.shields.io/github/license/adupaix/Historical_YFT_condition)](https://github.com/adupaix/Historical_YFT_condition/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/416344484.svg)](https://zenodo.org/badge/latestdoi/416344484)
[![Latest Release](https://img.shields.io/github/release/adupaix/Historical_YFT_condition)](https://github.com/adupaix/Historical_YFT_condition/releases)

This repository contains all the scripts to run the analysis performed in the manuscript in revision as a Note to *Marine Ecology Progress Series (MEPS)*:

__Dupaix A., Dagorn L., Deneubourg J-L., Duparc A., Guillou A., Capello M.__ (in revision). No evidence from long-term analysis of yellowfin tuna condition that Drifting Fish Aggregating Devices act as ecological traps.

Pre-print available here: [https://hal.archives-ouvertes.fr/hal-03690665](https://hal.archives-ouvertes.fr/hal-03690665)

Please contact me if you have any question (amael.dupaix@ird.fr).

# Run

The scripts run with [R version 4.0.3](https://www.r-project.org/) (R Core Team, 2020).

The [conda](https://docs.conda.io/projects/conda/en/latest/) environment to run the model is provided. To create, type : `conda env create -f r-tuna-condition.yml`

Templates of the config files are provided: `cfg/job_template.pbs` (to run the bootstrap on a cluster) and `cfg/cfg_template.R`. Config files used to obtain the results of the study are available in `cfg/paper_cfgs.zip`.

Outputs obtained with `.pbs` files are then used to run `main_post_cluster.R`

# Data availability

Data used in the study is available in Guillou et al. (2021).

# References

Guillou Aurelie, Bodin Nathalie, Chassot Emmanuel, Duparc Antoine, Fily Theotime, Sabarros Philippe, Depetris Mathieu, Amande Monin Justin, Lucas Juliette, Diaha Constance, Floch Laurent, Barde Julien, Pascual Alayon Pedro J, Baez José Carlos, Cauquil Pascal, Briand Karine, Lebranchu Julien (2022). Tunabio: biological traits of tropical tuna and bycatch species caught by purse seine fisheries in the Western Indian and Eastern Central Atlantic Oceans. SEANOE. https://doi.org/10.17882/73500

R Core Team. (2020). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing. https://www.R-project.org/
