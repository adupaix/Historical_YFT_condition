# Fish Aggregating Devices (FADs) as ecological traps: no evidence displayed by a long-term analysis of yellowfin tuna condition

This repository contains all the scripts to run the analysis performed in the manuscript: __Dupaix A., Dagorn L., Deneubourg J-L., Duparc A., Guillou A., Capello M. (in press). Fish Aggregating Devices (FADs) as ecological traps: no evidence displayed by a long-term analysis of yellowfin tuna condition.__

Please contact me if you have any question. 

Mail: amael.dupaix@ird.fr

# Run

The scripts run with [R version 4.0.3](https://www.r-project.org/) (R Core Team, 2020).

The [conda](https://docs.conda.io/projects/conda/en/latest/) environment to run the model is provided. To create, type : `conda env create -f r-tuna-condition.yml`

Templates of the config files are provided: `cfg/job_template.pbs` (to run the bootstrap on a cluster) and `cfg/cfg_template.R`

Outputs obtained with `.pbs` files are then used to run `main_post_cluster.R`

# Data availability

Data used in the study are available upon request to the IRD’s [Ob7](https://www.ob7.ird.fr/pages/datacall.html)–“Observatoire des Ecosystèmes Pélagiques Tropicaux Exploités”.

# References

R Core Team. (2020). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing. https://www.R-project.org/
