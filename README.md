# Long-term variations of yellowfin tuna physiological condition linked with FAD use

This repository contains all the scripts to run the analysis performed in our manuscript. Please contact me if you have any question.

Mail: amael.dupaix@ens-lyon.fr

# Run

The scripts run with [R version 4.0.3](https://www.r-project.org/)

The [conda](https://docs.conda.io/projects/conda/en/latest/) environment to run the model is provided. To create, type : `conda env create -f r-tuna-condition.yml`

Templates of the config files are provided: `cfg/job_template.pbs` (to run the bootstrap on a cluster) and `cfg/cfg_template.R`

Outputs obtained with `.pbs` files are then used to run `main_post_cluster`
