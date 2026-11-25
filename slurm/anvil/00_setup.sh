#!/bin/bash
## Run ONCE on a LOGIN node (needs internet):   bash slurm/00_setup.sh
## Installs any missing R packages into your personal R 4.4 library and downloads HeartSteps V1.
module load modtree/cpu
module load r/4.4.1
mkdir -p $SCRATCH/starem logs
Rscript -e 'lib <- Sys.getenv("R_LIBS_USER"); dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  need <- setdiff(c("data.table", "ranger", "FNN", "ggplot2"), rownames(installed.packages()))
  if (length(need)) install.packages(need, lib = lib, repos = "https://cloud.r-project.org")
  for (p in c("data.table", "ranger", "FNN", "ggplot2", "rpart")) cat(p, as.character(packageVersion(p)), "\n")'
[ -f $SCRATCH/starem/suggestions.csv ] || wget -q -O $SCRATCH/starem/suggestions.csv \
  https://raw.githubusercontent.com/klasnja/HeartStepsV1/main/data_files/suggestions.csv
ls -l $SCRATCH/starem/suggestions.csv
