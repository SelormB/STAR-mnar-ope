#!/bin/bash
## Alpine (CU Boulder, through RMACC): one-time setup of the R environment.
## Run inside a compile session:   acompile --ntasks=2 --time=01:00:00
## then:                           bash /projects/$USER/starem/00_alpine_setup.sh
module purge
module load miniforge
mkdir -p /projects/$USER/starem/logs /scratch/alpine/$USER/starem
conda create -y --prefix /projects/$USER/starem/envs/r44 -c conda-forge \
  r-base=4.4 r-data.table r-ranger r-fnn r-ggplot2 r-rpart
conda activate /projects/$USER/starem/envs/r44
Rscript -e 'library(data.table); library(ranger); library(FNN); library(ggplot2); library(rpart); cat("ok\n")'
## prep.rds and suggestions.csv must be in /scratch/alpine/$USER/starem before running 07 or 08:
## copy them from the Anvil run (prep step), or run the prep step on Alpine.
