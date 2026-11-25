#!/bin/bash
#SBATCH -A YOUR_ANVIL_ALLOCATION   # your ACCESS allocation on Anvil
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 02:00:00
#SBATCH -J starem_prep
#SBATCH -o logs/prep_%j.out
## Panel, descriptive tables, choice of K, true policy values (5000 people per policy),
## learned missingness model -> $SCRATCH/starem/prep.rds. Run once; everything else needs it.
module load modtree/cpu
module load r/4.4.1
cd ~/starem
export STAR_OUT=$SCRATCH/starem STAR_MODE=prep STAR_REPS_BIG=300
Rscript starem_anvil.R
