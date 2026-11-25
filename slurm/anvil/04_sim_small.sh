#!/bin/bash
#SBATCH -A YOUR_ANVIL_ALLOCATION   # your ACCESS allocation on Anvil
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 03:00:00
#SBATCH -J starem_sim
#SBATCH -o logs/sim_%A_%a.out
#SBATCH --array=1-280%100   # scenarios 1-7; reruns of timed-out tasks used --array=<ids> --time=12:00:00 on the command line
## The 7 scenarios with N = 37: 1000 replications each, 25 per task -> tasks 1-280.
## A finished task writes sim_XXX.rds; resubmitting skips finished tasks automatically.
module load modtree/cpu
module load r/4.4.1
cd ~/starem
export STAR_OUT=$SCRATCH/starem STAR_MODE=sim STAR_REPS_BIG=300
Rscript starem_anvil.R
