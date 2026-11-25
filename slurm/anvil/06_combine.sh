#!/bin/bash
#SBATCH -A YOUR_ANVIL_ALLOCATION   # your ACCESS allocation on Anvil
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 01:00:00
#SBATCH -J starem_combine
#SBATCH -o logs/combine_%j.out
## All simulation tables and figures from the finished tasks (can be run any time; it reports
## how many replications each scenario has so far).
module load modtree/cpu
module load r/4.4.1
cd ~/starem
export STAR_OUT=$SCRATCH/starem STAR_MODE=combine STAR_REPS_BIG=300
Rscript starem_anvil.R
