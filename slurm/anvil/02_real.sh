#!/bin/bash
#SBATCH -A YOUR_ANVIL_ALLOCATION   # your ACCESS allocation on Anvil
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -J starem_real
#SBATCH -o logs/real_%j.out
## HeartSteps analysis on the fine delta grid (13 x 2 tilts). Needs prep.rds.
module load modtree/cpu
module load r/4.4.1
cd ~/starem
export STAR_OUT=$SCRATCH/starem STAR_MODE=real STAR_REPS_BIG=300
Rscript starem_anvil.R
