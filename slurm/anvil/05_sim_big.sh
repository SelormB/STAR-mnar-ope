#!/bin/bash
#SBATCH -A YOUR_ANVIL_ALLOCATION   # your ACCESS allocation on Anvil
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 04:00:00
#SBATCH -J starem_big
#SBATCH -o logs/big_%A_%a.out
#SBATCH --array=281-430%50
## The hidden-state scenario with N = 150 participants: 300 replications, 2 per task -> tasks 281-430
## (change STAR_REPS_BIG below AND the --array range together; see pilot_timing.txt).
## 2 cores per task (more memory for the larger data). Submit after checking the balance.
module load modtree/cpu
module load r/4.4.1
cd ~/starem
export STAR_OUT=$SCRATCH/starem STAR_MODE=sim STAR_REPS_BIG=300
Rscript starem_anvil.R
