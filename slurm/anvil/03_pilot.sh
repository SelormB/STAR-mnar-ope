#!/bin/bash
#SBATCH -A YOUR_ANVIL_ALLOCATION   # your ACCESS allocation on Anvil
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH -J starem_pilot
#SBATCH -o logs/pilot_%j.out
## Timing pilot: 2 replications of the heaviest N = 37 scenario and of the N = 150 scenario.
## Writes $SCRATCH/starem/pilot_timing.txt with the projected SU and the time limits to use.
module load modtree/cpu
module load r/4.4.1
cd ~/starem
export STAR_OUT=$SCRATCH/starem STAR_MODE=sim STAR_PILOT=1
Rscript starem_anvil.R
