#!/bin/bash
#SBATCH --account=rmacc-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --job-name=starem_combine
#SBATCH --output=/projects/%u/starem/logs/combine_%j.out
## All 8 scenarios: tables and figures from every sim_XXX.rds in the scratch folder.
module purge
module load miniforge
conda activate /projects/$USER/starem/envs/r44
export STAR_OUT=/scratch/alpine/$USER/starem STAR_MODE=combine STAR_REPS_BIG=300
Rscript /projects/$USER/starem/starem_anvil.R
