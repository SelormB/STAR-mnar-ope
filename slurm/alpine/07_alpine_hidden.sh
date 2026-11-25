#!/bin/bash
#SBATCH --account=rmacc-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=14:00:00
#SBATCH --job-name=starem_hidden
#SBATCH --output=/projects/%u/starem/logs/hidden_%A_%a.out
#SBATCH --array=274,277,281-430%100
## STAR-EM hidden-activity-state scenarios on Alpine (no SU budget; fairshare only)
##   2 scenario-7 tasks (274, 277; N = 37, 25 reps each) that did not fit in the Anvil budget
##   tasks 281-430: scenario 8, N = 150, 300 replications (2 per task; peak memory ~4.5 GB)
## Same script and per-replication seeds as on Anvil; finished tasks are skipped automatically.
module purge
module load miniforge
conda activate /projects/$USER/starem/envs/r44
export STAR_OUT=/scratch/alpine/$USER/starem STAR_MODE=sim STAR_REPS_BIG=300
Rscript /projects/$USER/starem/starem_anvil.R
