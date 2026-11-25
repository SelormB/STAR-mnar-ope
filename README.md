# STAR: Sensitivity analysis for off-policy evaluation when mobile health rewards are missing not at random

Code, HPC job scripts, results and paper for **STAR**, a sensitivity-analysis framework for evaluating
mobile health decision rules (off-policy evaluation) from micro-randomized trials when the wearable
outcome is missing intermittently, in spells, and possibly not at random.

STAR has two estimators that share one sensitivity parameter, an exponential tilt
δ(A) = δ₀ + δ_A·A:

- **STAR-EM**: a joint model of the reward and its missingness (with an optional hidden activity state),
  fitted by a semiparametric EM algorithm, followed by multiple imputation and per-decision importance sampling.
- **STAR-DR**: a tilted doubly robust estimator (tilted pseudo-rewards plus fitted-Q evaluation).

The methods are tested in a plasmode simulation built from HeartSteps V1 (8 scenarios, 7,300 simulated
trials) and applied to HeartSteps V1 itself. The full write-up is in [`paper/STAR_paper.pdf`](paper/STAR_paper.pdf).

## Main findings

| Finding | Evidence |
| --- | --- |
| Ignoring MNAR biases policy values | MAR-DR off by −58 steps/day under reward-dependent MNAR; +130 to +250 when missingness depends on the prompt |
| Both STAR estimators work when δ is right | Within 30 steps of the truth in scenarios 1–6; STAR-EM coverage 0.93–0.96 |
| STAR-DR is the safer default | Only correction near unbiased under hidden activity states (+6, +9); ~8× cheaper than STAR-EM |
| The hidden-state (K = 2) extension did not help | Same bias as K = 1; hit the 20-iteration cap in 99–100% of fits under a nonzero tilt |
| HeartSteps policies cannot be told apart | No contrast differs from zero under MAR or anywhere on δ₀ ∈ [−1.5, 1.5] |

## Repository layout

```
R/starem_anvil.R          all methods; runs in four modes set by STAR_MODE (prep, real, sim, combine)
slurm/anvil/              SLURM scripts for Anvil (Purdue RCAC), the main run
  00_setup.sh             install R packages, download HeartSteps V1 (login node)
  01_prep.sh              data, true policy values, generator checks -> prep.rds, tasks.csv
  02_real.sh              HeartSteps analysis over the full tilt grid
  03_pilot.sh             timing pilot
  04_sim_small.sh         scenarios 1-7 (tasks 1-280, 25 replications each)
  05_sim_big.sh           scenario 8, N = 150 (tasks 281-430, 2 replications each)
  06_combine.sh           all simulation tables and figures
slurm/alpine/             SLURM scripts for Alpine (CU Boulder, via RMACC)
  00_alpine_setup.sh      conda R 4.4 environment
  07_alpine_hidden.sh     scenario 8 + scenario-7 tasks 274 and 277 (run on Alpine)
  08_alpine_combine.sh    final combine over all 430 task files
results/simulation/       final tables (CSV) and figures from the combine step
results/heartsteps/       real-data tables (CSV) and figures
paper/                    LaTeX source (Overleaf-ready), figures, compiled PDF
```

## Data

HeartSteps V1 is public: <https://github.com/klasnja/HeartStepsV1>. The setup script downloads
`suggestions.csv` from that repository; no data are stored here.

## Reproducing the results

The pipeline was run on two clusters. Each SLURM array task writes one `sim_XXX.rds`; resubmitting a script
skips finished tasks, and every replication has its own seed derived from `set.seed(100)` and its task number,
so results do not depend on which cluster ran a task.

**Anvil** (set your allocation in each script's `#SBATCH -A` line):

```bash
bash slurm/anvil/00_setup.sh          # login node
sbatch slurm/anvil/01_prep.sh
sbatch slurm/anvil/02_real.sh
sbatch slurm/anvil/03_pilot.sh        # check $SCRATCH/starem/pilot_timing.txt
sbatch slurm/anvil/04_sim_small.sh
sbatch slurm/anvil/05_sim_big.sh
sbatch slurm/anvil/06_combine.sh
```

**Alpine** (used here for scenario 8 and two unfinished scenario-7 tasks when the Anvil allocation ran low):

```bash
acompile --ntasks=2 --time=01:00:00
bash slurm/alpine/00_alpine_setup.sh  # then copy prep.rds and suggestions.csv to /scratch/alpine/$USER/starem
exit
module load slurm/alpine
sbatch slurm/alpine/07_alpine_hidden.sh
# copy the Anvil sim_*.rds files into the same folder, then:
sbatch slurm/alpine/08_alpine_combine.sh
```

Run notes: scenario 7 tasks (25 replications with STAR-EM at K = 1 and K = 2) took 5–7 hours each on Anvil,
longer than the pilot projected; timed-out tasks were resubmitted with a longer `--time` on the command line.
A scenario-8 task (2 replications, N = 150) took about 57 minutes on Alpine. The final combine read all
430 task files: 1,000 replications for scenarios 1–7 and 300 for scenario 8.

## Software

R 4.4 with data.table, ranger, rpart, FNN and ggplot2.

## Citation

Buaka, S. (2026). *STAR: A sensitivity-analysis framework for off-policy evaluation when mobile health
rewards are missing not at random.* University of Northern Colorado. See `CITATION.cff`.

## Author

Selorm Buaka, PhD Student, Applied Statistics and Research Methods, University of Northern Colorado ·
ORCID [0009-0005-3991-4859](https://orcid.org/0009-0005-3991-4859)
