#!/bin/bash
#
#SBATCH --job-name=sim_study
#SBATCH --output="./logs/sim_%A_%a.txt"
#
#SBATCH --ntasks=1
#SBATCH --time=25:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir="./"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=l.berchtold@umcutrecht.nl
#
# Array size = nrow(scenario_grid) * nsim
# Get the number before submitting:
#   Rscript -e "source('config/input_baseline.R'); cat(nrow(scenario_grid) * nsim, '\n')"
#SBATCH --array=1-1100

mkdir -p logs

CONFIG="config/input_mean_shift_variation.R"

Rscript "simulation/simulation_execution_hpc.R" $CONFIG