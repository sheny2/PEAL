#!/bin/bash
#SBATCH --output=logs/simulation.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=PEAL_Plot

module load R/4.3
R CMD BATCH run_sim_analysis.R