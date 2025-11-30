#!/bin/bash
#SBATCH --output=simulation.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --job-name=YS_PEAL2

module load R/4.3
R CMD BATCH Visual2.R