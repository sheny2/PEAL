#!/bin/bash
#SBATCH --output=simulation.out
#SBATCH --cpus-per-task=60
#SBATCH --mem=150G
#SBATCH --job-name=YS_PEAL_DiffRho

module load R/4.3
R CMD BATCH Simulation-RhoDiff-Full.R 