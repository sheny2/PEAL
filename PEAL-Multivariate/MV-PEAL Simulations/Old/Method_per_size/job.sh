#!/bin/bash
#SBATCH --output=simulation_%A_%a.out
#SBATCH --cpus-per-task=45
#SBATCH --mem=150G
#SBATCH --job-name=Method_Size
#SBATCH --array=1-8

module load R/4.3

# Run the R script; separate Rout per array index
R CMD BATCH --no-save --no-restore LargeSim.R LargeSim_${SLURM_ARRAY_TASK_ID}.Rout
