#!/bin/bash
#SBATCH --job-name=MV_PEAL_Sim
#SBATCH --output=logs/sim_%a.out
#SBATCH --error=logs/sim_%a.err
#SBATCH --array=1-80
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G

# Create dirs
mkdir -p results
mkdir -p logs

# Load R module 
module load R/4.3

echo "------------------------------------------------"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Date: $(date)"
echo "------------------------------------------------"

# Run the worker script
# The R script will automatically detect the 5 CPUs using 'parallelly'
Rscript run_sim_worker.R $SLURM_ARRAY_TASK_ID