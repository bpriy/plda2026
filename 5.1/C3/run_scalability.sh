#!/bin/bash
#SBATCH --job-name=scale_sim
#SBATCH --output=logs/scale_%A_%a.out
#SBATCH --error=logs/scale_%A_%a.err
#SBATCH --array=1-300
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH -p standard

echo "Running Scalability Simulation Task: ${SLURM_ARRAY_TASK_ID}"

# Load necessary modules 
module load goolf R

# Execute the R script
Rscript scalability_sim.R