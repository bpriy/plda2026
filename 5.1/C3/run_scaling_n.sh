#!/bin/bash
#SBATCH --job-name=scale_N
#SBATCH --output=logs/scale_N_%A_%a.out
#SBATCH --error=logs/scale_N_%A_%a.err
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=02:00:00
#SBATCH -p standard

echo "Running Large-N Scaling Task: ${SLURM_ARRAY_TASK_ID}"

# Load necessary modules
module load goolf R

# Execute the R script
Rscript scaling_n_sim.R