#!/bin/bash
#SBATCH --array=1-40
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=02:00:00
#SBATCH -p standard

module load goolf R

Rscript "master_sim.R" "$SETTING" "$SLURM_ARRAY_TASK_ID"