#!/bin/bash
#SBATCH --job-name=exact_adj_sim
#SBATCH --array=1-40
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=00:25:00
#SBATCH -p standard

# The scenario ID (1, 2, or 3) is passed as the first argument
SCENARIO=$1

if [ -z "$SCENARIO" ]; then
  echo "Error: Provide a scenario ID (1, 2, or 3)."
  exit 1
fi

module load goolf R
Rscript "s41_prop_times.R" $SCENARIO