#!/bin/bash
#SBATCH --job-name=survey_sim
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --array=1-50            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1          
#SBATCH --time=01:30:00            
#SBATCH --mem=10G  
#SBATCH --partition=standard

module load goolf R

Rscript run_survey.R