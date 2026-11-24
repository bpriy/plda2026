#!/bin/bash

# Define all 12 settings for the 5-step implementation
SETTINGS=(
  "Step1_Null" 
  "Step2_PosCorr" 
  "Step3_K5" "Step3_K10" "Step3_K20" 
  "Step4_E0.1" "Step4_E0.5" "Step4_E1.0" 
  "Step5A_Probit" "Step5B_Cloglog" "Step5C_Interact"
)

# Submit a 40-task array for each setting concurrently
for SETTING in "${SETTINGS[@]}"; do
  echo "Submitting job array for $SETTING..."
  sbatch --export=ALL,SETTING=$SETTING run_sim.slurm
done