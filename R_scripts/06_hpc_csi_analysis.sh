#!/bin/bash

#SBATCH --job-name=csi_analysis
#SBATCH --array=1-21
#SBATCH --cpus-per-task=4
#SBATCH --mem=18gb
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL

# Load in modules
module load r/4.4.0

# Run the task for each index in the job array
Rscript HPC/scripts/06_hpc_csi_analysis.R ${SLURM_ARRAY_TASK_ID}
