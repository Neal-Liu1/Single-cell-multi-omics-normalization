#!/bin/bash
#SBATCH --job-name=transform_simulated_data
#SBATCH --output=/vast/scratch/users/liu.ne/transformGamPoi_Output/results/simulation_results/logs/out/%A_%a.out
#SBATCH --error=/vast/scratch/users/liu.ne/transformGamPoi_Output/results/simulation_results/logs/errors/%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --array=1-475

module load R/4.4.1
Rscript transform_simulated_data.R $SLURM_ARRAY_TASK_ID