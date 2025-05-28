#!/bin/bash
#SBATCH --job-name=transform_simulated_data
#SBATCH --output=/vast/scratch/users/liu.ne/transformGamPoi_Output/results/downsampling_results/logs/out/%A_%a.out
#SBATCH --error=/vast/scratch/users/liu.ne/transformGamPoi_Output/results/downsampling_results/logs/errors/%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --time=09:00:00
#SBATCH --array=1-190

module load R/4.4.1
Rscript transform_downsampling_data.R $SLURM_ARRAY_TASK_ID