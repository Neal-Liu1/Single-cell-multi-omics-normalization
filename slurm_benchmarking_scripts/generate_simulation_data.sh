#!/bin/bash
#SBATCH --job-name=generate_simulation_data
#SBATCH --output=/vast/scratch/users/liu.ne/transformGamPoi_simulation_data/logs/output/%A_%a.out
#SBATCH --error=/vast/scratch/users/liu.ne/transformGamPoi_simulation_data/logs/errors/%A_%a.err
#SBATCH --array=1-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --time=01:00:00

module load R/4.4.1
Rscript generate_simulation_data.R $SLURM_ARRAY_TASK_ID