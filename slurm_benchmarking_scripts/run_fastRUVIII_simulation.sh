#!/bin/bash
#SBATCH --job-name=fastRUV_sim
#SBATCH --output=/vast/scratch/users/liu.ne/transformGamPoi_Output/results/simulation_results/logs/out/%A_%a.out
#SBATCH --error=/vast/scratch/users/liu.ne/transformGamPoi_Output/results/simulation_results/logs/errors/%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=012:00:00
#SBATCH --array=1-75

module load R/4.4.1
Rscript run_fastRUVIII_simulation.R $SLURM_ARRAY_TASK_ID