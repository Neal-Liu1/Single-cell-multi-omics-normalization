#!/bin/bash
#SBATCH --job-name=consistency_benchmark
#SBATCH --output=/vast/scratch/users/liu.ne/transformGamPoi_Output/results/consistency_results/logs/out/%A_%a.out
#SBATCH --error=/vast/scratch/users/liu.ne/transformGamPoi_Output/results/consistency_results/logs/errors/%A_%a.err
#SBATCH --array=179-187
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --time=07:00:00

module load R/4.4.1
Rscript consistency_benchmark.R $SLURM_ARRAY_TASK_ID
