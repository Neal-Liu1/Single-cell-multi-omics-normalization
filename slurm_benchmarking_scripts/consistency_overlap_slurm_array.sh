#!/bin/bash
#SBATCH --job-name=consistency_knn_overlap
#SBATCH --output=/vast/scratch/users/liu.ne/transformGamPoi_Output/results/consistency_results/logs/out/%A_%a.out
#SBATCH --error=/vast/scratch/users/liu.ne/transformGamPoi_Output/results/consistency_results/logs/errors/%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --array=1-1000

# Positional arguments
FILE_LIST=$1   # path to file_list.txt
OFFSET=$2      # offset to add to array index

# Compute actual file index
INDEX=$((SLURM_ARRAY_TASK_ID + OFFSET))

# Extract file from list
FILE=$(sed -n "${INDEX}p" "$FILE_LIST")

# Run your R script on the selected file
module load R/4.4.1
Rscript calculate_consistency_overlap.R "$FILE"