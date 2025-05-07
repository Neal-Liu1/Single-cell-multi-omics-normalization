library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
result_id <- args[1]
file_path <- paste0("/vast/scratch/users/liu.ne/transformGamPoi_Output/results/consistency_results/transformed_knns/",args[1])

cat("Processing:", file_path, "\n")

KNNs <- readRDS(file_path)
stopifnot(all(dim(KNNs[[1]]) == dim(KNNs[[2]])))
n_cells <- nrow(KNNs[[1]])


cons <- mean(sapply(seq_len(n_cells), function(cell_idx){
  length(intersect(KNNs[[1]][cell_idx,], KNNs[[2]][cell_idx,]))
}))

transformations <- c(
  "logp1", 
  "logp1_hvg", 
  "logp1_zscore", 
  "logp1_hvg_zscore",
  "logp_cpm", 
  "logp1_size_normed", 
  "acosh", 
  "logp_alpha",
  "pearson", 
  "pearson_clip", 
  "pearson_analytic", 
  "sctransform",
  "rand_quantile", 
  "pearson_clip_hvg", 
  "pearson_clip_zscore",
  "pearson_clip_hvg_zscore", 
  "normalisr_norm", 
  "raw_counts", 
  "scaled_raw_counts"
)

transformations_sorted <- transformations[order(nchar(transformations), decreasing = TRUE)]
transformation <- transformations_sorted[which.max(str_detect(result_id, fixed(transformations_sorted)))]

res <- tibble(
  mean_overlap = cons, 
  transformation_id = result_id, 
  dataset = str_extract(result_id, "^GSE\\d+"), 
  seed = 1,
  pca_dim = str_extract(result_id, "(?<=pc:)\\d+"), 
  knn = str_extract(result_id, "(?<=nn:)\\d+"), 
  transformation = transformation, 
  alpha = str_extract(result_id, "(?<=alpha:)\\d*\\.?\\d+"))

write_tsv(res, paste0("/vast/scratch/users/liu.ne/transformGamPoi_Output/results/consistency_results/knn_overlap_metrics/", result_id))

