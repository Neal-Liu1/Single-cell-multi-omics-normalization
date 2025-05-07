library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
file_path <- args[1]

cat("Processing:", file_path, "\n")

KNNs <- readRDS(file_path)
stopifnot(all(dim(KNNs[[1]]) == dim(KNNs[[2]])))
n_cells <- nrow(KNNs[[1]])


cons <- mean(sapply(seq_len(n_cells), function(cell_idx){
  length(intersect(KNNs[[1]][cell_idx,], KNNs[[2]][cell_idx,]))
}))

res <- tibble(
  mean_overlap = cons, 
  transformation_id = pa$input_id, 
  dataset = pa$dataset, 
  seed = pa$seed,
  pca_dim = pa$pca_dim, 
  knn = pa$knn, 
  transformation = pa$transformation, 
  alpha = pa$alpha)

write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))