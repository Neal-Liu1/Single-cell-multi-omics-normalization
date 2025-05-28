source("/home/users/allstaff/liu.ne/scMultiOmics-normalization/transformGamPoi_Paper_benchmark/benchmark/src/transformations/transformation_helper.R")

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
message(Sys.time(), " - Starting task for param_grid row: ", i)

param_grid <- readRDS("/home/users/allstaff/liu.ne/scMultiOmics-normalization/slurm_benchmarking_scripts/downsampling_params.rds")
row <- param_grid[i, ]
data <- row$data
trans <- row$trans
mode <- row$mode

output_dir <- if (mode == "full") {
  "/vast/scratch/users/liu.ne/transformGamPoi_Output/results/downsampling_results/full_knns"
} else {
  "/vast/scratch/users/liu.ne/transformGamPoi_Output/results/downsampling_results/reduced_knns"
}
data_dir <- "/vast/scratch/users/liu.ne/transformGamPoi_downsampling_data"

message(Sys.time(), " - Data: ", data, " | Transformation: ", trans, " | Mode: ", mode)

knns <- c(10, 20, 50)
pcs <- c(5, 10, 50)

alphas <- list(
  logp1 = FALSE,
  logp1_hvg = FALSE,
  logp1_zscore = FALSE,
  logp1_hvg_zscore = FALSE,
  logp_cpm = FALSE,
  logp1_size_normed = FALSE,
  acosh = TRUE,
  logp_alpha = TRUE,
  pearson = TRUE,
  pearson_clip = TRUE,
  pearson_clip_hvg = TRUE,
  pearson_clip_zscore = TRUE,
  pearson_clip_hvg_zscore = TRUE,
  pearson_analytic = TRUE,
  sctransform = TRUE,
  rand_quantile = TRUE,
  normalisr_norm = FALSE,
  sanity_map = FALSE,
  sanity_dists = FALSE,
  glmpca = FALSE,
  raw_counts = FALSE,
  scaled_raw_counts = FALSE
)


current_alphas <- alphas[[trans]]
if (!is.vector(current_alphas)) {
  current_alphas <- list(current_alphas)
}

######### Start Transformation #######

UMI <- readRDS(file.path(data_dir, data))
UMI <- if (mode == "full") UMI$full else UMI$reduced

message(Sys.time(), " - Loaded data. Filtering expressed cells and genes...")
expressed_cells <- matrixStats::colSums2(UMI) > 0
expressed_genes <- matrixStats::rowSums2(UMI) > 0
UMI <- UMI[expressed_genes, expressed_cells]

sf <- MatrixGenerics::colSums2(UMI)
sf <- sf / mean(sf)

for (alpha in current_alphas) {
  message(Sys.time(), " - Applying transformation: ", trans, " | alpha: ", paste(alpha, collapse = ", "))
  trans_dat <- all_transformations[[trans]](UMI, sf, alpha)
  
  for (pc in pcs) {
    for (nn in knns) {
      message(Sys.time(), " - Computing KNN | pc: ", pc, " | knn: ", nn)
      KNN <- make_knn_graph(trans, trans_dat, pc, nn)
      
      result_id <- paste0(data, "_", trans, "_alpha:", paste(alpha, collapse = ","), "_pc:", pc, "_nn:", nn)
      saveRDS(KNN, file.path(output_dir, result_id))
      message(Sys.time(), " - Saved result to: ", file.path(output_dir, result_id))
    }
  }
}

message(Sys.time(), " - Task complete for row: ", i)





