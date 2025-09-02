source(
  "/home/users/allstaff/liu.ne/scMultiOmics-normalization/transformGamPoi_Paper_benchmark/benchmark/src/transformations/transformation_helper.R"
)

output_dir <- "/vast/scratch/users/liu.ne/transformGamPoi_Output/results/simulation_results/transformed_knns"

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
cat("Starting task for param_grid row:", i, "\n")

param_grid <- readRDS("/home/users/allstaff/liu.ne/scMultiOmics-normalization/slurm_benchmarking_scripts/simulation_params.rds")

row <- param_grid[i, ]
data <- row$data
path <- row$path
trans <- row$trans

cat("Data:", data, "\n")

knns <- c(10, 20, 50)
pcs <- c(5, 10, 50)

alphas <- list(
  logp1 = FALSE,
  logp1_hvg = FALSE,
  logp1_zscore = FALSE,
  logp1_hvg_zscore = FALSE,
  logp_cpm = FALSE,
  logp1_size_normed = FALSE,
  acosh = c(0.05, 0, TRUE),
  logp_alpha = c(0.05, 0, TRUE),
  pearson = c(0.05, 0, TRUE),
  pearson_clip = c(0.05, 0, TRUE),
  pearson_analytic = c(0.05, 0, TRUE),
  sctransform = c(0.05, 0, TRUE),
  rand_quantile = c(0.05, 0, TRUE),
  pearson_clip_hvg = c(0.05, 0, TRUE),
  pearson_clip_zscore = c(0.05, 0, TRUE),
  pearson_clip_hvg_zscore = c(0.05, 0, TRUE),
  normalisr_norm = FALSE,
  sanity_map = FALSE,
  sanity_dists = FALSE,
  glmpca = c(FALSE, 0.05),
  raw_counts = FALSE,
  scaled_raw_counts = FALSE
)

current_alphas <- alphas[[trans]]
if (!is.vector(current_alphas)) {
  current_alphas <- list(current_alphas)
}

######### Start Transformation #######

UMI <- readRDS(path)$UMI

expressed_cells <- MatrixGenerics::colSums2(UMI) > 0
expressed_genes <- MatrixGenerics::rowSums2(UMI) > 0
UMI <- UMI[expressed_genes, expressed_cells]

sf <- MatrixGenerics::colSums2(UMI)
sf <- sf / mean(sf)

for (alpha in current_alphas) {
  trans_dat <- all_transformations[[trans]](UMI, sf, alpha)
  for (pc in pcs) {
    for (nn in knns) {
      cat(">> Transform:", trans, "| Alpha:", alpha, "| PC:", pc, "| KNN:", nn, "\n")
      time1 <- Sys.time()
      KNN <- make_knn_graph(trans, trans_dat, pc, nn)
      cat("Completed in", Sys.time() - time1 , "seconds\n")
      result_id <- paste0(data, "_", trans, "_alpha:", alpha, "_pc:", pc, "_nn:", nn)
      saveRDS(KNN, file.path(output_dir, result_id))
    }
  }
}


