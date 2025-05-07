source(
  "/home/users/allstaff/liu.ne/scMultiOmics-normalization/transformGamPoi_Paper_benchmark/benchmark/src/transformations/transformation_helper.R"
)

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
cat("Starting task for param_grid row:", i, "\n")

param_grid <- readRDS("consistency_params.rds")

row <- param_grid[i, ]
data <- row$data
trans <- row$trans
cat("Data:", data, "| Transformation:", trans, "\n")

datasets <- c(
  "GSE142647", "GSE178765", "GSE179831", "GSE164017", "GSE150068", 
  "GSE130931", "GSE163505", "GSE158941", "GSE179714", "GSE184806"
)
data_paths <- paste0("/vast/scratch/users/liu.ne/transformGamPoi_10x_data/", datasets) |> as.list()
names(data_paths) <- datasets

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
if (!is.vector(current_alphas)) current_alphas <- list(current_alphas)

# Load UMI matrix
cat("Loading UMI matrix for", data, "...\n")
UMI <- readRDS(data_paths[[data]])
expressed_cells <- matrixStats::colSums2(UMI) > 0
expressed_genes <- matrixStats::rowSums2(UMI) > 0
UMI <- UMI[expressed_genes, expressed_cells]

# Split genes
cat("Splitting genes for consistency benchmark...\n")
first_gene_half <- sample(seq_len(nrow(UMI)), round(nrow(UMI) / 2))
second_gene_half <- setdiff(seq_len(nrow(UMI)), first_gene_half)
UMI_1 <- UMI[first_gene_half, , drop = FALSE]
UMI_2 <- UMI[second_gene_half, , drop = FALSE]

# Size factors
sf_1 <- MatrixGenerics::colSums2(UMI_1)
sf_1 <- sf_1 / mean(sf_1)

sf_2 <- MatrixGenerics::colSums2(UMI_2)
sf_2 <- sf_2 / mean(sf_2)

# Ensure output directory exists
out_dir <- "/vast/scratch/users/liu.ne/transformGamPoi_Output/results/consistency_results/transformed_knns"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Apply methods and log progress
for (alpha in current_alphas) {
  for (pc in pcs) {
    for (nn in knns) {
      cat(">> Transform:", trans, "| Alpha:", alpha, "| PC:", pc, "| KNN:", nn, "\n")
      duration <- system.time({
        trans_dat1 <- all_transformations[[trans]](UMI_1, sf_1, alpha)
        trans_dat2 <- all_transformations[[trans]](UMI_2, sf_2, alpha)
        
        KNN_1 <- make_knn_graph(trans, trans_dat1, pc, nn)
        KNN_2 <- make_knn_graph(trans, trans_dat2, pc, nn)
      })
      cat("Completed in", round(duration[3], 2), "seconds\n")
      
      result_id <- paste0(data, "_", trans, "_alpha:", alpha, "_pc:", pc, "_nn:", nn)
      save_path <- file.path(out_dir, result_id)
      saveRDS(list(KNN_1 = KNN_1, KNN_2 = KNN_2), save_path)
      cat("Saved result to", save_path, "\n")
    }
  }
}

cat("Done with task", i, "\n")