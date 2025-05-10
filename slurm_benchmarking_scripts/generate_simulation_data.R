

args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])

options(error = function() { cat("Error occurred. Exiting.\n"); quit("no", 1) })
options(Ncpus = 1L)
options(dyngen_download_cache_dir = "/vast/scratch/users/liu.ne/transformGamPoi_simulation_data")

set.seed(seed)
cat("Starting simulations for seed", seed, "\n")

# ======================= scDesign2 Simulation =======================
cat("[1/5] Starting scDesign2 simulation\n")
start_time <- Sys.time()

library(scDesign2)

sce <- get_GSE130931_data()
colData(sce)$cluster_id <- quickCluster(sce, min.size = 20)
sce <- sce[, colSums2(assay(sce)) > 0]
sce <- sce[rowSums2(assay(sce)) > 0, ]
sce <- logNormCounts(sce)

mat <- counts(sce)
colnames(mat) <- as.character(sce$cluster_id)
cluster_names <- unique(colnames(mat))
fit <- fit_model_scDesign2(mat, cell_type_sel = cluster_names, sim_method = "copula", marginal = "nb")

UMI <- simulate_count_scDesign2(
  fit,
  n_cell_new = ncol(mat),
  cell_type_prop = table(colnames(mat))[cluster_names] / ncol(mat),
  sim_method = "copula"
)

ground_truth <- do.call(cbind, lapply(fit[cluster_names], function(fi) {
  mat <- matrix(NA, nrow = nrow(mat), ncol = fi$n_cell)
  mat[fi$gene_sel1, ] <- fi$marginal_param1[, 3]
  mat[fi$gene_sel2, ] <- fi$marginal_param2[, 3]
  mat[fi$gene_sel3, ] <- 1e-8
  mat
}))

sim <- SingleCellExperiment(
  assays = list(counts = UMI, LogExpressionMean = log10(ground_truth)),
  colData = data.frame(cluster = colnames(UMI))
)
sim$seq_depth <- colSums2(counts(sim))

saveRDS(list(
  ground_truth = assay(sim, "LogExpressionMean"),
  UMI = assay(sim, "counts")),
  file.path("/vast/scratch/users/liu.ne/transformGamPoi_simulation_data", paste0("scDesign_seed_", seed, ".rds"))
)

end_time <- Sys.time()
cat("Finished scDesign2 in", round(difftime(end_time, start_time, units = "secs")), "seconds\n")

# ======================= Dyngen Simulation =======================
cat("[2/5] Starting Dyngen simulation\n")
start_time <- Sys.time()

library(dyngen)
library(SingleCellExperiment)
library(Matrix)
source("/home/users/allstaff/liu.ne/scMultiOmics-normalization/transformGamPoi_Paper_benchmark/benchmark/src/consistency_benchmark/download_helper.R")

num_cells <- 5000
num_features <- 1000
backbone <- backbone_consecutive_bifurcating()
num_tfs <- nrow(backbone$module_info)
num_targets <- round((num_features - num_tfs) / 2)
num_hks <- num_features - num_targets - num_tfs

config <- initialise_model(
  backbone = backbone,
  num_tfs = num_tfs,
  num_targets = num_targets,
  num_hks = num_hks,
  num_cells = num_cells
)

res <- generate_dataset(config, format = "sce")
sim <- res$dataset

UMI <- counts(sim)
expressed_cells <- colSums2(UMI) > 10
expressed_genes <- rowSums2(UMI) > 0
sim <- sim[, expressed_cells]
sim <- sim[expressed_genes, ]

saveRDS(list(
  ground_truth = t(reducedDim(sim, "MDS")),
  UMI = assay(sim, "counts")),
  file.path("/vast/scratch/users/liu.ne/transformGamPoi_simulation_data", paste0("Dyngen_seed_", seed, ".rds"))
)

end_time <- Sys.time()
cat("Finished Dyngen in", round(difftime(end_time, start_time, units = "secs")), "seconds\n")


# ======================= Linear Walk Simulation =======================
cat("[3/5] Starting Linear Walk simulation\n")
start_time <- Sys.time()

library(scRNAseq)
library(scuttle)
library(scran)
library(matrixStats)

sce <- BaronPancreasData("human")
sce <- logNormCounts(sce)
colData(sce)$cluster_id <- quickCluster(sce)

reference_data_counts <- assay(sce)
reference_data_counts <- reference_data_counts[, colSums2(reference_data_counts) > 0]
reference_data_counts <- reference_data_counts[rowSums2(reference_data_counts) > 0, ]
n_genes <- nrow(reference_data_counts)
n_cells <- ncol(reference_data_counts)

delta_true <- matrix(NA, n_genes, n_cells)
parents <- rep(NA, n_cells)

branch_length <- 1200
branch_idx <- 0
start_point <- NULL
end_point <- NULL

for (idx in 0:(n_cells - 1)) {
  if (idx == 0) {
    parents[idx + 1] <- -1
    start_point <- rnorm(n_genes)
    end_point <- rnorm(n_genes)
  } else if (idx %% branch_length == 0) {
    start_id <- sample.int(idx, size = 1)
    parents[idx + 1] <- start_id - 1
    branch_idx <- 0
    start_point <- rnorm(n_genes, mean = delta_true[, parents[idx + 1]], sd = 0.1)
    end_point <- rnorm(n_genes, mean = 0, sd = 5)
  } else {
    branch_idx <- branch_idx + 1
    parents[idx + 1] <- idx - 1
  }
  delta <- rnorm(n_genes, mean = start_point + (end_point - start_point) * branch_idx / branch_length, sd = 0.1)
  delta_true[, idx + 1] <- delta
}

N_c <- colSums2(reference_data_counts)
mu_tilde_g <- log(rowSums2(reference_data_counts) / sum(reference_data_counts))
sig2_g <- rexp(n_genes, rate = 1 / 2)
lambda <- sqrt(sig2_g / rowVars(delta_true))
delta_true <- (delta_true - rowMeans2(delta_true)) * lambda
mu_g <- mu_tilde_g - sig2_g / 2

linearwalk_UMI <- matrix(
  rnbinom(n_genes * n_cells, mu = t(t(exp(mu_g + delta_true)) * N_c), size = 1 / 0.01),
  n_genes, n_cells)

saveRDS(list(
  ground_truth = delta_true,
  UMI = linearwalk_UMI),
  file.path("/vast/scratch/users/liu.ne/transformGamPoi_simulation_data", paste0("Linearwalk_seed_", seed, ".rds"))
)

end_time <- Sys.time()
cat("Finished Linear Walk in", round(difftime(end_time, start_time, units = "secs")), "seconds\n")


# ======================= Muscat Simulation =======================
cat("[4/5] Starting Muscat simulation\n")
start_time <- Sys.time()

library(muscat)
library(dplyr)
library(tidyr)
library(SingleCellExperiment)

data(example_sce)
sce_preped <- prepSim(example_sce, verbose = FALSE)

muscat_sim <- simData(
  sce_preped, rel_lfc = c(1, 0.5, 0.1, 0.05),
  nc = 5e3, nk = 4, p_dd = c(0.7, 0, 0.3, 0, 0, 0),
  lfc = 2, ng = 1e3, force = TRUE
)

tmp <- as_tibble(metadata(muscat_sim)$gene_info) %>%
  pivot_longer(starts_with("sim_mean"), names_sep = "\\.", names_to = c(".value", "group_id")) %>%
  dplyr::select(gene, cluster_id, group_id, sim_mean) %>%
  full_join(colData(muscat_sim) %>% as_tibble(rownames = "cell_id") %>%
              mutate(cell_id = factor(cell_id, levels = cell_id)),
            by = c("cluster_id", "group_id"))

ground_truth_mat <- tmp %>%
  arrange(cell_id) %>%
  dplyr::select(gene, cell_id, sim_mean) %>%
  pivot_wider(id_cols = gene, names_from = cell_id, values_from = sim_mean) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

assay(muscat_sim, "LogExpressionMean") <- log10(ground_truth_mat + 1e-4)

saveRDS(list(
  ground_truth = assay(muscat_sim, "LogExpressionMean"),
  UMI = assay(muscat_sim, "counts")),
  file.path("/vast/scratch/users/liu.ne/transformGamPoi_simulation_data", paste0("Muscat_seed_", seed, ".rds"))
)

end_time <- Sys.time()
cat("Finished Muscat in", round(difftime(end_time, start_time, units = "secs")), "seconds\n")


# ======================= Random Walk Simulation =======================
cat("[5/5] Starting Random Walk simulation\n")
start_time <- Sys.time()

# Reuse previous reference_data_counts
delta_true <- matrix(NA, n_genes, n_cells)
parents <- rep(NA, n_cells)
branch_length <- 13

for (idx in seq_len(n_cells)) {
  if (idx == 1) {
    delta <- rnorm(n_genes)
    parents[idx] <- 0
  } else if (idx %% branch_length == 0) {
    parents[idx] <- sample.int(idx - 1, size = 1)
    delta <- rnorm(n_genes, mean = delta_true[, parents[idx]], sd = 1)
  } else {
    parents[idx] <- idx - 1
    delta <- rnorm(n_genes, mean = delta_true[, parents[idx]], sd = 1)
  }
  delta_true[, idx] <- delta
}

lambda <- sqrt(sig2_g / rowVars(delta_true))
randomwalk_delta_true <- (delta_true - rowMeans2(delta_true)) * lambda
mu_g <- mu_tilde_g - sig2_g / 2

randomwalk_UMI <- matrix(
  rnbinom(n_genes * n_cells, mu = t(t(exp(mu_g + delta_true)) * N_c), size = 1 / 0.01),
  n_genes, n_cells
)

saveRDS(list(
  ground_truth = randomwalk_delta_true,
  UMI = randomwalk_UMI),
  file.path("/vast/scratch/users/liu.ne/transformGamPoi_simulation_data", paste0("Randomwalk_seed_", seed, ".rds"))
)

end_time <- Sys.time()
cat("Finished Random Walk in", round(difftime(end_time, start_time, units = "secs")), "seconds\n")
cat("All simulations complete for seed", seed, "\n")
