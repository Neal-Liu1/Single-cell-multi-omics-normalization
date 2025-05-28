library(tidyverse)

# Get SLURM array task ID
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
message("Starting SLURM task ID: ", i)

# Load result group
data_info <- readRDS("/home/users/allstaff/liu.ne/scMultiOmics-normalization/slurm_benchmarking_scripts/simulation_result_paths.rds")
result_list <- data_info[[i]]
message("Loaded result list with ", nrow(result_list), " entries.")

# Load ground truth
ground_truth_path <- paste0(
  "/vast/scratch/users/liu.ne/transformGamPoi_simulation_data/",
  str_extract(result_list[1, ]$filename, "^.+?\\.rds")
)
ground_truth <- readRDS(ground_truth_path)$ground_truth
message("Loaded ground truth from: ", ground_truth_path)

# Loop over each result
for (j in seq_len(nrow(result_list))) {
  result_info <- result_list[j, ]
  start_time <- Sys.time()
  
  message("Processing file ", j, "/", nrow(result_list), ": ", result_info$filename)
  
  # Extract alpha, pc, nn
  matches <- str_match(
    result_info$filename,
    "alpha:([0-9.]+|TRUE|FALSE)_pc:([0-9]+)_nn:([0-9]+)"
  )
  alpha <- matches[, 2]
  pc    <- as.numeric(matches[, 3])
  nn    <- as.numeric(matches[, 4])
  
  knn <- readRDS(result_info$path)
  stopifnot(ncol(ground_truth) == nrow(knn))
  
  ground_truth_knn <- BiocNeighbors::findAnnoy(
    t(ground_truth),
    k = nn,
    warn.ties = FALSE,
    get.distance = FALSE
  )$index
  
  knn_overlap <- mean(sapply(seq_len(nrow(knn)), function(cell_idx) {
    length(intersect(ground_truth_knn[cell_idx, ], knn[cell_idx, ]))
  }))
  
  ground_truth_knn_graph <- bluster::neighborsToKNNGraph(ground_truth_knn)
  ground_truth_clustering_obj <- igraph::cluster_walktrap(ground_truth_knn_graph)
  ground_truth_clustering <- factor(igraph::membership(ground_truth_clustering_obj))
  n_clusters <- length(levels(ground_truth_clustering))
  message("  Ground truth clusters: ", n_clusters)
  
  knn_graph <- bluster::neighborsToKNNGraph(knn)
  tmp_clustering_obj <- igraph::cluster_walktrap(knn_graph)
  cluster_assignment <- factor(igraph::membership(tmp_clustering_obj))
  n_clusters_counts <- length(levels(cluster_assignment))
  message("  Clusters found: ", n_clusters_counts)
  
  ARI <- aricode::ARI(ground_truth_clustering, cluster_assignment)
  AMI <- aricode::AMI(ground_truth_clustering, cluster_assignment)
  NMI <- aricode::NMI(ground_truth_clustering, cluster_assignment)
  
  res <- tibble(
    ARI = ARI,
    AMI = AMI,
    NMI = NMI,
    mean_knn_overlap = knn_overlap,
    n_clusters = n_clusters,
    n_clusters_counts = n_clusters_counts,
    simulator = result_info$simulation,
    seed = result_info$seed,
    pca_dim = pc,
    knn = nn,
    transformation = result_info$transformation,
    alpha = alpha
  )
  
  out_path <- file.path(
    "/vast/scratch/users/liu.ne/transformGamPoi_Output/results/simulation_results/knn_overlap_metrics",
    paste0(result_info$filename, ".tsv")
  )
  
  write_tsv(res, out_path)
  end_time <- Sys.time()
  duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
  message("Fnished ", result_info$filename, " in ", duration, " seconds.")
}


