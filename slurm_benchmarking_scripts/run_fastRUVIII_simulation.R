
library(tidyr)
library(Seurat)
library(BiocSingular)
library(uwot)
library(fastRUVIII)
library(BiocNeighbors)
library(bluster)
library(igraph)
library(aricode)
library(tibble)
library(readr)
source("/home/users/allstaff/liu.ne/scMultiOmics-normalization/transformGamPoi_Paper_benchmark/benchmark/src/transformations/transformation_helper.R")

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

simulation <- c("Dyngen", "scDesign", "Muscat", "Linearwalk", "Randomwalk")
cluster_res <- c(0.1, 0.5, 1)
ruv_params <- tidyr::crossing(seed = 1:5, simulation, cluster_res)
knns <- c(10, 20, 50)
pcs <- c(5, 10, 50)

row <- ruv_params[i,]

seed <- row$seed
simulation <- row$simulation
cluster_res <- row$cluster_res

ground_truth_path <- paste0(
  "/vast/scratch/users/liu.ne/transformGamPoi_simulation_data/", simulation,"_seed_",seed,".rds")
ground_truth_data <- readRDS(ground_truth_path)

matrix <- ground_truth_data$UMI %>% as("dgCMatrix")
if(is.null(colnames(matrix))){colnames(matrix) <- seq_len(ncol(matrix))}
if(simulation %in% c("scDesign")){colnames(matrix) <- seq_len(ncol(matrix))}

hvgs <- Seurat::VST(matrix, nselect = 0.2 * nrow(matrix))

pca <- BiocSingular::runSVD(
  t(matrix[hvgs$variable,]), center= T, scale = F, 
  BSPARAM = BiocSingular::bsparam(), 
  k = 20)

#umap <- uwot::umap(pca$u, n_neighbors = 20, n_components = 3)

neighbors <- Seurat::FindNeighbors(as(pca$u,'dgCMatrix'), k.param = 20)
clusters <- Seurat::FindClusters(neighbors$snn, resolution = cluster_res)

libsize <- log2(colSums(matrix))
metadata <- data.frame(libsize = libsize, cluster = unlist(clusters))
k <- c(1,2)

prpc <- fastRUVIII::CreatePRPC(
  matrix, 
  bio_vars = 'cluster', 
  uv_vars = 'libsize', 
  separate_bins_by_biology = T,
  group_by_vars = NA,
  metadata = metadata)

ncgs <- fastRUVIII::FindNCG(
  matrix, 
  unwanted_variables = 'libsize', 
  bio_variables = 'cluster', 
  metadata = metadata, 
  no.ncg = 0.1 * nrow(matrix))

ground_truth <- ground_truth_data$ground_truth

for(k in k){
  norm_data <- fastRUVIII::fastRUVIII(
    matrix, 
    k= k, 
    replicates = prpc, 
    control_genes = ncgs)$newY %>% t()
  
  for(pc in pcs){
    for(nn in knns){
      
      knn <- make_knn_graph("fastRUVIII", dat =  norm_data, pca_dim = pc, k_nearest_neighbors = nn)
      
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
        simulator = simulation, 
        seed = seed, 
        pca_dim = pc,
        knn = nn,
        transformation = "fastRUVIII",
        alpha = FALSE
      )
      
      out_path <- file.path(
        "/vast/scratch/users/liu.ne/transformGamPoi_Output/results/simulation_results/knn_overlap_metrics",
        paste0(simulation, "_seed_",seed,"_fastRUVIII_k", k,"_clust_res:",cluster_res,"_pc:",pc,"_nn:",nn,".tsv")
      )
      
      write_tsv(res, out_path)
    }
  }}



