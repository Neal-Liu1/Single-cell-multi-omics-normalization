
# Compute ARIs
setGeneric('ComputeARIs',
           function(
    obj, labels, hclust_method = 'complete', 
    distance_type  = 'euclidean', 
    sample_fraction = 0.3, 
    num_cross_validation = 5)
           {standardGeneric('ComputeARIs')})

setMethod('ComputeARIs', 
          signature = c(obj = 'BenchmarkMetrics', 
                        labels = 'character'),
          function(
    obj, labels, hclust_method, distance_type, 
    sample_fraction, 
    num_cross_validation){
            label_name <- labels
            labels <- obj@Metadata[[labels]]
            ARIs <- vector("list", length = length(obj@PCs))
            names(ARIs) <- names(obj@PCs)
            for (pc_name in names(obj@PCs)) {
              pc <- obj@PCs[[pc_name]]
              ari_values <- numeric(num_cross_validation)
              for (i in 1:num_cross_validation) {
                sample_indices <- sample(seq_len(nrow(pc)), 
                                         size = ceiling(sample_fraction * nrow(pc)))
                sampled_pc <- pc[sample_indices, , drop = FALSE]
                sampled_labels <- labels[sample_indices]
                
                clusters <- cutree(fastcluster::hclust(d = dist(sampled_pc, method = distance_type), 
                                                       method = hclust_method),
                                   k = length(unique(sampled_labels)))
                
                ari_values[i] <- mclust::adjustedRandIndex(clusters, sampled_labels)}
              ARIs[[pc_name]] <- mean(ari_values)
            }
            obj@ARI[[label_name]] <- ARIs
            return(obj)
          })


# plot ARIs
setGeneric('PlotARIs',
           function(obj, variable, title = "ARIs of different methods")
           {standardGeneric('PlotARIs')})

setMethod('PlotARIs',
          signature= 'BenchmarkMetrics',
          function(obj, variable, title){
            merged_ARIs_df <- data.frame(names = names(obj@ARI[[variable]]), ARI = unlist(obj@ARI[[variable]]))
            merged_ARIs_df$names <- factor(merged_ARIs_df$names, levels = obj@Algorithm)
            ggplot(merged_ARIs_df, aes(x= names , y = ARI, fill = names))+
              geom_bar(stat = "identity", colour = 'black') +
              theme_minimal() +
              labs(title = title,
                   x = NULL,
                   y = "ARI scores")+
              theme(legend.position = 'none',
                    panel.border=element_rect(colour = "grey87", fill=NA, size=0.7),
                    aspect.ratio = 1/1.1,
                    axis.line = element_line(colour = "grey45", linewidth = 0.8),
                    panel.grid.major = element_line(color = "grey96"),
                    axis.text.x = element_text(size = 10,angle = 45,hjust = 1))
          })



# HVG conservation
setGeneric('PlotHVG_Conservation', 
           function(obj, batch_variable = 'batch', flavour = 'vst', title = 'HVG conservation plot')
           {standardGeneric('PlotHVG_Conservation')})

# Take a metrics obj with populated adjusted data slots, 
# compute the percentage of HVGs between the raw data and all the adjusted datasets,
# outputting a bar graph of percentage for each adjusted dataset.
# if raw data has multiple batches, use common HVGs across all batches instead.

setMethod('PlotHVG_Conservation', 
          signature = c(obj = 'BenchmarkMetrics'),
          function(obj, batch_variable, flavour, title)
          {
            # Split raw data by batch, then find common HVGs across all batches
            batch_vector <- obj@Metadata[[batch_variable]]
            unique_batches <- unique(batch_vector)
            batch_data_list <- lapply(unique_batches, function(batch){
              obj@Raw_data[, batch_vector == batch]})
            
            batch_hvgs <- lapply(batch_data_list, function(x)
            {Seurat::FindVariableFeatures(
              x,
              selection.method = flavour)[['variable']]
            }) %>% as.data.frame()
            
            common_hvgs <- rowSums(batch_hvgs) == length(unique_batches)
            
            # Find HVGs for each adjusted dataset
            normalized_hvgs <- lapply(obj@Adj_data, function(x)
            {Seurat::FindVariableFeatures(x, selection.method = flavour)[['variable']]
            })
            
            # Use bitwise AND to find conserved HVGs between adjusted and raw data
            hvg_conserv_scores <- lapply(normalized_hvgs, function(x)
            {mean(x & common_hvgs)}) 
            
            # Plot
            df <- data.frame(method = names(hvg_conserv_scores),
                             values = unlist(hvg_conserv_scores))
            df$method <- factor(df$method, levels = obj@Algorithm)
            
            return(
              ggplot(df, aes(x= method, y = values, fill = method))+
                geom_bar(stat = "identity", colour = 'black') +
                theme_minimal() +
                labs(title = title, x = NULL, y = "% HVG conservation")+
                theme(legend.position = 'none',
                      panel.border=element_rect(colour = "grey87", fill=NA, size=0.7),
                      aspect.ratio = 1/1.1,
                      axis.line = element_line(colour = "grey45", linewidth = 0.8),
                      panel.grid.major = element_line(color = "grey96"),
                      axis.text.x = element_text(size = 10,angle = 45,hjust = 1))+
                coord_cartesian(ylim = c(0,1)))
            
          })



# Plot umap
plot_UMAP <- function(matrix, metadata_vector, title = 'UMAP', aspect_ratio = 1/1, run_umap = T, label_is_continuous = F, 
                      continuous_var_upper_lim = NULL, alpha = 1){
  # taking a matrix and a vector of metadata, plot UMAP of the matrix colored by the groups in the metadata vector
  
  if(run_umap){
    umap_result = umap(t(matrix))
    df = data.frame(umap_result$layout)}
  if(!run_umap){df <- as.data.frame(matrix)}
  colnames(df) = c("UMAP1", "UMAP2")
  
  if(!is.null(continuous_var_upper_lim)){
    if(class(continuous_var_upper_lim) != 'numeric'){stop("You didn't enter a numerical value for the continous variable uppe limit. Please only enter numbers.")}
    else(metadata_vector <- lapply(metadata_vector,
                                   function(x) ifelse(x > continuous_var_upper_lim,
                                                      continuous_var_upper_lim,
                                                      x)) %>% unlist())}
  
  df$metadata = metadata_vector 
  
  if(!label_is_continuous){
    centroids <- aggregate(cbind(UMAP1, UMAP2) ~ metadata, df, mean)
    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = metadata)) +
      geom_point(size = 0.07, alpha = alpha) +
      ggtitle(title) +
      theme_minimal() +
      theme(axis.line = element_line(colour = "grey50", linewidth = 0.9),
            panel.border = element_blank(),  #element_rect(colour = "grey90", fill=NA, size=0.7),
            panel.grid.major = element_blank(),  #element_line(color = "grey96"),
            panel.grid.minor = element_blank(),
            aspect.ratio = aspect_ratio,
            legend.position = "none")+
      geom_text(data = centroids, aes(label = metadata), size = 3, color = "black", hjust = 0.5, vjust = 0.5)
  }
  
  if(label_is_continuous){return(
    ggplot(df, aes(x = UMAP1, y = UMAP2, color = metadata)) +
      geom_point(size = 0.07, alpha = alpha) +
      scale_fill_viridis() +
      scale_color_viridis() +
      ggtitle(title) +
      theme_minimal() +
      theme(axis.line = element_line(colour = "grey50", linewidth = 0.9),
            panel.border = element_blank(),  #element_rect(colour = "grey90", fill=NA, size=0.7),
            panel.grid.major = element_blank(),  #element_line(color = "grey96"),
            panel.grid.minor = element_blank(),
            aspect.ratio = aspect_ratio,
            legend.position = "none")
  )}
  
  if(!label_is_continuous){
    xdens <- axis_canvas(p, axis = "x")+
      geom_density(df, mapping = aes(x = UMAP1, fill = metadata_vector), color= 'grey55', alpha = 0.50, size = 0.2) +
      theme(legend.position = "none")
    
    ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
      geom_density(df, mapping = aes(x = UMAP2, fill = metadata_vector), color= 'grey55', alpha = 0.50, size = 0.2) +
      theme(legend.position = "none")+
      coord_flip()
    
    
    p1 <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
    p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
    pList <- ggdraw(p2)
  }
  
  return(pList)
}




# FindNCG function

# NOT ROBUST TO LOWLY EXPRESSED GENES
# try doing intersection instead of ranking
# for biology, calculate the corr within batch not globally (try for libsize as well)

setGeneric(name = 'find_corr',
           function(data, variable_vector){standardGeneric('find_corr')})

setMethod('find_corr',
          signature = c(variable_vector = 'numeric'),
          function(data, variable_vector){
            results <- Rfast::correls(y=variable_vector, x=t(data), type = 'spearman') %>% as.data.frame()
            return(as.vector(results$correlation))
          })

setMethod('find_corr',
          signature = c(variable_vector = 'character'),
          function(data, variable_vector){
            #results <- ftest(data, variable_vector)
            results <- parallel::mclapply(1:nrow(data),
                                          function(x){summary(aov(as.vector(data[x,])~variable_vector))[[1]][4][1,1]})
            return((unlist(results)))
          })

setGeneric(name = 'FindNCG',
           function(object,
                    unwanted_variables,
                    bio_variables,
                    no.ncg = 1000,
                    apply_log = T,
                    sample_fraction = 0.1)
           {standardGeneric('FindNCG')})


setMethod('FindNCG', 
          signature = c(object = 'Seurat',
                        unwanted_variables = 'character',
                        bio_variables = 'character'),
          function(object,
                   unwanted_variables,
                   bio_variables,
                   no.ncg,
                   apply_log,
                   sample_fraction){
            
            sample_num <- sample_fraction * ncol(object)
            sample_ <- sample(ncol(object), sample_num)
            message(paste0('Sampling ', sample_num,' cells from your data'))
            
            # Hardcoding main data matrix to assays > RNA > layers > counts
            if(apply_log){data <- log2(object@assays$RNA@layers$counts+1)}
            else{data <- object@assays$RNA@layers$counts}
            
            data <- data[,sample_] %>% as.matrix()
            metadata <- object@meta.data[sample_,]
            
            message('Calculating Spearman correlation for continuous variables & F score for categorical variables')
            corr_data <- parallel::mclapply(c(unwanted_variables, bio_variables), 
                                            function(x){find_corr(data, metadata[[x]])})
            
            names(corr_data) <- c(unwanted_variables, bio_variables)
            
            # Flip bio variables 
            message('Finding genes highly affected by unwanted variables but not affected by biological variables')
            for(name in names(corr_data)){
              if(name %in% bio_variables){
                corr_data[[name]] <- corr_data[[name]] * (-1)
              }
            }
            corr_data <- do.call(cbind, corr_data) %>% as.data.frame()
            rownames(corr_data) <- rownames(object[['RNA']])
            ranks <- apply(corr_data, 2, rank) %>% as.data.frame()
            ranks$avg_expr <- rowMeans2(data)
            prod_ranks <- apply(ranks, 1, prod)
            final_gene_ranks <- prod_ranks[base::order(-prod_ranks)]
            message('FindNCG Completed!')
            
            return(names(final_gene_ranks)[1:no.ncg])
            
          })




# The benchmark metrics object 

setClass('BenchmarkMetrics', slots = list(Algorithm = 'character',
                                          Raw_data = 'matrix',
                                          Metadata = 'data.frame',
                                          Adj_data = 'list',
                                          PCs = 'list',
                                          UMAPs = 'list',
                                          Latent_dims = 'list',
                                          RunningTime = 'list',
                                          Silhouette = 'list',
                                          ARI = 'list',
                                          LISI = 'list'))




# Plot silhouette for multiple datasets

setGeneric("ComputeMultipleSilhouette", 
           function(obj,
                    metrics_obj, 
                    reductions, 
                    silhouette_format = 'per_cluster', 
                    labels) 
             standardGeneric("ComputeMultipleSilhouette"))

setMethod('ComputeMultipleSilhouette', 
          signature = c(obj = 'Seurat', 
                        metrics_obj = 'BenchmarkMetrics',
                        reductions = 'character',
                        labels = 'character'), 
          function(obj, 
                   metrics_obj, 
                   reductions, 
                   silhouette_format, 
                   labels){
            
            reductions_list <- list()
            for (name in reductions){
              reductions_list[[name]] <- Embeddings(obj, reduction = name)
            }
            
            label_vector <- obj[[labels]][,1]
            silhouettes <- parallel::mclapply(
              reductions_list,
              function(x){
                compute_silhouette(matrix = x,
                                   label_vector = label_vector,
                                   result_format = silhouette_format)})
            
            for (x in names(silhouettes)) {
              metrics_obj@Silhouette[[x]] <- silhouettes[[x]]
            }
            
            return(metrics_obj)
          })


setMethod('ComputeMultipleSilhouette', 
          signature = c(obj = 'list', 
                        metrics_obj = 'BenchmarkMetrics',
                        labels = 'character'), 
          function(obj, 
                   metrics_obj, 
                   reductions, 
                   silhouette_format, 
                   labels){
            # input a list of matrices of dimensionally reduced data (rows are 
            # samples and cols are components), calculate silhouette and store in metrics obj.
            
            reductions_list <- obj
            silhouettes <- parallel::mclapply(
              reductions_list,
              function(x){
                compute_silhouette(matrix = x,
                                   label_vector = labels,
                                   result_format = silhouette_format)})
            
            for (x in names(silhouettes)) {
              metrics_obj@Silhouette[[x]] <- silhouettes[[x]]
            }
            
            return(metrics_obj)
          })


setGeneric('PlotMultipleSilhouette',
           function(metrics_obj, 
                    subset = 'all', 
                    plot_type = 'violin',
                    title = NULL,
                    aspect_ratio = 1.2){
             standardGeneric('PlotMultipleSilhouette')
           }
)

setMethod('PlotMultipleSilhouette', 
          signature = c(metrics_obj = 'BenchmarkMetrics'),
          function(metrics_obj, 
                   subset, 
                   plot_type,
                   title,
                   aspect_ratio){
            
            if(all(subset == 'all')){
              merged_data <- do.call(rbind, metrics_obj@Silhouette)
              merged_data$method <- rep(names(metrics_obj@Silhouette), 
                                        each = length(unique(merged_data$labels)))
            }
            else if(sum(subset %in% names(metrics_obj@Silhouette)) == length(subset)){
              merged_data <- do.call(rbind, metrics_obj@Silhouette[subset])
              merged_data$method <- rep(names(metrics_obj@Silhouette[subset]), 
                                        each = length(unique(merged_data$labels)))}
            
            else{stop('Some or all the names you entered for the subset param do not exist. 
                      Please make sure you are entering a vector of strings 
                      corresponding to the names of the reductions.')}

            merged_data$method <- factor(merged_data$method, levels = metrics_obj@Algorithm)
            if(plot_type == 'boxplot'){
              plot <- plot_boxplot_categorical(merged_data$silhouette_score,
                                               merged_data$method,
                                               names = c('silhouette score', 'method')) + 
                theme(axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
                      legend.position = 'none')
            }
            
            if(plot_type == 'violin'){
              plot <- ggplot(merged_data, aes(x=method, y=silhouette_score, fill=method))+
                geom_violin()+
                labs(x = 'Method', y = 'Silhouette', fill='Method')+
                ggtitle(title)+
                theme_minimal() +
                theme(panel.border=element_rect(colour = "grey80", fill=NA, size=0.8),
                      axis.line = element_line(colour = "grey75", linewidth = 1.1),
                      panel.grid.major = element_line(color = "grey96"),
                      aspect.ratio = 1/aspect_ratio,
                      axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
                      legend.position = 'None')
              plot <- plot + geom_boxplot(width=0.1, alpha=0.5, outlier.shape =NA)}
            return(plot)
          })









setGeneric('ComputeMultipleLISI',
           function(obj, reductions , variables, metadata, metrics_obj)
           {standardGeneric('ComputeMultipleLISI')})


# Takes a list of reduction components (eg PCs) as dataframes or matrices (cells as rows and features as cols), 
# and a dataframe of metadata and a vector of metadata variables,
# Compute the LISI scores for every reduction for every variable
# Returns a list of dataframes
setMethod('ComputeMultipleLISI', 
          signature = c(obj = 'list',
                        reductions = NULL,
                        variables = 'character',
                        metadata = 'data.frame',
                        metrics_obj = NULL),
          function(obj, reductions , variables, metadata, metrics_obj){
            LISI_scores <- lapply(obj, function(x){compute_lisi(x,
                                                                meta_data = metadata,
                                                                label_colnames = variables)})
            return(LISI_scores)
          })


# Takes a metrics obj with reductions in the Reductions slot (reductions are list of dataframes of reduction components (eg PCs) (cells as rows and features as cols), 
# and matadata as dataframe in the metadata slot,
# Compute the LISI scores for every reduction for every variable
# Returns the metrics obj with LISI scores in the LISI slots.
setMethod('ComputeMultipleLISI', 
          signature = c(obj = 'BenchmarkMetrics',
                        reductions = NULL,
                        variables = 'character',
                        metadata = NULL,
                        metrics_obj = NULL),
          function(obj, reductions , variables, metadata, metrics_obj){
            if(length(obj@Metadata) == 0)
            {stop("Your Metrics obj doesn't have any metadata. Please add metadata first.")}
            if(length(obj@PCs) == 0)
            {stop("Your Metrics obj doesn't have any PCs. Please compute PCA first.")}
            
            reduction_matrix_list <- obj@PCs
            metadata <- obj@Metadata
            
            LISI_scores <- lapply(reduction_matrix_list,
                                  function(x){
                                    compute_lisi(
                                      x,
                                      meta_data = metadata,
                                      label_colnames = variables)})
            obj@LISI <- LISI_scores
            return(obj)
          })


setGeneric('PlotMultipleLISI',
           function(obj, variable, reductions, aspect_ratio = 1.3, 
                    title = NULL, levels = NULL)
           {standardGeneric('PlotMultipleLISI')})

setMethod('PlotMultipleLISI',
          signature = c(obj = 'BenchmarkMetrics',
                        reductions = NULL,
                        variable = 'character'),
          function(obj, variable, reductions, aspect_ratio, title, levels){
            
            if(!(variable %in% names(obj@LISI[[1]])))
            {stop('LISI scores have not yet been computed for this variable you entered.')}
            
            scores <- unlist(lapply(obj@LISI, function(x){x[[variable]]}))
            names <- rep(names(obj@LISI), each = ncol(obj@Raw_data))
            data <- data.frame(LISI = scores, method = names)
            
            if(!is.null(levels))
            {data$method <- factor(data$method, levels = levels)}
            
            ggplot(data, aes(x=method, y = LISI, fill=method))+
              geom_boxplot(outlier.shape = NA)+
              labs(x = NULL) +
              ggtitle(title)+
              theme_minimal()+
              theme(legend.position = 'none',
                    panel.border=element_rect(colour = "grey87", fill=NA, size=0.7),
                    aspect.ratio = 1/aspect_ratio,
                    axis.line = element_line(colour = "grey50", linewidth = 0.7),
                    panel.grid.major = element_line(color = "grey96"),
                    axis.text.x = element_text(size = 10,angle = 45,hjust = 1))
          })




# Linear and vector correlations as S4 methods

setGeneric('PlotCorrelations', 
           function(obj,
                    variable,
                    reductions = NULL,
                    num_pcs = 10,
                    title = 'Correlation plot')
           {standardGeneric('PlotCorrelations')})

setMethod('PlotCorrelations',
          signature = c(obj = 'BenchmarkMetrics',
                        variable = 'character',
                        reductions = NULL),
          function(obj,
                   variable,
                   reductions,
                   num_pcs,
                   title){
            
            var = obj@Metadata[[variable]]
            reductions = obj@PCs
            
            if(class(var) %in% c('factor', 'character')){
              dummies <- fastDummies::dummy_cols(var)[,-1]
              cancor_scores <- sapply(reductions, function(m) {lapply(1:num_pcs, function(y) {cca <- stats::cancor(x= m[,1:y, drop=F], 
                                                                                                                   y= dummies) 
              1 - prod(1 - cca$cor^2)})}) %>% unlist()
              PCs <- rep(paste0('PC1:', 1:num_pcs), length(reductions))
              Datasets <- rep(names(reductions), each = num_pcs)
              Datasets <- factor(Datasets, levels = obj@Algorithm)
              pc_correlations <- data.frame(PCs, cancor_scores, Datasets)
              pc_correlations$PCs <- factor(pc_correlations$PCs, levels = paste0('PC1:', 1:num_pcs))
              
              return(
                ggplot(pc_correlations, aes(x = PCs, y = cancor_scores, color = Datasets, group = Datasets)) +
                  geom_line(size=0.5,alpha=0.8) +
                  geom_point(alpha=0.8) +
                  labs(x = 'Principal Components', y = 'Correlation', color = 'Dataset') +
                  ylim(0,1) +
                  theme_minimal() +
                  theme(axis.line = element_line(colour = "grey83", linewidth = 1.1),
                        panel.border = element_rect(colour = "grey90", fill=NA, size=0.7),
                        panel.grid.major = element_line(color = "grey96"),
                        aspect.ratio = 1/1.2) +
                  ggtitle(title)
              )}
            
            if(class(var) %in% 'numeric'){
              R_squared <- sapply(reductions, function(matrix, var, num_pcs) {
                sapply(1:num_pcs, function(y) {
                  lm_model <- summary(lm(var ~ matrix[,1:y]))$r.squared})
              }, num_pcs = num_pcs, var = var) %>% as.vector()
              PCs <- rep(paste0('PC1:', 1:num_pcs), length(reductions))
              Datasets <- rep(names(reductions), each = num_pcs) %>% factor(levels = obj@Algorithm)
              pc_correlations <- data.frame(PCs, R_squared, Datasets)
              pc_correlations$PCs <- factor(pc_correlations$PCs, levels = paste0('PC1:', 1:num_pcs))
              
              return(
                ggplot(pc_correlations, aes(x = PCs, y = R_squared, color = Datasets, group = Datasets)) +
                  geom_line(size=0.5,alpha=0.8) +
                  geom_point(alpha=0.8) +
                  labs(x = 'Principal Components', y = 'R-squared', color = 'Dataset') +
                  ylim(0,1) +
                  ggtitle(title)+
                  theme_minimal() +
                  theme(axis.line = element_line(colour = "grey83", linewidth=1.1),
                        panel.border = element_rect(colour = "grey90", fill=NA, size=0.7),
                        panel.grid.major = element_line(color = "grey96"),
                        aspect.ratio = 1/1.2)
              )}
          })


# Plot multiple correlation plots in a panel

setGeneric('PlotMultipleCorrelations', function(obj, variables, reductions = 'all', num_pcs = 10, titles = NULL)
{standardGeneric('PlotMultipleCorrelations')})

setMethod('PlotMultipleCorrelations',
          signature = c(obj='BenchmarkMetrics',
                        variables = 'character'),
          function(obj, variables, reductions, num_pcs, titles){
            
            theme <- theme(legend.position = 'none',
                           axis.line = element_line(colour = "grey33", linewidth=0.5),
                           panel.border = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.text.x = element_text(angle = 45,hjust = 1),
                           axis.title.x  = element_blank(),
                           axis.title.y = element_blank(),
                           plot.title = element_text(hjust = 0.5),
                           axis.line.y = element_blank(),
                           axis.text.y = element_blank())
            
            plots <- lapply(variables, function(x)
            {PlotCorrelations(obj, variable = x, num_pcs = num_pcs, title = x) +
                theme})
            
            if(!(is.null(titles))){
              plots <- lapply(1:length(plots), 
                              function(x){plots[[x]]$labels$title <- titles[[x]]; plots[[x]]})}
            
            plots[[1]] <- plots[[1]] +
              theme(legend.position = 'none',
                    axis.line = element_line(colour = "grey33", linewidth=0.5),
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_text(angle = 45,hjust = 1),
                    plot.title = element_text(hjust = 0.5),
                    axis.line.y = element_line(),
                    axis.text.y = element_text())
            
            return(ggpubr::ggarrange(plotlist = c(plots),
                                     ncol = length(plots),
                                     widths = c(rep(1,length(plots))), 
                                     common.legend = T, 
                                     legend = 'right',
                                     align = 'v')
            )
          }
)

# Plot runtimes stored in the metrics object
PlotRuntime <- function(obj, title = 'Runtime'){
  df <- data.frame(time = unlist(obj@RunningTime),
                   method = names(obj@RunningTime))
  df$method <- factor(df$method, levels = obj@Algorithm)
  return(
    ggplot(df, aes(x=method, y= time, fill = time)) +
      geom_bar(stat = "identity", colour = 'black') +
      theme_minimal() +
      labs(title = title, x = NULL, y = "Runtime (seconds)")+
      theme(legend.position = 'none',
            panel.border=element_rect(colour = "grey87", fill=NA, size=0.7),
            aspect.ratio = 1/1.1,
            axis.line = element_line(colour = "grey45", linewidth = 0.8),
            panel.grid.major = element_line(color = "grey96"),
            axis.text.x = element_text(size = 10,angle = 45,hjust = 1)))
}


