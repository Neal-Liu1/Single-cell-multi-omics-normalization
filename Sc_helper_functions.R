

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
                                          Dataset = 'list',
                                          Params = 'list',
                                          RunningTime = 'numeric',
                                          Silhouette = 'list',
                                          ARI = 'data.frame',
                                          LISI = 'data.frame'))



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
                    plot_type = 'violin'){
             standardGeneric('PlotMultipleSilhouette')
           }
)

setMethod('PlotMultipleSilhouette', 
          signature = c(metrics_obj = 'BenchmarkMetrics'),
          function(metrics_obj, 
                   subset, 
                   plot_type){
            
            if(all(subset == 'all')){
              merged_data <- do.call(rbind, metrics_obj@Silhouette)
              merged_data$method <- rep(names(metrics_obj@Silhouette), 
                                        each = length(unique(merged_data$labels)))
            }
            else if(sum(subset %in% names(metrics_obj@Silhouette)) == length(subset)){
              merged_data <- do.call(rbind, metrics_obj@Silhouette[subset])
              merged_data$method <- rep(names(metrics_obj@Silhouette[subset]), 
                                        each = length(unique(merged_data$labels)))
              
            }
            else{stop('Some or all the names you entered for the subset param do not exist. 
                      Please make sure you are entering a vector of strings 
                      corresponding to the names of the reductions.')}
            
            if(plot_type == 'boxplot'){
              plot <- plot_boxplot_categorical(merged_data$silhouette_score,
                                               merged_data$method,
                                               names = c('silhouette score', 'method')) + 
                theme(axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
                      legend.position = 'none')
            }
            
            if(plot_type == 'violin'){
              plot <- plot_violin(merged_data$silhouette_score, 
                                  merged_data$method, 
                                  names = c('silhouette score', 'method'), 
                                  overlay_type = 'boxplot') + 
                theme(axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
                      legend.position = 'None')
            }
            return(plot)
          })















