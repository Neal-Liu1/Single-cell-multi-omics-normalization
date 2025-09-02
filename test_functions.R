
FindCorrectedMultimodalNeighbours <- function(
    seurat_obj, 
    assays, 
    uv_variables, 
    npcs,
    graph.name = 'harmony_wsnn',
    normalization_methods = c('LogNormalize', 'CLR')){
  
  if(length(assays) != length(normalization_methods))
  {stop('The number of assays and normalization methods must be the same (and in the same order)')}
  
  for(i in 1:length(assays)){
    Seurat::DefaultAssay(seurat_obj) <- assays[[i]]
    seurat_obj <- Seurat::NormalizeData(
      seurat_obj,
      normalization.method = normalization_methods[[i]],
      margin = 2) %>%
      Seurat::FindVariableFeatures() %>% 
      Seurat::ScaleData() %>% 
      Seurat::RunPCA(
        reduction.name = paste0('pca_', assays[[i]]),
        npcs = npcs[[i]])
    message('Normalization & PCA completed, removing batch effect with harmony')
    seurat_obj <- harmony::RunHarmony(
      seurat_obj, group.by.vars = uv_variables,
      reduction.use = paste0('pca_', assays[[i]]), 
      reduction.save = paste0('harmony_',assays[[i]]))
  }
  
  reductions <- paste0('harmony_', assays)
  dims <- list(1:ncol(seurat_obj[[reductions[[1]]]]), 
               1:ncol(seurat_obj[[reductions[[2]]]]))

  message('Computing multi-modal neighbours with')
  seurat_obj <- Seurat::FindMultiModalNeighbors(
    seurat_obj, 
    reduction.list = reductions, 
    dims.list = dims,
    snn.graph.name = graph.name)
  
  return(seurat_obj)
}


residop = function(A, B){
  decomp = qr(B)
  qr.resid(decomp, A)
}

#' Runs an expression, times it, and reports the time taken
#' @param expr An expression to run
#' @param description A character scalar, describing the operation as a noun
time_eval <- function(
    expr,
    description    
){
  message(paste0('Performing ', description))
  elapsed = system.time(expr)[["elapsed"]]
  message(paste0(
    'Finished performing ',
    description,
    ' in ',
    round(elapsed, digits = 2),
    ' seconds.'
  ))
}

tological <-function (ctl, n) 
{
  ctl2 = rep(FALSE, n)
  ctl2[ctl] = TRUE
  return(ctl2)
}

Sparse_RUV_III <- function (
    Y,
    Yrep,
    M,
    ctl,
    k = NULL,
    eta = NULL,
    Ynord = NULL,
    eigVec = NULL,
    include.intercept = TRUE,
    average = FALSE,
    return.info = FALSE,
    inputcheck = FALSE) 
{ 
  require(ruv)
  require(Matrix)
  m <- nrow(Y)
  n <- ncol(Y)
  ctl <- tological(ctl, n)
  message('check the inputs finished')
  ############# RUV-I
  time_eval({
    Y <- RUV1(Y, eta, ctl, include.intercept = include.intercept)
    Yrep <- RUV1(Yrep, eta, ctl, include.intercept = include.intercept)
  }, "RUV1 on Y and Yrep separately")
  if (ncol(M) >= m)
    newY = Y
  else if (is.null(k)) {
    ycyctinv = solve(Y[, ctl] %*% t(Y[, ctl]))
    newY = M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*% ycyctinv) %*% Y
    fullalpha = NULL
  }
  else if (k == 0) {
    newY = Y
    fullalpha = NULL
  }
  else {
    if (is.null(Ynord) & is.null(eigVec)) {
      ############# residual operators
      message('Y0 and eigVec are not provided')
      
      time_eval({
        Y0 <- residop(Yrep, M)
      }, "residual operator on Yrep")
      
      ############# eigen vectors
      time_eval({
        eigenvector <- BiocSingular::runSVD(
          x = Y0,
          k = k,
          BSPARAM = BiocSingular::FastAutoParam(),
          center = FALSE,
          scale = FALSE
        )$u
        if (!return.info){
          rm(Y0)
          gc()
        }
      }, "eigendecomposition on Y0")
      ############# fullalpha
      time_eval({
        fullalpha <- t(eigenvector[, 1:k, drop = FALSE]) %*% Yrep
      }, "calculation of fullalpha ram")
    }
    if (!is.null(Ynord) & is.null(eigVec)) {
      ############# eigen vectors
      message('Ynord is provided')
      time_eval({
        eigenvector <- BiocSingular::runSVD(
          x = Ynord %*% t(Ynord),
          k = k,
          BSPARAM = BiocSingular::bsparam(),
          center = TRUE,
          scale = FALSE
        )$u
      }, "eigendecomposition on Ynord")
      ############# fullalpha
      time_eval({
        fullalpha <- t(eigenvector[, 1:k, drop = FALSE]) %*% Yrep
      }, "obtaining fullalpha")
    }
    if (is.null(Ynord) & !is.null(eigVec)) {
      message('eigVec is provided')
      time_eval({
        fullalpha <- t(eigVec[, 1:k, drop = FALSE]) %*% Yrep
      }, "obtaining fullalpha")
    }
    ############# alpha
    time_eval({
      alpha <- fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
      ac <- alpha[, ctl, drop = FALSE]
    }, "obtaining alpha")
    ############# Standardize
    
    # THIS MIGHT BE A PROBLEM
    time_eval({
      Y_stand <- scale(Y[ , ctl], center=TRUE, scale=FALSE)
    }, 'standardization of the data')

    
    ############# W
    time_eval({
      W <- Y_stand %*% t(ac) %*% solve(ac %*% t(ac))
      rm(Y_stand)
      gc()
    }, 'calculationg of W')
    #####################################  data adjustment
    time_eval({
      newY <- Y - (W %*% alpha)
    }, "adjustment of Y")
  }
  if (average)
    newY = ((1 / apply(M, 2, sum)) * t(M)) %*% newY
  if (!return.info)
    return(newY)
  if (is.null(Ynord) & is.null(eigVec))
    return(list(
      newY = newY,
      eigenvector = eigenvector,
      W = W,
      Ynord = Y0
    ))
  if (is.null(eigVec))
    return(list(
      newY = newY,
      eigenvector = eigenvector,
      W = W
    ))
  else
    return(list(
      newY = newY,
      W = W))
}


