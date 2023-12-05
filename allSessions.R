# Data and functions shared across all sessions in the same R process
# Utility functions not depending on output and input
# Libraries
suppressMessages(library(shiny))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(gplots))
suppressMessages(library(MASS))
suppressMessages(library(RColorBrewer))
suppressMessages(library(DT))
suppressMessages(library(visNetwork))
suppressMessages(library(shinyjs))
suppressMessages(library(sRACIPE))
suppressMessages(library(data.table))
suppressMessages(library(shinyBS))
suppressMessages(library(SingleCellExperiment))
suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(FNN)
  library(ggfortify)
  library(varhandle)
  library(MASS)
  library(metR)
  library(ggnewscale)
  library(gridExtra)
  library(patchwork)
})

plotColor <- c("#5E4FA2", "#4F61AA", "#4173B3", "#3386BC", "#4198B6",
               "#51ABAE", "#62BEA6", "#77C8A4", "#8ED1A4", "#A4DAA4",
               "#B8E2A1", "#CBEA9D", "#DEF199", "#EAF69F", "#F2FAAC",
               "#FAFDB8", "#FEFAB6", "#FEF0A5", "#FEE695", "#FDD985",
               "#FDC978", "#FDB96A", "#FCA75E", "#F99254", "#F67D4A",
               "#F26943", "#E85A47", "#DE4B4B", "#D33C4E", "#C1284A",
               "#AF1446", "#9E0142")
col2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
          "#D55E00", "#CC79A7")

options(shiny.sanitize.errors = FALSE)

# Utility functions
.sracipePlotDensity = function(plotData,binCount=40, col=plotColor, 
                               label1 = "x",
                               label2 = "y",imagetitle = NULL, ...){
  
  colnames(plotData[,1:2]) <- c(label1, label2)
  h1 <- hist(as.numeric(plotData[,1]), breaks=binCount, plot=FALSE)
  h2 <- hist(as.numeric(plotData[,2]), breaks=binCount, plot=FALSE)
  top <- max(h1$counts, h2$counts)
  kernelDensity <- kde2d(plotData[,1], plotData[,2], n=binCount)
  
  oldpar <- par()
  par(mgp=c(2,1,0),mar=c(3,3,1,1))
  layout(matrix(c(2,0,1,3),2,2,byrow=TRUE),c(3,1), c(1,3))
  image(kernelDensity, cex = 1, cex.axis = 1, cex.lab = 1, 
        xlab= label1, ylab=label2, col=plotColor) 
  title(imagetitle, line = -2, adj = 0)
  par(mar=c(0,2,1,0))
  barplot(h1$counts, axes=FALSE, ylim=c(0, top), space=0, col='red')
  
  par(mar=c(2,0,0.5,1))
  barplot(h2$counts, axes=FALSE, xlim=c(0, top), space=0, col='red', horiz=TRUE)
  
}

.sracipePlotHeatmapInt <- function(plotData, nSamples = 500, 
                                   plotColor = plotColor, ...) {
  if(is.null(plotColor)){
    plotColor <- c("#5E4FA2", "#4F61AA", "#4173B3", "#3386BC", "#4198B6",
                   "#51ABAE", "#62BEA6", "#77C8A4", "#8ED1A4", "#A4DAA4",
                   "#B8E2A1", "#CBEA9D", "#DEF199", "#EAF69F", "#F2FAAC",
                   "#FAFDB8", "#FEFAB6", "#FEF0A5", "#FEE695", "#FDD985",
                   "#FDC978", "#FDB96A", "#FCA75E", "#F99254", "#F67D4A",
                   "#F26943", "#E85A47", "#DE4B4B", "#D33C4E", "#C1284A",
                   "#AF1446", "#9E0142")
  }
  hmap <- heatmaply(plotData, col=plotColor)
}



.sracipePlotHeatmap <- function(plotData, nSamples = 500, col = plotColor, 
                                nClusters = 2, 
                                assignedClusters = NULL, ...) {
  plotData <- plotData[,(1:min(nSamples,ncol(plotData)))]
  
  if(is.null(assignedClusters)){
    distfun=function(x) as.dist((1-cor(t(x), method = "s"))/2)
    #ref_cor <- cor((plotData), method = "spearman")
    #distance <- as.dist((1 - ref_cor) / 2)
    distance <- distfun(t(plotData))
    clusters <- hclust(distance, method = "ward.D2")
    clustersRow <- hclust(distfun(plotData), method = "ward.D2")
    assignedClusters <- cutree(clusters, nClusters)
  }
  clustNames <- unique(assignedClusters)
  nClusters <- length(clustNames)
  clustColors <- numeric(length(assignedClusters))
  for(tmp1 in seq_len(length(clustColors))){
    clustColors[tmp1] <- which(clustNames == assignedClusters[tmp1] )
  }
  col2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
            "#D55E00", "#CC79A7")
  clustColors <- col2[clustColors]
  names(clustColors) <- assignedClusters
  lmat <- rbind( c(5,4), c(5,1), c(3,2) )
  lhei <- c(0.75,.25, 4)
  lwid <- c(0.5, 4)
  gplots::heatmap.2(plotData,  lmat=lmat, 
                    lhei=lhei, lwid=lwid,
                    col=plotColor, 
                    Colv = as.dendrogram(clusters),
                    Rowv = as.dendrogram(clustersRow),
                    hclustfun = function(x) hclust(x,method = 'ward.D2'), 
                    distfun=function(x) as.dist((1-cor(t(x), method = "s"))/20), 
                    trace="none",
                    ColSideColors = clustColors, labCol = FALSE,key.title = NA, 
                    keysize = 0.25,
                    key.xlab = NA,key.ylab = NA,
                    xlab="Models",
                    margins = c(2, 5), ...
  )
}

.sracipePlotPca <- function(plotData, pca = NULL, title = NULL, ...){
  if(is.null(pca)) {pca = prcomp(t(plotData), scale. = FALSE)}
  pcaData <- data.frame(x=pca$x[,1],y=pca$x[,2])
  .sracipePlotDensity(pcaData, 
                      label1 = paste0("PC1(",100*summary(pca)$importance[2,1],"%)"),
                      label2 = paste0("PC2(",100*summary(pca)$importance[2,2],"%)"),
                      imagetitle = title, ...)
}

.sracipePlotPcaPoints <- function(plotData, pca = NULL, plotColor = plotColor, ...){
  if(is.null(pca)) {pca = prcomp(t(plotData), scale. = FALSE)}
  pcaData <- data.frame(x=pca$x[,1],y=pca$x[,2])
  ggplot2::ggplot(pcaData) +
    geom_point(aes(x = pcaData[,1], y =pcaData[,2])) +
    xlab(paste0("PC1(",100*summary(pca)$importance[2,1],"%)")) +
    ylab(paste0("PC2(",100*summary(pca)$importance[2,2],"%)")) +
    theme_bw() +
    NULL
}

.sracipePlotTS <- function(plotData, genes = NULL, ...){
  sexprs <- stack(as.data.frame(plotData))
  colnames(sexprs) <- c("geneExp", "Gene")
  sexprs$time <- rep(as.numeric(rownames(plotData)), ncol(plotData))
  theme_set(theme_bw(base_size = 18))
  if(is.null(genes)) genes <- colnames(plotData)
  sexprs <- sexprs[which(sexprs$Gene %in% genes),,drop = FALSE ]
  ggplot2::qplot(time, geneExp, data = sexprs, group = Gene, colour = Gene,
                 geom = "line", ylab = "Gene Expression", xlab = "time" )
}



# Create a little question mark link that shows a help popup on hover
# https://github.com/daattali/ddpcr/blob/master/inst/shiny/ui/helpers.R

helpPopup <- function(content, title = NULL) {
  a(href = "#",
    class = "popover-link",
    `data-toggle` = "popover",
    `data-title` = title,
    `data-content` = content,
    `data-html` = "true",
    `data-trigger` = "hover",
    icon("question-circle")
  )
}

# https://github.com/daattali/advanced-shiny/blob/master/busy-indicator/helpers.R

withBusyIndicatorCSS <- "
.btn-loading-container {
margin-left: 10px;
font-size: 1.2em;
}
.btn-done-indicator {
color: green;
}
.btn-err {
margin-top: 10px;
color: red;
}
"

withBusyIndicatorUI <- function(button) {
  id <- button[['attribs']][['id']]
  div(
    shinyjs::useShinyjs(),
    singleton(tags$head(
      tags$style(withBusyIndicatorCSS)
    )),
    `data-for-btn` = id,
    button,
    span(
      class = "btn-loading-container",
      shinyjs::hidden(
        icon("spinner", class = "btn-loading-indicator fa-spin"),
        icon("check", class = "btn-done-indicator")
      )
    ),
    shinyjs::hidden(
      div(class = "btn-err",
          div(icon("exclamation-circle"),
              tags$b("Error: "),
              span(class = "btn-err-msg")
          )
      )
    )
  )
}

# Call this function from the server with the button id that is clicked and the
# expression to run when the button is clicked
withBusyIndicatorServer <- function(buttonId, expr) {
  # UX stuff: show the "busy" message, hide the other messages, disable the button
  loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", buttonId)
  doneEl <- sprintf("[data-for-btn=%s] .btn-done-indicator", buttonId)
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  shinyjs::disable(buttonId)
  shinyjs::show(selector = loadingEl)
  shinyjs::hide(selector = doneEl)
  shinyjs::hide(selector = errEl)
  on.exit({
    shinyjs::enable(buttonId)
    shinyjs::hide(selector = loadingEl)
  })
  
  # Try to run the code when the button is clicked and show an error message if
  # an error occurs or a success message if it completes
  tryCatch({
    value <- expr
    shinyjs::show(selector = doneEl)
    shinyjs::delay(2000, shinyjs::hide(selector = doneEl, anim = TRUE, 
                                       animType = "fade", time = 0.5))
    value
  }, error = function(err) { errorFunc(err, buttonId) })
}

# When an error happens after a button click, show the error
errorFunc <- function(err, buttonId) {
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
  shinyjs::html(html = errMessage, selector = errElMsg)
  shinyjs::show(selector = errEl, anim = TRUE, animType = "fade")
}


.sracipePlotCircuit <- function(.object, plotToFile = FALSE){
  topology <- sracipeCircuit(.object)
  node_list <-
    unique(c(topology[, 1], topology[, 2]))
  nodes <-
    data.frame(
      id = node_list,
      label = node_list,
      font.size = 50,
      value = c(rep(1, length(node_list)))
    )
  edge_col <- data.frame(c(1, 2), c("blue", "darkred"))
  arrow_type <- data.frame(c(1, 2), c("arrow", "circle"))
  colnames(arrow_type) <- c("type", "color")
  colnames(edge_col) <- c("type", "color")
  edges <-
    data.frame(
      from = c(topology[, 1]),
      to = c(topology[, 2]),
      arrows.to.type	= arrow_type$color[c(as.numeric(topology[, 3]))],
      width = 3,
      color = edge_col$color[c(as.numeric(topology[, 3]))]
    )
  network <-
    visNetwork::visNetwork(nodes, edges, height = "1000px", 
                           width = "100%", background = "#F8F8F8") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = FALSE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
  return(network)
  
}



############################### STICCC UTILITIES ###############################
.distfun=function(x) {
  as.dist((1-cor(t(x), method = 's'))/20)
}


.getClusters <- function(k = 2,                  # Number of clusters
                        matrix,                 # Input matrix - rows = samples, cols = features
                        features = NA,          # List of matrix colnames
                        clustfun = "ward.D2"    # Clustering function 
) {
  # Check that features are properly entered
  if(is.na(features)) {
    features <- colnames(matrix)
  } else if(any(!features %in% colnames(matrix))) {
    print("Error: specified features not in matrix colnames")
  }
  
  # Clustering
  hclust_dend <- hclust(.distfun(matrix[,which(colnames(matrix) %in% features)]), method = clustfun)
  clusters <- cutree(hclust_dend, k=k)
  return(clusters)
  
}




.prepMetadata <- function(exprMat,       # expression matrix - rows = features, cols = samples
                         cluster = T,   # Logical - whether to perform clustering
                         ...            # Parameters for getClusters function
) {
  metadata <- as.data.frame(t(exprMat))
  if(cluster) {
    clusters <- .getClusters(matrix = metadata, ...)
    metadata$Cluster <- clusters
  }
  metadata$SampleID <- rownames(metadata)
  return(metadata)
  
}



.computeGrid <- function(sce,
                        grid.length = 15
) {
  # Get posMat
  posMat <- reducedDim(sce, "PCA")
  
  # Define grid points, create metadata structure for later
  ## Create uniform grid
  # Set grid area (2-dimensional?) w uniform points
  if(exists("sce@metadata$params$xMin")) {
    x.min <- sce@metadata$params$xMin
    x.max <- sce@metadata$params$xMax
    y.min <- sce@metadata$params$yMin
    y.max <- sce@metadata$params$yMax
    
  } else {
    xMin <- floor(min(posMat[,1]))
    xMax <- ceiling(max(posMat[,1]))
    yMin <- floor(min(posMat[,2]))
    yMax <- ceiling(max(posMat[,2]))
    x.min <- xMin - 1
    x.max <- xMax + 1
    y.min <- yMin - 1
    y.max <- yMax + 1
  }
  
  
  x.points <- seq(x.min, x.max, length.out = grid.length)
  y.points <- seq(y.min, y.max, length.out = grid.length)
  grid.dist <- sqrt((x.points[2] - x.points[1])^2 + (y.points[2] - y.points[1])^2) / 2
  grid.x.dist <- x.points[2] - x.points[1]
  grid.y.dist <- y.points[2] - y.points[1]
  
  grid.df <- as.data.frame(expand_grid(x.points, y.points))
  grid.df$GridPoint <- 1:nrow(grid.df)
  grid.df$dx <- NA
  grid.df$dy <- NA
  
  sce@metadata$params$xMin <- x.min
  sce@metadata$params$xMax <- x.max
  sce@metadata$params$yMin <- y.min
  sce@metadata$params$yMax <- y.max
  
  sce@metadata$grid.df <- grid.df
  sce@metadata$params$grid.dist <- grid.dist
  sce@metadata$params$grid.x.dist <- grid.x.dist
  sce@metadata$params$grid.y.dist <- grid.y.dist
  
  return(sce)
}




.computeDist <- function(sce,
                        usePCs=TRUE,
                        numPCs=10) {
  
  if(usePCs) {
    PCAmat <- reducedDim(sce,"PCA")
    sce@metadata$dist_mat <- as.matrix(dist(PCAmat[,1:min(numPCs, ncol(PCAmat))], method = "euclidean"))
    sce@metadata$max_dist <- max(sce@metadata$dist_mat)
  } else {
    print("Error: this function is in development. Try using usePCs=TRUE for now")
  }
  
  return(sce)
}


.listToNumeric <- function(x) {
  return(as.numeric(unname(unlist(x))))
}

.DCComputeVector <- function(sce, query_point, useGinv=F, v2=T, bias=F) {
  ## Prepare results
  rs_list <- list("X"=NA,"Y"=NA,"fewNeighbors"=FALSE,"numNeighbors"=NA,"usedGinv"=FALSE,"nonInvertible"=FALSE,"selfCorNA"=FALSE,"selfCor"=NA,"b_vec"=NA,"b_vec_in"=NA,"det"=NA)
  
  ## extract variables from sce
  # Create sign vector
  topo <- sce@metadata$topo
  sign_vector <- revalue(factor(topo$Type), c('1'='1','2'='-1','3'='1','4'='-1','5'='1','6'='-1'), warn_missing = FALSE)
  sign_vector <- as.numeric(levels(sign_vector))[sign_vector]
  
  # Compute sampling radius
  sampling_radius <- sce@metadata$params$sample_radius * sce@metadata$max_dist
  
  # get posMat (position matrix - this should be either identical or closely correlated w/ the data used to compute dist_mat)
  if(is.null(sce@metadata$params$nPCs)) {
    nPCs <- nrow(sce)
    sce@metadata$params$nPCs <- nPCs
    
  } else {
    nPCs <- sce@metadata$params$nPCs
  }
  posMat <- reducedDim(sce, "PCA")[colnames(sce),1:nPCs]
  posList <- colnames(posMat)
  rownames(posMat) <- colnames(sce)
  
  # get exprMat (normalized expression matrix - high-dimensional data (rows=features, cols=samples))
  exprMat <- assay(sce,"normcounts")
  
  
  ## Get query point data
  query_data <- posMat[query_point,]
  rs_list[["X"]] <- query_data[1]
  rs_list[["Y"]] <- query_data[2]
  
  ## Find neighbors
  neighbors <- which(sce@metadata$dist_mat[query_point,colnames(sce)] <= sampling_radius)
  #neighbors <- neighbors[which(neighbors %in% colnames(sce))]
  
  ## Skip sample if number of neighbors too small
  if(length(neighbors) <= (sce@metadata$params$minNeighbors+1)) {
    rs_list[["fewNeighbors"]] <- TRUE
    rs_list[["numNeighbors"]] <- length(neighbors)-1
    return(rs_list)
  }
  
  ## Select neighbors
  subset_models <- as.data.frame(posMat[neighbors,])
  subset_models <- subset_models[-which(rownames(subset_models) %in% query_point),]
  rs_list[["numNeighbors"]] <- nrow(subset_models)
  
  ## Multiple regression
  # Create relative positional matrix of neighbor samples
  subset_models[,1:ncol(subset_models)] <- apply(subset_models[,1:ncol(subset_models)], 2, .listToNumeric)
  deriv_df <- sweep(subset_models, 2, as.numeric(query_data), "-")
  
  # get relative expression matrix & transpose (neighboring cells only)
  geneex_mat <- as.matrix(deriv_df[,posList])
  geneex_mat_transposed <- t(geneex_mat)
  
  # Multiply relative positional matrix by its transpose
  mat_product <- geneex_mat_transposed %*% geneex_mat
  rs_list[["det"]] <- det(mat_product)
  
  # Take the inverse if possible - if not, skip this sample
  if(useGinv) {
    tol = 10^-6
    if(rcond(mat_product) < tol){
      inverse_mat = ginv(mat_product)
      rs_list[["usedGinv"]] <- TRUE
      
    } else{
      inverse_mat <- tryCatch({
        solve(mat_product)
      }, error = function(e) {
        "non-invertible"
      })
    }
  } else {
    inverse_mat <- tryCatch({
      solve(mat_product)
    }, error = function(e) {
      "non-invertible"
    })
  }
  
  
  if(inverse_mat[1] == "non-invertible") {
    rs_list[["nonInvertible"]] <- TRUE
    return(rs_list)
  }
  
  # Multiply the inverted matrix by the transposed positional matrix
  x_mat <- inverse_mat %*% geneex_mat_transposed
  
  
  ## Calculate delayed correlation
  # Self correlation - skip this sample if NA
  src_activity <- as.numeric(exprMat[topo$Source,query_point] * sign_vector)
  self_target_activity <- as.numeric(exprMat[topo$Target,query_point]) 
  self_cor <- cor(self_target_activity, src_activity)
  if(is.na(self_cor)) {
    rs_list[["selfCorNA"]] <- TRUE
    return(rs_list)
  } 
  rs_list[["selfCor"]] <- self_cor
  
  # For each neighbor:
  deriv_df$DelayedCorr <- NA
  for(model in rownames(deriv_df)) {
    # Correlate neighbor target activity with query source activity * sign_vector
    neighbor_target_activity <- as.numeric(exprMat[topo$Target, model]) 
    if(sd(neighbor_target_activity) > 0 & sd(src_activity)) {
      neighbor_cor <- cor(neighbor_target_activity,src_activity) 
      deriv_df[model,"DelayedCorr"] <- (neighbor_cor - self_cor)
    } else {
      #print(paste("Issue with correlation in models ",query_point, " and ", model))
      deriv_df[model,"DelayedCorr"] <- NA
    }
    
    ## QQ REMOVE THIS
    #plot(src_activity, neighbor_target_activity)
  }
  
  # Remove NA values caused by zero variance vectors
  na_models <- which(is.na(deriv_df$DelayedCorr))
  if(length(na_models) > 0) {
    deriv_df <- deriv_df[-na_models,]
    subset_models <- subset_models[-na_models,]
    x_mat <- x_mat[,-na_models]
  }
  
  
  # Compute product of earlier matrix product and delayed corr vector
  corr_vec <- as.matrix(deriv_df[,"DelayedCorr"])
  b_vec <- x_mat %*% corr_vec
  rs_list[["b_vec"]] <- b_vec
  
  
  if(bias) {
    ## Split data into train/test 80:20
    train_idx <- sample(nrow(deriv_df), round(nrow(deriv_df)*0.8))
    df.train <- deriv_df[train_idx,]
    df.test <- deriv_df[-train_idx,]
    
    alpha.guess <- 0
    beta.guess <- rep(0,nPCs)
    
    h <- 100
    error <- rep(NA,h)
    lambda <- seq(0,1,length.out = h)
    
    fit <- lm(DelayedCorr~.-1, data=df.train)
    yhat.test <- predict(fit, newdata = df.test)
    yhat.guess <- alpha.guess + as.matrix(df.test[,-ncol(df.test)]) %*% beta.guess
    
    error <- sapply(lambda, function(j) sum(((1-j)*yhat.test + j*yhat.guess - df.test$DelayedCorr)^2))
    rs_list[["error"]] <- error
  }
  
  
  
  ## Compute v2 as well, if specified
  if(v2) {
    # For each neighbor:
    deriv_df$DelayedCorr_in <- NA
    for(model in rownames(deriv_df)) {
      # Correlate neighbor source activity * sign_vector with query target activity 
      neighbor_src_activity <- as.numeric(exprMat[topo$Source, model]) * sign_vector
      neighbor_cor <- cor(neighbor_src_activity,self_target_activity)
      deriv_df[model,"DelayedCorr_in"] <- (neighbor_cor - self_cor)
    }
    
    # Remove NA values caused by zero variance vectors
    na_models <- which(is.na(deriv_df$DelayedCorr))
    if(length(na_models) > 0) {
      deriv_df <- deriv_df[-na_models,]
      subset_models <- subset_models[-na_models,]
      x_mat <- x_mat[,-na_models]
    }
    
    
    # Compute product of earlier matrix product and delayed corr vector
    corr_vec <- as.matrix(deriv_df[,"DelayedCorr_in"])
    b_vec_in <- x_mat %*% corr_vec
    rs_list[["b_vec_in"]] <- b_vec_in
    
  }
  
  
  
  return(rs_list)
}




.DCComputeTrajectorySCE_2022 <- function(sce, v2=TRUE, useGinv=FALSE, bias=F, subset=NA) {
  
  # Create sign vector
  topo <- sce@metadata$topo
  sign_vector <- revalue(factor(topo$Type), c('1'='1','2'='-1','3'='1','4'='-1','5'='1','6'='-1'), warn_missing = FALSE)
  sign_vector <- as.numeric(levels(sign_vector))[sign_vector]
  
  # Compute sampling radius
  sampling_radius <- sce@metadata$params$sample_radius * sce@metadata$max_dist
  
  # get posMat (position matrix - this should be either identical or closely correlated w/ the data used to compute dist_mat)
  if(is.null(sce@metadata$params$nPCs)) {
    nPCs <- nrow(sce)
    sce@metadata$params$nPCs <- nPCs
    
  } else {
    nPCs <- sce@metadata$params$nPCs
  }
  posMat <- reducedDim(sce, "PCA")[colnames(sce),1:nPCs]
  #posMat <- posMat[,1:nPCs]
  posList <- colnames(posMat)
  rownames(posMat) <- colnames(sce)
  
  # get exprMat (normalized expression matrix - high-dimensional data (rows=features, cols=samples))
  exprMat <- assay(sce,"normcounts")
  
  # Begin velocity calculation
  if(sce@metadata$params$verbose) {
    print("Computing inferred velocity")
  }
  
  sample_list <- colData(sce)$SampleID
  colData(sce)$numNeighbors <- NA
  colData(sce)$selfCor <- NA
  colData(sce)$X <- NA
  colData(sce)$Y <- NA
  colData(sce)$dX <- NA
  colData(sce)$dY <- NA
  colData(sce)$fewNeighbors <- FALSE
  colData(sce)$nonInvertible <- FALSE
  colData(sce)$selfCorNA <- FALSE
  colData(sce)$usedGinv <- FALSE
  sce@metadata$det_list <- list()
  
  
  
  if(bias) {
    error <- matrix(NA, length(sample_list), 100)
  }
  
  numSamples <- length(sample_list)
  if(!is.na(subset)) {
    #sample_list <- sample_list[which(sample_list %in% subset)] # subset by IDs, not index!
    numSamples <- subset
  }
  
  weighted_vector_list <- vector(mode = "list", length = length(sample_list))
  weighted_vector_list_in <- vector(mode = "list", length = length(sample_list))
  det_list <- list()
  pb = txtProgressBar(min = 0, max = numSamples, 
                      initial = 0, style = 3)
  stepi <- 0
  for(query_point in sample_list) {
    
    if(stepi > numSamples) {
      break
    }
    
    ## Compute vector
    rs_list <- .DCComputeVector(sce, query_point, v2 = v2, useGinv = useGinv, bias = bias)
    
    
    #"X"=NA,"Y"=NA,"fewNeighbors"=FALSE,"numNeighbors"=NA,"usedGinv"=FALSE,"nonInvertible"=FALSE,"selfCorNA"=FALSE,"selfCor"=NA,"b_vec"=NA
    # Store data
    colData(sce)[query_point,"X"] <- rs_list[["X"]]
    colData(sce)[query_point,"Y"] <- rs_list[["Y"]]
    colData(sce)[query_point,"fewNeighbors"] <- rs_list[["fewNeighbors"]]
    colData(sce)[query_point,"numNeighbors"] <- rs_list[["numNeighbors"]]
    colData(sce)[query_point,"usedGinv"] <- rs_list[["usedGinv"]]
    colData(sce)[query_point,"nonInvertible"] <- rs_list[["nonInvertible"]]
    colData(sce)[query_point,"selfCorNA"] <- rs_list[["selfCorNA"]]
    colData(sce)[query_point,"selfCor"] <- rs_list[["selfCor"]]
    sce@metadata$det_list[[query_point]] <- rs_list[["det_list"]]
    
    if(length(rs_list[["b_vec"]]) == 1) {
      next
    } 
    
    
    ## Update progress bar
    stepi <- stepi+1
    setTxtProgressBar(pb, stepi)
    
    
    if(bias) {
      error[stepi,] <- rs_list[["error"]]
    }
    
    
    ## Save vector to master list
    saved_vector <- as.data.frame(t(rs_list[["b_vec"]]))
    colnames(saved_vector) <- paste0("d",posList)
    weighted_vector_list[[stepi]] <- saved_vector
    names(weighted_vector_list)[[stepi]] <- query_point
    
    ## Update colData
    colData(sce)[query_point,"dX"] <- saved_vector[1]
    colData(sce)[query_point,"dY"] <- saved_vector[2]
    
    
    ## Add in v2 if necessary
    if(v2) {
      ## Save vector to master list
      saved_vector_in <- as.data.frame(t(rs_list[["b_vec_in"]]))
      colnames(saved_vector_in) <- paste0("d",posList)
      weighted_vector_list_in[[stepi]] <- saved_vector_in
      names(weighted_vector_list_in)[[stepi]] <- query_point
      
      ## Update colData
      colData(sce)[query_point,"dX_in"] <- saved_vector_in[1]
      colData(sce)[query_point,"dY_in"] <- saved_vector_in[2]
    }
    
  }
  
  ## Consolidate vectors and write to file
  weighted_vectors <- do.call(rbind, weighted_vector_list)
  sce@metadata$vectors <- weighted_vectors
  
  if(v2) {
    weighted_vectors_in <- do.call(rbind, weighted_vector_list_in)
    sce@metadata$vectors_in <- weighted_vectors_in
  }
  
  
  return(sce)
  
}




.DCComputeGridVectors <- function(sce,
                                 scalingFactor = NA,
                                 unitVectors = TRUE,
                                 inVectors = FALSE,
                                 combine = FALSE,
                                 how = NA) {
  if(is.na(scalingFactor)) {
    scalingFactor <- sce@metadata$params$gridPlotScalingFactor
  }
  
  # Get gridpoints from metadata
  plot_df <- as.data.frame(colData(sce))
  grid.df <- sce@metadata$grid.df
  grid.df$numNeighbors <- NA
  grid.dist <- sce@metadata$params$grid.dist
  x.dist <- sce@metadata$params$grid.x.dist
  y.dist <- sce@metadata$params$grid.y.dist
  arrow.max <- grid.dist * 0.5
  
  ## Iterate thru each gridpoint & calculate weighted k-means average vector
  for(point in grid.df$GridPoint) {
    
    
    ## Find neighbors
    query <- grid.df[which(grid.df$GridPoint == point),c(1:2)]
    colnames(query) <- c("X","Y")
    neighbors <- which(plot_df$X <= (query$X + x.dist) & 
                         plot_df$X >= (query$X - x.dist) &
                         plot_df$Y <= (query$Y + y.dist) & 
                         plot_df$Y >= (query$Y - y.dist))
    # Decide which vector (or sum) to be plotted
    if(inVectors) {
      subset_data <- plot_df[neighbors,c("X","Y","dX_in","dY_in")]
    } else {
      
      if(combine) {
        
        if(how == "v1+v2") {
          subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
          subset_data$dX = subset_data$dX + plot_df$dX_in[neighbors]
          subset_data$dY = subset_data$dY + plot_df$dY_in[neighbors]
        } else if(how == "v1-v2") {
          subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
          subset_data$dX = subset_data$dX - plot_df$dX_in[neighbors]
          subset_data$dY = subset_data$dY - plot_df$dY_in[neighbors]
        } else if(how == "avg+") {
          subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
          subset_data$dX = (subset_data$dX + plot_df$dX_in[neighbors]) / 2
          subset_data$dY = (subset_data$dY + plot_df$dY_in[neighbors]) / 2
        } else if(how == "avg-") {
          subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
          subset_data$dX = (subset_data$dX - plot_df$dX_in[neighbors]) / 2
          subset_data$dY = (subset_data$dY - plot_df$dY_in[neighbors]) / 2
        } else if(how == "v1") {
          subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
          subset_data$dX = subset_data$dX
          subset_data$dY = subset_data$dY
        } else if(how == "v2") {
          subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
          subset_data$dX = plot_df$dX_in[neighbors]
          subset_data$dY = plot_df$dY_in[neighbors]
        } else {
          print("Error: must specify `how`, if `combine` is set to TRUE")
        }
        
        
      } else {
        subset_data <- plot_df[neighbors,c("X","Y","dX","dY")]
      }
      
    }
    
    subset_data <- na.omit(subset_data)
    grid.df[point,"numNeighbors"] <- nrow(subset_data)
    
    if(nrow(subset_data) >= 3) {
      
    } else if(nrow(subset_data) < 3) {
      ## Add to grid.df as 0
      grid.df[which(grid.df$GridPoint == point),"dx"] <- 0
      grid.df[which(grid.df$GridPoint == point),"dy"] <- 0
      next
    }
    
    ## Calculate distance & proximity (1/dist) to each neighbor
    neighbor.dists <- as.matrix(dist(rbind(query,subset_data[,c(1:2)]), 
                                     method = "euclidean"))
    neighbor.dists <- neighbor.dists[which(rownames(neighbor.dists) == point),-which(colnames(neighbor.dists) == point)]
    subset_data$Dist <- neighbor.dists
    subset_data$Proximity <- 1 / subset_data$Dist
    
    
    
    ## Calculate weighted vector (closer = more weight)
    avg_dx <- .Avg_IDW(subset_data$dX, subset_data$Proximity)
    avg_dy <- .Avg_IDW(subset_data$dY, subset_data$Proximity)
    #avg_dx <- sum(subset_data$dX * subset_data$Proximity) / sum(subset_data$Proximity)
    #avg_dy <- sum(subset_data$dY * subset_data$Proximity) / sum(subset_data$Proximity)
    
    ## Add to grid.df
    grid.df[which(grid.df$GridPoint == point),"dx"] <- avg_dx
    grid.df[which(grid.df$GridPoint == point),"dy"] <- avg_dy
  }
  
  
  if(unitVectors) {
    scaleR <- function(x) {x / sqrt(sum(x^2))}
    scaled_vecs <-  as.data.frame(t(apply(grid.df[,c("dx","dy")], 1, scaleR)))
    grid.df$dx = scaled_vecs[,1] * grid.dist * 0.7
    grid.df$dy = scaled_vecs[,2] * grid.dist * 0.7
    grid.df$Magnitude <- sqrt(grid.df$dx^2 + grid.df$dy^2)
    
  } else {
    ## Scale all vectors to the target mean
    grid.df$Magnitude <- sqrt(grid.df$dx^2 + grid.df$dy^2)
    #mean.magnitude <- mean(grid.df$Magnitude[which(grid.df$Magnitude > 0)])
    #scaling.factor <- grid.dist / mean.magnitude
    
    #grid.df$dx <- grid.df$dx * scaling.factor * scalingFactor
    #grid.df$dy <- grid.df$dy * scaling.factor * scalingFactor
    #grid.df$Magnitude <- grid.df$Magnitude * scaling.factor * scalingFactor
  }
  
  
  ## Write to object
  sce@metadata$grid.df <- grid.df
  
  return(sce)
}

.Avg_IDW <- function(vals, weights) {
  
  return(sum(vals * weights) / sum(weights))
  # avg_dx <- sum(subset_data$dX * subset_data$Proximity) / sum(subset_data$Proximity)
  # avg_dy <- sum(subset_data$dY * subset_data$Proximity) / sum(subset_data$Proximity)
}



.DCPlotGrid <- function(sce,
                       colorVar,
                       plotLoadings = F,
                       loadingFactor = 3.5,
                       colorPalette = NA,
                       outputDir = NA,
                       scalingFactor = 1,
                       minMagnitude = 0.01,
                       minimal = F,
                       arrowheadSize = 0.3,
                       arrowSize = 1,
                       legend=T,
                       pca = NA) {
  
  topoName <- sce@metadata$topoName
  expName <- sce@metadata$experimentName
  
  # Get 2D plotting coordinates of samples from first 2 cols of position matrix
  if(all(is.na(pca))) {
    plot_df <- as.data.frame(colData(sce))  
  } else {
    plot_df <- pca
    plot_df <- as.data.frame(cbind(plot_df, colData(sce)[,colorVar]))
    colnames(plot_df)[1:2] <- c("X","Y")
    colnames(plot_df)[ncol(plot_df)] <- colorVar
  }
  
  if(!colorVar %in% colnames(plot_df)) {
    print("Error: colorVar not found in colnames of colData")
  }
  
  xMin <- sce@metadata$params$xMin
  xMax <- sce@metadata$params$xMax
  yMin <- sce@metadata$params$yMin
  yMax <- sce@metadata$params$yMax
  
  # Get proportion of variance for PCs 1 and 2 for axis labels
  if(!is.null(sce@metadata$pca_summary)) {
    pc1_weight <- round(100*sce@metadata$pca_summary$importance[2,1],2)
    pc2_weight <- round(100*sce@metadata$pca_summary$importance[2,2],2)
    plot_xlab <- paste("PC1 (",pc1_weight,"%)",sep="")
    plot_ylab <- paste("PC2 (",pc2_weight,"%)",sep="")  
  } else {
    plot_xlab <- "PC1"
    plot_ylab <- "PC2"
  }
  
  
  # get loadings if applicable
  if(plotLoadings) {
    loadingDF <- as.data.frame(sce@metadata$pca_data$rotation)
    loadingDF$gene <- rownames(loadingDF)
    loadingDF$X <- 0
    loadingDF$Y <- 0
  }
  
  
  # Add relevant metadata to plotting df
  grid.df <- sce@metadata$grid.df
  x.points <- unique(sce@metadata$grid.df$x.points)
  y.points <- unique(sce@metadata$grid.df$y.points)
  plot_df$colorVar <- plot_df[,colorVar]
  
  if(colorVar == "Cluster") {
    plot_df$colorVar <- factor(plot_df$colorVar)
  }
  
  ## Plot grid points
  image <- ggplot(grid.df[which(grid.df$Magnitude > minMagnitude),], aes(x=x.points,y=y.points)) +
    geom_point(data=plot_df[,],mapping=aes(x=X,y=Y, color=colorVar)) + 
    geom_point() + 
    xlab(plot_xlab) +
    ylab(plot_ylab) +
    scale_size(range=c(1.0, 3)) +
    xlim(xMin,xMax) +
    ylim(yMin,yMax) +
    #scale_color_manual(values=c(cbPalette[2:8])) +
    geom_segment(aes(xend=x.points+dx*scalingFactor, yend=y.points+dy*scalingFactor), arrow = arrow(length = unit(arrowheadSize,"cm")), linewidth=arrowSize) +
    guides(alpha="none", size="none", color=guide_legend(title = colorVar, override.aes = list(size = 5))) +
    theme(axis.text = element_text(size=28), axis.title = element_text(size=36)) 
  if(all(!check.numeric(plot_df$colorVar)) | colorVar == "Cluster") {
    if(all(is.na(colorPalette))) {
      colorPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
    }
    image <- image + scale_color_manual(values=c(colorPalette))
  } else if(length(colorPalette) > 1) {
    image <- image + scale_color_manual(values=c(colorPalette))
  } else if(all(check.numeric(plot_df$colorVar)) | colorVar == "Magnitude") {
    image <- image + scale_color_gradient(low="grey",high="red")
  }
  if(legend==F) {
    image <- image+guides(color="none")
  }
  if(plotLoadings) {
    loadingLabelFactor <- loadingFactor+0.5
    image <- image + geom_segment(data = loadingDF[,],aes(xend=PC1*loadingFactor, yend=PC2*loadingFactor, x = X, y = Y), arrow = arrow(length = unit(0.6,"cm"))) +
      geom_text(data = loadingDF[,], aes(x = PC1*loadingLabelFactor, y = PC2*loadingLabelFactor, label=gene), size=8, fontface="bold")
  }
  
  if(minimal) {
    image <- image + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
  }
  
  return(image)
}



.DCPlotVectors <- function(sce, 
                          colorVar,
                          plotLoadings = F,
                          loadingFactor = 3.5,
                          colorPalette = NA,
                          scalingFactor = NA,
                          pca = NA,
                          type = "v1") {
  
  
  topoName <- sce@metadata$topoName
  expName <- sce@metadata$experimentName
  
  
  # Get 2D plotting coordinates of samples from first 2 cols of position matrix
  if(all(is.na(pca))) {
    plot_df <- as.data.frame(colData(sce))  
  } else {
    plot_df <- pca
    plot_df <- as.data.frame(cbind(plot_df, colData(sce)[,colorVar]))
    colnames(plot_df)[1:2] <- c("X","Y")
    colnames(plot_df)[ncol(plot_df)] <- colorVar
    
    
    if(type == "v1") {
      plot_df$dX <- colData(sce)$dX[which(colData(sce)$SampleID %in% rownames(plot_df))]
      plot_df$dY <- colData(sce)$dY[which(colData(sce)$SampleID %in% rownames(plot_df))]
    }
    if(type == "v2") {
      plot_df$dX <- colData(sce)$dX_in[which(colData(sce)$SampleID %in% rownames(plot_df))]
      plot_df$dY <- colData(sce)$dY_in[which(colData(sce)$SampleID %in% rownames(plot_df))]
    }
    if(type == "avg+") {
      plot_df$dX <- colData(sce)$dX[which(colData(sce)$SampleID %in% rownames(plot_df))] +
        colData(sce)$dX_in[which(colData(sce)$SampleID %in% rownames(plot_df))] / 2
      plot_df$dY <- colData(sce)$dY[which(colData(sce)$SampleID %in% rownames(plot_df))] + 
        colData(sce)$dY_in[which(colData(sce)$SampleID %in% rownames(plot_df))] / 2
      
    }
    if(type == "avg-") {
      plot_df$dX <- colData(sce)$dX[which(colData(sce)$SampleID %in% rownames(plot_df))] -
        colData(sce)$dX_in[which(colData(sce)$SampleID %in% rownames(plot_df))] / 2
      plot_df$dY <- colData(sce)$dY[which(colData(sce)$SampleID %in% rownames(plot_df))] - 
        colData(sce)$dY_in[which(colData(sce)$SampleID %in% rownames(plot_df))] / 2
    }
    
    
  }
  if(!colorVar %in% colnames(plot_df)) {
    print("Error: colorVar not found in colnames of colData")
  }
  
  xMin <- sce@metadata$params$xMin
  xMax <- sce@metadata$params$xMax
  yMin <- sce@metadata$params$yMin
  yMax <- sce@metadata$params$yMax
  
  # Get proportion of variance for PCs 1 and 2 for axis labels
  pc1_weight <- round(100*sce@metadata$pca_summary$importance[2,1],2)
  pc2_weight <- round(100*sce@metadata$pca_summary$importance[2,2],2)
  plot_xlab <- paste("PC1 (",pc1_weight,"%)",sep="")
  plot_ylab <- paste("PC2 (",pc2_weight,"%)",sep="")
  
  # get loadings if applicable
  if(plotLoadings) {
    loadingDF <- as.data.frame(sce@metadata$pca_data$rotation)
    loadingDF$gene <- rownames(loadingDF)
    loadingDF$X <- 0
    loadingDF$Y <- 0
  }
  
  
  # Add relevant metadata to plotting df
  plot_df$colorVar <- plot_df[,colorVar]
  
  if(colorVar == "Cluster") {
    plot_df$colorVar <- factor(plot_df$colorVar)
  }
  
  
  
  image <- ggplot(plot_df, aes(x=X,y=Y)) +
    geom_point(mapping=aes(size=1, color=colorVar)) + 
    xlab(plot_xlab) +
    ylab(plot_ylab) +
    scale_size(range=c(1.75, 3)) +
    guides(alpha="none", size="none", color=guide_legend(title = colorVar, override.aes = list(size = 5))) +
    xlim(xMin,xMax) +
    ylim(yMin,yMax) +
    ggtitle(paste0("PCA on ",expName," cells"))
  if(all(!check.numeric(plot_df$colorVar)) | colorVar == "Cluster") {
    if(all(is.na(colorPalette))) {
      colorPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
    }
    image <- image + scale_color_manual(values=c(colorPalette))
  } else if(all(check.numeric(plot_df$colorVar)) | colorVar == "Magnitude") {
    image <- image + scale_color_gradient(low="grey",high="red")
  }
  if(plotLoadings) {
    loadingLabelFactor <- loadingFactor+0.5
    image <- image + geom_segment(data = loadingDF[,],aes(xend=PC1*loadingFactor, yend=PC2*loadingFactor, x = X, y = Y), arrow = arrow(length = unit(0.6,"cm"))) +
      geom_text(data = loadingDF[,], aes(x = PC1*loadingLabelFactor, y = PC2*loadingLabelFactor, label=gene), size=8, fontface="bold")
  }
  
  
  # Return plot with vectors
  image <- image + geom_segment(aes(xend=X+dX*scalingFactor, yend=Y+dY*scalingFactor), 
                             arrow = arrow(length = unit(0.2,"cm")))
  return(image)
  
}

.plotTrajectories <- function(vic,
                              plotLoadings = T,
                              scalingFactor = 30,
                              loadingFactor = 3.5,
                              gridSmoothing = T
                              ) {
  
  
  # Extract pca
  pca <- reducedDim(vic, "PCA")
  
  # invert v2 for interpretability
  # Multiply in vectors by -1
  colData(vic)$dX_in <- -1 * colData(vic)$dX_in
  colData(vic)$dY_in <- -1 * colData(vic)$dY_in
  
  # Plot results
  minMagnitude <- 0.001
  arrowheadSize <- 0.4
  
  types <- c("v1","v2","avg+","avg-")
  titles <- c("V1 (Outward)", "V2 (Inward)", "Net Flow (v1+v2)", "Reversibility (v1-v2)")
  plotList <- list()
  
  for(typeNo in seq_along(types)) {
    type <- types[typeNo]
    if(gridSmoothing) {
      
      vic <- .computeGrid(vic)
      vic <- .DCComputeGridVectors(vic, inVectors = F, combine = T, unitVectors = F, how=type)
      
      plot <- .DCPlotGrid(sce = vic,
                 plotLoadings = plotLoadings,
                 loadingFactor = loadingFactor,
                 minMagnitude = minMagnitude, 
                 scalingFactor = scalingFactor, 
                 colorVar = "Cluster", 
                 minimal = F,
                 arrowheadSize = arrowheadSize,
                 pca = pca)
      plotList[[typeNo]] <- plot
      
    } else {
      
      plot <- .DCPlotVectors(sce = vic,
                          plotLoadings = plotLoadings,
                          loadingFactor = loadingFactor,
                          scalingFactor = scalingFactor, 
                          colorVar = "Cluster",
                          pca = pca,
                          type = type)
      plotList[[typeNo]] <- plot
      
    }
    
    
  }
  
  
  ##### PLOTTING OPTION 1 #####
  # Redundant axis titles, plot titles, and legends
  # Now that we have a list of all four plots, we can arrange them
  # grid.arrange(
  #   grobs = lapply(1:4, function(i) {
  #     gridExtra::arrangeGrob(plotList[[i]], #+ ggplot2::ggtitle(titles[i]),
  #                            top = titles[i])
  #   }),
  #   ncol = 2
  # )
  
  ##### PLOTTING OPTION 2 #####
  # Broken in that it just removes axis titles, but doesn't replace them at the higher level
  # # Remove axis titles and legends from individual plots
  # for (i in 1:length(plotList)) {
  #   plotList[[i]] <- plotList[[i]] + theme(
  #     legend.position = "none",  # Remove legend
  #     axis.title.x = element_blank(),  # Remove x-axis title
  #     axis.title.y = element_blank()   # Remove y-axis title
  #   )
  # }
  # 
  # # Combine the plots into a 2x2 grid and add shared legends and axis titles
  # combined_plot <- (plotList[[1]] + plotList[[2]]) / (plotList[[3]] + plotList[[4]]) + 
  #   plot_layout(guides = 'collect') +
  #   plot_annotation(
  #     title = 'Circuit Transition Predictions',
  #     caption = 'TEST Caption',
  #     theme = theme(
  #       plot.title = element_text(size = 14),
  #       plot.caption = element_text(size = 8)
  #     )
  #   )
  # 
  # # Print the combined plot
  # print(combined_plot)
  
  
  ##### PLOTTING OPTION 3 #####
  library(gridExtra)
  library(cowplot)
  library(ggplot2)
  
  
  # Extract the legend from the fourth plot
  legend <- get_legend(plotList[[4]])
  
  # Remove legends and axis titles from all plots and add individual titles
  for (i in 1:4) {
    plotList[[i]] <- plotList[[i]] +
      ggplot2::theme(legend.position = "none",
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank()) +
      ggtitle(titles[i])
  }
  
  # Align the plots in a 2x2 grid without the legend
  plot_grid_combined <- plot_grid(plotlist = plotList, ncol = 2, align = 'hv')
  
  # Combine the grid of plots with the legend
  final_plot <- plot_grid(plot_grid_combined, legend, ncol = 2, rel_widths = c(1, 0.2))
  
  # Create labels for the common axis titles
  pc1_weight <- round(100*vic@metadata$pca_summary$importance[2,1],2)
  pc2_weight <- round(100*vic@metadata$pca_summary$importance[2,2],2)
  plot_xlab <- paste("PC1 (",pc1_weight,"%)",sep="")
  plot_ylab <- paste("PC2 (",pc2_weight,"%)",sep="")  
  x_title <- cowplot::draw_label(plot_xlab, fontface = 'bold')
  y_title <- cowplot::draw_label(plot_ylab, fontface = 'bold', angle = 90)
  
  # Arrange the common axis titles and the final plot
  final_layout <- cowplot::plot_grid(
    cowplot::plot_grid(NULL, y_title, NULL, nrow = 3, rel_heights = c(1, 10, 1)),
    cowplot::plot_grid(NULL, final_plot, x_title, nrow = 3, rel_heights = c(1, 10, 1)),
    ncol = 3,
    rel_widths = c(1, 10, 1)
  )
  
  # Print the final layout
  print(final_layout)
  
}
