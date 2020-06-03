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
    visOptions(manipulation = FALSE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
  return(network)
  
}

