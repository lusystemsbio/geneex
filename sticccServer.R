###########################################
# STICCC Server
###########################################


sticRacipe <- RacipeSE()
sticVals = reactiveValues()
sracipeCircuit(sticRacipe) <- demoCircuit
sticVals$rSet <- sticRacipe


sticVals$simExp <- NULL
sticVals$refClust <- NULL
sticVals$annealFlag <- FALSE
sticVals$filterModels <- NULL
sticVals$stochData <- NULL

output$sticCircuit <- renderVisNetwork({
  .sracipePlotCircuit(sticVals$rSet, plotToFile = FALSE)
})



observeEvent(input$computeTrajectories, {
  withBusyIndicatorServer("computeTrajectories",{
    shinyjs::show("sticText")
    shinyjs::show("sticPCA")
    
    sticRacipe <- sticVals$rsRacipe
    
    output$sticText <- renderUI({HTML(
      "Transition vectors overlaid on principal component analysis. Panels include V1 (outward), V2 (inward), Net Flow, and Reversibility.")})
    
    
    output$sticPCA <- renderPlot({
      if(input$computeTrajectories == 0) return()
      withProgress(message = 'Plotting pca', value = 0.25, {
        
        
        ## Go from racipe object, to vic object and compute trajectories from it
        # Get racipe data (done)
        # transpose for compatibility with VIC constructor
        normData <- assay(sticRacipe())
        
        # create SCE object
        vic <- SingleCellExperiment(assays = SimpleList(normcounts=normData))
        vic@metadata$experimentName <- "Geneex_STICCC"
        vic@metadata$topoName <- "Geneex_STICCC_Circuit"
        vic@metadata$topo <- sracipeCircuit(sticRacipe())
        vic@metadata$params <- list(sample_radius=isolate(input$samplingRadius), plotScalingFactor=1, gridPlotScalingFactor=1, minNeighbors=15, verbose=T)
        
        # add metadata
        metadata <- .prepMetadata(normData, cluster = T)
        rownames(metadata) <- metadata$SampleID
        colData(vic) <- DataFrame(metadata)
        colnames(vic) <- colData(vic)$SampleID
        
        # add PCA to SCE object
        print("racipeVals$pca:")
        print(racipeVals$pca)
        sticPCA <- racipeVals$pca
        reducedDim(vic, "PCA") <- sticPCA$x
        vic@metadata$params$nPCs <- ncol(sticPCA$x)
        vic@metadata$pca_data <- sticPCA[1:4]
        vic@metadata$pca_summary <- summary(sticPCA)
        
        # set borders?
        borders <- list(xmin=floor(min(sticPCA$x[,1]))-0.5,
                        xmax=ceiling(min(sticPCA$x[,1]))+0.5,
                        ymin=floor(min(sticPCA$x[,2]))-0.5,
                        ymax=ceiling(min(sticPCA$x[,2]))+0.5)
        vic@metadata$params$xMin <- borders[['xmin']]
        vic@metadata$params$xMax <- borders[['xmax']]
        vic@metadata$params$yMin <- borders[['ymin']]
        vic@metadata$params$yMax <- borders[['ymax']]
        
        # compute grid based on PCA
        vic <- .computeGrid(vic)
        
        # compute pairwise distance between points
        vic <- .computeDist(vic, numPCs=vic@metadata$params$nPCs)

        # Compute trajectories
        vic <- .DCComputeTrajectorySCE_2022(vic, v2=T)
        
        # Plot resulting vectors with a function that only takes vic as input
        .plotTrajectories(vic = vic,
                          plotLoadings = isolate(input$plotLoadings),
                          scalingFactor = 3,
                          loadingFactor = 3.5,
                          gridSmoothing = isolate(input$plotLoadings))
        

        

      })
    })
  })
})







