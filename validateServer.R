###########################################
# Validate Server
###########################################
validateVars = reactiveValues()
validateVars$simExp <- NULL
validateVars$refClust <- NULL

output$downloadSampleValidateRefData <- downloadHandler(
  filename <- function() {
    paste("sampleValidateRefData", "csv", sep=".")
  },
  
  content <- function(file) {
    file.copy("data/validateDataRef.csv", file)
  },
  contentType = "text/csv"
)

output$downloadSampleValidateSimData <- downloadHandler(
  filename <- function() {
    paste("sampleValidateSimData", "csv", sep=".")
  },
  
  content <- function(file) {
    file.copy("data/validateDataSim.csv", file)
  },
  contentType = "text/csv"
)


# Update the network if [Load Circuit] button is clicked from [Load File] tab
observeEvent(input$fileRefExp,{
  refExp <-  input$fileRefExp
  if(is.null(refExp)) {
    return(NULL)
  }
  refExp <- as.matrix(read.table(
    refExp$datapath, sep = ",", header = TRUE,row.names = 1))
  validateVars$refExp <- reactive(refExp)
})

observeEvent(input$fileSimExp,{
  simExp <-  input$fileSimExp
  if(is.null(simExp)) {
    return(NULL)
  }
  simExp <- as.matrix(read.table(
    simExp$datapath, sep = ",", header = TRUE,row.names = 1))
  validateVars$simExp <- reactive(simExp)
})

# Update the network if [Load Circuit] button is clicked from [Load File] tab
observeEvent(input$uploadRefClust,{
  refClust <-  input$fileRefClust
  if(is.null(refClust)) {
    return(NULL)
  }
  refClust <- as.character(read.table(refExp$datapath,sep="", header =FALSE,
                                      stringsAsFactors = FALSE))
  validateVars$refClust <- reactive(refClust)
})


observeEvent(input$compareValidate,{
  withBusyIndicatorServer("compareValidate",{
    validateVars$pValue <- reactive(isolate(input$validatePValue))
    validateVars$nClust <- reactive(isolate(input$validateNClust))
    validateVars$validateNPermut <- reactive(isolate(input$stepSizeExplorer ))
    refExp <-  as.matrix(validateVars$refExp())
    simExp <- as.matrix(validateVars$simExp())
    similarity <- sracipeHeatmapSimilarity(refExp,simExp, 
                                           returnData = TRUE,
                                           pValue = (input$validatePValue),
                                           nClusters = (input$validateNClust),
                                           permutations = (input$validateNPermut))
    
    plotData <- data.frame(t(similarity$simulated.refCor))
    plotData$row <- as.character(seq(1:ncol(similarity$simulated.refCor)))
    plotData <- reshape2::melt(plotData, value.name = "correlation")
    plotData$refName <- rep(colnames(similarity$simulated.refCor),
                            nrow(similarity$simulated.refCor))
    plotData$simName <- rep(rownames(similarity$simulated.refCor),
                            ncol(similarity$simulated.refCor))
    ggplot(data=plotData,
           aes(x=variable, y=row, fill=correlation)) + geom_tile() + 
      xlab("Simulation Data Samples") +
      ylab("Ref Data Samples") +
      scale_x_discrete(labels= plotData$simName) +
      scale_y_discrete(labels= plotData$refName) +
      theme_bw()
    
    output$validateRefHeatmap <- renderPlot({
      if(is.null(validateVars$refExp())) return()
      plotData <- validateVars$refExp()
      gplots::heatmap.2(similarity$dataReference, trace = "none",
                        dendrogram = "none", Colv=FALSE, col = plotColor,
                        ColSideColors = refClusters,
                        main = "Reference Data",distfun=function(x) as.dist(1-cor(t(x), method = "s")))
      
      
    })
    simClusters <- as.character(col2[(1+similarity$simClusters)])
    refClusters <- as.character(col2[(1+similarity$refClusters)])
    
    output$validateSimHeatmap <- renderPlot({
      if(is.null(validateVars$simExp())) return()
      gplots::heatmap.2(similarity$dataSimulation, trace = "none", 
                        dendrogram = "none", Colv=FALSE, col = plotColor,
                        ColSideColors = simClusters,
                        main = "Simulated Data",distfun=function(x) as.dist(1-cor(t(x), method = "s")))
    })
    output$validateRefSim <- renderPlot({
      if(is.null(similarity)) return()
      
      
      gplots::heatmap.2(as.matrix(similarity$simulated.refCor), Rowv = FALSE,
                        Colv = FALSE, ColSideColors = refClusters,
                        RowSideColors = simClusters ,trace = "none",
                        ylab =  "Simulation Data Samples",
                        xlab = "Ref Data Samples", col = plotColor, 
                        dendrogram = 'none',
                        main = "Correlation between reference
                        and simulated samples" ,distfun=function(x) as.dist(1-cor(t(x), method = "s")))
    })
    
    output$validateRefClustTable <- renderTable({ if(is.null(similarity)) return()
      tableData <- data.frame(similarity$ref.cluster.freq*100)
      colnames(tableData) <- c("Cluster", "Percentage in Reference Data")
      tableData <- tableData[-1,]
      return(tableData)
    })
    output$validateSimClustTable <- renderTable({if(is.null(similarity)) return()
      tableData <- data.frame(similarity$simulated.cluster.freq*100)
      colnames(tableData) <- c("Cluster", "Percentage in Simulated Data")
      tableData <- tableData[-1,]
      return(tableData)
    })
    output$validateKL <- renderText({
      return(paste("KL distance between the simulated and reference data: ", 
                   similarity$KL,"."))
      
    })
  })
  
})

