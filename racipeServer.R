
###########################################
# RACIPE
###########################################

rSetRacipe <- RacipeSE()
racipeVals <- reactiveValues()
sracipeCircuit(rSetRacipe) <- demoCircuit
racipeVals$rSet <- rSetRacipe
# Use a flag 
racipeVals$annealFlag <- FALSE
racipeVals$filterModels <- NULL
racipeVals$stochData <- NULL

output$racipeCircuit <- renderVisNetwork({
  .sracipePlotCircuit(racipeVals$rSet, plotToFile = FALSE)
})



observeEvent(input$simulateRacipe, {
  withBusyIndicatorServer("simulateRacipe",{
    shinyjs::show("racipeDeterministicText")
    shinyjs::show("racipeHeatmap")
    shinyjs::show("racipePca")
    shinyjs::show("parametricAnalysisRacipe")
    shinyjs::show("stochasticRacipe")
    shinyjs::show("paDescrpt")
    shinyjs::show("downloadRacipeData")
    shinyjs::show("downloadRacipeDataType")

    shinyjs::show("validateRacipe")
    
    rsRacipe <- racipeVals$rSet
    sracipeConfig(rsRacipe)$simParams["integrateStepSize"] <- 
      isolate(input$stepSizeRacipe)
    sracipeConfig(rsRacipe)$simParams["simulationTime"] <-  
      isolate(input$simTimeRacipe)
    sracipeConfig(rsRacipe)$stochParams["numModels"] <-
      isolate(input$numModels)
    rsRacipe <- sracipeSimulate(rsRacipe, plots = FALSE,genIC = TRUE, 
                                genParams = TRUE, integrate = TRUE, 
                                numModels = isolate(input$numModels),
                                integrateStepSize = isolate(input$stepSizeRacipe),
                                simulationTime = isolate(input$simTimeRacipe)
    )
    rsRacipe <- sracipeNormalize(rsRacipe)
    racipeVals$rsRacipe <- reactive(rsRacipe)
    sticVals$rsRacipe <- reactive(rsRacipe)

    output$racipeDeterministicText <- renderUI({HTML(
      "Hierarchical clustering and principal component 
    analysis of deterministic simulations ")})
    
    
    output$racipePca <- renderPlot({
      if(input$simulateRacipe == 0) return()
      withProgress(message = 'Plotting pca', value = 0.25, {
        plotData <- assay(rsRacipe)
        pca = prcomp(t(plotData), scale. = FALSE)
        racipeVals$pca <- pca
        pcaData <- data.frame(x=pca$x[,1],y=pca$x[,2])
        .sracipePlotDensity(pcaData, 
                            label1 = paste0("PC1(",100*summary(pca)$importance[2,1],"%)"),
                            label2 = paste0("PC2(",100*summary(pca)$importance[2,2],"%)"),
                            title = title)
        
      })
    })
    output$racipeHeatmap <- renderPlot({
      if(input$simulateRacipe == 0) return()
      withProgress(message = 'Hierarchical clustering', value = 0.25, {
        plotData <- assay(rsRacipe)
        .sracipePlotHeatmap(plotData, nClusters = input$racipeNClusters)
      })
    })
  })
  
  
  observeEvent(input$validateRacipe,{
    validateVars$simExp <- isolate(assay(racipeVals$rsRacipe(),1))
    output$validateSimHeatmap <- renderPlot({     return()})
    output$validateRefSim <- renderPlot({     return()})
    output$validateRefClustTable <- renderTable({  return() })
    output$validateSimClustTable <- renderTable({return()})
    output$validateKL <- renderText({return()})
    output$validateSimHeatmap <- renderPlot({
      if(is.null(validateVars$simExp)) return()
      gplots::heatmap.2(validateVars$simExp, trace = "none", 
                        Colv=TRUE, col = plotColor,
                        main = "Simulated Data",distfun=function(x) as.dist(1-cor(t(x), method = "s"))/2)
    })
  })
  ###########################################
  # Parametric Analysis
  ###########################################
  observeEvent(input$parametricAnalysisRacipe, {
    shinyjs::show("filteredOutputRacipe")
    shinyjs::show("filterSliderRacipe")
    shinyjs::show("filteredOutputRacipe2")
    shinyjs::show("filterSliderRacipe2")
    shinyjs::show("filteredOutputRacipe3")
    shinyjs::show("filterSliderRacipe3")
    
    
    shinyjs::show("racipeHeatmapFiltered")
    shinyjs::show("racipePcaFiltered")
    shinyjs::show("validateRacipePar")
    shinyjs::show("downloadRacipeParData")
    shinyjs::show("downloadRacipeParDataType")
    shinyjs::show("validateRacipeStoch")
    
    if(!is.null(rsRacipe)){
      parameterNames <- sracipeGenParamNames(rsRacipe)
      parameters <- as.matrix(sracipeParams(rsRacipe))
    }
    output$filteredOutputRacipe <- renderUI({
      selectInput("selectedParameter", "Parameter",
                  parameterNames,
                  selected = parameterNames[1]
      )
    })
    
    output$filteredOutputRacipe2 <- renderUI({

      selectInput("selectedParameter2", "Parameter",
                  parameterNames,
                  selected = parameterNames[2]
      )
    })
    
    output$filteredOutputRacipe3 <- renderUI({

      selectInput("selectedParameter3", "Parameter",
                  parameterNames,
                  selected = parameterNames[3]
      )
    })
    
    dataSimulation <- data.frame(t(assay(rsRacipe,1))) #sracipeNormalize(rsRacipe)
    pca = prcomp((dataSimulation), scale. = FALSE, center = FALSE)
    racipeVals$filterModels <- reactive(((parameters[,input$selectedParameter] >= 
                                            0.01*(input$parameterInput[1]*max(parameters[,input$selectedParameter]))) & 
                                           (parameters[,input$selectedParameter] < 
                                              0.01*(input$parameterInput[2]*max(parameters[,input$selectedParameter]))) &
                                           (parameters[,input$selectedParameter2] >= 0.01*(input$parameterInput2[1]*
                                                                                             max(parameters[,input$selectedParameter2]))) & 
                                           (parameters[,input$selectedParameter2] < 
                                              0.01*(input$parameterInput2[2]*max(parameters[,input$selectedParameter2]))) &
                                           (parameters[,input$selectedParameter3] >= 0.01*(input$parameterInput3[1]*
                                                                                             max(parameters[,input$selectedParameter3]))) & 
                                           (parameters[,input$selectedParameter3] < 0.01*(input$parameterInput3[2]*
                                                                                            max(parameters[,input$selectedParameter3])))
    ))
    filtered <- reactive({
      if (is.null(dataSimulation)) {
        return(NULL)
      }
      dataSimulation[racipeVals$filterModels(),]
      
    })
    
    output$filterSliderRacipe <- renderUI({
      if(is.null(parameters)) return(NULL)
      sliderInput("parameterInput", "Parameter Range", min = 0,
                  max = 100, value = c(0,100))
    })
    
    output$filterSliderRacipe2 <- renderUI({
      if(is.null(parameters)) return(NULL)
      sliderInput("parameterInput2", "Parameter Range", min = 0,
                  max = 100, value = c(0,100))
    })
    
    
    output$filterSliderRacipe3 <- renderUI({
      if(is.null(parameters)) return(NULL)
      sliderInput("parameterInput3", "Parameter Range", min = 0,
                  max = 100, value = c(0,100))
    })
    
    
    output$racipePcaFiltered <- renderPlot({
      if(is.null(filtered())) return(NULL)
      pcaData <- filtered()
      pcaData <- scale(pcaData, pca$center, pca$scale) %*% pca$rotation
      pcaData <- data.frame(PC1=pcaData[,1],PC2=pcaData[,2])
      .sracipePlotDensity(pcaData,label1 = paste0("PC1(",100*summary(pca)$importance[2,1],"%)"),
                          label2 = paste0("PC2(",100*summary(pca)$importance[2,2],"%)"))
    })
    
    output$racipeHeatmapFiltered <- renderPlot({
      if(is.null(filtered())) return(NULL)
      withProgress(message = 'Hierarchichal clustering', value = 0.25, {
        plotData <- filtered()
        .sracipePlotHeatmap(t(plotData),nClusters = input$racipeNClusters)
      })
    })
    observeEvent(input$validateRacipePar,{
      validateVars$simExp <- t(as.matrix(isolate(filtered())))
      output$validateSimHeatmap <- renderPlot({     return()})
      output$validateRefSim <- renderPlot({     return()})
      output$validateRefClustTable <- renderTable({  return() })
      output$validateSimClustTable <- renderTable({return()})
      output$validateKL <- renderText({return()})
      output$validateSimHeatmap <- renderPlot({
        if(is.null(validateVars$simExp)) return()
        gplots::heatmap.2(validateVars$simExp, trace = "none", 
                          Colv=TRUE, col = plotColor,
                          main = "Simulated Data",distfun=function(x) as.dist(1-cor(t(x), method = "s"))/2)
      })
    })
    
  })
  
  
  
  ########################### Download and Upload ########################
  output$downloadRacipeData <- downloadHandler(
    filename = function(){
      if(input$downloadRacipeDataType == "csv"){
        return(paste0(Sys.time(),"_",annotation(racipeVals$rsRacipe()),".zip"))}
      else return(paste0(Sys.time(),"_",annotation(racipeVals$rsRacipe()),".RDS"))
    },
    
    content = function(file){
      if(input$downloadRacipeDataType == "csv"){
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL;
        
        fileName <- paste(annotation(racipeVals$rsRacipe()),"_GE",".csv",
                          sep = "")
        write.csv(assay(racipeVals$rsRacipe()),fileName)
        files <- c(fileName,files)
        
        fileName <- paste(annotation(racipeVals$rsRacipe()),"_IC",".csv",
                          sep = "")
        write.csv(sracipeIC(racipeVals$rsRacipe()),fileName, row.names = TRUE)
        files <- c(fileName,files)
        
        fileName <- paste(annotation(racipeVals$rsRacipe()),"_network",".csv",
                          sep = "")
        write.csv(sracipeCircuit(racipeVals$rsRacipe()),fileName, row.names = FALSE, quote=FALSE)
        files <- c(fileName,files)
        
        fileName <- paste(annotation(racipeVals$rsRacipe()),"_params",".csv",
                          sep = "")
        write.csv(sracipeParams(racipeVals$rsRacipe()),fileName, row.names = FALSE)
        files <- c(fileName,files)
        
        #create the zip file
        zip(file,files, flags = "-r9Xj")
      }
      else{
        
        saveRDS(racipeVals$rsRacipe(), file)
        
      }
      
      
    })
  output$downloadRacipeParData <- downloadHandler(
    filename = function(){
      if(input$downloadRacipeParDataType == "csv"){
        return(paste0(Sys.time(),"_",annotation(racipeVals$rsRacipe()),".zip"))}
      else return(paste0(Sys.time(),"_",annotation(racipeVals$rsRacipe()),".RDS"))
    },
    
    content = function(file){
      dwnRacipe <- isolate(racipeVals$rsRacipe())
      dwnRacipe <- RacipeSE(dwnRacipe[,which(isolate(racipeVals$filterModels()))])
      if(input$downloadRacipeDataType == "csv"){
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL;
        
        fileName <- paste(annotation(dwnRacipe),"_GE",".csv",
                          sep = "")
        
        write.csv(assay(dwnRacipe),fileName)
        files <- c(fileName,files)
        fileName <- paste(annotation(dwnRacipe),"_IC",".csv",
                          sep = "")
        write.csv(sracipeIC(dwnRacipe),fileName, row.names = TRUE)
        files <- c(fileName,files)
        
        fileName <- paste(annotation(racipeVals$rsRacipe()),"_network",".csv",
                          sep = "")
        write.csv(sracipeCircuit(dwnRacipe),fileName, row.names = FALSE, quote=FALSE)
        files <- c(fileName,files)
        
        fileName <- paste(annotation(dwnRacipe),"_params",".csv",
                          sep = "")
        write.csv(sracipeParams(dwnRacipe),fileName, row.names = FALSE)
        files <- c(fileName,files)
        
        #create the zip file
        zip(file,files, flags = "-r9Xj")
      }
      else{
        
        saveRDS(dwnRacipe, file)
        
      }
      
    })
    })
###########################################
# sRACIPE
###########################################
observeEvent(input$stochasticRacipe, {
  
  shinyjs::show("sRacipeOption")
  shinyjs::show("sRacipeNoise")
  shinyjs::show("simulateSRacipe")
  shinyjs::show("sRacipeHeatmap")
  shinyjs::show("sRacipeHeatmapDet")
  shinyjs::show("sRacipePca")
  shinyjs::show("sRacipePcaDet")
  shinyjs::show("saveSRacipeData")
  shinyjs::show("downloadSRacipeData")
  shinyjs::show("downloadSRacipeDataType")
  observeEvent(input$simulateSRacipe, {
    shinyjs::show("validateRacipeStoch")
    
    withBusyIndicatorServer("simulateSRacipe",{
      output$CN <- renderText("")
      output$Anneal <- renderText("")

      isolate(
        if(input$sRacipeOption == "constantNoise")
        {
          output$CN <- renderText("Constant Noise Plots")
          rsSRacipe <- racipeVals$rsRacipe()
          detParam <- FALSE
          if(is.null(sracipeParams(rsSRacipe))) detParam <- TRUE
          withProgress(message = 'Simulating', value = 0.25, {
            shinyBS::createAlert(session, anchorId = "racipeAlert",
                                 alertId =  "racipeProcessing", title = "Processing",
                                 content = "Please wait...", append = FALSE)
            
            rsSRacipe <- sracipeSimulate(
              rsSRacipe, plots = FALSE, anneal = FALSE,
              nNoise = 1, numModels = isolate(input$numModels),
              integrateStepSize = isolate(input$stepSizeRacipe),
              simulationTime = isolate(input$simTimeRacipe),
              # paramRange = isolate(input$parameterRange),
              initialNoise =  50*(0.5)^(20-isolate(input$sRacipeNoise)),
              genParams = detParam, genIC = detParam)
            rsSRacipe <- sracipeNormalize(rsSRacipe)
            
            racipeVals$rsSRacipe <- reactive(rsSRacipe)
          })

          assayDataTmp <- assays(rsSRacipe)
          metadataTmp <- metadata(rsSRacipe)
          plotData <- assayDataTmp[[1]]
          pca1 = prcomp(t(plotData), scale. = FALSE) 

          output$sRacipePcaDet <-renderPlot({
            if(is.null(plotData)) return(NULL)
            .sracipePlotPca(plotData = plotData, pca = pca1,
                            title = "Deterministic")
          })
          
          output$sRacipeHeatmapDet <- renderPlot({
            if(input$simulateSRacipe == 0) return()
            .sracipePlotHeatmap(plotData, nClusters = input$racipeNClusters)
          })
          plotDataSt <- assayDataTmp[[2]]
          racipeVals$stochData <- reactive(assayDataTmp[[2]])

          stochasticPca <- (scale(t(plotDataSt), pca1$center, pca1$scale) %*%
                              pca1$rotation)
          
          
          output$sRacipePca <-renderPlot({
            if(is.null(plotDataSt)) return(NULL)
            
            pcaData <- data.frame(x=stochasticPca[,1],y=stochasticPca[,2])
            
            pca1$x <- stochasticPca
            .sracipePlotPca(plotData = stochasticPca, pca = pca1,
                            title = "Stochastic" )
          })
          output$sRacipeHeatmap <- renderPlot({
            if(input$simulateSRacipe == 0) return()
            .sracipePlotHeatmap(plotDataSt, nClusters = input$racipeNClusters)
          })
          
        }
      )
      
      if(input$sRacipeOption == "annealing")
      {
        output$Anneal <- renderText("Annealing Simulation Data")
        
        if(!racipeVals$annealFlag){
          rsSRacipeAnneal <- racipeVals$rsRacipe()
          detParam <- FALSE
          if(is.null(sracipeParams(rsSRacipeAnneal))) detParam <- TRUE
          withProgress(message = 'Simulating', value = 0.25, {
            shinyBS::createAlert(
              session, anchorId = "racipeAlert",
              alertId =  "racipeProcessing", title = "Processing",
              content = "Please wait...", append = FALSE)
            rsSRacipeAnneal <- sracipeSimulate(
              rsSRacipeAnneal, annealing = TRUE, nNoise = 20, 
              numModels = isolate(input$numModels),
              integrateStepSize = isolate(input$stepSizeRacipe),
              simulationTime = isolate(input$simTimeRacipe),
              initialNoise =  50/sqrt(length(names(rsSRacipeAnneal))),
              noiseScalingFactor = 0.5,
              genIC = detParam, genParams = detParam)
            rsSRacipeAnneal <- sracipeNormalize(rsSRacipeAnneal)
            
            racipeVals$annealFlag <- TRUE
            racipeVals$rsSRacipeAnneal <- reactive(rsSRacipeAnneal)
          })
        }
        
        rsSRacipeAnneal <- isolate(racipeVals$rsSRacipeAnneal())
        assayDataTmp <- assays(rsSRacipeAnneal)
        metadataTmp <- metadata(rsSRacipeAnneal)
        plotData <- assayDataTmp[[1]]
        pca1 = prcomp(t(plotData), scale. = FALSE)

        output$sRacipePcaDet <-renderPlot({
          if(is.null(plotData)) return(NULL)
          .sracipePlotPca(plotData = plotData, pca = pca1, 
                          title = "Deterministic")
        })
        
        output$sRacipeHeatmapDet <- renderPlot({
          if(input$simulateSRacipe == 0) return()
          .sracipePlotHeatmap(plotData, nClusters = input$racipeNClusters)
        })
        plotDataSt <- assayDataTmp[[(22 - input$sRacipeNoise)]]
        racipeVals$stochData <- reactive(plotDataSt)

        stochasticPca <- (scale(t(plotDataSt), pca1$center, pca1$scale) %*%
                            pca1$rotation)
        
        
        output$sRacipePca <-renderPlot({
          if(is.null(plotDataSt)) return(NULL)
          
          pcaData <- data.frame(x=stochasticPca[,1],y=stochasticPca[,2])
          pca1$x <- stochasticPca
          .sracipePlotPca(plotData = stochasticPca, pca = pca1,
                          title = "Stochastic")
        })
        output$sRacipeHeatmap <- renderPlot({
          if(input$simulateSRacipe == 0) return()
          .sracipePlotHeatmap(plotDataSt, nClusters = input$racipeNClusters)
        })
        shinyBS::closeAlert(session, alertId = "racipeProcessing")
        #    })
        
      }
    })
  })
  
  output$downloadSRacipeData <- downloadHandler(
    filename = function(){
      if(input$downloadSRacipeDataType == "csv"){
        return(paste0(Sys.time(),"_",annotation(racipeVals$rSet),".zip"))}
      else return(paste0(Sys.time(),"_",annotation(racipeVals$rSet),".RDS"))
    },
    
    content = function(file){
      if(input$downloadSRacipeDataType == "csv"){
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL
        tmpRSet <- RacipeSE()
        if(input$sRacipeOption == "annealing")
          tmpRSet <- racipeVals$rsSRacipeAnneal()
        if(input$sRacipeOption == "constantNoise")
        {
          tmpRSet <- racipeVals$rsSRacipe()}
        fileName <- paste(annotation(tmpRSet),"_GE",".csv",sep = "")
        write.csv(assays(tmpRSet),fileName)
        files <- c(fileName,files)
        
        fileName <- paste(annotation(tmpRSet),"_IC",".csv",sep = "")
        write.csv(sracipeIC(tmpRSet),fileName, row.names = TRUE)
        files <- c(fileName,files)
        
        fileName <- paste(annotation(tmpRSet),"_network",".csv",sep = "")
        write.csv(sracipeCircuit(tmpRSet),fileName, row.names = FALSE, quote = FALSE)
        files <- c(fileName,files)
        
        fileName <- paste(annotation(tmpRSet),"_params",".csv",sep = "")
        write.csv(sracipeParams(tmpRSet),fileName, row.names = FALSE)
        files <- c(fileName,files)
        
        #create the zip file
        zip(file,files, flags = "-r9Xj")
      }
      else{
        tmpRSet <- RacipeSE()
        if(input$sRacipeOption == "annealing")
          tmpRSet <- racipeVals$rsSRacipeAnneal()
        if(input$sRacipeOption == "constantNoise")
        {
          tmpRSet <- racipeVals$rsSRacipe()}
        saveRDS(tmpRSet, file)
        
      }
      
      
    })
  
  
})

observeEvent(input$validateRacipeStoch,{
  validateVars$simExp <- isolate(racipeVals$stochData())
  output$validateSimHeatmap <- renderPlot({     return()})
  #  output$validateRefHeatmap <- renderPlot({     return()})
  output$validateRefSim <- renderPlot({     return()})
  output$validateRefClustTable <- renderTable({  return() })
  output$validateSimClustTable <- renderTable({return()})
  output$validateKL <- renderText({return()})
  output$validateSimHeatmap <- renderPlot({
    if(is.null(validateVars$simExp)) return()
    gplots::heatmap.2(validateVars$simExp, trace = "none", 
                      Colv=TRUE, col = plotColor,
                      main = "Simulated Data",distfun=function(x) as.dist(1-cor(t(x), method = "s"))/2)
  })
})