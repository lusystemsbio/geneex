###########################################
# Model Explorer
###########################################
gvVals = reactiveValues()
demoCircuit = data.table(
  Source = c("A", "B"),
  Target = c("B", "A"),
  Interaction = c(2, 2))
rSetTmp <- sracipeSimulate(circuit = demoCircuit, integrate = FALSE, 
                           numModels = 1)
annotation(rSetTmp) <- "Circuit1"
gvVals$rSet <- reactive(rSetTmp)
gvVals$parameterNamesME <- reactive(sracipeGenParamNames(rSetTmp))
gvVals$icNamesME <- reactive(names(rSetTmp))
gvVals$parametersME <- reactive(sracipeParams(rSetTmp))
gvVals$icME <- reactive(sracipeIC(rSetTmp))
gvVals$rSetBifur <- reactive(rSetTmp)
######################### Network Input ######################################

output$circuitGv <- renderVisNetwork({
  if(is.null(gvVals$rSet())) { return (NULL)}
  rSet <- gvVals$rSet()
  .sracipePlotCircuit(rSet, plotToFile = FALSE)
})
# Array to hold reactive variables like network topology



gvVals$simTimeExplorer <- reactive({
    return(input$simTimeExplorer)
})
# Use maximum and minimum limits on the initial step size
gvVals$stepSizeExplorer <- reactive({
  if(input$stepSizeExplorer > 1) return(0.99)
  if(input$stepSizeExplorer < .0001) return(0.0001)
  else return(input$stepSizeExplorer)
})

gvVals$initialNoise <- reactive({

    return(input$noiseLevel)
})

gvVals$nNoise <- reactive({

  n = 0L
  if (input$noiseLevel>0) n = 1L
  n
})
##################### Parameter Modification #######################

output$modelParams <- renderUI({
     if(is.null(gvVals$parametersME())) return(NULL)
  selectInput("selectedParameterME", "Parameter",
              gvVals$parameterNamesME(),
              selected = NULL)
})

output$modelIc <- renderUI({
  if(is.null(gvVals$icME())) return(NULL)
  selectInput("selectedIcME", "Initial Condition",
              gvVals$icNamesME(),
              selected = NULL)
})

output$modelPlotGene <- renderUI({
  if(is.null(gvVals$icME())) return(NULL)
  selectInput("selectedPlotGeneME", "Genes to plot",
              gvVals$icNamesME(),
              selected = gvVals$icNamesME(), multiple = TRUE)
})

output$modelPlotGeneBif <- renderUI({
  if(is.null(gvVals$icME())) return(NULL)
  selectInput("selectedPlotGeneBifME", "Genes to plot",
              gvVals$icNamesME(),
              selected = gvVals$icNamesME(), multiple = TRUE)
})


output$newModelParamValue <- renderUI({
  if(is.null(gvVals$parametersME())) return(NULL)
  textInput("parameterValue", "Value",
            placeholder = gvVals$parametersME()[input$selectedParameterME],
            value = gvVals$parametersME()[input$selectedParameterME])
  
})

output$newModelIcValue <- renderUI({
  if(is.null(gvVals$icME())) return(NULL)
  textInput("icValue", "Value",
            placeholder = gvVals$icME()[input$selectedIcME,],
            value = gvVals$icME()[input$selectedIcME,])
  
})


observeEvent(input$updateGvParam,{
  params <- gvVals$parametersME()
  params[input$selectedParameterME] <- as.numeric(input$parameterValue)
  gvVals$parametersME <- reactive(params)
})

observeEvent(input$updateGvIc,{
  ics <- gvVals$icME()
  ics[input$selectedIcME,] <- as.numeric(input$icValue)
  gvVals$icME <- reactive(ics)
})


observeEvent(input$simulateGv, {
  withBusyIndicatorServer("simulateGv",{
  shinyjs::show("bifurcationExplorer")
  shinyjs::show("downloadMEData")
  shinyjs::show("downloadMEDataType")
  shinyjs::show("GvTS")
  rs <- gvVals$rSet()
  sracipeParams(rs) <- gvVals$parametersME()
  sracipeIC(rs) <- gvVals$icME()
   rs <- sracipeSimulate(rs, timeSeries = TRUE, plots = FALSE, genIC = FALSE, 
                         genParams = FALSE, integrate = TRUE,
                        simulationTime  = isolate(gvVals$simTimeExplorer()),
                       integrateStepSize = isolate(gvVals$stepSizeExplorer()),
                       initialNoise = isolate(gvVals$initialNoise()), 
                       nNoise = isolate(gvVals$nNoise()) 
                       )
   tsPlot <- reactive(if(input$simulateGv == 0) return()
                           else{
                             shinyjs::show("logPlotGV")
                             shinyjs::show("modelPlotGene")
                           plotData <- t(metadata(gvVals$rSet())$timeSeries)
                           p <- .sracipePlotTS(plotData, input$selectedPlotGeneME)
                           if(input$logPlotGV)
                             p <- p + scale_y_log10()
                           return(p)
                           }
                           )
   
   output$GvTS <- renderPlot({
     tsPlot()

   })
   
  gvVals$rSet <- reactive(rs)

  })

  output$downloadMEData <- downloadHandler(
      filename = function(){
        if(input$downloadMEDataType == "csv"){
        return(paste0(annotation(gvVals$rSet()),".zip"))}
        else return(paste0(annotation(gvVals$rSet()),".RDS"))
      },

      content = function(file){
        if(input$downloadMEDataType == "csv"){
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL;

        fileName <- paste(annotation(gvVals$rSet()),"_GE",".csv",sep = "")
        write.csv(metadata(gvVals$rSet())$timeSeries,fileName)
        files <- c(fileName,files)

        fileName <- paste(annotation(gvVals$rSet()),"_network",".csv",sep = "")
        write.csv(sracipeCircuit(gvVals$rSet()),fileName, row.names = FALSE, quote = FALSE)
        files <- c(fileName,files)

        fileName <- paste(annotation(gvVals$rSet()),"_params",".csv",sep = "")
        write.csv(sracipeParams(gvVals$rSet()),fileName, row.names = FALSE)
        files <- c(fileName,files)

        #create the zip file
        zip(file,files, flags = "-r9Xj")

        
        
        }
        else{
            saveRDS(gvVals$rSet(), file)
      }
    })
  ##################### Bifurcation Diagram #######################
  observeEvent(input$bifurcationExplorer,{
    shinyjs::show("modelParamsBif")
    shinyjs::show("modelParamBifMin")
    shinyjs::show("modelParamBifMax")
    shinyjs::show("modelNumBifurs")
    shinyjs::show("bifurcationME")
    shinyjs::show("modifiedBifME")
    
    
  })
  # Limit the maximum number of models to 5000.

  output$modelParamsBif <- renderUI({
    if(is.null(gvVals$rSet())) return(NULL)
    selectInput("selectedParameterMEBif", "Parameter",
                gvVals$parameterNamesME(),
                selected = 1)
  })

  output$modelParamBifMin <- renderUI({
    if(is.null(gvVals$parametersME())) return(NULL)
    textInput("parameterValueBifMin", "Min value", 
              value = 0.9*gvVals$parametersME()[input$selectedParameterMEBif], 
       placeholder = 0.9*gvVals$parametersME()[input$selectedParameterMEBif])
  })

  output$modelParamBifMax <- renderUI({
    if(is.null(gvVals$parametersME())) return(NULL)
    textInput("parameterValueBifMax", "Max value", 
              value = 1.1*gvVals$parametersME()[input$selectedParameterMEBif], 
       placeholder = 1.1*gvVals$parametersME()[input$selectedParameterMEBif])
  })


  observeEvent(input$bifurcationME, {
    withBusyIndicatorServer("bifurcationME",{
    shinyjs::show("downloadBifData")
    modelNumBifurs <-   isolate(input$modelNumBifurs)

    newParametersME <- 
      isolate(gvVals$parametersME()[rep(seq_len(nrow(gvVals$parametersME())), 
                                                      modelNumBifurs),])
    modPar <- isolate(seq(from = as.numeric(input$parameterValueBifMin), 
                  to = as.numeric(input$parameterValueBifMax),
                  length.out = modelNumBifurs))
    newParametersME[input$selectedParameterMEBif] <- modPar
    rs <- RacipeSE()
    
    sracipeCircuit(rs) <- sracipeCircuit(gvVals$rSet())
    rs <- sracipeSimulate(rs, genIC = TRUE, genParams = TRUE, integrate = FALSE,
                          initialNoise = 0, nNoise = 0, 
                          numModels = modelNumBifurs)
    sracipeParams(rs) <- newParametersME
    rs <- sracipeSimulate(rs, genIC = TRUE, genParams = FALSE, integrate = TRUE,
                          initialNoise = 0, nNoise = 0, 
                          numModels = modelNumBifurs,
                          simulationTime  = isolate(gvVals$simTimeExplorer()),
                          integrateStepSize = isolate(gvVals$stepSizeExplorer())
                          
                          )
    gvVals$rSetBifur <- reactive(rs)
    gvVals$modelNumBifurs <- reactive(modelNumBifurs)
    gvVals$newParametersME <- reactive(newParametersME)
    output$modifiedBifME <- renderPlot({
      if(input$bifurcationME == 0) return()
      shinyjs::show("logPlotGVBif")
      shinyjs::show("modelPlotGeneBif")
      library(ggplot2)
      sexprs <- assay(gvVals$rSetBifur(),1)
      sexprs <- reshape2::melt(t(sexprs))
      colnames(sexprs) <- c("bifurParameter","Gene","geneExp")
      test <- modPar#sracipeParams(gvVals$rSetBifur())[,input$selectedParameterMEBif]
      test <- rep(test,times = dim(gvVals$rSetBifur())[1])

      sexprs$bifurParameter <- test
      theme_set(theme_bw(base_size = 18))
      sexprs <- sexprs[which(sexprs$Gene %in% input$selectedPlotGeneBifME),,drop = FALSE ]
      p <- ggplot2::ggplot(sexprs) +
        geom_point(aes(x=bifurParameter,y=geneExp,color = Gene)) +
        xlab(isolate(input$selectedParameterMEBif)) +
        ylab("Gene Expression") +
        
        NULL
      if(input$logPlotGVBif)
        p <- p + scale_y_log10()
      return(p)
    })
    })
  }
  )
  output$downloadBifData <- downloadHandler(
    filename <- paste0(Sys.time(),"_",annotation(gvVals$rSetBifur()),".RDS" ),
    content = function(con) {
      saveRDS(gvVals$rSetBifur(), con)
    }
  )
})




