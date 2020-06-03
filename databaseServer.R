databaseVals <- reactiveValues()
databaseVals$rs <- NULL
responseMsg1 <- ""
responseMsg1 <-  eventReactive(input$shinySelectNetworkDb,
                               "No network with matching name found!")
responseMsg2 <- ""
responseMsg2 <-  eventReactive(input$shinySelectGeneDb, 
                               "No network with this gene found!")

output$shinySelectedNetworkGene <- renderText(responseMsg1())
output$shinySelectedNetwork <- renderText(responseMsg2())
databaseVals <- reactiveValues()
databaseVals$databaseTableAll <- readRDS("data/GeneExDatabase_10112019.RDS")
databaseVals$databaseTable <- reactive(databaseVals$databaseTableAll)

observeEvent(input$biologicalDB,{
  databaseVals$databaseTable <- reactive(databaseVals$databaseTableAll[Biological == "TRUE"])
})

observeEvent(input$syntheticDB,{
  databaseVals$databaseTable <- reactive(databaseVals$databaseTableAll[Biological == "FALSE"])
})

observeEvent(input$allDB,{
  databaseVals$databaseTable <- reactive(databaseVals$databaseTableAll)
})

databaseVals$databaseTable <- reactive(readRDS("data/GeneExDatabase_10112019.RDS"))


output$databaseTable <- DT::renderDT({
  databaseTable <- databaseVals$databaseTable()
  selectedRow <- reactive(input$databaseTable_row_last_clicked)
  if(is.null(selectedRow())){return(databaseTable[,1:10])}
  observeEvent(input$databaseTable_row_last_clicked,{
    selectedRow <- reactive(input$databaseTable_row_last_clicked)

    shinyjs::show("loadNetworkDatabase")
    shinyjs::hide("downloadDbData")
    shinyjs::hide("downloadDbDataType")
    output$msg <- renderText("")
    networkName <- databaseTable[selectedRow(),"Name"]
    circuit <- data.frame(databaseTable[selectedRow(),"Circuit"]$Circuit  )
    rs <- RacipeSE()
    sracipeCircuit(rs) <- circuit
    
    annotation(rs) <- networkName
    databaseVals$rSet <- reactive(rs)
    output$tableDbNetwork <- DT::renderDT({
      circuit
    }, options = list(
      pageLength = 100))
    
    output$plotDbNetwork <- renderVisNetwork({
      .sracipePlotCircuit(rs, plotToFile = FALSE)
    })
    databaseVals$rs <- reactive(rs)
  })
   
  return(databaseTable[,1:10])  
}, selection = 'single', editable = FALSE, rownames= FALSE, 
options = list(pageLength = 5)
)

observeEvent(input$loadNetworkDatabase, {
  
  circuitTmp <- sracipeCircuit(isolate(databaseVals$rs()))
  colnames(circuitTmp) <- c("Source", "Target", "Interaction")
  circuitVariables$circuit <- data.frame(circuitTmp)
  circuitVariables$deletedRows <- NULL 
  circuitVariables$deletedRowIndices <- list()
  
  rSet <- RacipeSE()
  sracipeCircuit(rSet) <- circuitVariables$circuit
  annotation(rSet) <- paste(annotation(isolate(databaseVals$rs())), Sys.Date(),basename(tempfile()), sep = "_")
  
  circuitVariables$rSet <- rSet
  
  racipeVals$rSet <- rSet
  racipeVals$filterModels <- NULL
  racipeVals$stochData <- NULL
  
  gvVals$rSet <- reactive(sracipeSimulate(rSet,integrate = FALSE,
                                          numModels = 1))
  gvVals$parameterNamesME <- reactive(sracipeGenParamNames(gvVals$rSet()))
  gvVals$parametersME <- reactive(sracipeParams(gvVals$rSet()))
  gvVals$icNamesME <- reactive(names(gvVals$rSet()))
  gvVals$icME <- reactive(sracipeIC(gvVals$rSet()))
  gvVals$rSetBifur <- reactive(gvVals$rSet())
  
  validateVars$simExp <- NULL
  validateVars$refClust <- NULL
  
  hideLoadCkt()
  
  
})