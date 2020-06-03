
circuitVariables <- reactiveValues()
circuitVariables$circuit <- data.frame(
  Source = c("A", "B"),
  Target = c("B", "A"),
  Interaction = c(2, 2), stringsAsFactors = FALSE
)
circuitVariables$deletedRows <- NULL 
circuitVariables$deletedRowIndices <- list()
circuitVariables$rSet <- NULL
observeEvent(input$circuitFile,{
  circuitData <- input$circuitFile
  
  circuitVariables$circuit <- data.table::fread(
    circuitData$datapath, sep="auto", header ="auto", stringsAsFactors = FALSE)
  circuitVariables$deletedRows <- NULL
  circuitVariables$deletedRowIndices = list()
  
})

output$downloadSampleCircuit <- downloadHandler(
  filename = "circuit.csv",
  content = function(con) {
 data = data.frame(
      Source = c("A", "B"),
      Target = c("B", "A"),
      Interaction = c(2, 2))
    write.csv(data, con, row.names = FALSE, quote = FALSE)
  }
)
output$downloadCircuit <- downloadHandler(
  filename = function() {
    if(is.null(circuitVariables$rSet)) return("circuit.csv")
    paste(annotation(circuitVariables$rSet), '.csv', sep='')
  },
  content = function(con) {
    if(is.null(circuitVariables$rSet)) { data = data.frame(
      Source = c("A", "B"),
      Target = c("B", "A"),
      Interaction = c(2, 2), stringsAsFactors = FALSE)}
    else {data = sracipeCircuit(circuitVariables$rSet)}
    write.csv(data, con, row.names = FALSE, quote = FALSE)
  }
)

observeEvent(input$toggleAddInteraction, {

  shinyjs::toggle("newIntSrc") 
  shinyjs::toggle("newIntTgt") 
  shinyjs::toggle("newIntType")  
  shinyjs::toggle("addInteraction") 
})

observeEvent(input$circuitTable_cell_edit, {
  circuitVariables$circuit[input$circuitTable_cell_edit$row,
                        (input$circuitTable_cell_edit$col)] <- 
    input$circuitTable_cell_edit$value
})

observeEvent(input$addInteraction, {
  rowNum <- (1 + nrow(circuitVariables$circuit))
  dataRow <- data.frame(Source = input$newIntSrc, Target = input$newIntTgt, 
                        Interaction = as.integer(input$newIntType),
                        stringsAsFactors = FALSE)
  circuitVariables$circuit <- addRowAt(
    circuitVariables$circuit, dataRow, 1)
  
})

observeEvent(input$deletePressed, {
  rowNum <- parseDeleteEvent(input$deletePressed)
  dataRow <- circuitVariables$circuit[rowNum,]
  
  # Put the deleted row into a data frame so we can undo
  # Last item deleted is in position 1
  circuitVariables$deletedRows <- rbind(dataRow, circuitVariables$deletedRows)
  circuitVariables$deletedRowIndices <- 
    append(circuitVariables$deletedRowIndices, rowNum, after = 0)
  
  # Delete the row from the data frame
  circuitVariables$circuit <- circuitVariables$circuit[-rowNum,]
})

observeEvent(input$undoDelete, {
  if(nrow(circuitVariables$deletedRows) > 0) {
    row <- circuitVariables$deletedRows[1, ]
    circuitVariables$circuit <- addRowAt(
      circuitVariables$circuit, row, circuitVariables$deletedRowIndices[[1]])
    
    # Remove row
    circuitVariables$deletedRows <- circuitVariables$deletedRows[-1,]
    # Remove index
    circuitVariables$deletedRowIndices <- circuitVariables$deletedRowIndices[-1]
  }
})

# Disable the undo button if we have not deleted anything
output$undoUI <- renderUI({
  if(
    !is.null(circuitVariables$deletedRows) && 
    nrow(circuitVariables$deletedRows) > 0) {
    actionButton(
      'undoDelete', label = 'Undo Delete', icon('undo'), 
      style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  } else {
    actionButton('undo', label = 'Undo delete', icon('undo'), disabled = TRUE)
  }
})

output$circuitTable <- DT::renderDataTable(
  # Add the delete button column
  deleteButtonColumn(circuitVariables$circuit, 'delete_button') 
)

observeEvent(input$updateCircuit,{
    rSet <-  RacipeSE()
    
    sracipeCircuit(rSet) <- circuitVariables$circuit
    annotation(rSet) <- paste(input$circuitName, Sys.Date(),
                              basename(tempfile()), sep = "_")
    circuitVariables$rSet <- rSet
    racipeVals$rSet <- rSet
    gvVals$rSet <- reactive(sracipeSimulate(rSet,integrate = FALSE,
                                            numModels = 1))
    gvVals$parameterNamesME <- reactive(sracipeGenParamNames(gvVals$rSet()))
    gvVals$parametersME <- reactive(sracipeParams(gvVals$rSet()))
    
    
    racipeVals$filterModels <- NULL
    racipeVals$stochData <- NULL
    
    gvVals$icNamesME <- reactive(names(gvVals$rSet()))
    gvVals$icME <- reactive(sracipeIC(gvVals$rSet()))
    gvVals$rSetBifur <- reactive(gvVals$rSet())
    
    validateVars$simExp <- NULL
    validateVars$refClust <- NULL
    
    hideLoadCkt()
})

observeEvent(input$updateCircuitTable,{
  rSet <-  RacipeSE()
  sracipeCircuit(rSet) <- circuitVariables$circuit
  annotation(rSet) <- paste(input$circuitName, Sys.Date(),basename(tempfile()),
                            sep = "_")
  circuitVariables$rSet <- rSet
  racipeVals$rSet <- rSet
  gvVals$rSet <- reactive(sracipeSimulate(rSet,integrate = FALSE,
                                          numModels = 1))
  gvVals$parameterNamesME <- reactive(sracipeGenParamNames(gvVals$rSet()))
  gvVals$parametersME <- reactive(sracipeParams(gvVals$rSet()))
  hideLoadCkt()
})

output$circuit <- renderVisNetwork({
  if(is.null(circuitVariables$rSet)) { return (NULL)}
  .sracipePlotCircuit(circuitVariables$rSet, plotToFile = FALSE)
})

#' A column of delete buttons for each row in the data frame 
#' for the first column
#'
#' @param df data frame
#' @param id id prefix to add to each actionButton. 
#' The buttons will be id'd as id_INDEX.
#' @return A DT::datatable with escaping turned off that has the delete buttons
#'  in the first column and \code{df} in the other
deleteButtonColumn <- function(df, id, ...) {
  # function to create one action button as string
  f <- function(i) {
    # https://shiny.rstudio.com/articles/communicating-with-js.html
    as.character(
      actionButton(
        paste(id, i, sep="_"), label = NULL, icon = icon('trash'),
        onclick = 'Shiny.setInputValue(\"deletePressed\",  this.id, 
        {priority: "event"})'))
  }
  
  deleteCol <- unlist(lapply(seq_len(nrow(df)), f))
  
  # Return a data table
  DT::datatable(cbind(df,Delete = deleteCol),
                # Need to disable escaping for html as string to work
                escape = FALSE,
                options = list(pageLength = 100, 
                               # Disable sorting for the delete column
                               columnDefs = list(list(targets = 3, sortable = FALSE))
                ),editable = FALSE, rownames= FALSE)
}

#' Extracts the row id number from the id string
#' @param idstr the id string formated as id_INDEX
#' @return INDEX from the id string id_INDEX
parseDeleteEvent <- function(idstr) {
  res <- as.integer(sub(".*_([0-9]+)", "\\1", idstr))
  if (! is.na(res)) res
}

# Using datatable functions from https://github.com/stefaneng/Shiny-DeleteRowsDT
#' Adds a row at a specified index
#'
#' @param df a data frame
#' @param row a row with the same columns as \code{df}
#' @param i the index we want to add row at.
#' @return the data frame with \code{row} added to \code{df} at index \code{i}
addRowAt <- function(df, row, i) {
  # Slow but easy to understand
  if (i > 1) {
    rbind(df[1:(i - 1), ], row, df[-(1:(i - 1)), ])
  } else {
    rbind(row, df)
  }
  
}
