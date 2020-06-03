# Module UI function
shinyLoadNetworkUI <- function(id, label = "loadNetworkUI") {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      column(3, offset = 0,
             tags$h4("Circuit Information"),
             textInput(ns("filenameTopo"), "Circuit Name", "Circuit1"),
             
             tabsetPanel(
               tabPanel("Enter Text",
                        fluidRow(
                          
                          textAreaInput( inputId =  ns("uiTopology"),
                          label = "Enter the
                                  circuit interactions",
                                         value =  ""),
                          
                          actionButton(ns("updateTopologyfromText"),
                                       "Load Circuit"),
                          htmlOutput(ns("networkTextFormat"))
                        )),
               tabPanel("Upload File",
                        fluidRow(
    #         checkboxInput(ns("headerTopology"), "Header", TRUE),
                          fileInput(ns("file"),"Upload circuit file"),
                          actionButton(ns("updateTopologyfromFile"),
                                       "Load Circuit"),
                          downloadButton(ns('downloadSampleNet'),
                                         'Download a sample circuit')
                        )
               )),
             
             downloadButton(ns('downloadRSet'), 'Download circuit')
      ),
      column(3, offset = 0,
             DTOutput(ns("tb"))
             
      ),
      column(6, offset = 0,
             (visNetworkOutput(ns("network")))
      ))
  )
}

# Module server function
shinyLoadNetwork <- function(input, output, session, stringsAsFactors) {
  
  f_tpo <- reactive({
    if(input$updateTopologyfromFile ==0) {
      if(input$updateTopologyfromText ==0) { return ()}
    }
    
    f_tpo1 <-  eventReactive(input$updateTopologyfromText, {
      tmp <- read.table(text=input$uiTopology,
                        col.names=c('Source','Target', 'Interaction'),
                        sep = ",", stringsAsFactors = FALSE) 
      return(tmp)
    })
    
    f_tpo2 <-  eventReactive(input$updateTopologyfromFile, {
      data <- input$file
      if(is.null(data)){return()}
      tmp <- read.table(data$datapath,sep="",
                        header =TRUE, stringsAsFactors = FALSE)
      return(tmp)
    })
    
    
    if(input$updateTopologyfromText){
      return(f_tpo1())
    }
    
    if(input$updateTopologyfromFile){
      updateTextAreaInput(session, "uiTopology", value = as.matrix(f_tpo2())  )
      return(f_tpo2())}
    
    
  })
  output$networkTextFormat <- renderUI({HTML(
    "Enter the circuit as comma separated values in a single line with no spaces.
    Use 1 for activation and 2 for inhibition. For example,
    src1,tgt1,1,src2,tgt2,2")})
  
  
  rs <- reactive({
    rs <-  new("racipeSet")
    network(rs) <- f_tpo()
    annotation(rs) <- input$filenameTopo
    return(rs)
  })
  
  
  
  
  output$tb <-renderDT({
    
    return(f_tpo())
  }, selection = 'none', editable = FALSE, rownames= FALSE
  )
  
  output$network <- renderVisNetwork({
    if(input$updateTopologyfromFile ==0) {
      if(input$updateTopologyfromText ==0) { return ()}
    }
    plotRSet(rs(),"network")
  })
  
  output$downloadRSet <- downloadHandler(
    filename = function() {
      paste(annotation(rs()), Sys.Date(), '.RDS', sep='')
    },
    content = function(con) {
      if(is.null(rs())) return(NULL)
      saveRDS(rs(), con)
    }
  )
  
  output$downloadSampleNet <- downloadHandler(
    filename = function() {
      paste("demo", '.csv', sep='')
    },
    content = function(con) {
      demo <- data.frame(Source = c("A", "B"), Target = c("B", "A"), 
                         Interaction = c("2", "2"))
      write.csv(demo, con)
    }
  )
  
  # Return the reactive that yields the data frame
  return(f_tpo)
}
