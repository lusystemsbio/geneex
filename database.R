database <-
  tabPanel("Database",
           useShinyjs(),
           fluidRow(
             column(1,offset = 1,
           actionButton(
             "biologicalDB", "Biological",
             style="color: #fff; background-color: #337ab7;
                        border-color: #2e6da4", title = "Add this interaction
          to the circuit.
                      ")),
           column(1,offset = 1,
           actionButton(
             "syntheticDB", "Synthetic",
             style="color: #fff; background-color: #337ab7;
                        border-color: #2e6da4", title = "Add this interaction
          to the circuit.
                      ")),

           column(1,offset = 1,
           actionButton(
             "allDB", "All",
             style="color: #fff; background-color: #337ab7;
                        border-color: #2e6da4", title = "Add this interaction
          to the circuit.
                      "))
           ),

        br(),

           fluidRow(
             DTOutput("databaseTable")
           ),
        hr(),
bsAlert("dbAlert"),
hr(),

br(),
fluidRow(
  column(4,offset = 0,
        DTOutput("tableDbNetwork")
        ),
  column(7,offset = 0,
         hidden(actionButton("loadNetworkDatabase", "Load Circuit", 
style="color: #fff;background-color: #32CD32; border-color: #2e6da4")),
         visNetworkOutput("plotDbNetwork")
  )

  ),
shinyjs::hidden(downloadButton('downloadDbData', 'Download Data')),
shinyjs::hidden(radioButtons("downloadDbDataType", "Format",
                             c("RDS" = "RDS","CSV" = "csv") , 
                             selected = "RDS",
                             inline = TRUE)),

hr(),
hr()
  )
