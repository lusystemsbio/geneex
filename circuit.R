############## Network Input ##################################################

circuit <-  tabPanel(
  "Circuit", 
  tags$em(
    "The network/circuit to be simulated can be specified here.
There are multiple ways to do so. As an example, we have
loaded the canonical toggle switch circuit as the default circuit. 
It can be modified either by loading the circuit from a file containing 
three columns, 'Source', 'Target', 'Interaction' (1 - activation, 2 - inhibition
) separated by space and each row
specifying an interaction. One can also 'add', 'delete', 'undo delete'
interactions
in the table below. Clicking the 'Load Circuit' button will load the circuit 
and display it in the right panel. Database tab can be used to specify a circuit 
if one of the circuits from the database is to be loaded. Use 'Circuit Name' to
provide a specific name to your circuit.
'Download Circuit' will download the loaded circuit as a text file."),
  hr(),
  fluidRow(
    column(5,offset = 1,
           fluidRow(
             textInput("circuitName", "Circuit Name (optional)", value = "Circuit1"),
             br(),
             
             fileInput(("circuitFile"), #label = "Choose Circuit File",
                       div("Choose Circuit File (optional)",
                           br(),
                           downloadLink("downloadSampleCircuit",
                                        "Example circuit file")
                       ),
                       accept = c(
                         "text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv",".txt","tsv") 
             ),
             shinyBS::tipify(downloadButton(('downloadCircuit'), 
                                            'Download Circuit'),
                             "Download the circuit as a text file."),
             br(),
             br(),
             
             
             actionButton( inputId = "updateCircuit", label = "Load Circuit",
                           icon("thumbs-up"), 
style="color: #fff; background-color: #32CD32;
           border-color: #2e6da4",
                           title =
                             "If the interactions look good
           you can load the circuit for further
           analysis. It will show the circuit figure below this button.
           "
             )
             
           ),
           br(),
           
           fluidRow( column(10,offset = 0,
                            visNetwork::visNetworkOutput("circuit")
           ))
    ),
    column(5,offset = 0,
           
           
           fluidRow(
             DT::dataTableOutput("circuitTable"),
             shinyBS::bsTooltip(
               "circuitTable", "Current Circuit interactions",
               placement = "bottom", trigger = "hover", options = NULL)
           ),
           br(),
           
           fluidRow(
             column(1,offset = 0,
                    shinyBS::tipify(uiOutput('undoUI'),
                                    "Put back the deleted interaction")
             )
           ),
           br(),
           fluidRow(
             column(1,offset = 0,
                    
                    actionButton(
                      "toggleAddInteraction", "Add Interaction",
                      icon("plus-circle"),
                      style="color: #fff; background-color: #337ab7;
                        border-color: #2e6da4", title = "Add an interaction
          to the circuit.
                      ")
             )
           ),
           fluidRow(
             column(3,offset = 0,
                    
                    shinyjs::hidden(
                      textInput("newIntSrc", "Source",value = "srcGene"))
             ),
             column(3,offset = 0,
                    
                    shinyjs::hidden( textInput(
                      "newIntTgt", "Target",value = "tgtGene"))
             ),
             column(3,offset = 0,
                    
                    shinyjs::hidden(numericInput(
                      "newIntType", "Interaction",min = 1,max = 2, value = 1))
             ),        
             column(3,offset = 0,
                    
                    shinyjs::hidden(actionButton(
                      "addInteraction", "Add",
                      style="color: #fff; background-color: #337ab7;
                        border-color: #2e6da4", title = "Add this interaction
          to the circuit.
                      ")
                    ))
             
           ),       
           fluidRow(
             column(3,offset = 6,
                    actionButton( inputId = "updateCircuitTable", label = "Update Circuit",
                                  icon("thumbs-up"), style="color: #fff; background-color: #32CD32;
           border-color: #2e6da4",
                                  title =
                                    "If the interactions look good
           you can load the circuit for further
           analysis. It will show the circuit figure below this button.
           "
                    )
             )
           ),
           br()
    )
  ),
  hr(),
  hr()
)
