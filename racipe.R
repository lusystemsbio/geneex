source("modules.R")
racipe <-
  tabPanel("RACIPE",
           useShinyjs(),

############## RACIPE SIMULATIONS ##########################
  fluidRow(
    column(5, offset=0,
         fluidRow(
           column(9, offset=1,
           ( numericInput(inputId = "numModels", "Number of Models",  min = 1, 
                          max = 5000, value = 500.0)))),
         fluidRow(
           column(9, offset=1,
                  numericInput(inputId = "racipeNClusters", 
                               "Expected Number of Clusters",  min = 2, step = 1,
                               max = 10, value = 2))),
           fluidRow(
             column(9, offset=1,
           ( numericInput(inputId = "simTimeRacipe", "Simulation Time", 
                          min = 1, max = 5000, value = 50.0)))),
           fluidRow(
             column(9, offset=1,
           ( numericInput(inputId = "stepSizeRacipe", "Integration Time Step",
                          min = 0.001, max = 0.9, value = 0.05)),
           withBusyIndicatorUI(actionButton("simulateRacipe", "Simulate", 
                                               class = 'btn-primary',
    style="color: #fff; background-color: #32CD32; border-color: #2e6da4",
    title = "Simulate the circuit with given parameters"))))),
  column(6, offset=0,
  visNetworkOutput("racipeCircuit"))
),
# bsAlert("racipeAlert"),
       fluidRow(
       htmlOutput("racipeDeterministicText")),

       fluidRow(
column(12, offset = 1,
         shinyjs::hidden(plotOutput("racipePca", width = "500px", height = "500px")))),
         fluidRow(
         shinyjs::hidden(plotOutput("racipeHeatmap"))),
fluidRow(
column(3, offset=8,
       shinyjs::hidden(actionButton("validateRacipe", "Load for Validation", 
                           style="color: #fff; background-color: #32CD32; 
                           border-color: #2e6da4", title = "Use this 
                           simulated data for comparison with experimental
                           data in the validation tab")))
),
################### RACIPE DATA DOWNLOAD #########################

fluidRow(
  column(3, offset=1,
         fluidRow(
           hidden(downloadButton('downloadRacipeData', 'Download Data'))),
         fluidRow(
           hidden(radioButtons("downloadRacipeDataType", "Format",
                               c("RDS" = "RDS","CSV" = "csv") , 
                               selected = "csv",
                               inline = TRUE))))

),

           hr(),
################# PARAMETERIC ANALYSIS ###########################
       fluidRow(
         column(5, offset=4,
                shinyjs::hidden(  actionButton("parametricAnalysisRacipe", 
                                      "Parametric Analysis", 
style="color: #fff;background-color: #337ab7; border-color: #2e6da4", 
title = "Display the distributions when the parameter values are limited to 
a given percentage of the original distribution.")))
       ),
       fluidRow(
         column(3,
                shinyjs::hidden( uiOutput("filteredOutputRacipe")),
                shinyjs::hidden(uiOutput("filterSliderRacipe"))
),
column(3,
       hidden(uiOutput("filteredOutputRacipe2")),
       shinyjs::hidden( uiOutput("filterSliderRacipe2"))),
column(3,

       shinyjs::hidden(  uiOutput("filteredOutputRacipe3")),
       shinyjs::hidden(uiOutput("filterSliderRacipe3"))
)),


           fluidRow(
             shinyjs::hidden(    plotOutput("racipePcaFiltered", width = "500px", height = "500px")),
          shinyjs::hidden(  plotOutput("racipeHeatmapFiltered"))
           ),
fluidRow(
  column(3, offset=1,
         fluidRow(
           hidden(downloadButton('downloadRacipeParData', 'Download Data'))),
         fluidRow(
           hidden(radioButtons("downloadRacipeParDataType", "Format",
                               c("RDS" = "RDS","CSV" = "csv") , 
                               selected = "csv",
                               inline = TRUE)))
         )
  
),
fluidRow(
  column(3, offset=8,
         shinyjs::hidden(actionButton("validateRacipePar", "Load for Validation", 
                                      style="color: #fff; background-color: #32CD32; 
                           border-color: #2e6da4", title = "Use this 
                           simulated data for comparison with experimental
                           data in the validation tab")))
),
hr(),
################# STOCHASTIC RACIPE ##############################
fluidRow(
  column(5, offset=4,
         shinyjs::hidden(actionButton("stochasticRacipe", "Stochastic RACIPE", 
                       style="color: #fff; background-color: #337ab7; 
                       border-color: #2e6da4", title = "Displays stochastic
                simulation options")))
),
fluidRow(
  column(3,
         shinyjs::hidden(radioButtons("sRacipeOption", "Stochastic Simulation Type:",
                      c("Constant Noise" = "constantNoise",
                        "Annealing" = "annealing")))
  ),
  column(3, offset=0,
         shinyjs::hidden(  sliderInput("sRacipeNoise", "Noise Level",step = 1,  min = 1, 
                              max = 20, value = 0))),

  column(5, offset=0,
         br(),
         withBusyIndicatorUI(     shinyjs::hidden(  actionButton("simulateSRacipe", 
      "Perform Stochastic Simulations",  class = 'btn-primary',
                               style="color: #fff; background-color: #32CD32; 
                           border-color: #2e6da4"))))
),
fluidRow(
  column(6,
         shinyjs::hidden(plotOutput("sRacipePcaDet", width = "500px", height = "500px"))
  ),
column(6,
       shinyjs::hidden(plotOutput("sRacipePca", width = "500px", height = "500px"))
)),
fluidRow(
  column(6,
hidden(plotOutput("sRacipeHeatmapDet"))),
column(6,
       shinyjs::hidden(plotOutput("sRacipeHeatmap")))),

fluidRow(
  column(3, offset=8,
         shinyjs::hidden(actionButton("validateRacipeStoch", "Load for Validation", 
                                      style="color: #fff; background-color: #32CD32; 
                           border-color: #2e6da4", title = "Use this 
                           simulated data for comparison with experimental
                           data in the validation tab")))
),
################### SRACIPE DATA DOWNLOAD #########################

fluidRow(
  column(3, offset=1,
         fluidRow(
           shinyjs::hidden(downloadButton('downloadSRacipeData', 'Download Data'))),
         fluidRow(
           shinyjs::hidden(radioButtons("downloadSRacipeDataType", "Format",
                               c("RDS" = "RDS","CSV" = "csv") , 
                               selected = "csv",
                               inline = TRUE)))),

  column(3, offset=0,
         shinyjs::hidden(actionButton("saveSRacipeData", "Upload Options")))

),
hr()
  )

