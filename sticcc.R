source("modules.R")
sticcc <-
  tabPanel("STICCC",
           tags$em(
             "PLACEHOLDER TEXT. THIS WILL BE A BRIEF INTRODUCTION TO THE FUNCTION OF STICCC AND HOW TO USE IT."),
           br(),
           useShinyjs(),
           
           ############## STICCC INPUTS ##########################
           fluidRow(
             column(5, offset=0,
                    fluidRow(
                      column(9, offset=1,
                             ( sliderInput(inputId = "samplingRadius", "Sampling Radius",  min = 0.01, 
                                           max = 1, value = 0.15)))),
                    fluidRow(
                      column(9, offset=1,
                             checkboxInput(inputId = "gridSmoothing", 
                                           "Grid Smoothing",  value = T))),
                    fluidRow(
                      column(9, offset=1,
                             ( checkboxInput(inputId = "plotLoadings", "Plot PCA Loadings", 
                                             value = T)),
                      withBusyIndicatorUI(actionButton("computeTrajectories", "Compute Trajectories", 
                                                       class = 'btn-primary',
                                                       style="color: #fff; background-color: #32CD32; border-color: #2e6da4",
                                                       title = "Simulate the circuit with given parameters")),
                      ))),
             column(6, offset=0,
                   visNetworkOutput("sticCircuit"))
           ),
           fluidRow(
             htmlOutput("sticText")),
           
           fluidRow(
             column(12, offset = 1,
                    shinyjs::hidden(plotOutput("sticPCA", width = "800px", height = "800px")))),

           
hr()
)