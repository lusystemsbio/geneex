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
                             ( sliderInput(inputId = "downsamplingProportion", "Down-sampling Proportion",  min = 0.01, 
                                           max = 1, value = 0.1)))),
                    fluidRow(
                      column(9, offset=1,
                             checkboxInput(inputId = "gridSmoothing", 
                                           "Grid Smoothing",  value = T))),
                    fluidRow(
                      column(9, offset=1,
                             ( checkboxInput(inputId = "plotLoadings", "Plot PCA Loadings", 
                                             value = F)),
                      )),
                    fluidRow(
                      column(9, offset=1,
                             withBusyIndicatorUI(actionButton("computeTrajectories", "Compute Trajectories", 
                                                              class = 'btn-primary',
                                                              style="color: #fff; background-color: #32CD32; border-color: #2e6da4",
                                                              title = "Simulate the circuit with given parameters")),
                      )),
                    fluidRow(
                      column(9, offset=1,
                             withBusyIndicatorUI(actionButton("plotTrajectories", "Plot Trajectories", 
                                                              class = 'btn-primary',
                                                              style="color: #fff; background-color: #32CD32; border-color: #2e6da4",
                                                              title = "Simulate the circuit with given parameters")),
                      )),
                    fluidRow(
                      column(9, offset=1,
                             numericInput(inputId = "scalingFactor", "Vector Scaling Factor (0.01-1000)",  min = 0.01, 
                                         max = 1000, value = 1)
                      )),
                    ),
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