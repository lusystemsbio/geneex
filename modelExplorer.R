source("modules.R")
modelExplorer <-  tabPanel(
  "Model Explorer",        tags$em("
'GeneVyuha' tab simulates a circuit with a specific parameter set. 
        Once the circuit is loaded from the 'Circuit' or 'Database' tab, a 
        random set of parameters is generated. If the circuit from the 
database has its own parameter set, that set is used. The parameter set contains
two paramters for each gene and three parameters for each interactions. 
The gene parameters are production (G_'gene') and degradation (k_'gene') rate 
whereas the parameters for each interaction ('source'_'target') are 'threshold'
(TH_), 'hill coefficient of cooperativity' and 'fold change'. 
These parameters can be modified using the dropdown box given below. 
        The default value of the parameter is dispalyed when a parameter is 
        selected. Edit the 'value' and click 'Update' to modify this value. 
        Repeat  these steps if other parameters are to be modified. Other 
        simulation paramters can also be modified."),
  br(),
  useShinyjs(),
  fluidRow(
    column(6, offset = 0,
           
           uiOutput("modelParams"),
           fluidRow(
             column(4, offset = 0,
                    style = "margin-top:-1em",
                    uiOutput("newModelParamValue")),
             column(3, offset = 0, style = "margin-top: 10px",
                    actionButton("updateGvParam",   "Update",  width = "100%",
                                 style="color: #fff;background-color: #337ab7; border-color: #2e6da4",
                                 title = "Modify the parameters using this button. 
  Enter the new value in the 'value' input and 'update'. ")
             )
           ),
           
           uiOutput("modelIc"),
           fluidRow(
             column(4, offset = 0,
                    style = "margin-top:-1em",
                    uiOutput("newModelIcValue")),
             column(3, offset = 0, style = "margin-top: 10px",
                    actionButton("updateGvIc",   "Update",  width = "100%",
                                 style="color: #fff;background-color: #337ab7; border-color: #2e6da4",
                                 title = "Modify the initial conditions using this button. 
  Enter the new value in the 'value' input and 'update'. ")
             )
           ),
           
           br(),
           shinyBS::tipify(numericInput(inputId = "simTimeExplorer", 
                                        "Simulation Time",  min = 1, 
                                        max = 5000, value = 50.0),"Simulation time after which the 
               gene expression values will be recorded. 
               Use a larger value if you 
expect longer transient dynamics."),
           shinyBS::tipify(numericInput(inputId = "stepSizeExplorer",
                                        "Integration Time Step", min = 0.001, max = 0.9,
                                        value = 0.1),
                           "Time step to be used in numerical integration methods."),
           
           numericInput("noiseLevel", 
                        "Noise",step = 0.1,  min = 0, max = 100, value = 0),
           
           br(),

           withBusyIndicatorUI(
             
             actionButton("simulateGv", "Simulate",
                          class = "btn-primary",
                          style="color: #fff; background-color: #32CD32; border-color: #2e6da4",
                          title = "Simulate the circuit using the above parameter values. 
    Change the circuit using the Circuit tab."))),
    column(5, offset = 0,
           visNetworkOutput("circuitGv")
    ),
    column(1, offset = 0 )
  ),
  hr(),
  fluidRow(
    column(10, offset = 0,
           shinyjs::hidden(  uiOutput("modelPlotGene")),
           shinyjs::hidden(  checkboxInput("logPlotGV","Logscale")),
           shinyjs::hidden( plotOutput("GvTS"))
    )
  ),
  

  fluidRow(
    column(3, offset=1,
           fluidRow(
             hidden(downloadButton('downloadMEData', 'Download Data'))),
           fluidRow(
             hidden(radioButtons("downloadMEDataType", "Format",
                                 c("RDS" = "RDS","CSV" = "csv") , 
                                 selected = "RDS",
                                 inline = TRUE)))),
    
  ),


  hr(),
  # ############## Bifurcations ##############################################
  
  fluidRow(
    
    column(4, offset=4,
           hidden(actionButton("bifurcationExplorer", "Parameter Perturbation", 
                               style="color: #fff; background-color: #337ab7; border-color: #2e6da4")))
  ),
  
  br(),
  fluidRow(
    column(3, offset=0,
           shinyjs::hidden(uiOutput("modelParamsBif")),
           shinyjs::hidden(uiOutput("modelParamBifMin")),
           shinyjs::hidden( uiOutput("modelParamBifMax")),
           shinyjs::hidden( numericInput("modelNumBifurs", "Simulation Points", 300,
                                         min = 50, max = 5000, step = 50)),
           
           withBusyIndicatorUI( shinyjs::hidden(  actionButton("bifurcationME", "Simulate",
                                                               class = 'btn-primary',
                                                               style="color: #fff; background-color: #32CD32; border-color: #2e6da4",
                                                               title = "Simulate circuit to generate the gene expressions for models with
varying selected parameter.")))
    ),
    column(9, offset=0,
           
           shinyjs::hidden(  uiOutput("modelPlotGeneBif")),
           shinyjs::hidden(  checkboxInput("logPlotGVBif","Logscale")),
           shinyjs::hidden(  plotOutput("modifiedBifME")),
           shinyjs::hidden(downloadButton('downloadBifData', 'Download Bifurcation Data'))
    )),
  
  hr(),
  hr()
)
