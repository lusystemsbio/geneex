validate <-  tabPanel("Validate",
                      tags$em("One can compare the simulated expression
with other simulated or experimental expression. The simulated expression can 
 be uploaded from a file or from the RACIPE tab. 
The reference data is clustered into the specified number of clusters. Then each 
 simulated data sample is compared 
with each sample of the cluster to find the cluster
 most similar to its expression pattern. Gene permutation is used
 to generate the null hypothesis. The percentage of
samples belonging to each cluster in the reference and simulated 
expressions as well as the overall Kullback–Leibler divergence between the 
two distributions is reported. 
Heatmaps of both simulated and reference expressions as well as sample-sample 
correaltion is also plotted. Cluster 0 is the null cluster and by default
                            we add one sample belonging to null cluster
                              to the reference samples"),
                      hr(),
                      fluidRow(
                        
                        column(3, offset = 0,
                               tags$h5("Upload Experimental Data (.csv)"),
                               fluidRow(
                                 column(8, offset = 1,
                                        fileInput(("fileRefExp"),
                                                  div(
                                                    downloadLink("downloadSampleValidateRefData",
                                                                 "Sample data")
                                                  ),
                                                  accept = c(
                                                    "text/csv",
                                                    "text/comma-separated-values,text/plain",
                                                    ".csv"))
                                 )
                                 
                                 
                               ),
                               tags$h5("Upload Simulated Data (.csv)"),
                               fluidRow(
                                 column(8, offset = 1,
                                        fileInput(("fileSimExp"),
                                                  div(
                                                    downloadLink("downloadSampleValidateSimData",
                                                                 "Sample data")
                                                  ),
                                                  accept = c(
                                                    "text/csv",
                                                    "text/comma-separated-values,text/plain",
                                                    ".csv"))
                                 )
                               ),
                               fluidRow(
                                 column(12, offset=1,
                                        fluidRow(       
                                          numericInput(
                                            inputId = "validatePValue", 
                                            "p-value",  min = 0, max = 1, 
                                            value = 0.05)
                                        ),
                                        fluidRow(
                                          numericInput(
                                            inputId = "validateNClust", "Number of Clusters",  min = 2, max = 8, value = 2)
                                        ),
                                        fluidRow(
                                          numericInput(
                                            inputId = "validateNPermut", 
                                            "Number of permutations",  
                                            min = 100, max = 10000, value = 1000)
                                        ),
                                        withBusyIndicatorUI(
                                          (actionButton(
                                            "compareValidate",  "Validate", 
                                            class = 'btn-primary',
                                            style="color: #fff; background-color: #32CD32; border-color: #2e6da4")))
                                 ))
                        ),
                        column(4, offset=0,
                               plotOutput("validateRefHeatmap")
                        ),
                        column(4,offset = 0,
                               plotOutput("validateSimHeatmap")
                        )
                        
                        
                      ),
                      hr(),
                      
                      textOutput('validateKL'),
                      hr(),
                      fluidRow(
                        
                        column(4, offset=0,
                               tableOutput('validateRefClustTable')),
                        column(4, offset=0,
                               tableOutput('validateSimClustTable'))
                      ),
                      plotOutput("validateRefSim"),
                      hr(),
                      hr()
)
