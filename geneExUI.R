
library(shinyjs)
source('circuit.R', local = TRUE)
source('modelExplorer.R', local = TRUE)
source('racipe.R', local = TRUE)
source('validate.R', local = TRUE)
source('database.R', local = TRUE)
source('sticcc.R', local = TRUE)
source('about.R', local = TRUE)
source('forum.R', local = TRUE)
source('allSessions.R', local = TRUE)

#geneExUI <- fluidPage(tags$div(tags$img(src="logo4.png", height = 40, width = 250, style="margin-right:200px;"),
#                            tags$img(src="JAX.png", height = 40, width = 170, style="float:center; margin-left:0px;margin-right:5px;"),
#                             tags$img(src="Rice_CTBP_2Color_Logo.jpg", height = 40, width = 180, style="float:center; margin-left:0px;margin-right:5px;")
geneExUI <- fluidPage(tags$div(tags$img(src="logo4.png", height = 40, width = 250, style="margin-right:200px;"),
                            tags$img(src="lulab.png", height = 40, width = 170, style="float:center; margin-left:0px;margin-right:5px;"),
                             tags$img(src="neu_black.svg", height = 40, width = 180, style="float:center; margin-left:0px;margin-right:5px;")
),
  navbarPage(title="",
windowTitle = "GeneEx",
                       tabPanel("Circuit",circuit),
                       tabPanel("GeneVyuha", modelExplorer),
                       tabPanel("RACIPE", racipe),
                       tabPanel("Validate",validate),
                       tabPanel("Database",database),
                       tabPanel("STICCC",sticcc),
                       tabPanel("About", about),
                       tabPanel("Forum", forum),
                       tags$head(
                          tags$style(type = 'text/css',
                                     HTML('.navbar { background-color: white;
                                          font-size:      1.5em;}
                                  .navbar-default .navbar-brand{color: #0085CA;
                                            font-size:      1.5em;}
                               .tab-panel{ background-color: white; color: #2C8FC8
                                            font-size:      1.5em;}
                               .navbar-default .navbar-nav > .active > a,
                           .navbar-default .navbar-nav > .active > a:focus,
                           .navbar-default .navbar-nav > .active > a:hover {
                                          color: #002D72;
                                          background-color: white;
                                          }')
                                    )
                          )
)
)
