
library(markdown)

about <-
  tabPanel("About",
           fluidPage(
             column(3,offset = 0, 
                    tags$style("color: #337ab7; 
                               font-size: 125%;font-family: sans-serif;
                               text-align: left; 
                               icon: thumbs-up}"),
                    selectInput("aboutSelect", label = "", c(
                      "Overview", "Background", "Theory","Circuit","GeneVyuha",
                      "RACIPE","Validate","Database"),
                      selected = 4, multiple = F,
                      selectize = F, width = '100%', size = 8)),
             column(8,offset = 0,
                    uiOutput("aboutMDFile"),
                    hr(),
                    tags$p(
                      "The development of the web-app is supported by a startup fund
                  from The Jackson Laboratory and by the National Institutes of 
                  Health under Award Number [P30CA034196] and [R35GM128717].
                  Work at the the Center for Theoretical 
                  Biological Physics sponsored by the National Science 
                  Foundation NSF Grant PHY-1427654."),
                    hr(),
                    tags$p(
                      "Please report any feature requests, bugs or concerns to 
                  Vivek Kohar, vivek.kohar@jax.org and Mingyang Lu, 
                  mingyang.lu@jax.org ")
             )
           )
  )
