library(shinyjs)
source('allSessions.R')
source('geneExServer.R', local = TRUE)
source('geneExUI.R', local = TRUE)
shinyApp(
  ui = geneExUI,
  server = geneExServer
)



