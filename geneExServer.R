geneExServer <- shinyServer(function(input, output, session) {
  source('circuitServer.R', local = TRUE)
  source('modelExplorerServer.R', local = TRUE)
  source('racipeServer.R', local = TRUE)
  source('validateServer.R', local = TRUE)
  source('databaseServer.R', local = TRUE)
  source('sticccServer.R', local = TRUE)
  source('aboutServer.R', local = TRUE)
  source('hideShow.R', local = TRUE)
})

  