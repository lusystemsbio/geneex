output$aboutMDFile <- renderUI({
  if(is.null(input$aboutSelect)) {
    
    file <- "md/aboutOverview.md"
  }  else {
  file <- switch(input$aboutSelect,
                 Overview = "md/aboutOverview.md",
                 Background = "md/aboutBackground.md",
                 Theory = "md/aboutTheory.md",
                 Circuit = "md/aboutCircuit.md",
                 GeneVyuha = "md/aboutGeneVyuha.md",
                 RACIPE = "md/aboutRacipe.md",
                 Validate = "md/aboutValidate.md",
                 Database = "md/aboutDatabase.md"
  )
  }
  includeMarkdown(file)
})