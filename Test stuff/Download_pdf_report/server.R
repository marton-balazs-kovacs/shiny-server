library(shiny)

function(input, output) {
  
  regFormula <- reactive({
    as.character(input$x)
  })
  
  output$table <- renderTable(
    as.data.frame(matrix(cbind(c("Did you pre-regsiter your study", "Did you report all independent variables"),
                                                         c(input$x, input$y)), nrow = 2))
    )

  

  output$downloadReport <- downloadHandler(
    filename = function() {
      paste('my-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
    },
    
    content = function(file) {
      src <- normalizePath('report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd', overwrite = TRUE)
      
      library(rmarkdown)
      out <- render('report.Rmd', switch(
        input$format,
        PDF = pdf_document(), HTML = html_document(), Word = word_document()
      ))
      file.rename(out, file)
    }
  )
  
}