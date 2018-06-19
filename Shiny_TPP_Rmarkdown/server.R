library(shiny)
library(knitr)
library(rmarkdown)

shinyServer(function(input, output){
  output$markdown <- renderUI({
    HTML(markdown::markdownToHTML(knit('TPP_Manuscript_test.Rmd', quiet = TRUE)))
  })
})

