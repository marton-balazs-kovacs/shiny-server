library(shiny)
library(knitr)
library(rmarkdown)

shinyUI(
  fluidPage(
    p("Loading the contents of this page might take a few seconds, please wait."),
    uiOutput('markdown')
  )

  )
