fluidPage(
  title = 'Download a PDF report',
  sidebarLayout(
    sidebarPanel(
      helpText(),
      selectInput('x', 'Did you pre-regsiter your study',
                  choices = c("yes", "no", "yes, but my dog ate it")),
      selectInput('y', 'Did you report all independent variables',
                  choices = c("yes", "no", "yes, but there were a few exploratory analyses that I did not report because they produced null results")),
      radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                   inline = TRUE),
      downloadButton('downloadReport')
    ),
    mainPanel(
      tableOutput('table')
    )
  )
)