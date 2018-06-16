library(shiny)
library(shinydashboard)


shinyUI(
  dashboardPage(
    dashboardHeader(title = "TPP Pilot"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Welcome", tabName = "tab_welcome"),
        menuItem("Results summary" , tabName = "res_summary"),
        menuSubItem("Main Confirmatory Analysis", tabName = "conf_analysis"),
        menuSubItem("Robustness Analysis", tabName = "robust_analysis"),
        menuSubItem("Exploratory analysis", tabName = "explor_analysis")
      )
    ),
    
    dashboardBody(
      tabItems(
        tabItem(tabName = "tab_welcome",
                
                fluidRow(
                  h1("Transparent Psi Project - Pilot data results")
                ),
                fluidRow(
                  p("You can follow the results of the Transparent Psi Project here in real time")
                ),
                fluidRow(
                  p("Use the menu items on the left to get the results of the different analyses")
                ),
                fluidRow(
                  p("the data is refreshed every 2 seconds")
                ),
                fluidRow(
                  box(plotOutput("plot_refresh"))
                )
        ),
        
        tabItem(tabName = "res_summary",
                
                fluidRow(
                  h1("Summary of the results")
                ),
                fluidRow(
                  textOutput("text_summary")
                ),
                fluidRow(
                  box(plotOutput("plot1a"))
                ),
                fluidRow(
                  p("To support any model, all three Bayes Factor values need to pass the threshold")
                )
        ),
        
        
        tabItem(tabName = "conf_analysis",
                
                fluidRow(
                  h1("Main Confirmatory Analysis results")
                ),
                fluidRow(
                  box(plotOutput("plot1b"))
                  ),
                fluidRow(
                  p("To support any model, all three Bayes Factor values need to pass the threshold")
                )
        ),
        
        tabItem(tabName = "robust_analysis",
                
                fluidRow(
                  h1("Result of the Bayesian Robustness Parameter Estimation Robustness Analysis")
                ),
                fluidRow(
                  box(plotOutput("plot2"))
                )
        ),
        
        tabItem(tabName = "explor_analysis",
                
                fluidRow(
                  h1("Histogram overlay of the expected and observed distribution of successful guess rate")
                ),
                fluidRow(
                  box(plotOutput("plot3"))
                )

        )
      )
      
      
      
    )
  )
)
