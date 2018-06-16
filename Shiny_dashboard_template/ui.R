library(shiny)
library(shinydashboard)


shinyUI(
  dashboardPage(
    dashboardHeader(title = "This is the header"),
    dashboardSidebar(
      menuItem("Dashboard"),
      menuSubItem("Dashboard Finance"),
      menuSubItem("Dashboard Sales"),
      menuItem("Details"),
      menuItem("Raw data")
    ),
    dashboardBody()
  )
)
