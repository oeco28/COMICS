#' COMICS Shiny Application
#'
#' This file loads the UI and server components for the COMICS app
#' and runs the Shiny application.
#'

# Load required libraries
library(shiny)
library(ggplot2)
library(tidyr)
library(bslib)
library(mvtnorm)
library(ICS)
library(moments)
library(ICSOutlier)
library(shinyjs)

# Source the UI and server components
ui <- COMICS:::ui
server <- COMICS:::server

# Run the application
shinyApp(ui = ui, server = server)
