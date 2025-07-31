#' UI Components for COMICS Shiny Application
#'
#' This file defines the user interface for the COMICS (Combined Outlier Method for
#' Identifying Candidate Signals) Shiny application.
#'
#' @importFrom shiny fluidPage titlePanel sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny fileInput numericInput selectInput downloadButton
#' @importFrom shiny h4 h3 verbatimTextOutput tabsetPanel tabPanel plotOutput
#' @importFrom bslib bs_theme
#' @importFrom shinyjs useShinyjs
#'
#' @noRd
ui <- fluidPage(
  useShinyjs(), # Enable shinyjs
  theme = bslib::bs_theme(bootswatch = "flatly"),
  titlePanel("Integrating selection Scans via C.O.M.I.C.S. (Standard Deviation Method)"),
  sidebarLayout(
    sidebarPanel(
      # File inputs
      h4("Data Files", style = "color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 5px;"),
      fileInput(inputId = "file1", label = "Data file",
                accept = c("text/csv", "text/tab-seperated-vlaues,text/plain", ".csv"),
                width = NULL, buttonLabel = "Browse",
                placeholder = "No file selected"),
      fileInput(inputId = "configuration", label = "Genome Configuration File",
                accept = c("text/csv", "text/tab-seperated-values,text/plain", ".csv"),
                width = NULL, buttonLabel = "Browse",
                placeholder = "No file selected"),
      
      # ICS Method Parameters
      h4("ICS Method Parameters", style = "color: #2c3e50; border-bottom: 2px solid #e74c3c; padding-bottom: 5px;"),
      selectInput("scatter_method", "Scatter Matrix Combination",
                  choices = list(
                    "Covariance + 4th Moment" = "cov_cov4",
                    "Covariance + Cubic Covariance" = "cov_cov3",
                    "Covariance + Correlation" = "cov_cor",
                    "Correlation + 4th Moment" = "cor_cov4"
                  ),
                  selected = "cov_cov4"),
      
      selectInput("distance_method", "Distance Calculation Method",
                  choices = list(
                    "Euclidean Distance" = "euclidean",
                    "Component-wise (Best Component)" = "component_wise",
                    "Multi-component (Top 2)" = "multi_component"
                  ),
                  selected = "euclidean"),
      
      # Standard Deviation-based threshold
      h4("Threshold Parameters", style = "color: #2c3e50; border-bottom: 2px solid #f39c12; padding-bottom: 5px;"),
      sliderInput("sd_multiplier",
                  label = HTML("Standard Deviation Multiplier<br/><small>(Threshold = μ + k×σ)</small>"),
                  value = 2.5, min = 1.0, max = 5.0, step = 0.1,
                  post = "σ"),
      
      # Analysis Parameters
      h4("Analysis Parameters", style = "color: #2c3e50; border-bottom: 2px solid #9b59b6; padding-bottom: 5px;"),
      selectInput("dataset", "Choose a dataset",
                  choices = c("ICS Distance")),
      numericInput("TestOfInterest", label = "Test of interest", value = 1, min = 0, max = 100),
      numericInput("Chromosomes", label = "Chromosome of interest", value = 1, min = 0, max = 100),
      numericInput("First.Index", label = "First Index", value = 1, min = 1, max = 100),
      numericInput("Second.Index", label = "Second Index", value = 1, min = 1, max = 100),
      
      # Download buttons
      h4("Export Results", style = "color: #2c3e50; border-bottom: 2px solid #27ae60; padding-bottom: 5px;"),
      downloadButton("downloadData", "Download ICS Output",
                     style = "background-color: #3498db; color: white; margin-bottom: 10px;"),
      br(),
      downloadButton("downloadPlot", "Download ICS Figure",
                     style = "background-color: #e74c3c; color: white;"),
      
      # Help text
      br(), br(),
      div(
        style = "background-color: #ecf0f1; padding: 10px; border-radius: 5px; margin-top: 10px;",
        h5("Method Information:", style = "color: #2c3e50;"),
        p("This app uses the optimized Standard Deviation-based ICS method for outlier detection.",
          style = "font-size: 12px;"),
        p("Data is automatically standardized before analysis.",
          style = "font-size: 12px; color: #7f8c8d;"),
        p("Threshold = Mean + (Multiplier × SD)",
          style = "font-size: 12px; font-weight: bold; color: #e74c3c;"),
        p("Higher multipliers = more conservative detection.",
          style = "font-size: 12px; color: #7f8c8d;")
      )
