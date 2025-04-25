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
ui <- shiny::fluidPage(
  shinyjs::useShinyjs(), # Enable shinyjs
  theme = bslib::bs_theme(bootswatch = "flatly"),
  shiny::titlePanel("Integrating selection Scans via C.O.M.I.C.S."),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::fileInput(inputId = "file1", label = "Data file",
                accept = c("text/csv", "text/tab-seperated-vlaues,text/plain", ".csv"),
                width = NULL, buttonLabel = "Browse",
                placeholder = "No file selected"),
      shiny::fileInput(inputId = "configuration", label = "Genome Configuration File",
                accept = c("text/csv", "text/tab-seperated-values,text/plain", ".csv"),
                width = NULL, buttonLabel = "Browse",
                placeholder = "No file selected"),
      shiny::numericInput("Cutoff", label = "Statistical Cutoff", value = 5, min = 0, max = 10),
      shiny::selectInput("dataset", "Choose a dataset",
                  choices = c("ICS Distance")),
      shiny::numericInput("TestOfInterest", label = "Test of interest", value = 1, min = 0, max = 100),
      shiny::numericInput("Chromosomes", label = "Chromosome of interest", value = 1, min = 0, max = 100),
      shiny::numericInput("First.Index", label = "First Index", value = 1, min = 1, max = 100),
      shiny::numericInput("Second.Index", label = "Second Index", value = 1, min = 1, max = 100),
      shiny::downloadButton("downloadData", "Download ICS output"),
      shiny::downloadButton("downloadPlot", "Download ICS figure")
    ),
    
    shiny::mainPanel(
      shiny::tabsetPanel(
        id = "plotTabs",
        shiny::tabPanel("Summary",
                 shiny::h4("Summary.ICS"),
                 shiny::verbatimTextOutput("summary.ICS"),
                 shiny::h3("Summary.SingleTest"),
                 shiny::verbatimTextOutput("summary.SingleTest"),
                 shiny::verbatimTextOutput("nText")
        ),
        shiny::tabPanel("ICS Plots",
                 shiny::plotOutput("ICS.hist"),
                 shiny::plotOutput("ICS.chromosome.hist")),
        shiny::tabPanel("Genome Scan", shiny::plotOutput("GenomeScan")),
        shiny::tabPanel("Test Scans",
                 shiny::plotOutput("GenomeTest"),
                 shiny::plotOutput("GenomeMelt"))
      )
    )
  )
)
