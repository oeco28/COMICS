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
  titlePanel("Integrating selection Scans via C.O.M.I.C.S."),
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = "file1", label = "Data file",
                accept = c("text/csv", "text/tab-seperated-vlaues,text/plain", ".csv"),
                width = NULL, buttonLabel = "Browse",
                placeholder = "No file selected"),
      fileInput(inputId = "configuration", label = "Genome Configuration File",
                accept = c("text/csv", "text/tab-seperated-values,text/plain", ".csv"),
                width = NULL, buttonLabel = "Browse",
                placeholder = "No file selected"),
      numericInput("Cutoff", label = "Statistical Cutoff", value = 5, min = 0, max = 10),
      selectInput("dataset", "Choose a dataset",
                  choices = c("ICS Distance")),
      numericInput("TestOfInterest", label = "Test of interest", value = 1, min = 0, max = 100),
      numericInput("Chromosomes", label = "Chromosome of interest", value = 1, min = 0, max = 100),
      numericInput("First.Index", label = "First Index", value = 1, min = 1, max = 100),
      numericInput("Second.Index", label = "Second Index", value = 1, min = 1, max = 100),
      downloadButton("downloadData", "Download ICS output"),
      downloadButton("downloadPlot", "Download ICS figure")
    ),
    
    mainPanel(
      tabsetPanel( # Added tabsetPanel
        id = "plotTabs", #id for tabset
        tabPanel("Summary",
                 h4("Summary.ICS"),
                 verbatimTextOutput("summary.ICS"),
                 h3("Summary.SingleTest"),
                 verbatimTextOutput("summary.SingleTest"),
                 verbatimTextOutput("nText")
                 ),
        tabPanel("ICS Plots",
                 plotOutput("ICS.hist"),
                 plotOutput("ICS.chromosome.hist")),
        tabPanel("Genome Scan", plotOutput("GenomeScan")),
        tabPanel("Test Scans",
                 plotOutput("GenomeTest"),
                 plotOutput("GenomeMelt"))
      )
    )
  )
)
