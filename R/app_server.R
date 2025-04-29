#' Server Logic for COMICS Shiny Application
#'
#' This file defines the server logic for the COMICS (Combined Outlier Method for
#' Identifying Candidate Signals) Shiny application.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#'
#' @importFrom shiny renderPlot renderPrint downloadHandler observe
#' @importFrom ggplot2 ggplot aes geom_point geom_histogram geom_vline scale_color_manual theme element_text theme_minimal labs ggtitle xlab ylab geom_hline guides ggsave
#' @importFrom tidyr pivot_longer all_of
#' @importFrom ICS ics2 ics.distances
#'
#' @noRd
#' This the most finished version up to 04/2025


server <- function(input, output, session) {
  
  # Reactive data and calculations
  data_objects <- reactive({
    # Check if files are uploaded
    file1 = input$file1
    configuration = input$configuration
    
    if (is.null(file1) || is.null(configuration)) {
      return(list(
        dataX = NULL,
        chr.length = NULL,
        X_for_ics = NULL,
        Z2 = NULL
      ))
    }
    
    # Read data
    dataX = read.table(file1$datapath, header = TRUE)
    chr.length = read.table(configuration$datapath, header = TRUE)
    
    # ICS Analysis preparation
    X_for_ics <- dataX
    X_for_ics[, 1:2] <- NULL # Remove the first two columns
    
    list(
      dataX = dataX,
      chr.length = chr.length,
      X_for_ics = X_for_ics,
      Z2 = NULL # Initialize Z2 here, will be calculated later if needed
    )
  })
  
  # Update max chromosomes
  observe({
    data_objs <- data_objects()
    chr.length <- data_objs$chr.length
    if(!is.null(chr.length)){
      updateNumericInput(session, "Chromosomes",
                         max = nrow(chr.length),
                         value = min(input$Chromosomes, nrow(chr.length)))
    }
  })
  
  # Reactive value for ICS analysis
  ics_analysis_results <- reactive({
    data_objs <- data_objects()
    dataX <- data_objs$dataX
    chr.length <- data_objs$chr.length
    X_for_ics <- data_objs$X_for_ics
    
    if (is.null(dataX) || is.null(chr.length) || is.null(X_for_ics)) {
      return(NULL)
    }
    
    index.first <- input$First.Index
    index.second <- input$Second.Index
    genome <- input$Chromosomes
    n <- input$Cutoff
    
    # Get user-defined parameters for ics.outlier
    level_test <- input$level_test
    m_dist <- input$m_dist
    
    # Perform ICS analysis
    Z <- ics2(X_for_ics)
    Z_dist <- ics.distances(Z, index = index.first:index.second)
    Z_frame <- data.frame(Z_dist)
    
    Z2 <- data.frame(
      Chr = dataX$Chrom,
      Midpoint = dataX$Midpoint,
      ICS.distance = Z_frame$Z_dist
    )
    
    Z2$Log.ICS <- log10(Z2$ICS.distance)
    
    # Calculate cutoff for visualization (still useful for the plot)
    cutoff.applied <- quantile(Z2$Log.ICS, probs = 1 - n / 100)
    
    # NEW: Use the ICSOutlier package's ics.outlier function with user-selected test
    # This replaces the previous quantile-based outlier detection
    # Now using user-defined parameters
    outlier_test <- input$outlier_test  # Get the selected test method
    icsOutlierAG <- ics.outlier(Z, test = outlier_test, level.dist = 0.00001, 
                              level.test = level_test, mDist = m_dist)
    Z2$Outlier <- icsOutlierAG@outliers
    
    # Prepare data for visualization
    Z2$Position <- NA
    if (any(Z2$Chr == 1)) {
      Z2$Position[Z2$Chr == 1] <- Z2$Midpoint[Z2$Chr == 1]
    }
    
    for (j in 2:genome) {
      if (any(Z2$Chr == j)) {
        prev_chr_length <- sum(chr.length$V2[1:(j - 1)])
        Z2$Position[Z2$Chr == j] <- Z2$Midpoint[Z2$Chr == j] + prev_chr_length
      }
    }
    
    # Color coding for plot visualization
    Z2$Color <- "grey20"
    for (i in seq(1, nrow(chr.length), 2)) {
      Z2$Color[Z2$Chr == i] <- "grey58"
    }
    # Here we're coloring outliers detected by the ics.outlier function
    Z2$Color[Z2$Outlier == 1] <- "darkred"
    Z2$Color <- factor(Z2$Color)
    
    list(Z2 = Z2, cutoff.applied = cutoff.applied, icsOutlierAG = icsOutlierAG)
  })
  
  # Reactive value for Test of Interest data
  test_of_interest_data <- reactive({
    data_objs <- data_objects()
    dataX <- data_objs$dataX
    if(is.null(dataX)){
      return(NULL)
    }
    TestX <- input$TestOfInterest
    TestXF <- dataX[, 2 + TestX]
    
    # Use proper tidyverse syntax
    data.1.name <- colnames(dataX)[3:ncol(dataX)]
    data.1.melt <- pivot_longer(dataX, cols = all_of(data.1.name),
                              names_to = "variable", values_to = "value")
    
    list(TestXF = TestXF, data.1.melt = data.1.melt)
  })
  
  # Output: Genome Test plot
  output$GenomeTest <- renderPlot({
    data_objs <- data_objects()
    dataX <- data_objs$dataX
    test_data <- test_of_interest_data()
    if(is.null(dataX) || is.null(test_data)) {
      return(NULL)
    }
    TestXF <- test_data$TestXF
    
    ggplot(data = dataX, aes(x = Midpoint, y = TestXF)) +
      geom_point(size = 1.0, aes(colour = Chrom)) +
      theme(
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = rel(1.1), hjust = 1, vjust = 1),
        axis.text.y = element_text(size = rel(1.3), color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      ) +
      ylab(expression(test_statistic)) +
      xlab(expression(Position(bp))) +
      guides(colour = FALSE) +
      ggtitle("Test of interest")
  })
  
  # Output: Genome Melt plot
  output$GenomeMelt <- renderPlot({
    data_objs <- data_objects()
    dataX <- data_objs$dataX
    test_data <- test_of_interest_data()
    if(is.null(dataX) || is.null(test_data)) {
      return(NULL)
    }
    data.1.melt <- test_data$data.1.melt
    
    ggplot(data = data.1.melt, aes(x = Midpoint, y = value)) +
      geom_point(size = 1.0, aes(colour = Chrom)) +
      theme(
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = rel(1.1), hjust = 1, vjust = 1),
        axis.text.y = element_text(size = rel(1.3), color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      ) +
      ylab(expression(test_statistic)) +
      xlab(expression(Position(bp))) +
      guides(colour = FALSE) +
      ggtitle("Individual Selection Scans") +
      facet_grid(variable ~ ., scales = "free")
  })
  
  # Output: ICS histogram
  output$ICS.hist <- renderPlot({
    ics_results <- ics_analysis_results()
    if (is.null(ics_results)) {
      return(NULL)
    }
    Z2 <- ics_results$Z2
    cutoff.applied <- ics_results$cutoff.applied
    
    # Get the current test being used for display in subtitle
    test_method <- input$outlier_test
    test_name <- ifelse(test_method == "anscombe", "Anscombe", "Jarque-Bera")
    
    # Highlight points classified as outliers by ics.outlier
    ggplot(data = Z2, aes(x = Log.ICS, fill = factor(Outlier))) +
      geom_histogram(binwidth = 0.05, alpha = 0.7) +
      scale_fill_manual(values = c("grey50", "darkred"), 
                        labels = c("Normal", "Outlier"),
                        name = "Classification") +
      geom_vline(xintercept = cutoff.applied, color = "blue", linetype = "dashed") +
      theme_minimal() +
      labs(
        title = "ICS Distance Distribution",
        subtitle = paste0("Red bars: outliers detected by ics.outlier (", test_name, " test)"),
        x = "Log10(ICS Distance)",
        y = "Count"
      )
  })
  
  # Output: Chromosome-specific ICS histogram
  output$ICS.chromosome.hist <- renderPlot({
    ics_results <- ics_analysis_results()
    if (is.null(ics_results)) {
      return(NULL)
    }
    Z2 <- ics_results$Z2
    cutoff.applied <- ics_results$cutoff.applied
    
    # Get the current test being used for display in subtitle
    test_method <- input$outlier_test
    test_name <- ifelse(test_method == "anscombe", "Anscombe", "Jarque-Bera")
    
    chr.interest <- input$Chromosomes
    Z2.chr <- subset(Z2, Chr == chr.interest)
    
    ggplot(data = Z2.chr, aes(x = Log.ICS, fill = factor(Outlier))) +
      geom_histogram(binwidth = 0.1, alpha = 0.7) +
      scale_fill_manual(values = c("grey50", "darkred"), 
                        labels = c("Normal", "Outlier"),
                        name = "Classification") +
      geom_vline(xintercept = cutoff.applied, color = "blue", linetype = "dashed") +
      theme_minimal() +
      labs(
        title = paste("ICS Distance Distribution - Chromosome", chr.interest),
        subtitle = paste0("Red bars: outliers detected by ics.outlier (", test_name, " test)"),
        x = "Log10(ICS Distance)",
        y = "Count"
      )
  })
  
  # Output: Genome scan plot
  output$GenomeScan <- renderPlot({
    ics_results <- ics_analysis_results()
    if (is.null(ics_results)) {
      return(NULL)
    }
    Z2 <- ics_results$Z2
    cutoff.applied <- ics_results$cutoff.applied
    
    ggplot(data = Z2, aes(x = Position, y = Log.ICS, color = Color)) +
      geom_point(size = 1.0) +
      scale_color_manual(values = c("darkred", "grey20", "grey58")) +
      theme(
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = rel(1.1), hjust = 1, vjust = 1),
        axis.text.y = element_text(size = rel(1.3), color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      ) +
      ylab(expression(Log[10] ~ ICS ~ Distance)) +
      xlab(expression(Position(bp))) +
      guides(colour = FALSE) +
      ggtitle("ICS genome scan") +
      # Show the quantile cutoff as a reference (blue dashed line)
      geom_hline(yintercept = cutoff.applied, color = "blue", linetype = "dashed")
  })
  
  # Output: ICS summary
  output$summary.ICS <- renderPrint({
    ics_results <- ics_analysis_results()
    if (is.null(ics_results)) {
      return(NULL)
    }
    Z2 <- ics_results$Z2
    summary(Z2$Log.ICS)
  })
  
  # Output: ICSOutlier summary
  output$summary.ICSOutlier <- renderPrint({
    ics_results <- ics_analysis_results()
    if (is.null(ics_results)) {
      return(NULL)
    }
    icsOutlierAG <- ics_results$icsOutlierAG
    summary(icsOutlierAG)
  })
  
  # Output: Single Test Summary
  output$summary.SingleTest <- renderPrint({
    data_objs <- data_objects()
    dataX <- data_objs$dataX
    if (is.null(dataX)) {
      return(NULL)
    }
    summary(dataX[, -1:-2])
  })
  
  # Output: Number of outliers text
  output$nText <- renderText({
    ics_results <- ics_analysis_results()
    if (is.null(ics_results)) {
      return(NULL)
    }
    Z2 <- ics_results$Z2
    paste("Number of Outliers:", sum(Z2$Outlier))
  })
  
  # Download handlers
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("ICS_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      ics_results <- ics_analysis_results()
      if (!is.null(ics_results)) {
        write.csv(ics_results$Z2, file, row.names = FALSE)
      }
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("ICS_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ics_results <- ics_analysis_results()
      if (!is.null(ics_results)) {
        Z2 <- ics_results$Z2
        cutoff.applied <- ics_results$cutoff.applied
        
        p <- ggplot(data = Z2, aes(x = Position, y = Log.ICS, color = Color)) +
          geom_point(size = 1.0) +
          scale_color_manual(values = c("darkred", "grey20", "grey58")) +
          theme(
            plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(hjust = 0.5, size = 14),
            axis.text.x = element_text(size = rel(1.1), hjust = 1, vjust = 1),
            axis.text.y = element_text(size = rel(1.3), color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()
          ) +
          ylab(expression(Log[10] ~ ICS ~ Distance)) +
          xlab(expression(Position(bp))) +
          guides(colour = FALSE) +
          ggtitle("ICS genome scan") +
          geom_hline(yintercept = cutoff.applied, color = "blue", linetype = "dashed")
        
        ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
      }
    }
  )
}
