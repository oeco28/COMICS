#' Server Logic for COMICS Shiny Application
#' UI Components for COMICS Shiny Application
#'
#' This file defines the server logic for the COMICS (Combined Outlier Method for
#' This file defines the user interface for the COMICS (Combined Outlier Method for
#' Identifying Candidate Signals) Shiny application.
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#'
#' @importFrom shiny renderPlot renderPrint downloadHandler observe
#' @importFrom ggplot2 ggplot aes geom_point geom_histogram geom_vline scale_color_manual theme element_text theme_minimal labs ggtitle xlab ylab geom_hline guides ggsave
#' @importFrom tidyr pivot_longer all_of
#' @importFrom ICS ics2 ics.distances ics
#' @importFrom tidyr pivot_longer all_of
#' @importFrom ICS ics2 ics.distances ics
#' @importFrom shiny fluidPage titlePanel sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny fileInput numericInput selectInput downloadButton
#' @importFrom shiny h4 h3 verbatimTextOutput tabsetPanel tabPanel plotOutput
#' @importFrom bslib bs_theme
#' @importFrom shinyjs useShinyjs
#'
#' @noRd
#' This the most finished version up to 06/2025 - Modified with SD-based ICS

server <- function(input, output, session) {
  
  # Helper function to safely extract ICS scores
  extract_ics_scores <- function(ics_obj) {
    # Try different methods to extract scores from ICS object
    possible_slots <- c("scores", "Scores", "components", "Components", "X", "data")
    for (slot_name in possible_slots) {
      tryCatch({
        result <- slot(ics_obj, slot_name)
        if (is.matrix(result) || is.data.frame(result)) {
          return(as.matrix(result))
        }
      }, error = function(e) {
        # Continue to next slot
      })
    }
    
    # Try matrix multiplication approach
    tryCatch({
      if (hasSlot(ics_obj, "W") && hasSlot(ics_obj, "X")) {
        W_matrix <- slot(ics_obj, "W")
        X_data <- slot(ics_obj, "X")
        result <- as.matrix(X_data) %*% as.matrix(W_matrix)
        return(result)
      }
    }, error = function(e) {
      # Continue
    })
    
    # Try coef() method
    tryCatch({
      result <- coef(ics_obj)
      if (is.matrix(result) || is.data.frame(result)) {
        return(as.matrix(result))
      }
    }, error = function(e) {
      # Continue
    })
    
    stop("Could not extract scores from ICS object")
  }
  

  
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
        X_standardized = NULL
      ))
    }
    
    # Read data
    dataX = read.table(file1$datapath, header = TRUE)
    chr.length = read.table(configuration$datapath, header = TRUE)
    
    # ICS Analysis preparation
    X_for_ics <- dataX
    X_for_ics[, 1:2] <- NULL # Remove the first two columns
    
    # Standardize the data (center and scale)
    X_standardized <- as.data.frame(scale(X_for_ics, center = TRUE, scale = TRUE))
    
    list(
      dataX = dataX,
      chr.length = chr.length,
      X_for_ics = X_for_ics,
      X_standardized = X_standardized
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
  
  # Reactive value for ICS analysis with multiple scatter combinations
  ics_analysis_results <- reactive({
    data_objs <- data_objects()
    dataX <- data_objs$dataX
    chr.length <- data_objs$chr.length
    X_standardized <- data_objs$X_standardized
    
    if (is.null(dataX) || is.null(chr.length) || is.null(X_standardized)) {
      return(NULL)
    }
    
    genome <- input$Chromosomes
    scatter_method <- input$scatter_method
    
    # Perform ICS analysis with selected scatter matrices
    tryCatch({
      if (scatter_method == "cov_cov4") {
        Z <- ics(X_standardized, S1 = cov, S2 = cov4)
      } else if (scatter_method == "cov_cov3") {
        Z <- ics(X_standardized, S1 = cov, S2 = function(x) cov(x^3))
      } else if (scatter_method == "cov_cor") {
        Z <- ics(X_standardized, S1 = cov, S2 = cor)
      } else if (scatter_method == "cor_cov4") {
        Z <- ics(X_standardized, S1 = cor, S2 = cov4)
      } else {
        # Default to cov_cov4
        Z <- ics(X_standardized, S1 = cov, S2 = cov4)
      }
      
      # Extract scores using the safe function
      scores <- extract_ics_scores(Z)
      
      # Calculate distances using the selected method
      if (input$distance_method == "euclidean") {
        distances <- sqrt(rowSums(scores^2))
      } else if (input$distance_method == "component_wise") {
        # Use the most variable component
        component_vars <- apply(scores, 2, var)
        best_component <- which.max(component_vars)
        distances <- abs(scores[, best_component])
      } else if (input$distance_method == "multi_component") {
        # Use top 2 components
        component_vars <- apply(scores, 2, var)
        n_components <- min(2, ncol(scores))
        top_components <- order(component_vars, decreasing = TRUE)[1:n_components]
        selected_components <- scores[, top_components, drop = FALSE]
        
        # Use Euclidean distance in selected component space
        component_center <- colMeans(selected_components)
        distances <- sqrt(rowSums((selected_components - matrix(component_center,
                                                               nrow = nrow(selected_components),
                                                               ncol = length(component_center),
                                                               byrow = TRUE))^2))
      } else {
        # Default to euclidean
        distances <- sqrt(rowSums(scores^2))
      }
      
      # Apply standard deviation-based threshold
      mean_dist <- mean(distances)
      sd_dist <- sd(distances)
      sd_multiplier <- input$sd_multiplier
      cutoff_value <- mean_dist + sd_multiplier * sd_dist
      
      # Create results dataframe
      Z2 <- data.frame(
        Chr = dataX$Chrom,
        Midpoint = dataX$Midpoint,
        ICS.distance = distances
      )
      
      Z2$Log.ICS <- log10(Z2$ICS.distance + 1e-10)  # Add small constant to avoid log(0)
      
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
      
      # Apply standard deviation-based cutoff
      Z2$Outlier <- as.numeric(Z2$ICS.distance > cutoff_value)
      
      Z2$Color <- "grey20"
      for (i in seq(1, nrow(chr.length), 2)) {
        Z2$Color[Z2$Chr == i] <- "grey58"
      }
      Z2$Color[Z2$Outlier == 1] <- "darkred"
      Z2$Color <- factor(Z2$Color)
      
      # Calculate additional statistics
      n_outliers <- sum(Z2$Outlier)
      outlier_percentage <- round(n_outliers / nrow(Z2) * 100, 2)
      
      list(
        Z2 = Z2,
        cutoff_applied = cutoff_value,
        log_cutoff = log10(cutoff_value + 1e-10),
        scores = scores,
        distances = distances,
        mean_dist = mean_dist,
        sd_dist = sd_dist,
        n_outliers = n_outliers,
        outlier_percentage = outlier_percentage,
        scatter_method = scatter_method,
        distance_method = input$distance_method,
        sd_multiplier = sd_multiplier
      )
      
    }, error = function(e) {
      return(list(error = paste("ICS analysis failed:", e$message)))
    })
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
    
    data.1.name <- colnames(dataX)[3:ncol(dataX)]
    data.1.melt <- pivot_longer(dataX, cols = all_of(data.1.name),
                                names_to = "variable", values_to = "value")
    
    list(TestXF = TestXF, data.1.melt = data.1.melt)
  })
  
  # Output: Genome Test plot
  output$GenomeTest <- renderPlot({
    data_objs <- data_objects()
    dataX <- data_objs$dataX;
    test_data <- test_of_interest_data()
    TestXF <- test_data$TestXF
    if(is.null(dataX)){
      return(NULL)
    }
    
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
    data.1.melt <- test_data$data.1.melt
    if(is.null(dataX)){
      return(NULL)
    }
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
    if (is.null(ics_results) || !is.null(ics_results$error)) {
      return(NULL)
    }
    
    distances <- ics_results$distances
    cutoff_applied <- ics_results$cutoff_applied
    mean_dist <- ics_results$mean_dist
    sd_dist <- ics_results$sd_dist
    
    hist_data <- data.frame(distances = distances)
    
    ggplot(data = hist_data, aes(x = distances)) +
      geom_histogram(binwidth = (max(distances) - min(distances))/30,
                     fill = "lightblue", alpha = 0.7) +
      geom_vline(xintercept = cutoff_applied, color = "red", linewidth = 1.2) +
      geom_vline(xintercept = mean_dist, color = "blue", linetype = "dashed") +
      geom_vline(xintercept = mean_dist + sd_dist, color = "green", linetype = "dotted") +
      theme_minimal() +
      labs(
        title = paste("ICS Distance Distribution (SD-based method)"),
        subtitle = paste("Method:", ics_results$scatter_method, "|",
                        ics_results$distance_method, "| Multiplier:", ics_results$sd_multiplier),
        x = "ICS Distance",
        y = "Count"
      ) +
      annotate("text", x = cutoff_applied, y = Inf,
               label = paste("Cutoff =", round(cutoff_applied, 3)),
               vjust = 2, color = "red")
  })
  
  # Output: Chromosome-specific ICS histogram
  output$ICS.chromosome.hist <- renderPlot({
    ics_results <- ics_analysis_results()
    if (is.null(ics_results) || !is.null(ics_results$error)) {
      return(NULL)
    }
    
    Z2 <- ics_results$Z2
    cutoff_applied <- ics_results$cutoff_applied
    chr.interest <- input$Chromosomes
    
    Z2.chr <- subset(Z2, Chr == chr.interest, select = c(ICS.distance))
    
    ggplot(data = Z2.chr, aes(x = ICS.distance)) +
      geom_histogram(binwidth = (max(Z2.chr$ICS.distance) - min(Z2.chr$ICS.distance))/20,
                     fill = "lightgreen", alpha = 0.7) +
      geom_vline(xintercept = cutoff_applied, color = "red", linewidth = 1.2) +
      theme_minimal() +
      labs(
        title = paste("ICS Distance Distribution - Chromosome", chr.interest),
        x = "ICS Distance",
        y = "Count"
      )
  })
  
  # Output: Genome scan plot
  output$GenomeScan <- renderPlot({
    ics_results <- ics_analysis_results()
    if (is.null(ics_results) || !is.null(ics_results$error)) {
      return(NULL)
    }
    
    Z2 <- ics_results$Z2
    log_cutoff <- ics_results$log_cutoff
    
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
      ggtitle(paste("ICS Genome Scan (SD-based:", ics_results$sd_multiplier, "σ)")) +
      geom_hline(yintercept = log_cutoff, color = "red", linetype = "dashed")
  })
  
  # Output: ICS summary
  output$summary.ICS <- renderPrint({
    ics_results <- ics_analysis_results()
    if (is.null(ics_results) || !is.null(ics_results$error)) {
      if (!is.null(ics_results$error)) {
        cat("Error in ICS analysis:", ics_results$error)
      }
      return(NULL)
    }
    
    cat("=== ICS Analysis Summary (Standard Deviation Method) ===\n")
    cat("Scatter matrices:", ics_results$scatter_method, "\n")
    cat("Distance method:", ics_results$distance_method, "\n")
    cat("SD multiplier:", ics_results$sd_multiplier, "\n\n")
    
    cat("Distance statistics:\n")
    cat("Mean:", round(ics_results$mean_dist, 4), "\n")
    cat("SD:", round(ics_results$sd_dist, 4), "\n")
    cat("Cutoff (μ + ", ics_results$sd_multiplier, "σ):", round(ics_results$cutoff_applied, 4), "\n\n")
    
    cat("Outlier detection:\n")
    cat("Number of outliers:", ics_results$n_outliers, "\n")
    cat("Percentage of outliers:", ics_results$outlier_percentage, "%\n\n")
    
    cat("Distance distribution:\n")
    print(summary(ics_results$distances))
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
  
  # Output: Component analysis
  output$component.analysis <- renderPrint({
    ics_results <- ics_analysis_results()
    if (is.null(ics_results) || !is.null(ics_results$error)) {
      return(NULL)
    }
    
    scores <- ics_results$scores
    if (is.null(scores)) return(NULL)
    
    component_vars <- apply(scores, 2, var)
    component_ranking <- order(component_vars, decreasing = TRUE)
    
    cat("=== ICS Component Analysis ===\n")
    cat("Component importance ranking:\n")
    for (i in 1:length(component_ranking)) {
      cat("Rank", i, ": Component", component_ranking[i],
          "- Variance:", round(component_vars[component_ranking[i]], 4), "\n")
    }
  })
  
  #output nText
  output$nText <- renderText({
    ics_results <- ics_analysis_results()
    if (is.null(ics_results) || !is.null(ics_results$error)) {
      return("No results available")
    }
    paste("Number of Outliers:", ics_results$n_outliers,
          "(", ics_results$outlier_percentage, "%)")
  })
  
  # Download handlers
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("ICS_SD_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      ics_results <- ics_analysis_results()
      if (!is.null(ics_results) && is.null(ics_results$error)) {
        # Add analysis parameters to the output
        output_data <- ics_results$Z2
        output_data$scatter_method <- ics_results$scatter_method
        output_data$distance_method <- ics_results$distance_method
        output_data$sd_multiplier <- ics_results$sd_multiplier
        output_data$cutoff_value <- ics_results$cutoff_applied
        write.csv(output_data, file, row.names = FALSE)
      }
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("ICS_SD_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ics_results <- ics_analysis_results()
      if (!is.null(ics_results) && is.null(ics_results$error)) {
        Z2 <- ics_results$Z2
        log_cutoff <- ics_results$log_cutoff
        
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
          ggtitle(paste("ICS Genome Scan (SD-based:", ics_results$sd_multiplier, "σ)")) +
          geom_hline(yintercept = log_cutoff, color = "red", linetype = "dashed")
        
        ggsave(file, plot = p, width = 12, height = 8, dpi = 300)
      }
    }
  )
}
