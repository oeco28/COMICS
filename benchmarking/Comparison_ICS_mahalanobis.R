library("data.table");
library("dplyr");
library("tidyverse");
library("Matrix");
library("MASS");
library("ggplot2");
library("snpStats");
library("fields");
library("RColorBrewer");
library("extrafont");
library("ICSOutlier");
library("gridExtra");

#############################
#Data preprocessing
#############################

setwd("")
 
load("TwoRefSim.rda")

data2 <- data.frame(scale(na.omit(TwoRefSim[,c(8:10)]),center = TRUE, scale = TRUE))

data2$s_high <- TwoRefSim[rownames(na.omit(TwoRefSim[,c(8:10)])),4]

########################
# From Minotaur
########################

S <- stats::cov(data2[,2:4], use="pairwise.complete.obs")
S_inv <- solve(S)
M <- colMeans(data2[,2:4],na.rm=TRUE)

M <- as.vector(unlist(M))
output <- list(S=S, S_inv=S_inv, M=M)

Mahalanobis <- function(dfv, column.nums=1:ncol(dfv), subset=1:nrow(dfv), S=NULL, M=NULL){
 
#### perform simple checks on data
   #dfv_check <- data_checks(dfv, column.nums, subset, S, M, check.na=FALSE, check.M=TRUE)

# extract variables from dfv and dfv_check
   diff <- as.matrix(dfv[,column.nums,drop=FALSE])
   S_inv <- S
   M <- M
 
   # calculate Mahalanobis distance
   for (i in 1:ncol(diff)) {
     diff[,i] <- diff[,i] - M[i]
   }
   distance <- Mod(sqrt(as.complex(rowSums((diff %*% S_inv) * diff, na.rm = TRUE))))
 
   return(distance)
 } # end Mahalanobis


maha_1 <- Mahalanobis(data2,column.nums=2:4,S=S, M=M)

Q95 <- quantile(maha_1, 0.95)
outliers <- which(maha_1>Q95)

data2$Outlier_maha_mino <- rep(0,dim(data2)[1])

data2$Outlier_maha_mino[which(maha_1>Q95)] <- 1

table(data2$Outlier_maha_mino)


#######################################
# OE Cornejo ICSOutlier identification
#######################################


Z <- ics2(data2[2:4],S1 = MeanCov, S2 = Mean3Cov4)
Z_comp <- comp.norm.test(Z)
Z_dist <- ics.distances(Z)
icsOutlierJB <- ics.outlier(Z, test = "jarque", level.dist = 0.000005,level.test = 0.01, mDist = mean(Z_dist))
summary(icsOutlierJB)
icsOutlierAG <- ics.outlier(Z, test = "anscombe", level.dist = 0.00001, level.test = 0.01, mDist = 100)
summary(icsOutlierAG)

data2$Outlier_ics <- icsOutlierAG@outliers

data3 <- data2[1:4]

n <- 0.01
data3$Log.ICS <- log10(Z_dist)
cutoff.applied <- quantile(data3$Log.ICS, probs = 1 - n / 100)
data3$Outlier <- 0
data3$Outlier[data3$Log.ICS >= cutoff.applied] <- 1

par(mfrow=c(2,1))
plot(icsOutlierJB)
plot(icsOutlierAG)

Z2 <- cbind(rownames(data2),attributes(icsOutlierAG)$outliers,attributes(icsOutlierAG)$ics.distances)
colnames(Z2) <- c("Position","outlier","ics.distance")


data2$true_selection <- rep(0,dim(data2)[1])
data2$true_selection[which(data2$s_high != 0)] <- 1

scaled_distances <- Z_dist / max(Z_dist)
cat("Scaled ICS Distances (first 10):\n", scaled_distances[1:10], "\n\n")

threshold_ics <- icsOutlierAG@ics.dist.cutoff
distance_from_threshold <- icsOutlierAG@ics.distances - threshold
sigmoid <- function(x, k = 0.5) { 1 / (1 + exp(-k * x)) }
pseudo_probabilities_threshold <- sigmoid(distance_from_threshold)
cat("Pseudo-probabilities (vs. threshold, first 10):\n", pseudo_probabilities_threshold[1:10], "\n\n")

ecdf_distances <- ecdf(icsOutlierAG@ics.distances)
pseudo_probabilities_ecdf <- 1 - ecdf_distances(icsOutlierAG@ics.distances)
cat("Pseudo-probabilities (based on ECDF, first 10):\n", pseudo_probabilities_ecdf[1:10], "\n")





#######################
#' Comparison of performance
########################





#' Function to calculate precision, recall, F1-score, specificity, sensitivity, and accuracy'
#''
#' @param y_true Vector of true labels (0 and 1).'
#' @param y_pred Vector of predicted labels (0 and 1).'
#''
#' @return A list containing the calculated metrics. Returns NA for'
#'         any metric where the calculation would result in division by zero.'
#'         Handles edge cases where TP+FP, TP+FN, TN+FP, or TN+FN are zero.'
calculate_metrics2 <- function(y_true, y_pred) {
  # Check for input validity
  if (length(y_true) != length(y_pred)) {
    stop("Shape mismatch: y_true and y_pred must have the same length.")
  }
  if (length(y_true) == 0) {
    return(list(
      precision = NA_real_,
      recall = NA_real_,
      f1_score = NA_real_,
      specificity = NA_real_,
      sensitivity = NA_real_,
      accuracy = NA_real_
    ))
  }
  if (!all(y_true %in% c(0, 1)) || !all(y_pred %in% c(0, 1))) {
    stop("Inputs must contain only binary values (0 and 1).")
  }

  # Calculate confusion matrix values
  true_positives <- sum((y_true == 1) & (y_pred == 1))
  true_negatives <- sum((y_true == 0) & (y_pred == 0))
  false_positives <- sum((y_true == 0) & (y_pred == 1))
  false_negatives <- sum((y_true == 1) & (y_pred == 0))

  # Calculate metrics, handling potential division by zero
  precision_denominator <- true_positives + false_positives
  precision <- ifelse(precision_denominator > 0, true_positives / precision_denominator, NA_real_)

  recall_denominator <- true_positives + false_negatives
  recall <- ifelse(recall_denominator > 0, true_positives / recall_denominator, NA_real_)

  f1_denominator <- 2 * true_positives + false_positives + false_negatives
  f1_score <- ifelse(f1_denominator > 0, (2 * true_positives) / f1_denominator, NA_real_)

  specificity_denominator <- true_negatives + false_positives
  specificity <- ifelse(specificity_denominator > 0, true_negatives / specificity_denominator, NA_real_)

  sensitivity_denominator <- true_positives + false_negatives
  sensitivity <- ifelse(sensitivity_denominator > 0, true_positives / sensitivity_denominator, NA_real_)

  accuracy_denominator <- length(y_true)
  accuracy <- ifelse(accuracy_denominator > 0, (true_positives + true_negatives) / accuracy_denominator, NA_real_)

  return(list(
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    specificity = specificity,
    sensitivity = sensitivity,
    accuracy = accuracy
  ))
}

#' Function to plot the metrics'
#''
#' @param metrics A list of calculated metrics (output of calculate_metrics).'
#' @param title A title for the plot.'
#''
#' @return A ggplot object.'
plot_metrics <- function(metrics, title = "Classification Metrics") {
  # Convert metrics to a data frame for plotting
  metrics_df <- data.frame(
    metric = names(metrics),
    value = unlist(metrics)
  )

  # Remove NA values for plotting
  metrics_df <- metrics_df %>%
    dplyr::filter(!is.na(value)) # Use dplyr::filter for clarity

  # Create the bar plot
  p <- ggplot2::ggplot(metrics_df, ggplot2::aes(x = metric, y = value, fill = metric)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.4f", value)), vjust = -0.3, size = 3) + # add the value labels
    ggplot2::ylim(0, 1) + # set y axis limit from 0 to 1
    ggplot2::labs(
      title = title,
      x = "Metric",
      y = "Value"
    ) +
    ggplot2::theme_minimal() +
    scale_fill_manual(values=c("#042a2b","#5EB1BF","#0BC9CD","#CDEDF6","#EF7B45","#D84727")) +
    ggplot2::theme(ggplot2::legend.position = "none") # Remove legend

  return(p)
}


results_maha1 <- calculate_metrics2(data2$true_selection, data2$Outlier_maha_mino)
    print("Mahalanobis Results:")
    print(results_maha1)
    plot_maha1 <- plot_metrics(results_maha1, "Mahalanobis Results")
    print(plot_maha1)



results_ics <- calculate_metrics2(data2$true_selection, data2$Outlier_ics)
    print("ICS Results:")
    print(results_ics)
    plot_ics <- plot_metrics(results_ics, "Invariant Coordinate Selection (ICS)")
    print(plot_ics)

grid.arrange(plot_maha1, plot_ics)
