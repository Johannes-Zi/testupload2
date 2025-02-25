library(ggplot2)  # Required for creating boxplots
library(RColorBrewer)  # Required for color palettes in boxplots
library(rstudioapi)  # Required for getting the active document context
library(dplyr)  # Required for data manipulation operations

# Import Performance_Overview files
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

# New location
input_directory_path_LeaveOneOut <- "c:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/leaveOneOut_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "c:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/standard_regression/performance_evaluation/Performance_Overview.txt"

# Create dataframes
df_LeaveOneOut <- read.table(input_directory_path_LeaveOneOut, header = TRUE, sep = "\t")
df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
print("LeaveOneOut")
print(head(df_LeaveOneOut))
print("Standard")
print(tail(df_standard))


#' Preprocesses the performance evaluation dataframe
#'
#' This function takes an input dataframe and performs preprocessing steps to adapt the dataframe for further analysis.
#' It extracts information from the Sample_Name column, creates a new column for segmentation preselection correlation,
#' extracts gene names, and renames the Sample_Name column to gene_name.
#'
#' @param input_df The input dataframe to be preprocessed
#'
#' @return The preprocessed dataframe
#'
#' @examples
#' input_df <- data.frame(Sample_Name = c("gene1_condition1_sample1", "gene2_condition2_sample2"))
#' preprocess_performance_evaluation_df(input_df)
preprocess_performance_evaluation_df <- function(input_df) {
  # Extracts the original filenames for each df entry based on the Sample_Name column
  filename_parts <- strsplit(as.character(input_df$Sample_Name), "_")  # Split the Sample_Name column by "_"
  #head(filename_parts)  # Print the first few filename parts

  # Extract the segmentation preselection correlation type for each df entry
  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))
  #head(segmentation_preselection_corellation)  # Print the first few segmentation preselection correlations

  # Create a new column in the dataframe for the segmentation preselection correlation type
  adapted_df <- input_df
  adapted_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

  # Extract the entry specific gene names
  extracted_gene_names <- sapply(filename_parts, function(x) paste(x[2], x[3], sep = "_"))
  adapted_df$Sample_Name <- extracted_gene_names

  # Rename the Sample_Name column to gene_name
  adapted_df <- rename(adapted_df, gene_name = Sample_Name) # nolint: object_usage_linter.

  return(adapted_df)
}


#' Create boxplots for performance evaluation
#'
#' This function creates boxplots to compare the performance of two datasets. 
#'
#' @param df_LeaveOneOut The first dataset (Leave-One-Out cross-validation).
#' @param df_standard The second dataset (Standard cross-validation).
#' @param target_column The name of the target column to be compared.
#' @param output_directory The directory where the boxplot image will be saved.
#'
#' @return The created boxplot as a ggplot object.
#'
#' @details The size of the plot is scaled based on the range of the target column values. 
#' This scaling ensures that the boxplots are visible and outliers are not excluded from displaying. 
#' The y-axis limits are calculated based on the lower and upper whiskers of both datasets, 
#' and the plot is scaled by multiplying the limits by 1.1. This scaling factor ensures that 
#' the boxplots are not cut off and outliers are visible.
#'
#' @examples
#' # Example usage:
#' create_boxplots(df_LeaveOneOut, df_standard, "target_column", "output_directory")
create_boxplots <- function(df_LeaveOneOut, df_standard, target_column, output_directory) {
  # Extract x axis names
  df_LeaveOneOut_name <- strsplit(deparse(substitute(df_LeaveOneOut)), "_")[[1]][2]
  df_standard_name <- strsplit(deparse(substitute(df_standard)), "_")[[1]][2]

  # Create new data frames that only contain the target column and the new CV_Type column
  df_LeaveOneOut_new <- data.frame(CV_Type = df_LeaveOneOut_name, target = df_LeaveOneOut[[target_column]])
  df_standard_new <- data.frame(CV_Type = df_standard_name, target = df_standard[[target_column]])

  # Compute lower and upper whiskers for both data frames
  ylim1 = boxplot.stats(df_LeaveOneOut[[target_column]])$stats[c(1, 5)]
  ylim2 = boxplot.stats(df_standard[[target_column]])$stats[c(1, 5)]

  # Use the minimum lower whisker and maximum upper whisker as the y limits
  ylim = c(min(ylim1[1], ylim2[1]), max(ylim1[2], ylim2[2]))

  # Combine the new data frames into one
  df_combined <- rbind(df_LeaveOneOut_new, df_standard_new)

  print(head(df_combined))

  df_combined$CV_Type <- recode(df_combined$CV_Type, "LeaveOneOut" = "LOOCV")
  df_filtered <- df_combined %>% filter(target < quantile(target, 0.99)) 

  print("Combined data frame")
  print(tail(df_combined))
  # Create boxplot
  boxplot <- ggplot(df_filtered, aes(y = CV_Type, x = target, fill = CV_Type)) + # Swap x and y
    geom_boxplot(outlier.shape = NA) +
    labs(title = "Cross-Validation Techniques MSE values distribution", x = target_column, y = "Cross Validation Type") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12, angle = 30)) +
    coord_cartesian(xlim = c(0.1, 2)) + 
    scale_fill_manual(values = c("LOOCV" = "#00aeba", "standard" = "#e7b800")) +
    theme(legend.position = "none") # Remove legend

  # Save the plot to the target directory
  ggsave(filename = paste(output_directory, "/", target_column, "_boxplot.svg", sep = ""), plot = boxplot, width = 7, height = 4)
  ggsave(filename = paste(output_directory, "/", target_column, "_boxplot.png", sep = ""), plot = boxplot, width = 7, height = 4)

  # Create violin plot
  violin_plot <- ggplot(df_filtered, aes(y = CV_Type, x = target, fill = CV_Type)) + # Swap x and y
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +  # Use violin plot with full range and consistent width
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 1, fill = "white") +  # Add boxplot inside violin plot for better visualization
    labs(title = "Cross-Validation techniques MSE values distribution", x = "Mean squared error", y = "Cross validation type") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 15, color = "black"),
          axis.text = element_text(color = "black", size = 14),
          axis.text.y = element_text(angle = 30)) +
    coord_cartesian(xlim = c(0, 3)) + 
    scale_fill_manual(values = c("LOOCV" = "#00aeba", "standard" = "#e7b800")) +
    theme(legend.position = "none") # Remove legend
  
  # Save the plot to the target directory
  ggsave(filename = paste(output_directory, "/", target_column, "_violin_plot.svg", sep = ""), plot = violin_plot, width = 7, height = 4)
  ggsave(filename = paste(output_directory, "/", target_column, "_violin_plot.png", sep = ""), plot = violin_plot, width = 7, height = 4)
      
  return(boxplot)
}


# Set the output path to the directory of the active document
output_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(output_path)

# Preprocess the LeaveOneOut CV dataframe
df_LeaveOneOut_preprocessed_1 <- preprocess_performance_evaluation_df(df_LeaveOneOut)
head(df_LeaveOneOut_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
df_LeaveOneOut_preprocessed_2 <- df_LeaveOneOut_preprocessed_1[df_LeaveOneOut_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
print("LeaveOneOut")
print(head(df_LeaveOneOut_preprocessed_2))

# Preprocess the standard CV dataframe
df_standard_preprocessed_1 <- preprocess_performance_evaluation_df(df_standard)
head(df_standard_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
df_standard_preprocessed_2 <- df_standard_preprocessed_1[df_standard_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
print("Standard")
print(head(df_standard_preprocessed_2))

# Call the function with the dataframes as arguments
create_boxplots(df_LeaveOneOut_preprocessed_2, df_standard_preprocessed_2, target_column = "MSE", output_directory = output_path)
