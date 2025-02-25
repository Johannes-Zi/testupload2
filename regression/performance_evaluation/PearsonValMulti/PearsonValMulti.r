library(ggplot2)
library(RColorBrewer)
library(rstudioapi)
library(dplyr)

# Import Performance_Overview files
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

# Create dataframes
df_LeaveOneOut <- read.table(input_directory_path_LeaveOneOut, header = TRUE, sep = "\t")
df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
head(df_LeaveOneOut)
head(df_standard)

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

  # Extract the segmentation preselection correlation
  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))
  #head(segmentation_preselection_corellation)  # Print the first few segmentation preselection correlations

  # Create a new column in the dataframe for the segmentation preselection correlation
  adapted_df <- input_df
  adapted_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

  # Extract the entry specific gene names
  extracted_gene_names <- sapply(filename_parts, function(x) paste(x[2], x[3], sep = "_"))
  adapted_df$Sample_Name <- extracted_gene_names

  # Rename the Sample_Name column to gene_name
  adapted_df <- rename(adapted_df, gene_name = Sample_Name) # nolint: object_usage_linter.

  return(adapted_df)
}


#' Create violin plots for performance evaluation
#'
#' This function creates violin plots to compare the performance of two datasets.
#'
#' @param df_LeaveOneOut The first dataset (Leave-One-Out cross-validation).
#' @param df_standard The second dataset (Standard cross-validation).
#' @param target_column The column name of the target variable.
#' @param output_directory The directory where the violin plots will be saved.
#'
#' @return The created violin plot.
#'
#' @examples
#' create_violin_plots(df_LeaveOneOut, df_standard, "target_variable", "output_directory")
create_violin_plots <- function(df_LeaveOneOut, df_standard, target_column, output_directory) {
  
  # Extract x axis names
  df_LeaveOneOut_name <- strsplit(deparse(substitute(df_LeaveOneOut)), "_")[[1]][2]
  df_standard_name <- strsplit(deparse(substitute(df_standard)), "_")[[1]][2]

  # Create new data frames that only contain the target column and the new CV_Type column
  df_LeaveOneOut_new <- data.frame(CV_Type = df_LeaveOneOut_name, target = df_LeaveOneOut[[target_column]])
  df_standard_new <- data.frame(CV_Type = df_standard_name, target = df_standard[[target_column]])
  print(head(df_LeaveOneOut_new))
  print(head(df_standard_new))

  # Combine the new data frames into one
  df_combined <- rbind(df_LeaveOneOut_new, df_standard_new)

  # Create violin plots
  violinplot <- ggplot(df_combined, aes(x = CV_Type, y = target, fill = CV_Type)) + # nolint: object_usage_linter.
    geom_violin(trim = FALSE, width = 0.5) +
    labs(title = paste(target_column, "across Cross-Validation Techniques", sep = " "), y = target_column, x = "Cross Validation Type") +
    theme_gray() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12, face = "bold")) +
    scale_fill_brewer(palette = "Set2") + # Apply Brewer color palette
    theme(legend.position = "none") # Remove legend

  # Save the plot to the target directory
  ggsave(filename = paste(output_directory, "/", target_column, "_violinplot.png", sep = ""), plot = violinplot)
  return(violinplot)
}


# Set the output path to the directory of the active document
output_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(output_path)

# Preprocess the LeaveOneOut CV dataframe
df_LeaveOneOut_preprocessed_1 <- preprocess_performance_evaluation_df(df_LeaveOneOut)
head(df_LeaveOneOut_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
df_LeaveOneOut_preprocessed_2 <- df_LeaveOneOut_preprocessed_1[df_LeaveOneOut_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
head(df_LeaveOneOut_preprocessed_2)

# Preprocess the standard CV dataframe
df_standard_preprocessed_1 <- preprocess_performance_evaluation_df(df_standard)
head(df_standard_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
df_standard_preprocessed_2 <- df_standard_preprocessed_1[df_standard_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
head(df_standard_preprocessed_2)

# Call the function with the dataframes as arguments
create_violin_plots(df_LeaveOneOut_preprocessed_2, df_standard_preprocessed_2, target_column = "Pearson", output_directory = output_path)
