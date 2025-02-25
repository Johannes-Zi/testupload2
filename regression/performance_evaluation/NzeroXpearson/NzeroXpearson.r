# Load required libraries
library(dplyr)  # For data manipulation
library(MASS)  # For kernel density estimation
library(ggplot2)  # For data visualization
library(viridis)  # For color palettes


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

  # Extract the segmentation preselection correlation type
  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))

  # Create a new column in the dataframe for the segmentation preselection correlation type
  adapted_df <- input_df
  adapted_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

  # Extract the entry specific gene identifiers
  extracted_gene_names <- sapply(filename_parts, function(x) paste(x[2], x[3], sep = "_"))
  adapted_df$Sample_Name <- extracted_gene_names

  # Rename the Sample_Name column to gene_name (identifiers)
  adapted_df <- rename(adapted_df, gene_name = Sample_Name) # nolint: object_usage_linter.

  return(adapted_df)
}


#' Extract relevant elastic net information from an .RData file
#'
#' This function loads an .RData file, extracts relevant information such as the filename,
#' the number of features, and the nzero value based on the elasticnet lambda.min model stored 
#' in the file.
#'
#' @param file_path The path to the .RData file.
#'
#' @return A data frame with the extracted information of the lamda.min model.
#'
#' @examples
#' process_file("/path/to/file.RData")
#'
process_file <- function(file_path) {
  # Load the .RData file - accessible as 'elasticnet_model'
  load(file_path)

  # Get the filename without extension
  filename <- tools::file_path_sans_ext(basename(file_path))

  # Get ident of lambda min based on the list with all lambdas
  lambda_min_value <- elasticnet_model$model$lambda.min # nolint: object_usage_linter.
  lambda_min_index <- which(elasticnet_model$model$lambda == lambda_min_value) # nolint: object_usage_linter.

  # Get the corresponding nzero value
  nzero_value <- elasticnet_model$model$nzero[lambda_min_index] # nolint: object_usage_linter.
  number_of_features <- elasticnet_model$model$glmnet.fit$dim[1] # nolint: object_usage_linter.

  # Create a data frame with the results
  result <- data.frame(filename = filename,
                       number_of_features = number_of_features,
                       nzero_value = nzero_value)
  return(result)
}


#' Process the output directory of regression in form of .RData files
#'
#' This function takes a directory path, imports .RData files from that directory, and processes each file to extract relevant information.
#' The function limits the number of files to import based on the import_limit parameter and uses the file_pattern parameter to filter the files.
#' The processed information is stored in a data frame with columns for the filename, number of features, and nzero value.
#'
#' @param directory_path The path to the directory containing the .RData files.
#' @param import_limit The maximum number of files to import.
#' @param file_pattern The pattern to match the .RData files.
#'
#' @return A data frame with the processed information from the .RData files.
#'
#' @examples
#' process_directory("/path/to/directory", 100, "\\Pattern.RData$")
#'
process_directory <- function(directory_path, import_limit, file_pattern) {
  # Get the list of .RData files in the directory
  file_list <- list.files(directory_path, pattern = file_pattern,
                          full.names = TRUE)
  print(paste("Number of detected files:", length(file_list)))
  print(paste("Used filename pattern", file_pattern))

  # Limit the number of files to import
  if (length(file_list) > import_limit) {
    file_list <- file_list[1:import_limit]
    print(paste("Number of imported files limited to", import_limit))
  }

  # Initialize an empty data frame to store the results
  result_df <- data.frame(filename = character(),
                          number_of_features = numeric(),
                          nzero_value = numeric(),
                          stringsAsFactors = FALSE)

  # Process each file and append the results to the data frame
  for (i in seq_along(file_list)) {
    file_path <- file_list[i]
    result <- process_file(file_path)
    result_df <- rbind(result_df, result)
    
    # Print the current number of processed files every 50 files
    if (i %% 50 == 0) {
      print(paste("Processed", i, "files"))
    }
  }

  return(result_df)
}


#' Preprocess Regression Models DataFrame
#'
#' This function takes an input dataframe and performs several operations to process it.
#' It creates a new dataframe using the filename as a key, merges it with the input dataframe,
#' adapts the filename column, extracts segmentation preselection correlation type,
#' extracts entry specific gene identifiers, renames the filename column to gene_name,
#' and returns the processed dataframe.
#'
#' @param input_df The input dataframe to be processed.
#'
#' @return The processed dataframe.
#'
#' @examples
#' input_df <- data.frame(filename = c("file1", "file2", "file3"),
#'                        value = c(1, 2, 3))
#' processed_df <- process_regression_models_df(input_df)
#' print(processed_df)
#'
#' @export
process_regression_models_df <- function(input_df) {
  # Create a new dataframe using the filename as key - extracting only the filename
  new_df <- data.frame(filename = unique(input_df$filename))

  # Merge the new dataframe with the input_df dataframe
  merged_df <- merge(new_df, input_df, by = "filename", all.x = TRUE)

  # Adapt the filename column
  filename_parts <- strsplit(as.character(merged_df$filename), "_")

  # Extract the segmentation preselection correlation type
  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(x[7], sep = "_"))

  merged_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

  # Extract the entry specific gene identifiers
  extracted_gene_names <- sapply(filename_parts, function(x) paste(x[5], x[6], sep = "_"))

  merged_df$filename <- extracted_gene_names

  # Rename the filename column to gene_name (identifiers)
  merged_df <- rename(merged_df, gene_name = filename)

  return(merged_df)
}


#' Calculate the density of points in a two-dimensional space
#'
#' This function calculates the density of points in a two-dimensional space using kernel density estimation.
#'
#' @param x A numeric vector representing the x-coordinates of the points.
#' @param y A numeric vector representing the y-coordinates of the points.
#' @param ... Additional arguments to be passed to the `kde2d` function.
#'
#' @return A numeric vector representing the density of points at each (x, y) coordinate.
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' density <- get_density(x, y)
#' plot(x, y, col = rgb(0, 0, 0, density), pch = 16)
#'
#' @seealso \code{\link{MASS::kde2d}}
#'
#' @export
get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}


#' Create a dotplot to visualize the relationship between the number of nonzero feature coefficients and Pearson values.
#'
#' This function takes a data frame `df` containing the necessary columns `nzero_value` and `Pearson`, and an `output_path` where the dotplot image will be saved.
#' It calculates the density for each point in the data, creates a dotplot using ggplot2, and saves the dotplot as a PNG image.
#'
#' @param df A data frame containing the necessary columns `nzero_value` and `Pearson`.
#' @param output_path The path where the dotplot image will be saved.
#' @return None
#' @examples
#' df <- data.frame(nzero_value = c(1, 2, 3), Pearson = c(0.5, 0.8, 0.9))
#' create_dotplot(df, "output/")
create_dotplot <- function(df, output_path) {
  # Calculate the density for each point in your data
  df$density <- get_density(df$nzero_value, df$Pearson, n = 1000)

  # Create the dotplot
  dotplot <- ggplot(df, aes(x = nzero_value, y = Pearson, color = density)) +
    geom_point(size = 2) +
    theme_bw() +
    #labs(title = "Utilized features vs. respective CV Pearson coefficients",
    labs(title = "",
       x = "Number of utilized features",
       y = "Cross validation Pearson coefficient") +
    scale_color_viridis() +
    scale_x_continuous(breaks = seq(min(df$nzero_value), max(df$nzero_value), by = 2)) +  # Increase the number of breaks on the x-axis
    coord_cartesian(xlim = c(0, 29)) +  # Limit the x-axis to 0 - 30
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Add title and adjust size
      axis.title = element_text(size = 16),
      axis.title.y = element_text(margin = margin(r = 10,)),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.text = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  # Save the dotplot
  ggsave(filename = paste(output_path, "/dotplot.png", sep = ""), plot = dotplot, width = 9, height = 6)
  # Save svg dotplot
  ggsave(filename = paste(output_path, "/dotplot.svg", sep = ""), plot = dotplot, width = 9, height = 6)
}


#'
#' Import and preprocess Performance_Evaluation data
#'
# Import Performance_Overview files as dataframes
input_directory_path_LeaveOneOut <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/leaveOneOut_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/standard_regression/performance_evaluation/Performance_Overview.txt"

# Create dataframes
df_LeaveOneOut <- read.table(input_directory_path_LeaveOneOut, header = TRUE, sep = "\t")
#df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
head(df_LeaveOneOut)
#head(df_standard)

# Preprocess the standard CV dataframe
#df_standard_preprocessed_1 <- preprocess_performance_evaluation_df(df_standard)
#head(df_standard_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
#df_standard_preprocessed_2 <- df_standard_preprocessed_1[df_standard_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
#head(df_standard_preprocessed_2)
# Preprocess the LeaveOneOut CV dataframe
df_LeaveOneOut_preprocessed_1 <- preprocess_performance_evaluation_df(df_standard)
head(df_LeaveOneOut_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
df_LeaveOneOut_preprocessed_2 <- df_LeaveOneOut_preprocessed_1[df_LeaveOneOut_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
head(df_LeaveOneOut_preprocessed_2)


#'
#' Import and process the regression models (regression output)
#'
# Define the input directory paths
# For regression standard CV
#input_directory_path_standard <- "C:/Users/johan/Desktop/temp_inputfile_for_MA_plots/NzeroXpearson/standard_CV/regression_output/regression_output/"
#cat("imported directory path:\n", input_directory_path_standard, "\n")
# Process the directory and load the results - only for Pearson based feature (segment) preselection input
#output_df_standard <- process_directory(input_directory_path_standard, 20000, "\\Pearson.RData$")
#head(output_df_standard) # Crashes sometimes at first iteration and runs only when data is already loaded in cache
#  For regression LeaveOneOut CV
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/temp_inputfile_for_MA_plots/RegFeaturesXnzero/regression_output/regression_output/"
cat("imported directory path:\n", input_directory_path_LeaveOneOut, "\n")
# Process the directory and load the results - only for Pearson based feature (segment) preselection input
output_df_LeaveOneOut <- process_directory(input_directory_path_LeaveOneOut, 20000, "\\Pearson.RData$")
head(output_df_LeaveOneOut) # Crashes sometimes at first iteration and runs only when data is already loaded in cache

# Second preprocessing step before merging the two dataframes
#processed_regression_models_df_standard <- process_regression_models_df(output_df_standard)
#head(processed_regression_models_df_standard)
processed_regression_models_df_LeaveOneOut <- process_regression_models_df(output_df_LeaveOneOut)
head(processed_regression_models_df_LeaveOneOut)

# Combine the two dataframes
# For regression standard CV
#combined_df_standard <- merge(df_standard_preprocessed_2, processed_regression_models_df_standard , by = c("gene_name", "segmentation_preselection_corellation"), all = TRUE)
# Filter out Spearman based entries
#combined_filtered_df_standard <- combined_df_standard[combined_df_standard$segmentation_preselection_corellation == "Pearson", ]
#head(combined_filtered_df_standard)
#  For regression LeaveOneOut CV
combined_df_LeaveOneOut <- merge(df_LeaveOneOut_preprocessed_2, processed_regression_models_df_LeaveOneOut , by = c("gene_name", "segmentation_preselection_corellation"), all = TRUE)
# Filter out Spearman based entries
combined_filtered_df_LeaveOneOut <- combined_df_LeaveOneOut[combined_df_LeaveOneOut$segmentation_preselection_corellation == "Pearson", ]
head(combined_filtered_df_LeaveOneOut)


#'
#' Create the dpotplot for standard CV dataframe
#'
# Set the output path to the directory of the active document
#output_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
#print(output_path)
# Created dotplot 
#create_dotplot(combined_filtered_df_standard, output_path)


#'
#' Exclude rows with missing values from LeaveOneOut CV dataframe and create dotplot
#'

# Filter the data frame to only show rows with missing values in the nzero_value or Pearson columns
# Filters out datapoints where the regression model was created due to error and thus resulted in missin entries after merge
missing_values_df_LeaveOneOut <- combined_filtered_df_LeaveOneOut[is.na(combined_filtered_df_LeaveOneOut$nzero_value) | is.na(combined_filtered_df_LeaveOneOut$Pearson), ]
head(missing_values_df_LeaveOneOut)
# Remove rows with missing values in the nzero_value or Pearson columns
combined_filtered2_df_LeaveOneOut <- combined_filtered_df_LeaveOneOut[!is.na(combined_filtered_df_LeaveOneOut$nzero_value) & !is.na(combined_filtered_df_LeaveOneOut$Pearson), ]
head(combined_filtered2_df_LeaveOneOut)

create_dotplot(combined_filtered2_df_LeaveOneOut, output_path)
