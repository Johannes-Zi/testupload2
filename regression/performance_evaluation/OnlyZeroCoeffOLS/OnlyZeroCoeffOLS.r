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
  lambda_min_index <- which(c == lambda_min_value) # nolint: object_usage_linter.

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


#' Import and preprocess Performance_Evaluation data
#' 
#' This function reads the input file, preprocesses the dataframe, and selects entries based on Pearson based feature selection.
#' 
#' @param input_file The path to the input file.
#' 
#' @return The preprocessed and selected dataframe.
#' 
#' @examples
#' df_preprocessed_2 <- import_and_preprocess_data(input_directory_path)
import_and_preprocess_performance_evaluation_data <- function(input_file) {
    # Read the input file into a dataframe
    df <- read.table(input_file, header = TRUE, sep = "\t")
    print(head(df))
    
    # Preprocess the dataframe
    df_preprocessed_1 <- preprocess_performance_evaluation_df(df)
    print(head(df_preprocessed_1))
    
    # Select entries based on Pearson based feature selection
    df_preprocessed_2 <- df_preprocessed_1[df_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
    print(head(df_preprocessed_2))
    
    return(df_preprocessed_2)
}


#' Import and preprocess regression models output data
#'
#' This function imports and preprocesses regression output data from a specified directory.
#' 
#' @param input_directory_path The path to the directory containing the regression output files.
#' @param import_limit The maximum number of files to import.
#' @param file_pattern The pattern to match the regression output files.
#'
#' @return A data frame containing the preprocessed regression output data.
#'
#' @examples
#' import_and_preprocess_regression_output_data("path/to/directory", 10, "*.csv")
import_and_preprocess_regression_output_data <- function(input_directory_path, import_limit, file_pattern) {
    cat("imported directory path:\n", input_directory_path, "\n")

    # Process the directory and load the results - only for Pearson based feature (segment) preselection input
    output_df <- process_directory(input_directory_path, import_limit, file_pattern)
    print(head(output_df))
    
    # Second preprocessing step before merging the dataframes
    processed_regression_models_df <- process_regression_models_df(output_df)
    print(head(processed_regression_models_df))
    return(processed_regression_models_df)
}


#' Combine performance evaluation and regression data
#'
#' This function merges the performance evaluation dataframe and the regression models dataframe based on the gene name and segmentation preselection correlation.
#' It filters out entries that have a segmentation preselection correlation of "Pearson".
#'
#' @param performance_evaluation_df The performance evaluation dataframe.
#' @param regression_models_df The regression models dataframe.
#'
#' @return The merged and filtered dataframe.
#'
#' @examples
#' performance_df <- data.frame(gene_name = c("Gene1", "Gene2", "Gene3"),
#'                              segmentation_preselection_corellation = c("Pearson", "Spearman", "Pearson"))
#' regression_df <- data.frame(gene_name = c("Gene1", "Gene2", "Gene3"),
#'                             regression_model = c("Model1", "Model2", "Model3"))
#' combine_performance_evaluation_and_regression_data(performance_df, regression_df)
#'
#' @export
combine_performance_evaluation_and_regression_data <- function(performance_evaluation_df, regression_models_df) {
    # Merge the performance evaluation and regression models dataframes
    combined_df <- merge(performance_evaluation_df, regression_models_df, by = c("gene_name", "segmentation_preselection_corellation"), all = TRUE)
    
    # Filter out Spearman based entries
    combined_filtered_df <- combined_df[combined_df$segmentation_preselection_corellation == "Pearson", ]
    
    return(combined_filtered_df)
}


#' Process a BED file and extract relevant information.
#'
#' This function reads a BED file, filters out the lines that start with "chr",
#' and extracts relevant information such as the gene name, number of features,
#' number of zero coefficients, and whether there are only zero coefficients.
#'
#' @param file_path The path to the BED file.
#'
#' @return A data frame containing the extracted information.
#'
#' @examples
#' process_bed_file("path/to/bed_file.bed")
#'
#' @export
process_bed_file <- function(file_path) {
    
    # Read the file line by line
    lines <- readLines(file_path)

    # Filter out the lines that start with "chr"
    chr_lines <- lines[grep("^chr", lines)]

    # Set custom headers for the data frame
    headers <- c("chr", "start", "end", "coeffs", "idk")

    # Convert the filtered lines back into a data frame
    temp_df <- read.delim(text = chr_lines, header = FALSE, col.names = headers)

    # Get the number of lines in the temp_df dataframe
    num_lines <- nrow(temp_df)

    # Get the nzero value from the temp_df dataframe
    nzero_value <- sum(temp_df$coeffs != 0)

    # Extract the gene name from the filename
    filename <- tools::file_path_sans_ext(basename(file_path))
    filename_parts <- strsplit(as.character(filename), "_")  # Split the Sample_Name column by "_"
    extracted_gene_name <- sapply(filename_parts, function(x) paste(x[4], x[5], sep = "_"))

    # Create a new row for the dataframe with the extracted information
    result_row <- data.frame(gene_name = extracted_gene_name,
                            number_of_features_ols = num_lines,
                            nzero_ols = nzero_value,
                            only_zero_coefficients_ols = sum(!(num_lines == nzero_value)))
    return(result_row)
}


#' Process OLS Data . bed files in a given directory.
#'
#' This function processes the data from .bed files in a given directory.
#' It imports a specified number of files based on a given file pattern,
#' and stores the results in a data frame.
#'
#' @param directory_path The path to the directory containing the .bed files.
#' @param import_limit The maximum number of files to import.
#' @param file_pattern The pattern to match the file names.
#'
#' @return A data frame containing the processed results, including the
#'         filename, number of features, number of zero coefficients, and
#'         a logical value indicating if there are only zero coefficients.
#'
#' @examples
#' process_ols_data("path/to/directory", 10, "*.bed")
#'
process_ols_data <- function(directory_path, import_limit, file_pattern) {
  # Get the list of .bed files in the directory
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
                           number_of_features_ols = numeric(),
                           nzero_ols = numeric(),
                           only_zero_coefficients_ols = numeric())

  # Process each file and append the results to the data frame
  for (i in seq_along(file_list)) {
    file_path <- file_list[i]
    result <- process_bed_file(file_path)
    result_df <- rbind(result_df, result)
    
    # Print the current number of processed files every 50 files
    if (i %% 50 == 0) {
      print(paste("Processed", i, "files"))
    }
  }

  return(result_df)
}


#' Combine OLS and elastic net data
#'
#' This function merges the performance evaluation and regression models dataframes.
#'
#' @param ols_df The dataframe containing the OLS regression model results.
#' @param elnet_df The dataframe containing the Elastic Net regression model results.
#'
#' @return The merged dataframe containing the combined data from both regression models.
#'
#' @examples
#' ols_data <- read.csv("ols_results.csv")
#' elnet_data <- read.csv("elnet_results.csv")
#' combined_data <- combine_ols_and_elnet_data(ols_data, elnet_data)
#' print(combined_data)
#'
combine_ols_and_elnet_data <- function(ols_df, elnet_df) {
  # Merge the performance evaluation and regression models dataframes
  combined_df <- merge(elnet_df, ols_df, by = c("gene_name"), all = TRUE)
    
  return(combined_df)
}


#' Main workflow 
# Set the input directory paths for standard and LeaveOneOut regression performance evaluation data
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"

# Import and preprocess performance evaluation data for standard regression
df_standard_preprocessed_2 <- import_and_preprocess_performance_evaluation_data(input_directory_path_standard)
# Import and preprocess performance evaluation data for LeaveOneOut regression
df_LeaveOneOut_preprocessed_2 <- import_and_preprocess_performance_evaluation_data(input_directory_path_LeaveOneOut)


# Set the input directory paths for standard and LeaveOneOut regression models
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/regression_output/regression_output/"
#input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/regression_output_example/"
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/regression_output/"

# Import and preprocess regression output data for standard regression
processed_regression_models_df_standard <- import_and_preprocess_regression_output_data(input_directory_path_standard, 50000, "\\Pearson.RData$")
# Import and preprocess regression output data for LeaveOneOut regression
processed_regression_models_df_LeaveOneOut <- import_and_preprocess_regression_output_data(input_directory_path_LeaveOneOut, 50000, "\\Pearson.RData$")

# Combine performance evaluation and regression data for standard regression
combined_elnet_df_standard <- combine_performance_evaluation_and_regression_data(df_standard_preprocessed_2, processed_regression_models_df_standard)
# Combine performance evaluation and regression data for LeaveOneOut regression
combined_elnet_df_LeaveOneOut <- combine_performance_evaluation_and_regression_data(df_LeaveOneOut_preprocessed_2, processed_regression_models_df_LeaveOneOut)

print(head(combined_elnet_df_standard))
print(head(combined_elnet_df_LeaveOneOut))

# Load OLS data in form of .bed files
ols_df_standard <- process_ols_data(input_directory_path_standard, 50000, "\\Pearson.bed$")
ols_df_LeaveOneOut <- process_ols_data(input_directory_path_LeaveOneOut, 50000, "\\Pearson.bed$")
print(head(ols_df_standard))
print(head(ols_df_LeaveOneOut))

final_df_standard <- combine_ols_and_elnet_data(ols_df_standard, combined_elnet_df_standard)
final_df_LeaveOneOut <- combine_ols_and_elnet_data(ols_df_LeaveOneOut, combined_elnet_df_LeaveOneOut)
print(head(final_df_standard))
print(head(final_df_LeaveOneOut))



# Save final_df as a TSV file
#write.table(final_df_standard, file = "file.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Select standrd entries that are potentially problematic
selected_weird_entries_standard <- subset(final_df_standard, only_zero_coefficients_ols == 1 & pVal < 0.05 & qVal < 0.05)
print(nrow(selected_weird_entries_standard))
# Save weird OLS Models as file
write.table(selected_weird_entries_standard, file = "selected_weird_entries_standard.csv", sep = ",", quote = FALSE, row.names = FALSE)


# NOT POSSIBLE - NO pVal column in LeaveOneOut
# # Select leaveOneOut entries that are potentially problematic
# selected_weird_entries_leaveOneOut <- subset(final_df_LeaveOneOut, only_zero_coefficients_ols == 1 & pVal < 0.05)
# print(nrow(selected_weird_entries_leaveOneOut))

# # Save weird OLS Models as file
# write.table(selected_weird_entries_leaveOneOut, file = "selected_weird_entries_leaveOneOut.csv", sep = ",", quote = FALSE, row.names = FALSE)
