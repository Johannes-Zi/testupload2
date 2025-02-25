# Load necessary libraries
library(ggplot2)  # For creating plots and visualizations
library(RColorBrewer)  # For color palettes in plots
library(rstudioapi)  # For interacting with RStudio IDE
library(dplyr)  # For data manipulation and transformation
library(MASS)  # For kernel density estimation
library(viridis)  # For color scales in plots

# Import Performance_Overview files
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"
# Create dataframes
df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
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

  # Extract the segmentation preselection correlation type
  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))
  #head(segmentation_preselection_corellation)  # Print the first few segmentation preselection correlations

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


#' Create a dotplot of Pearson vs PearsonVar
#'
#' This function creates a dotplot of Pearson correlation coefficient vs Pearson Variance (PearsonVar) using the provided dataframe.
#' It filters out outliers based on the interquartile range (IQR) (plots in a range IQR +- 1.5 * IQR)
#' and calculates the density of points using kernel density estimation.
#' The resulting dotplot is saved as a PNG file in the specified output directory.
#'
#' @param df_standard The dataframe containing the Pearson and PearsonVar columns.
#' @param output_directory The directory where the dotplot PNG file will be saved.
#'
#' @return The dotplot object.
#'
#' @examples
#' df <- data.frame(Pearson = c(0.5, 0.6, 0.7), PearsonVar = c(0.1, 0.2, 0.3))
#' create_dotplot(df, output_directory = "C:/output")
#'
#' @seealso \code{\link{ggplot2::geom_point}}, \code{\link{MASS::kde2d}}, \code{\link{ggplot2::ggsave}}
#'
#' @export
create_dotplot <- function(df_standard, output_directory) {
    # Extract x axis names
    df_standard_name <- strsplit(deparse(substitute(df_standard)), "_")[[1]][2] # nolint
    
    # Create new data frames that only contain the target columns
    df_standard_new <- data.frame(Pearson = df_standard$Pearson, PearsonVar = df_standard$PearsonVar)

    # Calculate IQR for Pearson and PearsonVar
    Q1_Pearson <- quantile(df_standard_new$Pearson, 0.25)
    Q3_Pearson <- quantile(df_standard_new$Pearson, 0.75)
    IQR_Pearson <- IQR(df_standard_new$Pearson)

    Q1_PearsonVar <- quantile(df_standard_new$PearsonVar, 0.25)
    Q3_PearsonVar <- quantile(df_standard_new$PearsonVar, 0.75)
    IQR_PearsonVar <- IQR(df_standard_new$PearsonVar)

    # Filter out the outliers
    df_standard_new <- df_standard_new[!(df_standard_new$Pearson < (Q1_Pearson - 1.5 * IQR_Pearson) | 
                                                                             df_standard_new$Pearson > (Q3_Pearson + 1.5 * IQR_Pearson) |
                                                                             df_standard_new$PearsonVar < (Q1_PearsonVar - 1.5 * IQR_PearsonVar) |
                                                                             df_standard_new$PearsonVar > (Q3_PearsonVar + 1.5 * IQR_PearsonVar)), ]
        
    # Calculate density for each point
    # n specifies the number of points to evaluate the density with/ defines number of grid cells
    df_standard_new$density <- get_density(df_standard_new$Pearson, df_standard_new$PearsonVar, n = 1000)

    # Create dotplot
    dotplot <- ggplot(df_standard_new, aes(x = Pearson, y = PearsonVar, color = density)) + # nolint
        geom_point() +
        scale_color_viridis() +
        labs(title = "Pearson vs PearsonVar", x = "Pearson", y = "PearsonVar") +
        theme_gray() +
        theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                    axis.title = element_text(size = 14, face = "bold"),
                    axis.text = element_text(size = 12, face = "bold"))
    
    # Save the plot to the target directory
    ggsave(filename = paste(output_directory, "/dotplot.png", sep = ""), plot = dotplot)
    
    # Return the dotplot
    return(dotplot)
}

# Set the output path to the directory of the active document
output_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(output_path)

# Preprocess the standard CV dataframe
df_standard_preprocessed_1 <- preprocess_performance_evaluation_df(df_standard)
head(df_standard_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
df_standard_preprocessed_2 <- df_standard_preprocessed_1[df_standard_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
head(df_standard_preprocessed_2)

# Call the function with the dataframes as arguments
create_dotplot(df_standard_preprocessed_2, output_directory = output_path)
 