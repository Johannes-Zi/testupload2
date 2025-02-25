# Load required libraries
library(MASS)
library(ggplot2)
library(viridis)
library(glmnet)
library(ggExtra)


# Function to process each file
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

# Define the input directory path
#input_directory_path <- "C:/Users/johan/Desktop/standard_regression/regression_output/regression_output/"
input_directory_path <- "C:/Users/johan/Desktop/temp_inputfile_for_MA_plots/RegFeaturesXnzero/regression_output/regression_output/"


cat("imported directory path:\n", input_directory_path, "\n")

# Process the directory and get the results - Limit 9k - not more than 9k files
output_df <- process_directory(input_directory_path, 9050, "\\Pearson.RData$")
#output_df <- process_directory(input_directory_path, 20000, "\\Spearman.RData$")

# Print the resulting data frame
print(head(output_df))

p <- ggplot(output_df, aes(x = number_of_features, y = nzero_value)) +
  geom_count(aes(color = ..n..)) +
  scale_color_viridis_c() +
  labs(x = "number of input features", y = "number of features utilized", color = "Number of models", size = "") +
  ggtitle(expression(bold("Partion of model input features utilized for predictions"))) +
  theme_bw() +
  theme(
    # Make all text black
    text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "bottom" 
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = seq(0, 45, by = 10), limits = c(0, 45)) +
  scale_y_continuous(breaks = seq(0, 45, by = 10), limits = c(0, 45)) +
  coord_cartesian(xlim = c(0, 45), ylim = c(0, 45)) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  geom_jitter(data = subset(output_df, number_of_features > 45 | nzero_value > 45), 
              aes(x = pmin(number_of_features, 45), y = pmin(nzero_value, 45), shape = "Outliers"), 
              color = "#000000", fill = "#beefec", size = 2, width = 0, height = 0.5) +
  scale_shape_manual(name = "", values = c("Outliers" = 23))+
  guides(size = guide_legend(order = 2), shape = guide_legend(order = 1))


# Print cwd
print(getwd())

# directory of the script
script_dir <- dirname(sys.frame(1)$ofile)

# Create plot filepath with script directory
plot_filepath <- file.path(script_dir, "RegFeaturesXnzero")

# Save the plot as a SVG file
ggsave(paste(plot_filepath, ".svg", sep = ""), plot = p, width = 7, height = 6)

# Save the plot as a PNG file
ggsave(paste(plot_filepath, ".png", sep = ""), plot = p, width = 7, height = 6, dpi = 250)



density_plot_1 <- ggplot(output_df, aes(x = number_of_features)) +
  geom_density(fill = NA, color = "black", size = 1) +
  labs(title = "Density Distribution of Number of Features",
       x = "Number of Features",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  scale_x_continuous(limits = c(0, 45))


# Save the plot as a PNG file
ggsave(paste(plot_filepath, "_density_1.png", sep = ""), plot = density_plot_1, width = 11, height = 2, dpi = 200)
# Save the plot as a SVG file
ggsave(paste(plot_filepath, "_density_1.svg", sep = ""), plot = density_plot_1, width = 11, height = 2)

density_plot_2 <- ggplot(output_df, aes(x = nzero_value)) +
  geom_density(fill = NA, color = "black", size = 1.1) +
  labs(title = "Density Distribution of Number of Features",
       x = "Number of Features",
       y = "Density") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  ) +
  scale_x_continuous(limits = c(0, 45))


# Save the plot as a PNG file
ggsave(paste(plot_filepath, "_density_2.png", sep = ""), plot = density_plot_2, width = 11, height = 2, dpi = 200)
# Save the plot as a SVG file
ggsave(paste(plot_filepath, "_density_2.svg", sep = ""), plot = density_plot_2, width = 8, height = 2)
