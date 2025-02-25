library(ggplot2)
library(RColorBrewer)
library(rstudioapi)

# Import Performance_Overview files
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

# Create dataframes
df_LeaveOneOut <- read.table(input_directory_path_LeaveOneOut, header = TRUE, sep = "\t")
df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
head(df_LeaveOneOut)
head(df_standard)

create_boxplots <- function(df_LeaveOneOut, df_standard, target_column, output_directory) {
  
  # Extract x achsis names
  df_LeaveOneOut_name <- tail(strsplit(deparse(substitute(df_LeaveOneOut)), "_")[[1]], 1)
  df_standard_name <- tail(strsplit(deparse(substitute(df_standard)), "_")[[1]], 1)
  
  # Create new data frames that only contain the target column and the new CV_Type column
  df_LeaveOneOut_new <- data.frame(CV_Type = df_LeaveOneOut_name, target = df_LeaveOneOut[[target_column]])
  df_standard_new <- data.frame(CV_Type = df_standard_name, target = df_standard[[target_column]])
  
  # Compute lower and upper whiskers for both dataframes
  ylim1 = boxplot.stats(df_LeaveOneOut[[target_column]])$stats[c(1, 5)]
  ylim2 = boxplot.stats(df_standard[[target_column]])$stats[c(1, 5)]

  # Use the minimum lower whisker and maximum upper whisker as the y limits
  ylim = c(min(ylim1[1], ylim2[1]), max(ylim1[2], ylim2[2]))

  # Combine the new data frames into one
  df_combined <- rbind(df_LeaveOneOut_new, df_standard_new)
  
  # Create boxplot
  boxplot <- ggplot(df_combined, aes(x = CV_Type, y = target, fill = CV_Type)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = paste(target_column, "Comparison", sep = " "), y = target_column, x = "CV Type") +
    theme_gray() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) +
    coord_cartesian(ylim = ylim*1.1) + # Scale y limits based on ylim
    scale_fill_brewer(palette = "Set2") + # Apply Brewer color palette
    theme(legend.position = "none") # Remove legend

  # Save the plot to the target directory
  ggsave(filename = paste(output_directory, "/", target_column, "_violinplot.png", sep = ""), plot = boxplot)

  return(boxplot)
}

create_violin_plots <- function(df_LeaveOneOut, df_standard, target_column, output_directory) {
  
  # Extract x achsis names
  df_LeaveOneOut_name <- tail(strsplit(deparse(substitute(df_LeaveOneOut)), "_")[[1]], 1)
  df_standard_name <- tail(strsplit(deparse(substitute(df_standard)), "_")[[1]], 1)
  
  # Create new data frames that only contain the target column and the new CV_Type column
  df_LeaveOneOut_new <- data.frame(CV_Type = df_LeaveOneOut_name, target = df_LeaveOneOut[[target_column]])
  df_standard_new <- data.frame(CV_Type = df_standard_name, target = df_standard[[target_column]])

  # Combine the new data frames into one
  df_combined <- rbind(df_LeaveOneOut_new, df_standard_new)
  
  # Create violin plots
  violinplot <- ggplot(df_combined, aes(x = CV_Type, y = target, fill = CV_Type)) +
    geom_violin(trim = FALSE) +
    labs(title = paste(target_column, "Comparison", sep = " "), y = target_column, x = "CV Type") +
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

output_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(output_path)

# Call the function with the dataframes as arguments
create_violin_plots(df_LeaveOneOut, df_standard, target_column = "Pearson", output_directory = output_path)

create_violin_plots(df_LeaveOneOut, df_standard, target_column = "Spearman", output_directory = output_path)

create_boxplots(df_LeaveOneOut, df_standard, target_column = "MSE", output_directory = output_path)

#create_violin_plots(df_LeaveOneOut, df_standard, target_column = "pVal")

#create_violin_plots(df_LeaveOneOut, df_standard, target_column = "qVal")

