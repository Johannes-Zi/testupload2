# This script iterates over the elastic net models, that were trained during the
# regression step of STICHIT. It count the number of features (segments, that
# are putative regulatory elements) that were used for the model and count how
# many of those features have coefficients, that are 0 (and thus meaningless).

# Define Input parameters
input_file_path <- paste(
  "C:/Users/johan/Desktop/standard_regression/regression_output/regression_output/",
  "Elasticnet_Regression_Model_Segmentation_ENSG00000001167_10_Spearman.RData",
  sep = ""
)

# Echo path of imported file
cat(sprintf("File path for import:\n%s\n", input_file_path))

# Load the .RData file - accesible as 'elasticnet_model'
#load(file.choose())
load(input_file_path)

# Print structure of the loaded object
#str(elasticnet_model)
# Print keys that are accessible in a specif layer of the data structure
#names(elasticnet_model$model)

# insight into the elasticnet_model object
lambda_min_value <- elasticnet_model$model$lambda.min
cat("lambda min of the model:", lambda_min_value)

# Get ident of lambda min based on the list with all lambdas
lambda_min_index <- which(elasticnet_model$model$lambda == lambda_min_value)
cat("lambda min index:", lambda_min_index)

# Get the corresponding nzero value
nzero_value <- elasticnet_model$model$nzero[lambda_min_index]
number_of_features <- elasticnet_model$model$glmnet.fit$dim[1]

# Print the calculated values after each step
print("overview for lambda.min")
cat("Number of features used in the model:", number_of_features, "\n")
cat("Number of features with coefficients unequal to 0:", nzero_value, "\n")

